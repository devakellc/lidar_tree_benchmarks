#!/usr/bin/env Rscript
# #I3 model-runner contract. Invoke a GPU model arm as a SUBPROCESS against a
# frozen clip and read its detections CSV back as a contract-valid
# data.frame(x,y,z). Backends: a Python venv (TreeisoNet, #M7) and a Docker
# container (SAT #M6 / FF3D #M8 via run_docker_arm), both with the same crash
# discipline. (#19)
#
# Discipline: non-zero exit OR missing output OR a wrong-schema CSV -> NULL, so
# the equal-set guard drops the cell (never a fake 0-row from a stale/partial or
# malformed file). Only a VALID `x y z` CSV with no data rows -> a 0-row frame
# (legit ran-but-empty -> recall 0).
bs <- Find(file.exists, c(
  file.path("scripts", "bootstrap.R"),
  file.path("..", "..", "scripts", "bootstrap.R"),
  file.path(getwd(), "scripts", "bootstrap.R")))
if (!length(bs)) stop("bootstrap.R not found", call. = FALSE)
source(bs[1]); rm(bs)
source(.find("model_bench_lib.R"))

.log_python_failure <- function(label, out, st) {
  if (is.null(label) || !nzchar(label)) return(invisible(NULL))
  message(sprintf("  %s: GPU driver failed (status %s)", label,
                  if (is.null(st)) "?" else st))
  if (length(out)) message(paste(tail(out, 3), collapse = "\n"))
}

# Read an arm's detections CSV under the scorer contract: a valid `x y z` table
# -> data.frame(x,y,z) (a header with no data rows is a legit 0-row frame); a
# wrong-schema or unreadable file -> NULL. Callers map NULL -> skip the cell.
.read_detection_csv <- function(out_csv) {
  d <- tryCatch(read.table(out_csv, header = TRUE), error = function(e) NULL)
  if (is.null(d) || !identical(names(d), c("x", "y", "z"))) return(NULL)
  data.frame(x = as.numeric(d$x), y = as.numeric(d$y), z = as.numeric(d$z))
}

run_python_arm <- function(venv_python, script, input, out_csv,
                           extra = character(), timeout = 1800,
                           label = NULL) {
  if (file.exists(out_csv)) unlink(out_csv)          # never read a stale file
  # shQuote every arg: system2 pastes into a shell, so unquoted paths with
  # spaces/parens (e.g. the upstream "...(GPU4GB).json" config) would break sh.
  args <- shQuote(c(script, input, out_csv, extra))
  out <- tryCatch(suppressWarnings(system2(venv_python, args, stdout = TRUE,
                                           stderr = TRUE, timeout = timeout)),
                  error = function(e) NULL)
  st <- if (is.null(out)) 1L else attr(out, "status")   # non-NULL only when != 0
  if (!is.null(st) && st != 0) {
    .log_python_failure(label, out, st)
    return(NULL)                                     # crash -> skip cell
  }
  if (!file.exists(out_csv)) return(NULL)
  det <- .read_detection_csv(out_csv)                # NULL on wrong-schema/unreadable
  if (is.null(det)) return(NULL)
  assert_detection_contract(det)
  det
}

run_python_crown_arm <- function(venv_python, script, input, out_csv,
                                 extra = character(), timeout = 1800,
                                 label = NULL) {
  if (file.exists(out_csv)) unlink(out_csv)          # never read a stale file
  args <- shQuote(c(script, input, out_csv, extra))
  out <- tryCatch(suppressWarnings(system2(venv_python, args, stdout = TRUE,
                                           stderr = TRUE, timeout = timeout)),
                  error = function(e) NULL)
  st <- if (is.null(out)) 1L else attr(out, "status")
  if (!is.null(st) && st != 0) {
    .log_python_failure(label, out, st)
    return(NULL)
  }
  if (!file.exists(out_csv)) return(NULL)
  d <- tryCatch(read.table(out_csv, header = TRUE), error = function(e) NULL)
  if (is.null(d) || !identical(names(d), c("x", "y", "z", "crown_id")))
    return(NULL)
  pts <- data.frame(X = as.numeric(d$x), Y = as.numeric(d$y),
                    Z = as.numeric(d$z), crown_id = as.integer(d$crown_id))
  pts[is.finite(pts$X) & is.finite(pts$Y) & is.finite(pts$Z) &
        !is.na(pts$crown_id), , drop = FALSE]
}

# #19 Docker backend. Same contract as run_python_arm, but the arm runs inside a
# CUDA container (SAT #M6 / FF3D #M8). Identity bind-mounts (-v abs:abs) make
# in-container paths == host paths, so the model command receives real paths;
# `cmd` is required because Docker args after the image override an image CMD
# unless the image defines an ENTRYPOINT. `mounts` adds weight/config dirs,
# `gpus=NULL` drops --gpus. The runner OWNS the contract: non-zero exit or
# missing output -> NULL; the reader (default = the x y z CSV parse; pass
# read_instances_ply/read_instances_laz for a labeled-cloud output) is wrapped
# so a throwing or NULL-returning reader -> NULL (skip cell), and any non-NULL
# result is asserted against the detection contract.
run_docker_arm <- function(image, input, out_csv, extra = character(),
                           cmd = character(), mounts = NULL, gpus = "all",
                           docker = "docker", timeout = 1800, label = NULL,
                           reader = NULL) {
  if (is.null(cmd) || !length(cmd))
    stop("run_docker_arm: cmd must be supplied so input/out do not override an image CMD",
         call. = FALSE)
  if (file.exists(out_csv)) unlink(out_csv)          # never read a stale file
  in_abs  <- normalizePath(input,   mustWork = FALSE)
  out_abs <- normalizePath(out_csv, mustWork = FALSE)
  mdirs <- unique(c(normalizePath(dirname(in_abs),  mustWork = FALSE),
                    normalizePath(dirname(out_abs), mustWork = FALSE),
                    if (length(mounts)) normalizePath(mounts, mustWork = FALSE)))
  vol <- as.vector(rbind("-v", paste0(mdirs, ":", mdirs)))   # identity mounts
  gpu <- if (!is.null(gpus) && nzchar(gpus)) c("--gpus", gpus) else character()
  args <- c("run", "--rm", gpu, vol, image, cmd, in_abs, out_abs, extra)
  out <- tryCatch(suppressWarnings(system2(docker, shQuote(args), stdout = TRUE,
                                           stderr = TRUE, timeout = timeout)),
                  error = function(e) NULL)
  st <- if (is.null(out)) 1L else attr(out, "status")
  if (!is.null(st) && st != 0) {
    .log_python_failure(label, out, st)
    return(NULL)                                     # container crash -> skip cell
  }
  if (!file.exists(out_abs)) return(NULL)            # ran but produced nothing
  rd  <- if (is.null(reader)) .read_detection_csv else reader
  det <- tryCatch(rd(out_abs), error = function(e) NULL)
  if (is.null(det)) return(NULL)                     # schema failure -> skip cell
  assert_detection_contract(det)
  det
}
