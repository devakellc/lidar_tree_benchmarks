#!/usr/bin/env Rscript
# #I3 model-runner contract. Invoke a GPU model arm as a SUBPROCESS against a
# frozen clip and read its detections CSV back as a contract-valid
# data.frame(x,y,z). Backend = a Python venv (TreeisoNet, #M7); a Docker backend
# is added with #M6 (issue: GPU runner Docker + LAZ->PLY).
#
# Discipline: non-zero exit OR missing output OR a wrong-schema CSV -> NULL, so
# the equal-set guard drops the cell (never a fake 0-row from a stale/partial or
# malformed file). Only a VALID `x y z` CSV with no data rows -> a 0-row frame
# (legit ran-but-empty -> recall 0).
.script_path <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && length(ofile) && nzchar(ofile)) return(ofile)
  hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(hit)) return(sub("^--file=", "", hit[1]))
  NA_character_
}
.sp <- .script_path()
.ROOT <- if (!is.na(.sp)) normalizePath(file.path(dirname(.sp), ".."),
                                        mustWork = FALSE) else getwd()
.find <- function(rel) Find(file.exists, c(file.path(.ROOT, "scripts", rel),
                                           file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("model_bench_lib.R"))

run_python_arm <- function(venv_python, script, input, out_csv,
                           extra = character(), timeout = 1800) {
  if (file.exists(out_csv)) unlink(out_csv)          # never read a stale file
  # shQuote every arg: system2 pastes into a shell, so unquoted paths with
  # spaces/parens (e.g. the upstream "...(GPU4GB).json" config) would break sh.
  args <- shQuote(c(script, input, out_csv, extra))
  out <- tryCatch(suppressWarnings(system2(venv_python, args, stdout = TRUE,
                                           stderr = TRUE, timeout = timeout)),
                  error = function(e) NULL)
  st <- if (is.null(out)) 1L else attr(out, "status")   # non-NULL only when != 0
  if (!is.null(st) && st != 0) return(NULL)          # crash -> skip cell
  if (!file.exists(out_csv)) return(NULL)
  d <- tryCatch(read.table(out_csv, header = TRUE), error = function(e) NULL)
  if (is.null(d) || !identical(names(d), c("x", "y", "z")))
    return(NULL)                                     # wrong schema -> contract failure
  det <- data.frame(x = as.numeric(d$x), y = as.numeric(d$y), z = as.numeric(d$z))
  assert_detection_contract(det)
  det
}

run_python_crown_arm <- function(venv_python, script, input, out_csv,
                                 extra = character(), timeout = 1800) {
  if (file.exists(out_csv)) unlink(out_csv)          # never read a stale file
  args <- shQuote(c(script, input, out_csv, extra))
  out <- tryCatch(suppressWarnings(system2(venv_python, args, stdout = TRUE,
                                           stderr = TRUE, timeout = timeout)),
                  error = function(e) NULL)
  st <- if (is.null(out)) 1L else attr(out, "status")
  if (!is.null(st) && st != 0) return(NULL)
  if (!file.exists(out_csv)) return(NULL)
  d <- tryCatch(read.table(out_csv, header = TRUE), error = function(e) NULL)
  if (is.null(d) || !identical(names(d), c("x", "y", "z", "crown_id")))
    return(NULL)
  pts <- data.frame(X = as.numeric(d$x), Y = as.numeric(d$y),
                    Z = as.numeric(d$z), crown_id = as.integer(d$crown_id))
  pts[is.finite(pts$X) & is.finite(pts$Y) & is.finite(pts$Z) &
        !is.na(pts$crown_id), , drop = FALSE]
}
