# Repo-root resolution for Rscript entry points. Resolves scripts/ and gpu/
# paths regardless of CLAUDE_JOB_DIR / getwd().
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
