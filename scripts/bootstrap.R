# Load repo_paths.R (.ROOT + .find) into the caller's environment.
.load_repo <- function(envir = parent.frame()) {
  cands <- c(
    file.path("scripts", "repo_paths.R"),
    file.path("..", "..", "scripts", "repo_paths.R"),
    file.path(getwd(), "scripts", "repo_paths.R"))
  hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(hit))
    cands <- c(file.path(dirname(sub("^--file=", "", hit[1])), "repo_paths.R"),
               cands)
  rp <- Find(file.exists, cands)
  if (!length(rp)) stop("repo_paths.R not found", call. = FALSE)
  source(rp[1], local = envir)
}
.load_repo()
