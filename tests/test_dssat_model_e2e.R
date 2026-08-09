#!/usr/bin/env Rscript
# One-point real DSSAT smoke test for machines with a verified installation.

root <- normalizePath(getwd(), mustWork = TRUE)
fixture <- file.path(root, "tests", "fixtures", "dssat_model_smoke")

find_dssat_executable <- function() {
  configured <- trimws(Sys.getenv("DSSAT_EXE", ""))
  exe_name <- if (.Platform$OS.type == "windows") "DSCSM048.EXE" else "dscsm048"
  candidates <- character()
  if (nzchar(configured)) candidates <- c(candidates, path.expand(configured))
  candidates <- c(candidates, file.path(dirname(root), "DSSAT48", exe_name))
  on_path <- Sys.which(exe_name)
  if (nzchar(on_path)) candidates <- c(candidates, on_path)
  candidates <- unique(candidates)
  found <- candidates[file.exists(candidates)]
  if (length(found)) normalizePath(found[[1]], mustWork = TRUE) else ""
}

exe <- find_dssat_executable()
if (!nzchar(exe)) {
  message("SKIP: verified DSSAT executable not installed; set DSSAT_EXE")
  quit(status = 0)
}

profiles <- file.path(dirname(exe), c("DSSATPRO.V48", "DSSATPRO.L48"))
if (!any(file.exists(profiles))) {
  message("SKIP: DSSAT support profile missing next to ", exe)
  quit(status = 0)
}

work <- tempfile("dssat_model_smoke_")
dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
copied <- file.copy(list.files(fixture, full.names = TRUE), work)
stopifnot(all(copied))

old_wd <- setwd(work)
on.exit(setwd(old_wd), add = TRUE)
output <- system2(exe, c("B", "DSSBatch.V48"), stdout = TRUE, stderr = TRUE)
status <- attr(output, "status")
if (is.null(status)) status <- 0L
if (status != 0L) {
  stop("DSSAT exited ", status, "\n", paste(output, collapse = "\n"))
}

if (!file.exists("summary.csv")) stop("DSSAT completed without summary.csv")
header <- strsplit(readLines("summary.csv", n = 1L), ",", fixed = TRUE)[[1]]
# DSSAT 4.8 writes one trailing field beyond the header in summary.csv.
summary <- utils::read.csv(
  "summary.csv",
  header = FALSE,
  skip = 1L,
  col.names = c(header, "EXTRA"),
  check.names = FALSE
)
stopifnot(nrow(summary) == 1L)
stopifnot(trimws(summary$CR[[1]]) == "MZ")
stopifnot(as.integer(summary$HWAM[[1]]) == 6292L)
stopifnot(as.integer(summary$CWAM[[1]]) == 15971L)
message("PASS: one-point real DSSAT model smoke test")
