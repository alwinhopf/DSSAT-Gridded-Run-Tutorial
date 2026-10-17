# A receipt is produced only after the driver publishes results and provenance.
# The orchestrator must still reject partial or subsequently damaged results.
completed_receipt_path <- function(receipt) {
  tryCatch({
    value <- jsonlite::fromJSON(receipt)
    path <- value$result_path
    if (!isTRUE(value$schema == 1L) || length(path) != 1L ||
        !file.exists(paste0(path, ".manifest.json")) || !scenario_is_complete(path))
      return(NULL)
    path
  }, error = function(e) NULL)
}

run_driver <- function(rscript, script, config, receipt, log_file) {
  # system2(env=...) is not portable for Rscript on Windows. Each sweep worker
  # is a separate R process, so scoped inherited variables are safe here.
  keys <- c("DSSAT_CONFIG_FILE", "DSSAT_RESULT_RECEIPT")
  previous <- Sys.getenv(keys, unset=NA_character_)
  on.exit({
    Sys.unsetenv(keys[is.na(previous)])
    keep <- !is.na(previous)
    if (any(keep)) do.call(Sys.setenv, as.list(previous[keep]))
  }, add=TRUE)
  Sys.setenv(DSSAT_CONFIG_FILE=config, DSSAT_RESULT_RECEIPT=receipt)
  system2(rscript, args=shQuote(script), stdout=log_file, stderr=log_file)
}
