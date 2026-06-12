source("renv/activate.R")

# --- Prefer CRAN binaries on Windows/macOS (runs after renv activation) -------
# System-library packages — openssl, sf, terra, curl, ... — ship as prebuilt
# binaries on Windows/macOS but FAIL to compile from source (no system dev
# headers / build toolchain). renv::restore() and renv::install() both honour
# getOption("pkgType"), so forcing "binary" here makes the FRESH-INSTALL path
# and the pipeline's ensure_packages() pull binaries. This is what prevents the
# openssl -> httr -> daymetr -> dssatutils source-build cascade on a new laptop.
local({
  if (.Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin") {
    options(pkgType = "binary", install.packages.check.source = "no")
  }
  # Guarantee a real, non-interactive CRAN mirror. The default "@CRAN@" sentinel
  # pops a mirror-chooser that stalls/aborts under Rscript and breaks restore().
  cran <- getOption("repos")["CRAN"]   # single bracket -> NA (no error) if absent
  if (is.na(cran) || !nzchar(cran) || identical(unname(cran), "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
})
