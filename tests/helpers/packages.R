# Source opt-in mirrors Python. Installed packages remain the default.
if (identical(Sys.getenv("DSSAT_USE_WORKSPACE_SIBLINGS"), "1")) {
  root <- Sys.getenv("CROPMODEL_ROOT", "")
  if (!nzchar(root)) root <- if (dir.exists(".ci-deps/cropmodel")) ".ci-deps/cropmodel" else "../cropmodel"
  path <- file.path(root, "packages", "cropmodel_data", "R", "cropmodeldata")
  if (!file.exists(file.path(path, "DESCRIPTION"))) stop("Requested local package is missing: ", path)
  if (!requireNamespace("pkgload", quietly=TRUE)) stop("pkgload is required for local source opt-in")
  pkgload::load_all(path, quiet=TRUE)
  message("Testing explicit local cropmodeldata source: ", normalizePath(path))
}
