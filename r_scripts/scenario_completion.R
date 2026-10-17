complete_point_ids <- function(data, treatments) {
  if (is.null(data) || !nrow(data) || !all(c("point_id", "treatment") %in% names(data))) return(character())
  ids <- unique(as.character(data$point_id))
  ids[vapply(ids, function(id) {
    values <- suppressWarnings(as.numeric(data$treatment[as.character(data$point_id) == id]))
    values <- values[is.finite(values) & values == round(values)]
    all(as.integer(treatments) %in% values)
  }, logical(1))]
}

start_completion <- function(result_path) {
  path <- paste0(result_path,".status.json")
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  value <- tryCatch(jsonlite::fromJSON(path,simplifyVector=FALSE),error=function(e) list())
  value$status <- "running";value$pid <- Sys.getpid()
  tmp <- tempfile(".completion-",tmpdir=dirname(path))
  jsonlite::write_json(value,tmp,auto_unbox=TRUE,pretty=TRUE,null="null")
  if(!file.rename(tmp,path)) stop("Could not publish running scenario status")
}

resume_identity <- function(identity) {
  keys <- names(identity$resolved_inputs)
  identity$resolved_inputs <- identity$resolved_inputs[!grepl("^(weather_wth|soil_sol|soil_mapping)",keys)]
  identity
}

point_input_signatures <- function(grid,id_col,soil_dir,weather_dir) {
  wgs <- sf::st_transform(grid,4326);coords <- sf::st_coordinates(wgs)
  out <- list()
  digest <- function(path) if(file.exists(path)) unname(tools::md5sum(path)) else NULL
  for(i in seq_len(nrow(wgs))) {
    id <- as.character(wgs[[id_col]][i])
    out[[id]] <- list(soil=digest(file.path(soil_dir,paste0(id,".SOL"))),
      weather=digest(file.path(weather_dir,paste0(id,".WTH"))),longitude=coords[i,1],latitude=coords[i,2])
  }
  out
}

restore_aggregate_points <- function(result_path, run_dir, identity, treatments, expected_ids, input_signatures=NULL) {
  manifest <- paste0(result_path, ".manifest.json")
  if (!file.exists(result_path) || !file.exists(manifest)) return(character())
  previous <- tryCatch(jsonlite::fromJSON(manifest, simplifyVector = FALSE), error = function(e) NULL)
  current <- jsonlite::fromJSON(jsonlite::toJSON(identity, auto_unbox=TRUE, null="null"), simplifyVector=FALSE)
  reusable <- NULL
  if (!identical(previous, current)) {
    status <- tryCatch(jsonlite::fromJSON(paste0(result_path,".status.json"),simplifyVector=FALSE),error=function(e) NULL)
    normalized <- jsonlite::fromJSON(jsonlite::toJSON(resume_identity(identity),auto_unbox=TRUE,null="null"),simplifyVector=FALSE)
    if(!identical(status$resume_identity,normalized)) return(character())
    signatures <- jsonlite::fromJSON(jsonlite::toJSON(input_signatures,auto_unbox=TRUE,null="null"),simplifyVector=FALSE)
    reusable <- names(signatures)[vapply(names(signatures),function(id) {
      !is.null(signatures[[id]]$soil) && !is.null(signatures[[id]]$weather) && identical(signatures[[id]],status$point_input_signatures[[id]])
    },logical(1))]
  }
  data <- tryCatch(read.csv(result_path, colClasses = c(point_id="character")), error=function(e) NULL)
  complete <- intersect(complete_point_ids(data, treatments), as.character(expected_ids))
  if(!is.null(reusable)) complete <- intersect(complete,reusable)
  for (id in complete) {
    folder <- file.path(run_dir, id); dir.create(folder, recursive=TRUE, showWarnings=FALSE)
    write.csv(data[as.character(data$point_id)==id,,drop=FALSE], file.path(folder,paste0("results_",id,".csv")), row.names=FALSE)
  }
  complete
}

publish_completion <- function(result_path, data, expected_ids, treatments, excluded_ids=character(),identity=NULL,input_signatures=NULL) {
  expected <- unique(as.character(expected_ids))
  complete <- intersect(complete_point_ids(data, treatments), expected)
  excluded <- setdiff(intersect(as.character(excluded_ids), expected), complete)
  failures <- setdiff(expected, union(complete, excluded))
  value <- list(schema=1L, status=if(length(failures)) "partial" else "complete",
                expected_ids=as.list(sort(expected)), completed_ids=as.list(sort(complete)),
                excluded_ids=as.list(sort(excluded)), failed_ids=as.list(sort(failures)),
                treatments=as.list(as.integer(treatments)))
  if(!is.null(identity)) {
    value$resume_identity <- resume_identity(identity)
    value$point_input_signatures <- input_signatures
  }
  path <- paste0(result_path,".status.json")
  tmp <- tempfile(".completion-",tmpdir=dirname(path))
  jsonlite::write_json(value,tmp,auto_unbox=TRUE,pretty=TRUE,null="null")
  if (!file.rename(tmp,path)) stop("Could not publish scenario completion status")
  value
}

validate_cached_point_results <- function(run_dir, expected_ids, treatments) {
  for(id in as.character(expected_ids)) {
    path <- file.path(run_dir,id,paste0("results_",id,".csv"))
    if(!file.exists(path)) next
    data <- tryCatch(read.csv(path,colClasses=c(point_id="character")),error=function(e) NULL)
    if(!id %in% complete_point_ids(data,treatments)) {
      archive <- tempfile(paste0("incomplete_results_",id,"_"),tmpdir=dirname(path),fileext=".csv")
      if(!file.rename(path,archive)) stop("Could not preserve incomplete point result for retry")
    }
  }
}

scenario_is_complete <- function(result_path) {
  tryCatch({
    value <- jsonlite::fromJSON(paste0(result_path,".status.json"))
    data <- read.csv(result_path,colClasses=c(point_id="character"))
    actual <- complete_point_ids(data,value$treatments)
    isTRUE(value$schema==1L) && identical(value$status,"complete") && length(value$expected_ids)>0 &&
      !length(value$failed_ids) && all(value$completed_ids %in% actual) &&
      all(value$expected_ids %in% union(actual,value$excluded_ids))
  },error=function(e) FALSE)
}
