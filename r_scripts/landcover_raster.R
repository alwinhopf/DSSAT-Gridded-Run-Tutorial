suppressPackageStartupMessages({
  library(terra)
})

# landcover_raster.R
# ------------------
# PURPOSE
#   Clip a national landcover raster (NLCD or CDL) to a set of US states and
#   produce a binary cropland mask (1 = cropland, 0 = everything else).
#   The mask is the input for landcover_raster_to_gridpoints.R.
#
# RUN FROM THE PROJECT ROOT so that all relative paths below resolve correctly:
#   Rscript r_scripts/landcover_raster.R
# Or source it from dssat_main_pipeline.R (which sets the working directory first).
#
# -----------------------
# USER SETTINGS
# -----------------------

# Path to the raw landcover raster (relative to project root).
# Download NLCD: https://www.mrlc.gov/data     (Annual NLCD, GeoTIFF ~2 GB per year)
# Download CDL:  https://croplandcros.scinet.usda.gov/
input_raster <- file.path("data", "landcover", "Annual_NLCD_LndCov_2024_CU_C1V1",
                          "Annual_NLCD_LndCov_2024_CU_C1V1.tif")

# US state boundary shapefile (relative to project root).
# Download: see README section "Example boundary shapefile (US): TIGER/Line states"
boundary_vector <- file.path("shapefile", "tl_2024_us_state.shp")

# States to process — use full names or two-letter abbreviations (e.g., "MT" or "Montana").
state_names <- c(
  "Montana", "North Dakota", "South Dakota", "Wyoming",
  "Nebraska", "Kansas", "Colorado", "Oklahoma",
  "Minnesota", "Iowa", "Missouri", "Wisconsin", "Illinois", "Michigan",
  "Indiana", "Ohio", "Idaho", "Utah", "Washington", "Oregon", "Nevada"
)

# Landcover class codes to treat as "cropland" (all other classes become 0).
# NLCD class 82 = Cultivated Crops (the broadest cropland category).
# For CDL, replace with the crop-specific codes you need, e.g.:
#   c(1)   = corn
#   c(5)   = soybeans
#   c(1,5) = corn + soybeans
# CDL code list: https://www.nass.usda.gov/Research_and_Science/Cropland/sarsfaqs2.php
cropland_values <- c(82)   # NLCD: 82 = Cultivated Crops

# Where to save the outputs (relative to project root).
output_dir   <- file.path("data", "landcover", "derived")
output_mask  <- file.path(output_dir, "cropland_mask_region.tif")
# Write an individual mask file per state as well?
write_per_state <- TRUE

# -----------------------
# END USER SETTINGS
# -----------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

state_lookup <- data.frame(
  NAME = c(
    "Alabama","Alaska","Arizona","Arkansas","California","Colorado","Connecticut","Delaware",
    "District of Columbia","Florida","Georgia","Hawaii","Idaho","Illinois","Indiana","Iowa",
    "Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota",
    "Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire","New Jersey",
    "New Mexico","New York","North Carolina","North Dakota","Ohio","Oklahoma","Oregon",
    "Pennsylvania","Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah",
    "Vermont","Virginia","Washington","West Virginia","Wisconsin","Wyoming"
  ),
  STUSPS = c(
    "AL","AK","AZ","AR","CA","CO","CT","DE",
    "DC","FL","GA","HI","ID","IL","IN","IA",
    "KS","KY","LA","ME","MD","MA","MI","MN",
    "MS","MO","MT","NE","NV","NH","NJ",
    "NM","NY","NC","ND","OH","OK","OR",
    "PA","RI","SC","SD","TN","TX","UT",
    "VT","VA","WA","WV","WI","WY"
  ),
  stringsAsFactors = FALSE
)

normalize_states <- function(x) {
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  x
}

resolve_states <- function(state_names, lookup) {
  supplied <- normalize_states(state_names)
  if (length(supplied) == 0) return(NULL)
  is_abbr <- all(nchar(supplied) <= 2)
  if (is_abbr) {
    matched <- lookup$NAME[match(toupper(supplied), lookup$STUSPS)]
  } else {
    matched <- lookup$NAME[match(supplied, lookup$NAME)]
  }
  matched <- matched[!is.na(matched)]
  unique(matched)
}

cat("Loading landcover raster...\n")
lc <- rast(input_raster)

cat("Loading boundary vector...\n")
bnd <- vect(boundary_vector)

if (!all(c("NAME", "STUSPS") %in% names(bnd))) {
  stop("Boundary vector must contain NAME and STUSPS fields. Available fields: ",
       paste(names(bnd), collapse = ", "))
}

allowed_states <- resolve_states(state_names, state_lookup)
if (length(allowed_states) == 0) {
  stop("No valid states found in state_names.")
}

cat("Filtering boundary to requested states using NAME...\n")
bnd <- bnd[trimws(as.character(bnd$NAME)) %in% allowed_states, ]

if (nrow(bnd) == 0) {
  stop("No matching features found for the requested state_names. Check spelling/case.")
}

bnd_proj <- project(bnd, crs(lc))

cat("Cropping and masking landcover to region...\n")
lc_roi <- crop(lc, bnd_proj)
lc_roi <- mask(lc_roi, bnd_proj)

cat("Creating binary cropland mask...\n")
cropland_mask <- lc_roi %in% cropland_values

writeRaster(cropland_mask, output_mask, overwrite = TRUE)
cat("Saved cropland mask to:", output_mask, "\n")

if (isTRUE(write_per_state)) {
  cat("Writing per-state masks...\n")
  for (nm in allowed_states) {
    sel <- bnd[trimws(as.character(bnd$NAME)) == nm, ]
    if (nrow(sel) == 0) next
    sel_proj <- project(sel, crs(lc))
    lc_state <- crop(lc, sel_proj)
    lc_state <- mask(lc_state, sel_proj)
    m_state <- lc_state %in% cropland_values
    out_state <- file.path(output_dir, paste0("cropland_mask_", gsub("\\s+", "_", tolower(nm)), ".tif"))
    writeRaster(m_state, out_state, overwrite = TRUE)
    cat(" ", nm, "->", out_state, "\n")
  }
}

cat("Done.\n")