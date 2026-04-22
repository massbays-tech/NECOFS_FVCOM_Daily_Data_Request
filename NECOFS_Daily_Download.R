# ---- Install and Load Packages ----
pkgs <- c("ncdf4", "dplyr")

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Build Directories ----
# Always compute month folder like Data/NECOFS/2026-02
month_folder <- format(Sys.Date(), "%Y-%m")

# MASSBAYS FVCOM
out_dir_MB <- file.path("Data", "NECOFS", "MASSBAY_FORECAST", month_folder)
# Create directory if it doesn't exist
dir.create(out_dir_MB, showWarnings = FALSE, recursive = TRUE)

# NORTHEAST
out_dir_NE <- file.path("Data", "NECOFS", "NORTHEAST_FORECAST", month_folder)
# Create directory if it doesn't exist
dir.create(out_dir_NE, showWarnings = FALSE, recursive = TRUE)

# Create day stamp for saving files later
day_stamp <- format(Sys.Date(), "%Y-%m-%d")


# ---- Read site data and extract list of coordinates by site ----
site_coords_MA <- read.csv("Data/site_coords_MASSBAYS.csv")
site_coords_NE <- read.csv("Data/site_coords.csv")


# ---- Build custom function to determine which node is closest to a given coordinate point ----
nearest_node <- function(lon, lat, lon0, lat0) {
  dlon <- lon - lon0
  dlat <- lat - lat0
  
  # Scale lon distance by cos(latitude) so degrees behave roughly like meters
  # This is necessary because we need to calculate the node that is the shortest euclidean distance away. Isolating the closest lat
  # coordinate and closest lon coordinate separately won't work, because the mesh is not a regular grid.
  x <- dlon * cos(lat0 * pi/180)
  y <- dlat
  
  which.min(x^2 + y^2)
}

# ---- Build a small helper function to safely open remote netCDF URLs ----
# Remote OPeNDAP connections can occasionally fail due to server/network issues.
# This function retries nc_open() a few times before stopping.
safe_nc_open <- function(url, tries = 3, wait_sec = 10) {
  for (i in seq_len(tries)) {
    out <- try(nc_open(url), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    
    message("nc_open failed on attempt ", i, " of ", tries)
    if (i < tries) Sys.sleep(wait_sec)
  }
  stop("Could not open remote dataset after retries: ", url)
}

# ---- Build a helper function to get the current length of the time dimension ----
# Opening the time variable directly (for example with ?time) has proven slow and unstable
# for the very large NORTHEAST forecase dataset. 
# Instead, use the lightweight DDS metadata endpoint.
# The DDS is plain-text metadata that describes the structure of the remote dataset,
# including variable dimensions. In the DDS, the time variable appears like:
#   Float64 time[time = 4969];
#
# This function reads that metadata and extracts the current length of the time dimension
# without opening the time variable through ncdf4.
get_time_length <- function(base_URL) {
  dds_url <- paste0(base_URL, ".dds")
  dds_txt <- readLines(dds_url, warn = FALSE)
  
  # Match the actual line defining the numeric time variable
  time_line <- grep("Float64 time\\[time =", dds_txt, value = TRUE)
  
  if (length(time_line) == 0) {
    stop("Could not find time dimension in DDS metadata.")
  }
  
  nt <- as.integer(sub(".*time *= *([0-9]+).*", "\\1", time_line))
  
  if (is.na(nt)) {
    stop("Failed to parse time dimension length from DDS metadata.")
  }
  
  return(nt)
}


# ---- NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST ----
# Link to NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc THREDDS page:
# http://www.smast.umassd.edu:8080/thredds/catalog/models/fvcom/NECOFS/Forecasts/catalog.html?dataset=models/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc

# Set OpenDAP URL
URL <- "http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc"

# Open netCDF file to read data
nc <- safe_nc_open(URL)
  
# Get time range
# On any given day, the date range available for this data set is forecasted *two* days
# forward, plus one additional hourly time point at midnight the next day, and *three* days back.
# For example data requested the evening of March 5th will include data from midnight of
# March 2nd through March 7th, plus the midnight value on March 8th.
# This makes a total of 145 hourly values.
Times <- ncvar_get(nc, "Times")

# Extract all longitude and latitude values
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# # Extract the matrix of sigma layers (vertical layer through the water column)
# # There are 10 layers. Values in siglay matrix represent the normalized vertical position of each layer
# # (sigma coordinate), from ~0 at the surface to ~-1 at the bottom.
# siglay_matrix <- ncvar_get(nc, "siglay")

# Create empty columns for node index and not lat/lon coords
site_coords_MA$nearest_node <- NA_integer_
site_coords_MA$node_lat <- NA
site_coords_MA$node_lon <- NA

# Create an empty list of dfs the length of site_coords_MA (one df for each site)
hourly_temps <- vector("list", nrow(site_coords_MA))

# Iterate through site coordinates df identifying the closest node and corresponding coordinates
for (i in seq_len(nrow(site_coords_MA))) {
  # Determine nearest node index
  lon0 <- site_coords_MA$longitude[i]
  lat0 <- site_coords_MA$latitude[i]
  idx <- nearest_node(lon, lat, lon0, lat0)
  site_coords_MA$nearest_node[i] <- idx
  
  # Use index to determine lat and lon of nearest node
  site_coords_MA$node_lat[i] <- lat[idx]
  site_coords_MA$node_lon[i] <- lon[idx]
  
  # Get all temp measurements for one node at lowest siglay (siglay index = 10)
  ts <- ncvar_get(nc, "temp",
                  start = c(idx, 10, 1),  # (node index, sigma-layer index, first time index)
                  count = c(1, 1, -1))    # (one node, one sigma layer, all times)
  ts <- as.numeric(ts)
  
  # Build an hourly time series table for this site (one row per model timestamp) with temperature at the nearest node
  hourly_temps[[i]] <- data.frame(
    site_id = site_coords_MA$site.id[i],
    time = Times,
    temp_c = as.numeric(ts)
  )
}

# Create df of hourly temps
hourly_temps <- bind_rows(hourly_temps)

# Close .nc file
nc_close(nc)


# ---- Write CSV of daily hourly_temps MASSBAY ----
out_csv <- file.path(out_dir_MB, paste0("tempc_bot_hrly_MASSBAY_FVCOM_", day_stamp, ".csv"))
write.csv(hourly_temps, out_csv, row.names = FALSE)


# ---- NECOFS_FVCOM_OCEAN_NORTHEAST_FORECAST ----
# Link to NECOFS_FVCOM_OCEAN_NORTHEAST_FORECAST.nc THREDDS page:
# http://www.smast.umassd.edu:8080/thredds/catalog/models/fvcom/NECOFS/Forecasts/catalog.html?dataset=models/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_NORTHEAST_FORECAST.nc

# Set the base OpenDAP URL
base_URL <- "http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_NORTHEAST_FORECAST.nc"

# Read the number of times currently in the forecast data set
nt <- get_time_length(base_URL)

# On any given day, the date range available for this data set is a foretasted *five* days
# forward, plus one additional hourly time point at midnight the next day.
# The forecast also includes a seemingly variable number of retroactive days.
# To limit the amount of data requested, this script requests a total of 8 days prior to the
# date of the latest forecast, which results in *two* days back.
# For example data requested the evening of March 5th will include data from midnight of
# March 2nd through March 9th, plus the midnight value on March 10th.
# This makes a total of 193 hourly values.

# Set the starting day (time index) to be 193 hours prior to the furthest out forecast
last_idx_0 <- nt - 1        # 0-based last index
first_idx_0 <- max(0, last_idx_0 - 192) # 8 days prior to the furthest forecast

time_slice <- paste0("[", first_idx_0, ":1:", last_idx_0, "]")

# Request all lat & lon first, separately from the temperature request.
# These mesh geometry variables are lightweight compared with the forecast variables,
# and are needed to identify the nearest model node for each site.
geo_url <- paste0(base_URL, "?lon,lat")

# Open netCDF file to read geometry
nc_geo <- safe_nc_open(geo_url)

# Extract all longitude and latitude values and check ranges
lon <- ncvar_get(nc_geo, "lon")
lat <- ncvar_get(nc_geo, "lat")

# Close geometry file once lon/lat have been read
nc_close(nc_geo)

# Create empty columns for node index and node lat/lon coords
site_coords_NE$nearest_node <- NA_integer_
site_coords_NE$node_lat <- NA
site_coords_NE$node_lon <- NA

# Iterate through site coordinates df identifying the closest node and corresponding coordinates
for (i in seq_len(nrow(site_coords_NE))) {
  # Determine nearest node index
  lon0 <- site_coords_NE$longitude[i]
  lat0 <- site_coords_NE$latitude[i]
  idx <- nearest_node(lon, lat, lon0, lat0)
  site_coords_NE$nearest_node[i] <- idx

  # Use index to determine lat and lon of nearest node
  site_coords_NE$node_lat[i] <- lat[idx]
  site_coords_NE$node_lon[i] <- lon[idx]
}

# Create an empty list of dfs the length of site_coords_NE (one df for each site)
hourly_temps <- vector("list", nrow(site_coords_NE))

# Set up for loop:
# Establish total time steps and chunk size
nt <- last_idx_0 - first_idx_0 + 1
chunk_size <- 24 # 1 day

# Iterate through site coordinates df using the closest node and corresponding coordinates
for (i in seq_len(nrow(site_coords_NE))) {
  # Get nearest node index already calculated above
  idx <- site_coords_NE$nearest_node[i]

  # Convert R's 1-based index to OpenDAP's 0-based index for URL slicing
  idx_0 <- idx - 1

  # Build the data request URL for one site only.
  # Request limited range of Times and temp for:
  #   - the selected time window
  #   - the bottom sigma layer only ([44:1:44] in 0-based OpenDAP indexing)
  #   - the single nearest node for this site
  site_url <- paste0(
    base_URL,
    "?Times", time_slice,
    ",temp", time_slice, "[44:1:44][", idx_0, ":1:", idx_0, "]"
  )

  # Open netCDF file for this one site's request
  nc <- safe_nc_open(site_url)

  # Check time range
  Times <- ncvar_get(nc, "Times")

  # If Times is returned as a character matrix, collapse each entry into one timestamp string
  if (is.matrix(Times)) {
    Times <- apply(Times, 2, paste0, collapse = "")
  }

  temp_ts <- numeric(nt)

  # Read in chunks over time to report progress
  for (t0 in seq(1, nt, by = chunk_size)) {
    t1 <- min(t0 + chunk_size - 1, nt) # Set high value in chunk
    n_this <- t1 - t0 + 1 # Determine how many time values for the request

    # Get temp measurements
    vals <- ncvar_get(
      nc, "temp",
      start = c(1, 1, t0),    # (node index, sigma-layer index, time index)
      count = c(1, 1, n_this) # (one node, one sigma layer, number of time values)
    )

    temp_ts[t0:t1] <- as.numeric(vals)

    message(
      "Site ", site_coords_NE$site.id[i],
      " progress: ", t0, "-", t1, " / ", nt,
      " (", Times[t0], " to ", Times[t1], ")"
    )
  }

  # Close .nc file for this site
  nc_close(nc)

  # Build an hourly time series table for this site (one row per model timestamp) with temperature at the nearest node
  hourly_temps[[i]] <- data.frame(
    site_id = site_coords_NE$site.id[i],
    time = Times,
    temp = temp_ts
  )
}

# Create df of hourly temps
hourly_temps <- bind_rows(hourly_temps)

# ---- Write CSV of daily hourly_temps NORTHEAST ----
out_csv <- file.path(out_dir_NE, paste0("tempc_bot_hrly_NORTHEAST_NECOFS_", day_stamp, ".csv"))
write.csv(hourly_temps, out_csv, row.names = FALSE)
