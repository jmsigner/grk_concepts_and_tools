# Intro -----------------
# Author: Johannes Signer
# Date: 2025-06-27
# Verion: v1.0

# Aim: Calculate bioclimatic covariates for Germany
# See here for details: https://www.worldclim.org/data/bioclim.html

# Setup --------

# Packages we need for later
library(tidyverse) # to work with data
library(sf) # Work with vector data
library(mapview) # Visualize geographic data
library(geodata) # For data

library(gstat) # For variogram and interpolation
library(terra) # For raster data
library(mgcv) # For GAMS

# Computation outline ----------------------------------------------------------

# Steps:
# 1. Get weather data (DWD - Deutscher Wetterdienst; temperature and percipitation)
# 2. Interpolate monthly weather rasters
# 3. Calculate bioclimatic covariates

# 1. Get weather data ----------------------------------------------------------

# We have to go the CDC (Climate Data Center of the DWD)
# https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/monthly/kl/recent/

# There is a list of all stations
url <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/monthly/kl/recent/KL_Monatswerte_Beschreibung_Stationen.txt"


stations <- readLines(url)

# We are only interested in Lower Saxony & NRW (= Niedersachsen)

str_detect(stations, "Niedersachsen|Nordrhein-Westfalen") |> table()
stations_oi <- stations[str_detect(stations, "Niedersachsen|Nordrhein-Westfalen")]
head(stations)

# Station ID
ids <- str_extract(stations_oi, "^\\d{5}")

# Station from to
date <- str_extract_all(stations_oi, "\\d{8}") |>
    map(~ tibble(from = ymd(.x[1]), to = ymd(.x[2]))) |>
    bind_rows()

# Position in space
xy <- str_extract_all(stations_oi, "\\d{1,2}\\.\\d{1,4}") |>
    map(~ tibble(y = .x[1], x = .x[2])) |>
    bind_rows()

# Put it all together
sta <- bind_cols(id = ids, date, xy)

# Make stations spatial
sta <- st_as_sf(sta, coords = c("x", "y"), crs = 4326)

# Only take stations that cover 2024
sta1 <- sta |> filter(from < ymd("20231231"), to > ymd("20250101"))

# Have a look at them
ggplot(sta1) + geom_sf()

# Or interactively
mapview(sta1)

# Now let's download the actual data for all stations in sta1
sta1

# We have to go here
# https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/monthly/kl/recent/
#
# Then the files are named
# monatswerte_KL_<station id>_akt.zip
# Within this file we want:

# Strategy: 1) download file to tempfile, 2) read it and then 3) delete it again

# For the first weather station (id = 00044) this would be:
head(sta1)

base_url <- paste0("https://opendata.dwd.de/climate_environment/",
                   "CDC/observations_germany/climate/monthly/kl/recent")

url <- file.path(base_url, paste0("monatswerte_KL_", sta1$id[1] ,"_akt.zip"))

# Download file
temp <- tempfile()
download.file(url, temp)
names <- unzip(temp, list = TRUE)$Name
data <- read.csv2(unz(temp, names[grepl("^produkt_klima_monat_.*", names)]), na.strings = -999)
unlink(temp)

# We aree interested in the Start and end Date and MO_TT (temperature in 2m) and MX_RS (rainfall for the month)
data |> transmute(start = ymd(MESS_DATUM_BEGINN),
                  end = ymd(MESS_DATUM_ENDE),
                  mean_temp = MO_TT,
                  mean_rain = MX_RS)



# Lets package this into a function
get_monthly_data <- function(id) {
    base_url <- paste0("https://opendata.dwd.de/climate_environment/",
                       "CDC/observations_germany/climate/monthly/kl/recent")
    url <- file.path(base_url, paste0("monatswerte_KL_", id ,"_akt.zip"))

    # Download file
    temp <- tempfile()
    download.file(url, temp)
    names <- unzip(temp, list = TRUE)$Name
    data <- read.csv2(unz(temp, names[grepl("^produkt_klima_monat_.*", names)]),
                      na.strings = -999)
    unlink(temp)
    data |> transmute(
        start = ymd(MESS_DATUM_BEGINN),  end = ymd(MESS_DATUM_ENDE),
        mean_temp = as.numeric(MO_TT),  mean_rain = as.numeric(MX_RS))
}

# Now lets apply this for all stations in sta1
sta2 <- sta1 |> mutate(data = map(id, get_monthly_data))

# Lets extract the data for 2024
sta2$data[[3]] |> filter(year(start) == 2024) |>
    mutate(month = month(start, label = TRUE))

# Apply to all
sta2 <- sta2 |> mutate(data = map(data, ~ .x |> filter(year(start) == 2024) |>
    mutate(month = month(start, label = TRUE))))

# Lets focus on Temperature
sta2 <- sta2 |> select(id, data) |>
    unnest(cols = data)

sta2 |> filter(month == "Apr") |>
    ggplot(aes(col = mean_temp)) + geom_sf()

# Outline lower Saxony
ger <- gadm("Germany", level = 1, path = "data/") |> st_as_sf()
states <- filter(ger, NAME_1 %in% c("Niedersachsen", "Nordrhein-Westfalen"))

sta2 |> filter(month == "Jan") |>
    ggplot() +
    geom_sf(data = states) +
    geom_sf(aes(col = mean_rain))


# 2. Interpolation -------------------------------------------------------------

# First lets check if there is some spatial autocorrelation
# Are nearby weather stations more similar to each other than distant ones?

xy <- sta2 |> filter(month == "Jan") |>
    filter(!is.na(mean_temp))

out_rast <- rast(states, resolution = 0.05)
values(out_rast) <- 1
out_grid <- as.data.frame(out_rast, xy = TRUE) |>
    st_as_sf(coords = c("x", "y"), crs = 4326) |>
    st_filter(states)
ggplot(out_grid) + geom_sf(pch = ".")

# .. 0. Nearest neighbours ----
nn1 <- voronoi(vect(xy |> select(mean_temp)), bnd = vect(states))
nn1 |> st_as_sf() |>
    ggplot(aes(fill = mean_temp)) + geom_sf()

# Cut to study area
nn1 |> st_as_sf() |> st_intersection(states) |>
    ggplot(aes(fill = mean_temp)) + geom_sf()

# .. 1. IDW -------------

idw.1 <- gstat(formula = mean_temp ~ 1, locations = xy,
             nmax = nrow(xy), # maximum number of neighbours
             set = list(idp = 1)) # beta coefficient -> weighting factor

# idp: determines to which degree nearer observations are prefered.

idw.1 <- predict(idw.1, out_grid)

idw.1 <- rasterize(idw.1, out_rast, field = "var1.pred", fun = "mean")
plot(idw.1)

# If we set idp = 0 and nmax = nrow(xy); we get the mean
idw.2 <- gstat(formula = mean_temp ~ 1, locations = xy,
             nmax = nrow(xy), set = list(idp = 0))

idw.2 <- predict(idw.2, out_grid)
idw.2 <- rasterize(idw.2, out_rast, field = "var1.pred", fun = "mean")
plot(idw.2)
mean(xy$mean_temp)

# On the other hand, if we set idp = and nmax = 1, we get the same as in 0).
idw.2 <- gstat(formula = mean_temp ~ 1, locations = xy,
             nmax = 1, set = list(idp = 0))

idw.2 <- predict(idw.2, out_grid)
idw.2 <- rasterize(idw.2, out_rast, field = "var1.pred", fun = "mean")
plot(idw.2)

# .. 2. Ordinary Kriging ----

# First we have to fit a variogram to the data.
# A variogram is a method to measure how much spatial autocorrelation is in the
# data.
# See here:
# Parameters are:
# - nugget: Variability at the origin.
# - sill: Limit of that variance at infenitely large distances.
# - range: The distance after which autocorrelation fades out.
# Source: https://en.wikipedia.org/wiki/Variogram#Parameters and
# https://en.wikipedia.org/wiki/Variogram#/media/File:Schematic_variogram.svg

# Fit an empirical variogram
vg_emp <- variogram(mean_temp ~ 1, data = xy, cutoff = 250)
plot(vg_emp)

vg_ini <- vgm(psill = 0.8, model = "Sph",
              range = 120, nugget = 0.2)
plot(vg_emp, vg_ini, cutoff = 1000, cex = 1.5)

vg_fit <- fit.variogram(vg_emp, model = vg_ini)

plot(vg_emp, vg_fit)

vg_fit

# We can now use the variogram to interpolate the data with ordinary kriging
# (without covariates) or co-kriging (with covariates).

krige <- gstat(formula = mean_temp ~ 1, data = xy, model = vg_fit)
krige <- predict(krige, out_grid)
krige
krige.pred <- rasterize(krige, out_rast, field = "var1.pred", fun = "mean")
krige.var <- rasterize(krige, out_rast, field = "var1.var", fun = "mean")
plot(krige.pred)
plot(krige.var)

# .. 3. Generalised Additive Model (GAM) ----------------------------

# We could also youse a gam here

xy$x <- st_coordinates(xy)[, 1]
xy$y <- st_coordinates(xy)[, 2]

out_xy <- out_grid |> select(-lyr.1)
out_xy$x <- st_coordinates(out_xy)[, 1]
out_xy$y <- st_coordinates(out_xy)[, 2]
out_xy <- out_xy |> st_set_geometry(NULL)

g1 <- gam(mean_temp ~ s(x) + s(y), data = xy)
out_xy$gam <- predict(g1, newdata = out_xy)

ggplot(out_xy, aes(x, y, fill = gam)) + geom_raster()

# 3. Calculate bioclim variables ----------

# Lets pick a interpolation method. For sake of simplicity, I will use use the IDW
# here and interpolate the mean temperature for every month for our area of
# interest.

# This is the same for every month of the year
out_rast <- rast(states, resolution = 0.05)
values(out_rast) <- 1
out_grid <- as.data.frame(out_rast, xy = TRUE) |>
    st_as_sf(coords = c("x", "y"), crs = 4326) |>
    st_filter(states)

xy <- sta2 |> mutate(month1 = month(start)) # to have a numeric month
xy

# Iterate over the 12 months
results <- map(1:12, ~ {
    xy_temp <- filter(xy, month1 == .x) |>
        filter(!is.na(mean_temp))
    gstat(formula = mean_temp ~ 1, locations = xy_temp, nmax = nrow(xy),
          set = list(idp = 1)) |>
        predict(out_grid) |>
        rasterize(out_rast, field = "var1.pred", fun = "mean")
})

r_temp <- rast(results)
names(r_temp) <- 1:12
r_temp

plot(r_temp)

# Something likely went wrong with missing values!
xy <- sta2 |> mutate(month1 = month(start)) # to have a numeric month
range(xy$mean_temp, na.rm = TRUE)
xy$mean_temp <- ifelse(xy$mean_temp == -999, NA, xy$mean_temp)

# Iterate over the 12 months

results <- map(1:12, ~ {
    xy_temp <- filter(xy, month1 == .x) |>
        filter(!is.na(mean_temp))
    gstat(formula = mean_temp ~ 1, locations = xy_temp, nmax = nrow(xy),
          set = list(idp = 1)) |>
        predict(out_grid) |>
        rasterize(out_rast, field = "var1.pred", fun = "mean")
})

r_temp <- rast(results)
names(r_temp) <- 1:12
r_temp

plot(r_temp) # Looks better


# Now lets calculate some of bioclim variables:
# - Bio 1: Annual mean temp
# - Bio 4: Temp seasonality (sd of temp)
# - Bio 10: Mean temp of the warmest quarter

# Bioclim 1
mean(r_temp)

# or
app(r_temp, sd)

# Bioclim 4
app(r_temp, sd)

# Biocom 10
quarters <- list(1:3, 4:6, 7:9, 10:12)
map_dbl(quarters, ~ mean(r_temp[[.x]][], na.rm = TRUE)) # Q3 is the warmest

app(r_temp[[quarters[3]]], mean)

# 4. Points to improve --------------

# 1. Interpolation: add additional algorithms.
# 2. Interpolation: cross validate different algorithms, to find the best one.
# 3. Interpolation: perform an ensamble approach (e.g., https://www.paulamoraga.com/book-spatial/spatial-interpolation-methods.html#ensemble-approach).

# 4. Calculate all bioclimatic variables.
# 5. Possibly use a finer temporal and/or spatial resolution.
# 6. Calculate bioclimatic variables for all of Germany and different years.

# 5. Sources -------

# Some wources that might be useful

# - https://www.paulamoraga.com/book-spatial/
# - https://www.worldclim.org/data/bioclim.html
# - https://opendata.dwd.de/
# - https://r-spatial.org/book/
# - https://r.geocompx.org/

