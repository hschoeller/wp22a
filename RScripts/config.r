# Paths
OUT_DIR <- "/net/scratch/schoelleh96/WP2/WP2.2a/figs/"
LM_DIR <- "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/mods/"
GEO_DIR <- "/net/scratch/schoelleh96/WP2/WP2.2a/geo_data/"
ENS_DATA_DIR <- "/net/scratch/schoelleh96/WP2/WP2.2a/ens_data/"

# Geo data
CRS <- "+proj=laea +lat_0=60 +lon_0=-20 +x_0=0 +y_0=0
    +datum=WGS84 +units=m +no_defs"
LON_BOUND <- c(-80, 40)
LAT_BOUND <- c(30, 90)
GRAT_LAT <- seq(30, 90, 30)
GRAT_LON <- seq(-90, 45, 45)

# Data
TIME_STEPS <- 31013
N_COEFFS <- 30
YEAR_BOUND <- c(1940, 2024)
CP <- c("1948-09", "1957-08", "1979-03", "1998-08", "2009-01")

# Statistical tests
ALPHA <- 0.05
N_PERM <- 10000

# Plotting
CONT_SEQ_SCALE <- "batlow"
CAT_SCALE <- "batlowS"
COLORS_REGIMES <- c(
    "#4B0082", "red", "darkorange", "gold", "yellowgreen",
    "darkgreen", "blue", "grey"
)
