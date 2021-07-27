library(dynatopmodel)
library(raster)
library(ggplot2)
library(rgdal)
library(xts)
# If library not available, install using 
# install.packages('libraryname')

# Get current working directory
getwd()

# Set working directory
setwd('./Desktop/Learning/Dynatopmodel/')


# Read in DEM, landcover, channel, rainfall, flow
dem <- raster('./Data/DTM data/Brompton_dem.tif')
land_cover <- raster('./Data/Land cover Brompton/LC_Brompton.tif')
channels <- readOGR('./Data/MasterMap Water Network_3621266/Brompton_channels.shp')
rain <- as.xts(read.csv.zoo('./Data/prcp.csv'))
qobs <- as.xts(read.csv.zoo('./Data/qobs.csv'))

# Visualize
plot(dem)
plot(land_cover)
plot(channels)
plot(rain)
plot(qobs, col='blue')

# Check coordinates. Re-project if all layers do not have the 
# same projection
crs(dem)
crs(land_cover)
crs(channels)

# Construct raster of channel locations for discretisation
chans.rast <- build_chans(dem = dem, drn = channels,
                          chan.width = 3, buffer=10) # Average 
# channel width obtained from channels shapefiile
plot(chans.rast, col='blue')

# Generate network routing table
routing_table <- build_routing_table(dem=dem, 
                                     chans=chans.rast,
                                     breaks=5)
head(routing_table)


land_cover.5m <- disaggregate(land_cover, fact=5) # Increase
#land cover resolution from 25 to 5m

min_ext <- dem + land_cover.5m # Calculate minimum common extent
dem_crop <- crop(dem, min_ext)
landcover_crop <- crop(land_cover.5m, min_ext)

# Build layers to produce upslope contributing area and
# topographic wetness index from DEM
layers <- build_layers(dem_crop)

# Add land cover to layers stack
layers <- addLayer(layers, landcover_crop) # Now works

# Discretise
disc <- discretise(layers = layers, chans=chans.rast,
                cuts=c(atb=10), area.thresh=0.5/100)

# Calculate potential evapotranspiration
pe <- approx.pe.ts(start='2012-09-01', end='2012-10-30',
                   dt=0.25) # Time interval (dt) is 15 mins,
# expressed in hours as (15/60)

#Set model parameters
groups <- disc$groups

# Rainfall is hourly interval. Convert to 15 min
# Also select data for time period of interest
rain <- rain['2012-09-01 12:00::2012-10-30',]
rain <- as.xts(aggregate_xts(rain, dt=24))

qobs <- qobs['2012-09-01 12:00::2012-10-30',]
qobs <- as.xts(aggregate_xts(qobs, dt=24))

pe <- aggregate_xts(pe, dt=24)
  
# Set initial specific discharge value
qt0 <- as.numeric(qobs[1,])


m <- runif(2, min=0.0011, max=0.033)
td <- runif(2, min=0.01, max=100)
ln_t0 <- runif(2, min=3, max=16)
srz_max <- runif(2, min=0.01, max=0.3)
srz0 <- runif(2, min=0.5, max=1)
vchan <- runif(2, min=500, max=5000)
vof <- runif(2, min=10, max=150)

n_iter <- 2 # Number of iterations

pb <- (txtProgressBar(min = 0, max = length(n_iter), initial = 0, 
                      style=3))

# Run model with observed flows
for (i in 1:n_iter){
  
  setTxtProgressBar(pb, i)
  
  groups$m <-m[i]
  groups$td <- td[i]
  groups$ln_t0 <- ln_t0[i]
  groups$srz_max <- srz_max[i]
  groups$srz0 <- srz0[i]
  groups$vchan <- vchan[i]
  groups$vof <- vof[i]
  
  modsim <- run.dtm(groups=groups,
                    weights=disc$weights,
                    rain=rain,
                    pe=pe,
                    qobs=qobs,
                    qt0=qt0,
                    routing=routing_table,
                    dt=24)
  
  
  
  # store simQ in an extra object
  simQ <- modsim$qsim[c(1:length(modsim$qsim))] # extract only the rows that have qobs counterpart 
  if(i==1) {
    simDat <- simQ
  } else {
    simDat <- cbind(simDat, simQ)
  }
  
  # calculate and store efficiency
nse <- NSE(simQ, modsim$qobs, digits=2)
 
  # create a vector that stores all parameters
para_vec <- (c(i,nse, groups$m[1], groups$td[1], groups$ln_t0[1],
            groups$srz_max[1], groups$srz0[1], groups$vchan[1], 
            groups$vof[1], groups$sd_max[1]))
  
  # store parameters of all model runs in a dataframe
  if (i == 1) {
    para_sets <- para_vec    
  } else {
    para_sets <- rbind(para_sets, para_vec)
  }
}

# some cosmetics on the data 
para_sets <- data.frame(para_sets)
names(para_sets) <- c("run","NSE","m", "td", "ln_t0", "srz_max", "srz0", "vchan", "vof", "sd_max")
seq <- seq(1,i,1)
names(simDat) <- seq

View(para_sets)
View(simDat)


