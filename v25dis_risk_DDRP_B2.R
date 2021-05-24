#!/usr/bin/Rscript
.libPaths("/usr/lib64/R/library/")
options(echo=FALSE)

# Install pakcages if necessary and load them
pkgs <- c("doParallel", "dplyr", "foreach", "ggplot2", "ggthemes", 
          "lubridate", "mapdata", "mgsub", "optparse", "parallel",
          "purrr", "RColorBrewer", "rgdal", "raster", "readr", "sp", 
          "stringr", "tidyr", "tictoc", "tools", "toOrdinal")
ld_pkgs <- lapply(pkgs, library, character.only = TRUE) # load them

# Bring in states feature for summary maps (PNG files), extract the lower 
# 48 states, and project it. Requires these libraries: "mapdata" and "maptools."
cat("Downloading US states feature\n")
states <- map_data("state")

# Start timing the model run
tic("Total run time")

################################################################################
########                  Header Documentation Section             #############
#  DDRP for disease risk: Degree-Day, estab. Risk, and Pest event mapping system #
########  By Brittany Barker, Len Coop, Gericke Cook, Dan Upper, and ###########
########  Tyson Wepprich for APHIS PPQ and IPM needs ###########################
################################################################################

# Last modified by B. Barker on 5/24/21

# 11v. This is for development of a plant disease risk version of DDRP so use prism data such as:
#   A) Add Drystress looking similar to Heat and Cold (Chill) stress, using relhum=func(tmean,tdmean)
#   B) Add DDLW unique to Boxwood blight since BoxDD is a piecewise regression formula based on BB DD lookup table used in site vers.
#"PRISM_tddif_early_4kmD2_MTD_20210401=(PRISM_tmean_early_4kmD2_MTD_20210401-PRISM_tdmean_early_4kmD2_MTD_20210401)"
#"LW_20210401 = if(PRISM_ppt_early_4kmD2_MTD_20210401 > 8,1,if(PRISM_tddef_early_4kmD2_MTD_20210401 < 5,2,0))"
# the thresholds values can be calibrated and included as input params from the spp_param files
#

#### (1). PARAM HANDLING -----
#### * Process command line args ####
#### # Uses "optparse" package, a command line parser similar to Python's
option_list <- list(
  make_option(c("--spp"), action = "store", type = "character", 
              default = NA, help = "species parameter file name"),
  make_option(c("--forecast_data"), action = "store", type = "character", 
              default = NA, help = "weather data for forecast"),
  make_option(c("--start_year"),  action = "store", type = "character", 
              default = NA, help = "start year"),
  make_option(c("--start_doy"), action = "store", type = "integer", 
              default = NA, help = "start day of year"),
  make_option(c("--end_doy"),  action = "store", type = "integer",
              default = NA, help = "end day of year"),
  make_option(c("--keep_leap"), action = "store", type = "integer", 
              default = NA, help = "should leap day be kept? 0=no, 1=yes"),
  make_option(c("--region_param"), type = "character", action = "store", 
              default = NA, help = "study region: CONUS, EAST, WEST, or 
              state (2-letter abbr.)"),
  make_option(c("--exclusions_stressunits"), type = "integer", action = "store",
              default = NA, help = "0 = off, 1 = on"),
  make_option(c("--pems"), type = "integer", action = "store", 
              default = NA, help = "0 = off, 1 = on"),
  make_option(c("--mapA"), type = "integer", action = "store", 
              default = NA, help = "0 = off, 1 = on"),
  make_option(c("--mapE"), type = "integer", action = "store", 
              default = NA, help = "0 = off, 1 = on"),
  make_option(c("--mapL"), type = "integer", action = "store", 
              default = NA, help = "0 = off, 1 = on"),
  make_option(c("--mapP"), type = "integer", action = "store", 
              default = NA, help = "0 = off, 1 = on"),
  make_option(c("--out_dir"), type = "character", action = "store", 
              default = NA, help = "name of out directory"),
  make_option(c("--out_option"), type = "integer", action = "store", 
              default = NA, help = "sampling frequency: 1 = 30 days; 
              2 = 14 days; 3 = 10 days; 4 = 7 days; 5 = 2 days; 6 = 1 day")
)

# Read in commands 
# If command line isn't used, then opts list elements will be NA
opts <- parse_args(OptionParser(option_list = option_list))
if (!is.na(opts[1])) {
  spp <- opts$spp
  forecast_data <- opts$forecast_data
  start_year <- opts$start_year
  start_doy <- opts$start_doy
  end_doy <- opts$end_doy
  keep_leap <- opts$keep_leap
  region_param <- opts$region_param
  exclusions_stressunits <- opts$exclusions_stressunits
  pems <- opts$pems
  mapA <- opts$mapA
  mapE <- opts$mapE
  mapL <- opts$mapL
  mapP <- opts$mapP
  out_dir <- opts$out_dir
  out_option <- opts$out_option
  odd_gen_map <- opts$odd_gen_map
} else {
  #### * Default values for params, if not provided in command line ####
  spp           <- "BOXB" # Default species to use
  forecast_data <- "PRISM" # Forecast data to use (PRISM or NMME)
  start_year    <- "2020" # Year to use
  start_doy     <- 1 # Start day of year          
  end_doy       <- 365 # End day of year - need 365 if voltinism map 
  keep_leap     <- 1 # Should leap day be kept?
  region_param  <- "OR" # Region [CONUS, EAST, WEST, or state (2-letter abbr.)]
  exclusions_stressunits    <- 0 # Turn on/off climate stress unit exclusions
  pems          <- 1 # Turn on/off pest event maps
  mapA          <- 1 # Make maps for adult stage
  mapE          <- 0 # Make maps for egg stage
  mapL          <- 0 # Make maps for larval stage
  mapP          <- 0 # Make maps for pupal stage
  out_dir       <- "BOXB_test2" # Output dir
  out_option    <- 1 # Sampling frequency
}


#### (2). FUNCTIONS AND PLOT SETTINGS -----

# For generating log text (txt1 = line breaks, txt2 = log message)
LogWrap <- function(brk1, txt, brk2) { 
  cat(brk1, str_wrap(txt, width = 80), brk2, sep = "", 
      file = Model_rlogging, append = TRUE) 
}

# * DD Calculation Methods funcs ###
# Tmean (no upper threshold, Tmax+Tmin AVG w/horiz upper threshold, 
# and Single Triangle w/horiz. upper threshold

# Simple Mean Temp DD Calc method: ((tmean > LDT) * (tmean - LDT))
# same result as  max((tmax + tmin)/2 - LDT,0)  
# so no need for tmean PRISM data. 
SimpDD <- function(tmax,tmin,LDT) {
  return(max((tmax + tmin)/2 - LDT,0))
}

# Averaging DD Calc method (max + min/2 - tlow) but with horizontal (substitution) 
# upper threshold:
AvgDD <- function(tmax, tmin, LDT, UDT) {
  return(Cond(
    tmax < LDT, 0, Cond(tmin > UDT, 0, 
                        Cond(tmax > UDT, (UDT + tmin)/2 - LDT, 
                             Cond((tmax + tmin)/2-LDT < 0,0,(tmax + tmin)/2 - LDT)))))
}

# Single triangle with upper threshold (Sevachurian et al. 1977) - 
# also a good substitution for single sine method
TriDD <- function(tmax, tmin, LDT, UDT) {
  Tmp1=6*((tmax-LDT)*(tmax-LDT))/(tmax-tmin)
  Tmp2=6*((tmax-UDT)*(tmax-UDT))/(tmax-tmin)
  Cond(tmax < LDT,0,
       Cond(tmin >= UDT,UDT-LDT,
            Cond((tmax < UDT) & (tmin <= LDT), Tmp1/12,
                 Cond((tmin <= LDT) & (tmax >= UDT), (Tmp1-Tmp2)/12,
                      Cond((tmin > LDT) & (tmax >= UDT), 6*(tmax+tmin-2*LDT)/12 - (Tmp2/12),
                           Cond((tmin > LDT) & (tmax < UDT), 6*(tmax+tmin-2*LDT)/12,0))))))
} 

# * BOXB funcs ####
#### Relative Humidity - from tmean and tdmean
# RH: =100*(EXP((17.625*TD)/(243.04+TD))/EXP((17.625*T)/(243.04+T)))
# calc RH for todays avg Temp and Dewpoint Temp.  -verified formula
# use ppt and RH to estimate leaf wetness (LW)
#Later add: r.mapcalc "LW=if(ppt > $ppt_T3,4,if(RH > $RH_T3,3,if(ppt > $ppt_T2 && RH > $RH_T2,2,if(ppt > $ppt_T1 && RH > $RH_T1,1,0))))"
#thresholds for DDs; DONE: replace with lookup table or piecewise regression approximation of infection curve
calcrelhum <- function(tmean, tdmean){
  return(100*(exp((17.625*tdmean)/(243.04+tdmean))/exp((17.625*tmean)/(243.04+tmean))))
}

ppt_T1 <- 0.3
ppt_T2 <- 8
ppt_T3 <- 30
RH_T1 <- 68
RH_T2 <- 82
RH_T3 <- 88

calcLW <- function(ppt, relhum){
  return(Cond(ppt > ppt_T3,4,
              Cond(relhum > RH_T3,3,
                   Cond(ppt > ppt_T2 & relhum > RH_T2,2,
                        Cond(ppt > ppt_T1 & relhum > RH_T1,1,0)))))
}

calcDDs <- function(tmax, tmin){
  #divide 24 hour day into six 4hr bins for piecewise regression version of lookup table
  r2 <- tmin+(1*((tmax-tmin)/5))
  r3 <- tmin+(2*((tmax-tmin)/5))
  r4 <- tmin+(3*((tmax-tmin)/5))
  r5 <- tmin+(4*((tmax-tmin)/5))
  #for each bin, calculate risk using piecewise regression version of lookup table
  p4  <- Cond(tmin < 9,0,Cond(tmin < 19.4,1.068*tmin-9.08,Cond(tmin< 24.9,4.182*tmin-69.117,Cond(tmin<29,-8.25*tmin+238.3,0))))
  p8  <- Cond(r2 < 9,0,Cond(r2 < 19.4,1.068*r2-9.08,Cond(r2< 24.9,4.182*r2-69.117,Cond(r2<29,-8.25*r2+238.3,0))))
  p12 <- Cond(r3 < 9,0,Cond(r3 < 19.4,1.068*r3-9.08,Cond(r3< 24.9,4.182*r3-69.117,Cond(r3<29,-8.25*r3+238.3,0))))
  p16 <- Cond(r4 < 9,0,Cond(r4 < 19.4,1.068*r4-9.08,Cond(r4< 24.9,4.182*r4-69.117,Cond(r4<29,-8.25*r4+238.3,0))))
  p20 <- Cond(r5 < 9,0,Cond(r5 < 19.4,1.068*r5-9.08,Cond(r5< 24.9,4.182*r5-69.117,Cond(r5<29,-8.25*r5+238.3,0))))
  p24 <- Cond(tmax < 9,0,Cond(tmax <19.4,1.068*tmax-9.08,Cond(tmax< 24.9,4.182*tmax-69.117,Cond(tmax<29,-8.25*tmax+238.3,0))))
  #average the results equivalent to daily degree-days
  return (p4+p8+p12+p16+p20+p24)/6
}

# * Parallel processing funcs ####
#### RegCluster: register cluster for parallel computing ##
# The "RegCluster" function creates a set of copies of R running in parallel and 
# communicating over sockets (parallel socket clusters). The value may be 
# specified manually; however, here the value is estimated based on the number 
# of available cores on the computer or server DDRP is being run on. 
# Specifying too many clusters will overload the computer.
RegCluster <- function(value) {
  
  # Change the value if it is 1 (this could happen if number of cores is low)
  # Otherwise the whole process will not be run in parallel
  if (value == 1) {
    value <- 2
  }
  
  if (grepl("Windows", Sys.info()[1])) {
    cl <<- makePSOCKcluster(value)
  } else {
    cl <<- makeCluster(value)
  }
  
  # If run is being done on Hopper, need to specify the library for each worker
  if (Sys.info()["nodename"] == "hopper.science.oregonstate.edu") {
    clusterEvalQ(cl, .libPaths("/usr/local/lib64/R/library/"))
  }
  
  doParallel::registerDoParallel(cl)
  return(cl)
}

# * Raster processing funcs ####

# CombineMaps: merge raster tiles
# Merge tiles back together for CONUS or EAST
CombineMaps <- function(brick_files, type) {
  # Search for file type in brick file list
  fls_by_type <- brick_files[grep(pattern = paste0(type, "_tile"), x = brick_files, 
                                  fixed = TRUE)]
  mrgd <- Reduce(raster::merge, lapply(fls_by_type, brick))
  # Get datatype from brick - most are INT2S but a few FLT4S or INT1U
  datatype <- dataType(mrgd) 
  writeRaster(mrgd, filename = type, overwrite = TRUE, 
              datatype = datatype, format = "GTiff")
}

# Extract best PRISM file (highest quality for each date)
# Take .bil files from PRISM or NMME (= forecast_data) yearly directories 
# (= files). The function returns the best data for each day. If data are 
# from a leap year, then leap day data may be removed or kept (= keep_leap).
ExtractBestPRISM <- function(files, forecast_data, keep_leap) {
  numsplits <- str_count(string = files[1], pattern = "/")
  pfile <- str_split(string = files, pattern = coll("/"), numsplits) %>% 
    purrr::map(numsplits)
  qa <- str_split(string = pfile, pattern = coll("_"), 6) %>% purrr::map(3) %>% 
    unlist()
  
  # Extract dates - the dates will always be the largest number, so sort and
  # then take first element of vector
  dates <- regexpr(pattern = "[0-9]{8}", text = files)
  df <- data.frame(dates = regmatches(files, dates),
                   quality = substr(qa, start = 1, stop = 4),
                   rownum = 1:length(qa))
  
  # Added to Tyson's version of this function - ability to choose PRISM vs. 
  # NMME for future dates
  if (forecast_data == "PRISM") {
    df <- mutate(df, rank = ifelse(quality == "stab", 1, 
                                   ifelse(quality == "prov", 2, 
                                          ifelse(quality == "earl", 3, 
                                                 ifelse(quality == "10yr", 4, 
                                                        ifelse(quality == "30yr", 5, 
                                                               ifelse(quality == "nmme", 6, NA)))))))
  } else if (dat == "NMME") {
    df <- mutate(df, rank = ifelse(quality == "stab", 1, 
                                   ifelse(quality == "prov", 2, 
                                          ifelse(quality == "earl", 3, 
                                                 ifelse(quality == "nmme", 4, 
                                                        ifelse(quality == "10yr", 5, 
                                                               ifelse(quality == "30yr", 6, NA)))))))
  }
  
  # Sorting backwards matches data hierarchy
  # If PRISM, then stable > provisional > early > 10yr > 30yr > nmme
  # If NMME, then stable > provisional > early > nmme > 10yr > 30yr
  df2 <- df %>% group_by(dates) %>% 
    dplyr::filter(rank == min(rank)) %>%
    dplyr::filter(1:n() == 1)
  
  best <- files[df2$rownum]
  dates <- regexpr(pattern = "[0-9]{8}", text = best)
  
  fileorder <- order(regmatches(best, dates))
  files <- best[fileorder]
  
  # Remove leap day (2/29) if the start_year is not a leap year, or if it is 
  # a leap year but user wants it removed (keep_leap == 0). This doesn't apply 
  # if 30 yr average data or other non-specific year data are being used
  # (in this case, start_year is non-numeric)
  if (is.numeric(start_year)) {
    if (start_year %% 4 != 0 | start_year %% 4 == 0 & keep_leap == 0) {
      files <- files[!grepl(paste0(start_year, "0229"), files)]
    }
  }
  
  return(files)
}

# Function to substitute raster values with severe (-2) and moderate (-1) 
# climate stress exclusions (Excl2) in areas under climate stress. Useful for 
# Lifestage and NumGen rasts but could be applied to others as well.
# This function is modified from DDRP_v2 such that the calculation is done 
# within the daily loop. 
Rast_Subs_Excl <- function(brk, tile_num, type) {
  
  # Get AlLEXCL_brick
  if (region_param %in% c("CONUS", "EAST")) {
    AllEXCL_brick <- brick(paste0("All_Stress_Excl_tile", tile_num,".tif"))
  } else {
    AllEXCL_brick <- brick("All_Stress_Excl.tif")
  }
  
  # If PEM, then only take the last layer of allEXCL_brk (last date)
  # because PEMs are created only for the last date
  if (grepl("PEM", names(brk), ignore.case = FALSE)) {
    AllEXCL_brick <- AllEXCL_brick[[last(nlayers(AllEXCL_brick))]] 
  }
  
  # Identify pixels that have moderate (-1) or severe (-2) stress values.
  # The corresponding pixel in the input raster brick will be replaced with a
  # stress value. 
  lapply(1:nlayers(brk), function(lyr) {

    # For each layer in the brick, get the corresponding layer
    # in the All Stress Exclusion brick (AllEXCL_brick)
    # If PEM, then only take the last layer of allEXCL_brk (last date)
    # because PEMs are created only for the last date
    if (deparse(substitute(brk)) == "PEM") {
      brk_lyr <- last(nlayers(AllEXCL_brick))
    } else {
      brk_lyr <- brk[[lyr]]
      Excl <- AllEXCL_brick[[lyr]] 
    }
    
    # Replace pixels in brick layer with the stress values in
    # areas of overlap
    if (type == "Excl2") {
      Excl2 <- Excl < 0   # Both severe and moderate exclusion
      brk_lyr[Excl2] <- Excl[Excl2] 
    } else if (type == "Excl1") {
      Excl1 <- Excl == -2   # Only severe exclusion
      brk_lyr[Excl1] <- Excl[Excl1]
    }
    
    return(brk_lyr)
  })
  
}

# Save raster function
# r = raster, outnam = output file name; datatype = number of digits int the 
# output rasters (see "raster" library specificatoins); log_capt = caption
SaveRaster <- function(r, step, tile_num, outnam, datatype, log_capt) {
  
  if (step == "DailyLoop") {
    
    if (region_param %in% c("CONUS", "EAST")) {
      writeRaster(r, file = paste0(outnam, "_tile", tile_num),
                  format = "GTiff",  datatype = datatype, overwrite = TRUE)
      outnam2 <- paste0("\n\nSaving raster brick: ", 
                        outnam, "_tile_", tile_num, ".tif\n")
      
    } else {
      writeRaster(r, file = outnam, format = "GTiff", 
                  datatype = datatype, overwrite = TRUE)
      outnam2 <- paste0("\n\nSaving raster brick: ", outnam, ".tif\n")
    }
    
  } else if (step == "Processing") {
      writeRaster(r, file = outnam, format = "GTiff", 
                  datatype = datatype, overwrite = TRUE)
    if (nlayers(r) > 1) {
      outnam2 <- paste0("\n\nSaving raster brick: ", outnam, ".tif\n")
    } else {
      outnam2 <- paste0("\n\nSaving raster: ", outnam, ".tif\n")
    }
    
  }
  
  cat(outnam2, str_wrap(paste0(log_capt, "\n"), width = 80, exdent = 2), 
  sep = "", file = Model_rlogging, append = TRUE) 

}

### SplitRas: split raster into tiles
# The function spatially aggregates the original raster
# from https://stackoverflow.com/questions/29784829/
# r-raster-package-split-image-into-multiples
# It turns each aggregated cell into a polygon, then the extent of each polygon 
# is used to crop the original raster. The function returns a list with all the 
# pieces in case you want to keep them in the memory. The arguments are:
# raster = raster to be chopped            (raster object)
# ppside = pieces per side                 (integer)
# save   = write raster                    (TRUE or FALSE)
# plot   = do you want to plot the output? (TRUE or FALSE)
SplitRas <- function(raster, ppside, save, plot) {
  h        <- ceiling(ncol(raster)/ppside)
  v        <- ceiling(nrow(raster)/ppside)
  agg      <- aggregate(raster, fact = c(h, v))
  agg[]    <- 1:ncell(agg)
  agg_poly <- rasterToPolygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for (i in 1:ncell(agg)) {
    e1          <- extent(agg_poly[agg_poly$polis == i,])
    r_list[[i]] <- crop(raster, e1)
  }
  if (save == TRUE) {
    for (i in 1:length(r_list)) {
      ii <- formatC(i, width = 2, format = "d", flag = "0")
      writeRaster(r_list[[i]], filename = paste("SplitRas", ii, sep = ""),
                  format = "GTiff", datatype = "FLT4S", overwrite = TRUE)  
    }
  }
  if (plot == TRUE) {
    par(mfrow = c(ppside, ppside))
    for (i in 1:length(r_list)) {
      plot(r_list[[i]], axes = FALSE, legend = FALSE, bty = "n", box = FALSE)  
    }
  }
  return(r_list)
}

# * Daily Loop funcs ####

# If then else raster function [sim. to GRASS r.mapcalc if (x,a,b)]
Cond <- function(condition, trueValue, falseValue) {
  return(condition * trueValue + (!condition)*falseValue)
}


# Convert raster (TIF file) into a data frame for plotting
ConvDF <- function(rast) {
  spdf <- as(rast, "SpatialPixelsDataFrame")
  df <- as.data.frame(spdf)
  colnames(df) <- c("value", "x", "y")
  return(df)
}

# Classify data into bins so it can be visualized in categories (e.g., 1-10, 11-20, 21-30) 
Cut_bins <- function(df, breaks) {
  df$value_orig <- df$value # Keep old value so can sort factors against it
  # Round up max value to highest number divisible by 10
  df$value[df$value == max(df$value)] <- 10 * ceiling(max(df$value)/10)
  #df2 <- df %>% cut2(df$value, g = 10)
  # Cut values into bins and format results
  df2 <- df %>% mutate(value = cut_interval(df$value, n = breaks),
                       value = gsub("[()]|\\[|\\]", "", value),
                       bin1 = ceiling(as.numeric(str_split_fixed(value, ",", 2)[,1])),
                       bin2 = ceiling(as.numeric(str_split_fixed(value, ",", 2)[,2]))) 
  # Add a 1 to the first bin value unless it's the first bin, which is usually 0
  # The cut_interval and base R cut functions seems to round inappropriately so that
  # an 81.1 is assigned to 82 for example. Not sure how to change it.
  df3 <-  df2 %>% mutate(bin1 = ifelse(bin1 == min(bin1, na.rm = TRUE), bin1, bin1 + 1)) %>%
    mutate(value = paste(bin1, bin2, sep = "-"))
  return(df3)
}

### Function that converts a matrix (= m) to a raster, which involves specifying 
# the extent (= ext) from the template (= template), setting the coordinate system 
# and assigning it a spatial resolution (from the template). Used in daily loop.
Mat_to_rast <- function(m, ext, template) {
  rast <- raster(m, xmn = ext[1,1], xmx = ext[1,2], 
                    ymn = ext[2,1], ymx = ext[2,2])
  crs <- crs(template)
  res(rast) <- res(template)
  NAvalue(rast) <- NaN
  return(rast)
}

# * Summary map funcs ####

# Base features used for all summary (PNG) maps in "PlotMap" function
Base_map <- function(df) {
  p <- ggplot(states, aes(x = long, y = lat)) + 
    geom_raster(data = df, aes(x = x, y = y, fill = value)) + 
    geom_path(aes(group = group), color = "black", lwd = 0.4) +
    coord_quickmap(xlim = c(REGION@xmin, REGION@xmax), 
                   ylim = c(REGION@ymin, REGION@ymax), expand = FALSE) 
  
}

# Creates a gradient color ramp consisting of n colors
Colfunc <- function(x,y,n) {
  colorRampPalette(c(x,y))(n)
}

# Deal with 0 vs non-zero values when plotting chill and heat stress unit rasters
Stress_Val_Conv <- function(x) {
  if (all(x$value == 0)) {
    x2 <- x
  } else {
    x2 <- Cut_bins(x, 10)
  }
  # Need to fix bins if all data are < 10 
  if (all(x2$value_orig < 10)) {
    x2$value <- "0-10"
  }
  x2$value <- factor(x2$value, 
                     levels = unique(x2$value[order(x2$value_orig)]))
  return(x2)
}

# Create and write summary maps (PNG) of all rasters except for stress units
PlotMap <- function(r, d, titl, lgd, outfl) {
  
  dateForlog <- format(as.Date(d, "%Y%m%d"), "%m/%d/%Y")
  
  # Convert raster to data frame except Stage Count results already in df format
  if (grepl("StageCount", outfl)) {
    df <- r 
  } else {
    df <- ConvDF(r)
  }
  
  sp <- paste0(gsub(pattern = "_", replacement = " ", fullname),":")
  nam <- deparse(substitute(r)) # Get name of raster, capitalize first letter
  dat <- as.character(format(strptime(d,format="%Y%m%d"), format="%m/%d/%Y")) # format the date
  titl <- paste(titl,dat,sep=" ")
  titl_orig <- titl # Used for some log captions
  subtitl <- paste("Maps and modeling",
                   format(Sys.Date(), "%m/%d/%Y"), str_wrap("by Oregon State University IPPC 
    USPEST.ORG and USDA-APHIS-PPQ; climate data from OSU PRISM Climate Group",
                                                            width = 150))
  
  # Need to enforce a rule for wrapping title and subtitle on plot
  # The title and subtitle go off of page for some states (VT...any else?)
  if (asp > 1.55) {
    titl_width <- 45
    subtitl_width <- 55
  } else if (asp >= 1.4 & asp < 1.55) {
    titl_width <- 55
    subtitl_width <- 65
  } else {
    titl_width <- 65
    subtitl_width <- 75
  }
  
  # Code for plots will vary depending on the product type, as specified below
  #### DD accumulations ##
  if (grepl("DDtotal|Cum_Inf_Risk|Cum_DDs", outfl)) {
    
    # Log captions
    if (grepl("DDtotal", outfl)) {
      log_capt <- paste("- Number of accumulated degree-days on", dateForlog)
    } else if (grepl("Cum_Inf_Risk", outfl)) {
      log_capt <- paste("- Number of infection risk units on", dateForlog)
    } else if (grepl("Cum_DDs", outfl)) {
      log_capt <- paste("- Number of accumulated plant disease degree-days on", dateForlog)
    }
    if (exclusions_stressunits & grepl("Cum_Inf_Risk", outfl)) {
      log_capt <- paste("- No. of infection risk units with climate stress excl. on", dateForlog)
    }
        
    # Use pre-defined bins
    if (grepl("DDtotal", outfl)) {
      bins <- c("0", "1-250", "251-500", "501-1000", "1001-2000", "2001-3000", 
                "3001-4000", "4001-5000", "5001-6000", "6001-7000", ">7000")
      vals <- c(250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000)
    } else if (grepl("Cum_Inf_Risk", outfl)) {
      bins <- c("0", "1-5", "6-10", "11-20", "21-30", "31-40", "41-50", 
                "51-75", "76-125", "126-200",">200") 
      vals <- c(5, 10, 20, 30, 40, 50, 75, 125, 200)
    } else if (grepl("Cum_DDs", outfl)) {
      bins <- c("0", "1-400", "401-1000", "1001-2000", "2001-3000", "3001-4000", 
                "4001-5000", "5001-7000", "7001-9000", "9001-11000", ">11000")
      vals <- c(400, 1000, 2000, 3000, 4000, 5000, 7000, 9000, 11000)
     }
    
    # Assign bins to data and order factor levels
    df2 <- df %>% 
      mutate(value = case_when(value == -2 ~ "excl.-severe",
                              value == -1 ~ "excl.-moderate",
                              value == 0 ~ bins[1],
                              value > 0 & value <= vals[1] ~ bins[2],
                              value > vals[1] & value <= vals[2] ~ bins[3],
                              value > vals[2] & value <= vals[3] ~ bins[4],
                              value > vals[3] & value <= vals[4] ~ bins[5],
                              value > vals[4] & value <= vals[5] ~ bins[6],
                              value > vals[5] & value <= vals[6] ~ bins[7],
                              value > vals[6] & value <= vals[7] ~ bins[8],
                              value > vals[7] & value <= vals[8] ~ bins[9],
                              value > vals[8] & value <= vals[9] ~ bins[10],
                              value > vals[9] ~ bins[11]) )
    
    # Color palette
    cols <-  colorRampPalette(rev(brewer.pal(11, "Spectral")))(11)
    
    # Set factor levels
    if (grepl("Excl1", outfl)) {
      bins <- factor(c("excl.-severe", bins), levels = c("excl.-severe", bins))
      cols <- c("gray30", cols)
    } else if (grepl("Excl2", outfl)) {
      bins <- factor(c("excl.-severe", "excl.-moderate", bins), 
                     levels = c("excl.-severe", "excl.-moderate", bins))
      cols <- c("gray30", "gray70", cols)
    } else {
      bins <- factor(bins, levels = bins)
    }
    
    # Set data frame factors and set names for color palette
    df2$value <- factor(df2$value, levels = levels(bins))
    cols <- setNames(cols, levels(df2$value))
    
    # Make the plot
    p <- Base_map(df2) + 
      scale_fill_manual(values = cols, name = str_wrap(paste0(lgd), width = 15), 
                        drop = FALSE) +
      #scale_fill_manual(values = cols, name = str_wrap(paste0(lgd), width = 15)) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width), 
           subtitle = str_wrap(subtitl, width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme
    
    #### Lifestage ##
  } else if (grepl("Lifestage", outfl)) {
    
    # Log captions
    if (grepl("Lifestage_Excl1|Lifestage_Excl2", outfl)) { 
      log_capt <- paste("-", titl_orig, "with climate stress excl. on", dateForlog)
    } else {
      log_capt <- paste("-", titl_orig, "on", dateForlog)
    }
    
    # Slightly modify "stg_vals" dataframe for plotting purposes here
    stg_vals2 <- stg_vals %>% # Add "OW" back on to the OW stage
      mutate(stg_name = ifelse(stg_num == 1, paste("OW", stg_name), stg_name))
    #stg_vals2$stg_name <- factor(stg_vals2$stg_name, 
    #                             levels = unique(stg_vals2$stg_name[order(stg_vals2$stg_num)]))
    all_stages <- unique(stg_vals2$stg_name)
    
    # Replace stage numbers w/ stage names in results
    df$value <- as.character(df$value)
    df$value_orig <- df$value
    df <- left_join(df, stg_vals2, by = c("value" = "stg_num")) %>%
      mutate(stg_name = ifelse(value == "-2", "excl.-severe", 
                            ifelse(value == "-1", "excl.-moderate", stg_name)),
             value = stg_name) # Make value the stg_name column
    
    # Stage colors
    cols <- c("blue2", "orange2", "purple2", "sand", "teal")
    col_key <- data.frame(cols = c("blue2", "gold2", "purple2", "cyan2", "red2"), 
                          value = unique(stg_vals2$stg_name))

    #cols <-  colorRampPalette(rev(brewer.pal(5, "Spectral")))(5)
    
    # Set factor levels for values
    if (grepl("Excl1", outfl)) {
      col_key <- col_key %>% add_row(cols = "gray30", value = "excl.-severe")
      levls <- c("excl.-severe", unique(stg_vals2$stg_name))
    } else if (grepl("Excl2", outfl)) {
      col_key <- col_key %>% add_row(cols = c("gray30", "gray70"), 
                                     value = c("excl.-severe", "excl.-moderate"))
      levls <- c("excl.-severe", "excl.-moderate", unique(stg_vals2$stg_name))
    } else {
      levls <- unique(stg_vals2$stg_name)
    }
    df$value <- factor(df$value, levels = levls)

    # Make the plot
    cols <- setNames(as.character(col_key$cols), col_key$value) # Final colors
    p <- Base_map(df) + 
      scale_fill_manual(values = cols, 
                        name = str_wrap(paste0(lgd), width = 15), drop = FALSE) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width), 
           subtitle = str_wrap(subtitl, width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme
    
    #### No. of gens ##
  } else if (grepl("NumGen|NumGenExcl1|NumGenExcl2", outfl)) {
    
    # Log captions
    if (grepl("NumGen_Excl1|NumGen_Excl2", outfl)) { 
      log_capt <- paste("-", str_wrap("Number of gens. with climate stress excl. 
                                      on", width = 80), dateForlog)
    } else {
      log_capt <- paste("-", "Numbers of generations on", dateForlog)
    }
    
    # Order by original values so plots in numerical order
    df$value_orig <- df$value
    df <- mutate(df, value = ifelse(value == -1, "excl.-moderate", 
                                    ifelse(value == -2, "excl.-severe", value)))
    df$value <- factor(df$value, 
                       levels = unique(df$value[order(df$value_orig)]))
    
    # Make the color key for the legend, and keep only needed colors
    col_key <- suppressWarnings(suppressMessages(data.frame(
      cols = c("gray70", "gray30", "blue3", "firebrick4", "gold3", "darkgreen", 
               "magenta4", "darkorange3", "cyan4", "chartreuse4", "purple3",  "gold4", 
               "cornflowerblue", "palevioletred2", "seagreen4",  "lemonchiffon",
               "darkred", "deeppink3", "lightgreen", "sienna4", "maroon4", 
               "royalblue4",  "thistle4"),
      value = c("excl.-moderate", "excl.-severe", 0:20)) %>%
        semi_join(., df, by = "value")))
    col_key$value <- factor(col_key$value, levels = levels(df$value))
    cols <- setNames(as.character(col_key$cols), col_key$value) 
    
    # Make the plot
    p <- Base_map(df) + 
      scale_fill_manual(values = cols, name = str_wrap(paste0(lgd), width = 15)) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width), 
           subtitle = str_wrap(subtitl, width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme
    
    #### Stage count ##
  } else if (grepl("StageCount", outfl)) {
   # Caption for logging file
      if (grepl("StageCount_Excl1|StageCount_Excl2", outfl)) { 
        log_capt <- paste("-", 
          str_wrap("Stage count with climate stress excl. on", width = 80), dateForlog)
      } else {
        log_capt <- paste("-", "Stage count on", dateForlog)
      }    
      
    # Format data if there are climate stress exclusion values
      if (exclusions_stressunits) {
        df <- df %>% mutate(value = ifelse(value_orig == -2, "excl.-severe", 
                        ifelse(value_orig == -1, "excl.-moderate", value)))
        
        if (any(df$value_orig > 0)) {
          # Need to remove -2 and -1 values prior to binning values for plot
          df2 <- dplyr::filter(df, !value_orig < 0) 
          excl_df <- df %>% dplyr::filter(value_orig < 0) # Take only -2 and -1 
          # Put exlcusion values and stage count values back together
          df <- rbind(df2, excl_df) # Rename data frame
        
          # If clim. exclusions masks out all non-zero values, then just plot 
          # climate stress exclusions
        } else if (all(df$value_orig <= 0)) {
          df <- mutate(df, value = ifelse(value_orig == -2, "excl.-severe", 
                            ifelse(value_orig == -1, "excl.-moderate", value)))
            
        }
        
        # If stress values are missing in data, then add a row so the legend
        # still shows the stress category (Excl1 = excl.-severe, Excl2 = 
        # excl.-severe and excl.-moderate). Otherwise just recode stress values.
        if (grepl("StageCount_Excl1", outfl) & (!(-2 %in% df$value_orig))) {
          df <- df %>% 
            add_row(value = "excl.-severe", gen_stg = factor(-2))
        } else if (grepl("StageCount_Excl2", outfl)) {
          if (!(-2 %in% df$value_orig)) {
            df <- df %>% 
              add_row(value = "excl.-severe", gen_stg = factor(-2))
          }
          if (!(-1 %in% df$value_orig)) {
            df <- df %>% 
              add_row(value = "excl.-moderate",  gen_stg = factor(-1))
          }
        }
        df <- arrange(df, as.numeric(gen_stg))
      }
      
      # Define factor levels to order legend key properly
      sorted <- unique(as.numeric(df$gen_stg))  
      df$gen_stg <- factor(df$gen_stg, levels = sorted)
      
      # Make the color key for the legend 
      # Currently enough colors for 20 generations
      cols_df <- data.frame("cols" = 
        c(Colfunc("deepskyblue", "blue3", 4), # Gen 0 
          Colfunc("orangered", "firebrick4", 4), # Gen 1
          Colfunc("yellow", "gold3", 4), # Gen 2
          Colfunc("lightgreen", "darkgreen", 4), # Gen 3
          Colfunc("magenta", "magenta4", 4), # Gen 4
          Colfunc("tan1", "darkorange3", 4), # Gen 5
          Colfunc("cyan", "cyan4", 4), # Gen 6
          Colfunc("greenyellow", "chartreuse4", 4), # Gen 7
          Colfunc("mediumpurple1", "purple3", 4), # Gen 8
          Colfunc("lightgoldenrod", "gold4", 4), # Gen 9
          Colfunc("cadetblue1", "cornflowerblue", 4), # Gen 10
          Colfunc("mistyrose", "palevioletred2", 4), # Gen 11
          Colfunc("seagreen1", "seagreen4", 4), # Gen 12
          Colfunc("lemonchiffon", "gold", 4), # Gen 13
          Colfunc("red", "darkred", 4), # Gen 14
          Colfunc("lightpink", "deeppink3", 4), # Gen 15
          Colfunc("lightgreen", "darkgreen", 4), # Gen 16
          Colfunc("sienna", "sienna4", 4), # Gen 17
          Colfunc("maroon1", "maroon4", 4), # Gen 18
          Colfunc("royalblue1", "royalblue4", 4), # Gen 19
          Colfunc("thistle1", "thistle4", 4))) # Gen 20
        gens_df <- data.frame(gen = c(rep(1, 4), rep(2, 4),
          rep(3, 4), rep(4, 4), rep(5, 4), rep(6, 4), rep(7, 4), rep(8, 4), 
          rep(9, 4), rep(10, 4), rep(11, 4), rep(12, 4), rep(13, 4), rep(14, 4),
          rep(15, 4), rep(16, 4), rep(17, 4), rep(18, 4), rep(19, 4), 
          rep(20, 4)))
        gens_df$gen <- sapply(gens_df$gen, function(x) {
          paste(toOrdinal(x), "gen.") 
          })
        gens_df <- rbind(data.frame(gen = c(rep("OW gen.", 4))), gens_df)
        
        # Create the color key and named vector of the key
        #col_key <- cols_df %>% mutate(gen_name = paste0("G", gens_df$gen)) %>%
        col_key <- cols_df %>% mutate(gen = gens_df$gen) %>%
          mutate(stg_name = 
            rep_len(c("eggs", "larvae", "pupae", "adults"), nrow(gens_df))) %>%
          mutate(value = paste(gen, stg_name)) %>%
          semi_join(., df, by = "value") %>%
          dplyr::select(cols, value)
        
        # Order values by generationa nd life cycle stage order (i.e. egg, 
        # larvae, pupae, adult) or won't show up correctly in legend key
        df$value <- factor(df$value, 
                  levels = unique(df$value[order(df$gen_stg)]))
        
        # Add grayscale colors to legend colors if climate stress exclusions
        # Moderate stress exclusion
        if (any(df$value == "excl.-moderate")) {
          col_key <- rbind(data.frame("cols" = "gray70", 
                                       "value" = "excl.-moderate"), col_key)
        }
        # Severe stress exclusions
        if (any(df$value == "excl.-severe")) {
          col_key <- rbind(data.frame("cols" = "gray30", 
                                      "value" = "excl.-severe"), col_key)
        }
        
      # Make the plot
      cols <- setNames(as.character(col_key$cols), levels(df$value)) 
      p <- Base_map(df) + 
        scale_fill_manual(values = cols, 
                            name = str_wrap(paste0(lgd), width = 15)) +
        labs(title = str_wrap(paste(sp, titl), width = titl_width), 
              subtitle = str_wrap(subtitl, width = subtitl_width)) +
        theme_map(base_size = base_size) + 
        mytheme    
    # Need to adjust number of rows in legend for very small plots or 
    # the legend will go off the page
    if (asp < 0.5) {
      p <- p + guides(fill = guide_legend(nrow = 15))
    }
    
    #### Climate stress exclusion maps ##
  } else if (grepl("Heat_Stress_Excl|Cold_Stress_Excl|Dry_Stress_Excl|All_Stress_Excl", outfl)) {
    
    # Log captione
    log_capt <- paste("-", titl_orig, dateForlog)

    # Re-assign values (0, -1, -2) to their corresponding description
    df <- mutate(df, value = factor(ifelse(value == -2, "excl.-severe", 
                                           ifelse(value == -1, "excl.-moderate", 
                                                  ifelse(value == 0, "not excluded",NA)))))    
    df$value <- factor(df$value, levels = c("excl.-severe", "excl.-moderate",
                                            "not excluded")) # Order by values
    # Make the plot
    p <- Base_map(df) +        
      scale_fill_manual(
        values = c("excl.-severe" = "gray30","excl.-moderate" = "gray70",
                   "not excluded" = "green2"), 
        name = paste0(lgd), drop = FALSE) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width),
           subtitle = str_wrap(subtitl, width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme + 
      theme(legend.text = element_text(size = rel(1.5)), 
            legend.title = element_text(size = rel(1.4), face = "bold"))
    
    #### Pest Event Maps ##
  } else if (grepl("PEM", outfl)) {
    
    # Caption for log file
    log_capt <- paste("-", titl_orig) 
    
    # Optimize code in this section - probably more efficient ways to do
    # Make a start year in numeric format
    start_year.num <- as.numeric(start_year) 
    
    # Format the data value column
    df <- df %>% 
      dplyr::filter(!(value %in% c(0, 366))) # remove day 0 and day 366
    #df$value <- round(df$value)
    # TO DO: try to fix this so that don't end up with 6 weeks in Dec
    # Currently must change value for day 364 and 365 to day 363, or end up w/
    # 6 weeks in December for some years
    df <- df %>% mutate(value = ifelse(value >= 364, 363, value))
    df$value_orig <- df$value # Create for ordering recoded values later
    
    # Convert day of year to a date, convert date to a week of the year, and 
    # format to "month-day" format. First need to subtract 1 from all day of 
    # year values, because as.Date starts at day 1, not day 0
    df$value <- df$value - 1
    df$value <- as.Date(df$value, origin = 
                          as.Date(paste0(start_year.num, "-01-01"))) 
    
    # The resulting values may have dates before Jan-1 (Dec-30, Dec-31) because 
    # they occur in the same week as Jan-1. The ceiling_date function (dplyr)
    # rounds up to the next month so that they begin on Jan-1 
    # (e.g., 2019-12-31 == 2020-01-01)
    df$value <- as.character(cut.POSIXt(strptime(df$value, format = "%Y-%m-%d"), 
                                        breaks = "1 weeks"))
    badDates <- df %>% filter(grepl(as.character(start_year.num - 1), value))
    badDates$value <- as.character(ceiling_date(as.Date(badDates$value, 
                                                        format = "%Y-%m-%d"), "month"))
    
    # Now replace old data with data that have fixed date
    # Then add week of month column
    df <- df %>% filter(!(grepl(as.character(start_year.num - 1), value))) %>%
      bind_rows(., badDates)    
    df$week <- ceiling(mday(df$value)/7)
    
    # Reformat the dates to month-day (e.g., Jan-01, Jan-06, ...)
    # Clim. exc. values (-1 and -2) will become NA
    df$value <- format(strptime(df$value, format = "%Y-%m-%d"), 
                       format = "%b-%d")
    
    # Data frames w/ January may have 2 dates (Jan-01 + another) for week 1
    # This will result in the key having 2 dates w/ same color
    # Check to see if this is true, and if so, then add 1 week onto weeks for
    # January
    week1 <- df %>% filter(week == 1 & value_orig < 30) %>% distinct(value)
    
    if (nrow(week1) > 1) {
      df <- df %>% 
        mutate(week = ifelse(grepl("Jan-", value) & !grepl("Jan-01", value), 
                             week + 1, week))
    }
    
    # Generate a key for colors for every week of the year, allowing up to 
    # 5 weeks per month
    cols_df <- data.frame("cols" = 
                            c(Colfunc("deepskyblue", "blue3", 5),
                              Colfunc("red", "darkred", 5), 
                              Colfunc("yellow", "gold3", 5),
                              Colfunc("lightgreen", "darkgreen", 5), 
                              Colfunc("magenta", "magenta4", 5),
                              Colfunc("sienna1", "sienna4", 5),
                              Colfunc("cyan", "cyan4", 5),
                              Colfunc("greenyellow", "chartreuse4", 5),
                              Colfunc("mediumpurple1", "purple3", 5),
                              Colfunc("lightpink", "deeppink4", 5),
                              Colfunc("lightgoldenrod", "gold4", 5), 
                              Colfunc("cadetblue1", "cornflowerblue", 5))) # Colors
    weeks_df <- data.frame("mnth" = c(rep("Jan", 5), rep("Feb", 5), rep("Mar", 5), 
                                      rep("Apr", 5), rep("May", 5), rep("Jun", 5), 
                                      rep("Jul", 5), rep("Aug", 5), rep("Sep", 5), 
                                      rep("Oct", 5), rep("Nov", 5), rep("Dec", 5))) # 5 weeks per month
    weeks_df <- data.frame(weeks_df %>% group_by(mnth) %>% # Group by month
                             mutate(mnth_wk = row_number()) %>% # Assign unique row # to rep. week #
                             mutate(mnth_wk = paste(mnth, mnth_wk, sep = "_")))
    
    # Attach those data frames to make the key
    col_key <- cbind(cols_df, weeks_df)
    col_key$mnth <- as.character(col_key$mnth) # Add which month
    
    # Extract all unique weeks from data,and count no. of bins (months)
    dats.1 <- df %>% distinct(value) 
    dats.1$mnth <- str_split_fixed(dats.1$value, pattern = "-", 2)[,1]
    dats.1 <- dats.1 %>% arrange(., value) %>%
      group_by(mnth) %>% # Group by month
      # Assign a unique value to each row in a month - used for joining later
      # This will also be done for other data frames below 
      #mutate(week = ceiling(day(mnth_day) / 7)) %>%
      left_join(dplyr::select(df, value, week), by = "value") %>%
      mutate(mnth_day = format(
        as.Date(paste0(start_year, "-", value), 
                format = "%Y-%b-%d"), format = "%b-%d")) %>%
      mutate(mnth_wk = paste(mnth, week, sep = "_")) %>%
      distinct(., .keep_all = TRUE) %>%
      arrange(mnth_day)
      
    # Necessary for removing weeks that are not in the data in the col_key2
    mnth_ct <- data.frame(dats.1 %>% group_by(mnth) %>% 
                            dplyr::mutate(freq = n_distinct(value))) %>%
      group_by(mnth) #%>% # Group by month
    
    # Finally filter unncessary weeks out of generic color key (some months 
    # have only 4 weeks)
    col_key2 <- dplyr::semi_join(col_key, dats.1, by = "mnth") %>% 
      dplyr::left_join(., mnth_ct, by = c("mnth_wk")) %>%
      na.omit %>%
      dplyr::select(cols, value)
    
    # Format the dates dataframe (dats2) for joining to col_key2 (color key)
    dats.2 <- data.frame(dplyr::select(dats.1, value, mnth) %>% 
                          arrange(mnth, value))
    
    # Attach the colors to the value and format with needed colunms, etc.
    col_key2 <- left_join(col_key2, dplyr::select(dats.2, -mnth), by = "value")
    col_key2$year <- start_year
    col_key2$date <- paste0(col_key2$year, "-", col_key2$value)
    col_key2$date <- as.Date(col_key2$date, format = "%Y-%b-%d")
    col_key2 <- col_key2 %>% dplyr::select(value, cols)
    
    # If data have climate stress exclusion values, need to reformat data
    # because values are non-dates; recode them, then bind back to original
    # data. Finally add grayscale shades to the color key for legend.
    if (any(df$value_orig < 0)) {
      df2 <- filter(df, value_orig < 0)
      df2 <- mutate(df2, value = ifelse(df2$value_orig == -2, "excl.-sev.", 
                                        ifelse(df2$value_orig == -1, "excl.-mod.", 
                                               df2$value)))
      df3 <- rbind(filter(df, value_orig > 0), df2)
      # Order according to orig. vals (DOY)
      df3$value <- factor(df3$value, 
        levels = unique(df3$value[order(as.numeric(as.character(df3$value_orig)))]))
    } else {
      df3 <- df
      df3$value <- factor(df3$value, 
        levels = unique(df3$value[order(as.numeric(as.character(df3$value_orig)))]))
    }
    if (any(df$value_orig == -1)) {
      col_key2 <- rbind(data.frame("value" = "excl.-mod.", "cols" = "gray70"), col_key2)
    }
    if (any(df$value_orig == -2)) {
      col_key2 <- rbind(data.frame("value" = "excl.-sev.", "cols" = "gray30"), col_key2)
    }
    
    # Make legend colors df a vector and plot results
    cols <- setNames(as.character(col_key2$cols), col_key2$value)
    p <- Base_map(df3) + 
      scale_fill_manual(values = cols, name = str_wrap(paste0(lgd), width = 15)) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width), 
           subtitle = str_wrap(paste(subtitl), width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme
    
    # Need to adjust number of rows in legend for small aspect plots or 
    # the legend will go off the page
    if (asp < 0.6) {
      p <- p + guides(fill = guide_legend(nrow = 15))
    }
    
  }
  
  #### Save the plots ##
  # Save the plot, or else report that there was an error and skip
  # See "rmessages.txt" for error report
  # save the plot, or else report that there was an error and skip 
  # See "rmessages.txt" for error report
  tryCatch(
    {
      suppressMessages(ggsave(p, file = paste0(outfl, "_", d, ".png"), 
                              height = asp * 7, units = c('in'), dpi = 300))
      cat(paste0("\n\nSaving summary map: ", outfl, "_", d, ".png\n"),  
          str_wrap(paste0(log_capt, "\n"), width = 80, exdent = 2), sep = "",
          file = Model_rlogging, append = TRUE)    
      },
    error=function(e){
      LogWrap("\n", paste0("Could not create plot for ", outfl), "\n")
      print(noquote(paste0("Could not create plot for ", outfl)))
    } )
}


# Create summary maps (PNG) of heat stress and chill stress units, with max1 
# (Stress limit 1) and max2 (Stress limit 2) shown as "countour" lines
PlotMap_stress <- function(r, d, max1, max2, titl, lgd, outfl) {
  log_capt <- paste("- Number of accumulated", tolower(lgd), "on", 
                       format(as.Date(d, "%Y%m%d"), "%m/%d/%Y")) # Log file cap
  sp <- paste0(gsub(pattern = "_", replacement = " ", fullname), ":")
  dat <- as.character(format(strptime(d, format = "%Y%m%d"), 
                             format = "%m/%d/%Y"))
  titl <- paste(titl, dat, sep = " ")
  subtitl <- paste("Maps and modeling", format(Sys.Date(), "%m/%d/%Y"), 
                   str_wrap("by Oregon State University IPPC USPEST.ORG and USDA-APHIS-PPQ; 
             climate data from OSU PRISM Climate Group",  width = 150))  
  df <- ConvDF(r)
  df$value_orig <- df$value
  df2 <- Stress_Val_Conv(df) # Properly formats values
  
  # Need to wrap title and subtitle for narrow plots (e.g., RI)
  if (asp > 1.55) {
    titl_width <- 45
    subtitl_width <- 55
  } else {
    titl_width <- 55
    subtitl_width <- 75
  }
  
  # Create contours from raster values greater than limit 1 and limit 2, if 
  # raster values are greater than max1 and/or max2
  max1_c <- tryCatch(
    {
      # max_1 is object of class "SpatialLinesDataFrame"
      max1_c <- rasterToContour(r > max1) 
    },
    error = function(e) {
      max1_c <- 0 
    } )
  max2_c <- tryCatch(
    {
      # max_2 is object of class "SpatialLinesDataFrame"
      max2_c <- rasterToContour(r > max2) 
    },
    error = function(e) {
      max2_c <- 0 
    } )
  
  # If all values are 0 or are less than 10 (if stress limit < 10), then don't 
  # include contours (must include this code or it will throw an error)
  if (all(df$value == 0 | all(df$value < 10 & all(df$value < max1)))) {
    p <- Base_map(df2) +
      scale_fill_manual(values = c("#5E4FA2"), name = paste0(lgd)) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width), 
           subtitle = str_wrap(subtitl, width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme 
    
    # If any values are greater than 0, but not greater than max1, then don't 
    # plot either countour (must include or it will thrown an error)
  } else if (any(df$value > 0) & is.numeric(max1_c)) {
    p <- Base_map(df2) +
      scale_fill_brewer(palette = "Spectral", direction = -1, 
                        name = str_wrap(paste0(lgd), width = 15)) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width), 
           subtitle = str_wrap(subtitl, width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme
    
    # If any values are greater than limit 1 but less than limit 2, then 
    # plot countour line for just limit 1
  } else if (any(df$value > 0) & !is.numeric(max1_c) & is.numeric(max2_c)) { 
    # max1_c is class "SpatialLinesDataFrame" but max2_c is class "numeric" 
    # (max2_c = 0)
    p <- Base_map(df2) +
      scale_fill_brewer(palette = "Spectral", direction = -1, 
                        name = str_wrap(paste0(lgd), width = 15)) +
      scale_color_manual(name = "Stress Limits", 
                         values = c("Stress limit 1" = "magenta")) +
      geom_path(data = max1_c, aes(x = long, y = lat, group = group, 
                                   color = "Stress limit 1"), lwd = 0.15) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width), 
           subtitle = str_wrap(subtitl, width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme +
      guides(colour = guide_legend(override.aes = list(size = 1), order = 1)) 
    
    # If any values are greater than limit1 and limit2, then plot 
    # countour lines for both limit1 and limit2
  } else if (!is.numeric(max1_c) & !is.numeric(max2_c)) {
    p <- Base_map(df2) +
      scale_fill_brewer(palette = "Spectral", direction = -1, 
                        name = str_wrap(paste0(lgd), width = 15)) +
      scale_color_manual(name = "Stress Limits", 
                         values = c("Stress limit 1" = "magenta", 
                                    "Stress limit 2" = "mediumblue")) +
      geom_path(data = max1_c, aes(x = long, y = lat, group = group, 
                                   color = "Stress limit 1"), lwd = 0.15) +
      geom_path(data = max2_c, aes(x = long, y = lat, group = group, 
                                   color = "Stress limit 2"), lwd = 0.15) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width), 
           subtitle = str_wrap(subtitl, width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme +
      guides(colour = guide_legend(override.aes = list(size = 1), order = 1))
  }
  
  # Save the plot, or else report that there was an error and skip - 
  # See "rmessages.txt" for error report
  tryCatch(
    {
      ggsave(p, file = paste0(outfl, "_", d, ".png"), height = 7 * asp, 
             units = c('in'), dpi = 200) 
      cat(paste0("\n\nSaving summary map: ", outfl, "_", d, ".png\n"),  
          str_wrap(paste0(log_capt, "\n"), width = 80, exdent = 2), sep = "",
          file = Model_rlogging, append = TRUE)    
    },
    error = function(e) {
      LogWrap("\n", paste0("Could not create plot for ", outfl), "\n")
      print(noquote(paste0("Could not create plot for ", outfl)))
    } )
  
}


#### (3). DIRECTORY INIT ------

#### * Param inputs - species params; thresholds, weather, etc. ####
params_dir <- "/usr/local/dds/DDRP_B1/spp_params/"

#### * Weather inputs and outputs - climate data w/subdirs 4-digit year ####
# If outdir has 2 consec. numbers, assume webuser; otherwise just use base dir
if (grepl("16", start_year, perl = TRUE)) {
  base_dir <- "/mnt/ssd1/PRISM/"
} else {
  base_dir <- "/data/PRISM/"
}

prism_dir <- paste0(base_dir, start_year)

cat("\nBASE DIR: ", base_dir, "\n")
cat("\nWORKING DIR: ", prism_dir, "\n")

#### * Output directory, log file, and error message file ####
# MUST remove .tif files or script will crash during processing because it will 
# try to analyze previously processed results. 

#output_dir <- paste0("/home/httpd/html/CAPS/", output_dir)
output_dir <- paste0("/usr/local/dds/DDRP_B1/DDRP_results/", out_dir)

# If the directory already exists, then a backup directory will be created that
# contains the old run files. Old backup directories will be removed if present.
if (file.exists(output_dir)) {
  
  # Create a directory for the previous run to backup old run files to and 
  # delete old backup (if present)
  backup_dir <- paste0(output_dir, "/previous_run")
  unlink(backup_dir, recursive = TRUE)
  dir.create(backup_dir)
  
  # Copy dirs and files to backup dir
  flsToCopy <- list.files(output_dir, full.names = TRUE)
  flsToCopy <- str_subset(flsToCopy, pattern = backup_dir, negate = TRUE)
  file.copy(flsToCopy, backup_dir, recursive = TRUE)
  
  # Delete remaining folders and files now that they have been backed up
  # First ensure that flsToCopy actually contains output_dir in names to 
  # avoid potential fiascos of deleting files from other non-target directories
  if (all(grepl(output_dir, flsToCopy))) {
    lapply(flsToCopy, unlink, recursive = TRUE)
  }

  cat("\n", str_wrap(paste0("EXISTING OUTPUT DIR: ", output_dir, "\n", 
                            "moved old run files to", backup_dir, "\n"), 
                     width = 80), sep = "") 
  
  # If no files to backup, don't do anything
} else {
  dir.create(output_dir)
  cat("NEW OUTPUT DIR:", output_dir)
}

# Push out a rlogging file with all main messages in model
# Put all log, message, and metadata files in a separate folder
setwd(output_dir)
dir.create("Logs_metadata")
Model_rlogging <- sprintf("%s%s", "./", "/Logs_metadata/Model_rlogging.txt")

# Make header for logging file
cat(paste0(rep("#", 39), collapse = ""), "\n", 
    "### Log file for Plant Disease DDRP ###\n", 
    paste0(rep("#", 39), collapse = ""), "\n\n", sep = "", 
    file = Model_rlogging)

# Record PRISM and output dir
cat("BASE DIR: ", base_dir, "\n", file = Model_rlogging, append = TRUE)
cat("WORKING DIR: ", prism_dir, "\n", file = Model_rlogging, append = TRUE)

# Push out a message file with all R error messages
#msg <- file(paste0("/home/httpd/html/CAPS/", output_dir), open = "wt")
msg <- file(paste0(output_dir, "/Logs_metadata/rmessages.txt"), open = "wt")
sink(msg, type = "message")

#### (4). PARAMETER AND SETTINGS SETUP ----- 
LogWrap("", "PARAMETER AND SETTINGS SETUP: getting parameters to use in model", "\n")
cat("\n\nPARAMETER AND SETTINGS SETUP: getting parameters to use in model\n")

# Read from source param files in ./spp_params/SPP.params
param_file <- sprintf("%s%s", spp, ".params")
spp <- gsub(".params", "", param_file) # Get species abbr.
species_params <- sprintf("%s%s", params_dir, param_file) # Location of file

if (file.exists(species_params)) {
  LogWrap("", paste0("Species params: ", species_params), "\n")
  source(species_params) # Read in species parameters
  LogWrap("", paste0("Reading params for species: ", spp, "\nFullname: ", fullname), "\n")
  cat("\nReading params for species: ", spp, " Fullname: ", fullname, "\n")
} else {
  LogWrap("", paste0("Param file: ", species_params, "...not found; exiting program"), "\n")
  cat("\nParam file: ", species_params, "...not found; exiting program\n")
  q()  # No reason to keep going without any params
}

# Change year to numeric if it's a specific year
# If using climate normals, there may be letters in folder name
if (!grepl("[A-z]", start_year)) {
  start_year <- as.numeric(start_year)
} else {
  start_year <- 1975 # 
}

# Set up start and stop day of year depending on whether it's a leap year or
# not (need to modify the last day of the year if the year is a leap year 
# and user wants to include it)
# This does not apply to 30 yr climate data, which would not have a numeric
# start year
if (is.numeric(start_year)) {
  
  if (start_year %% 4 == 0 & keep_leap == 1) {
    LogWrap("", paste0(start_year, " is a leap year and leap day (2/29) will be 
                 included in the model"), "\n")
    
    # Need to add an extra day onto year if all 365 days are being included
    if (end_doy == 365) {
      end_doy <- 366
    }
    
  } else if (start_year %% 4 == 0 & keep_leap == 0) {
    LogWrap("", paste0(start_year, " is a leap year but leap day (2/29) will 
                        not be included in the model"), "\n")
  } else if (start_year %% 4 != 0) {
    LogWrap("", paste(start_year, "is not a leap year - ignoring 'keep_leap' 
                       parameter"), "\n")
  }
}

# Check for appropriate command line parameters
# Exit program if no pest event maps have been specified but users want pest
# event maps (pems = 1)
if (pems == 1 & !(1 %in% c(mapA, mapE, mapL, mapP))) {
  LogWrap("\n", "No pest event maps (mapA, mapE, mapL, mapP) specified; 
                     exiting program", "")
  cat("\n", str_wrap("No pest event maps (mapA, mapE, mapL, mapP) specified; 
          exiting program", width = 80), sep = "")
  q()
}

# Exit program if an incorrect sampling frequency has been specified
if (out_option %in% !c(1, 2, 3, 4, 5, 6)) {
  LogWrap("", "Out_option =", out_option, "is unacceptable; exiting program", "\n")
  cat("Out_option =", out_option, "is unacceptable; exiting program\n")
  q() 
}

# Exit if end day of year is inappropriate
if (end_doy > 366) {
  LogWrap("\n", paste("End day of year (end_doy) of", end_doy, "is 
                     unacceptable; exiting program"), "")
  cat("\n", str_wrap(paste("End day of year (end_doy) of", end_doy, "is 
                     unacceptable; exiting program"), width = 80), sep = "")
  q()
}

# Create a list of days to use for daily loop
sublist <- start_doy:end_doy

#### * Format threshold and DD params for Daily Loop ####

# Need to match length and order of stgorder, which is species specific

# Upper and lower thresholds
# OW stage will have the same threshold as actual stage (e.g., OWadult = adult)
# Need to match LDT or UDT value to stage, which first requires changing 
# stage abbr to the stage name ("E" = "egg")
stage_ldt_list <- list()
stage_udt_list <- list()
j <- 1

for (i in 1:length(stgorder)) {
  stg_nam <- stgorder[i]
  stg_nam <- mgsub(string = stg_nam, # Requires "mgsub" package
                   pattern = c("OE", "OL", "OP", "OA", "E", "L", "P", "A"), 
                   replacement = c("egg", "larvae", "pupae", "adult", "egg", 
                                   "larvae", "pupae", "adult"))
  stage_ldt_val <- get(paste0(stg_nam, "LDT")) # returns LDT value for stage
  stage_ldt_list[[j]] <- stage_ldt_val
  stage_udt_val <- get(paste0(stg_nam, "UDT")) # returns UDT value for stage
  stage_udt_list[[j]] <- stage_udt_val
  j <- j + 1
}

# DD parameters - OW stage has it's own DD, so change stage abbr. to stage 
# name or "OW" plus stage name
stage_dd_list <- list()
j <- 1

for (i in 1:length(stgorder)) {
  stg_nam <- stgorder[i]
  stg_nam <- mgsub(string = stg_nam, 
                   pattern = c("OE", "OL", "OP", "OA", "E", "L", "P", "A"), 
                   replacement = c("OWegg", "OWlarvae", "OWpupae", "OWadult", 
                                   "egg", "larvae", "pup", "adult"))
  stage_dd_val <- get(paste0(stg_nam, "DD")) # returns DD value for stage
  stage_dd_list[[j]] <- stage_dd_val
  j <- j + 1
}

# Put the values in the list into a numeric vector
stage_ldt <- as.numeric(do.call(cbind, stage_ldt_list))
stage_udt <- as.numeric(do.call(cbind, stage_udt_list))
stage_dd <- as.numeric(do.call(cbind, stage_dd_list))

#### (5). METADATA OUTPUT FILE -----

# Push out a metadata file with all inputs used in model
LogWrap("\n", "METADATA: creating metadata file for all inputs used in model", "\n")
cat("\nMETADATA: creating metadata file for all inputs used in model\n")

# Create the metadata file
setwd(output_dir)
metadata <- sprintf("%s%s", "./", "/Logs_metadata/metadata.txt")
cat("### Metadata ###\n", file = metadata)

# Document run date and time
cat("\nRun date and time:", strftime(Sys.time(), format = "%m/%d/%Y %H:%M"),
    file = metadata, append = TRUE)

# Document species information and method used to calculate degree-days
cat("\n\n### Model Species Parameters ###\n Species Abbrev:", spp, 
    "\n Full Name:", fullname, 
    "\n Pest of:", pestof,
    "\n Overwintering Stage:", owstage, 
    "\n Degree-day calculation method:", calctype, 
    file = metadata, append = TRUE)

# Document developmental threshold temperatures
cat("\n \n Developmental threshold temperatures",
    "\n Egg Lower Devel Threshold:", eggLDT, 
    "\n Egg Upper Devel Threshold:", eggUDT, 
    "\n Larvae Lower Devel Threshold:", larvaeLDT, 
    "\n Larvae Upper Devel Threshold:", larvaeUDT, 
    "\n Pupae Lower Devel Threshold:", pupaeLDT,
    "\n Pupae Upper Devel Threshold:", pupaeUDT, 
    "\n Adult Lower Devel Threshold:", adultLDT, 
    "\n Adult Upper Devel Threshold:", adultUDT, file =  metadata, append = T)

# Document stage durations
cat("\n\n Stage durations in degree-days (DDs)",
    "\n Egg DDs:", eggDD, 
    "\n Larvae DDs", larvaeDD, 
    "\n Pupae DDs:", pupDD, 
    "\n Adult DDs:", adultDD, 
    file = metadata, append = TRUE)

# Document climate stress exclusion parameter values, if applicable
if (exclusions_stressunits) {
  cat("\n \n Climate stress parameters",
      "\n Lower Cold Threshold:", coldstress_threshold, 
      "\n Upper Heat Threshold:", heatstress_threshold,
      "\n Upper Dry Threshold:", drystress_threshold,
      "\n Max Cold Units (lower bound):", coldstress_units_max1, 
      "\n Max Cold Units (upper bound):", coldstress_units_max2,
      "\n Max Heat Stress Units (lower bound):", heatstress_units_max1,
      "\n Max Heat Stress Units (upper bound):", heatstress_units_max2,
      # For BOXB
      "\n Max Dry Units (lower bound):", drystress_units_max1,
      "\n Max Dry Units (upper bound):", drystress_units_max2,
      file = metadata, append = TRUE)
}

# Document Pest Event Map parameter values, if applicable
if (pems) {
  cat("\n \n Pest Event Map parameters",
      "\n Number of generations to make Pest Event Maps (PEMs): ", PEMnumgens,
      "\n Egg Event DDs and Label: ", eggEventDD, " (", eggEventLabel,")", 
      "\n Larvae Event DDs and Label: ", larvaeEventDD, " (", larvaeEventLabel, ")",
      "\n Pupae Event DDs and Label: ", pupaeEventDD, " (", pupaeEventLabel, ")",
      "\n Adult Event DDs and Label: ", adultEventDD, " (", adultEventLabel, ")",
      sep = "", file = metadata, append = TRUE)
}

cat("\n\n### Model Input Parameters ###\n Start Year:", start_year, 
    "\n Weather data for forecasts:", forecast_data, 
    "\n Start day-of-year:", start_doy,
    "\n End day-of-year:", end_doy, 
    "\n Region:", region_param, 
    "\n Climate stress exclusion maps:", exclusions_stressunits, 
    "\n Pest Event Maps:", pems,    
    "\n Adult Event Maps:", mapA, 
    "\n Egg Event Maps:", mapE, 
    "\n Larvae Event Maps:", mapL, 
    "\n Pupae Event Maps:", mapP,
    "\n Output_Dir:", out_dir, 
    "\n Output option:", out_option,
    file = metadata, append = TRUE)

LogWrap("", "Done writing metadata file", "\n")
LogWrap("\n", paste0(forecast_data, " DATA PROCESSING"), "\n")
cat("\nDone writing metadata file\n\n", forecast_data, " DATA PROCESSING\n", sep = "")

#### (6). WEATHER DATA LOADING AND PROCESSING -----

# Weather inputs and outputs - PRISM climate data w/subdirs 4-digit year
# New feature - choose whether to use PRISM or NMME for weather forecasts 
# (forecast_data = PRISM, or forecast_data = NMME)

# Loop through each needed variable and create a list of needed files
vars <- c("tmin", "tmax", "tmean", "tdmean", "ppt")
fls_list <- c("tminfiles", "tmaxfiles", "tmeanfiles", "tdmeanfiles", "pptfiles")

for (i in seq_along(vars)) {
  # Create list of files
  fls <- list.files(
    path = prism_dir, 
    pattern = glob2rx(paste0("*PRISM_", vars[i], "*", start_year, "*.bil$*")), 
    all.files = FALSE, full.names = TRUE, recursive = TRUE
    )
  
  # Exit program if files are missing
  if (length(fls) == 0) {
    LogWrap("", paste0("Could not find ", vars[i], " files - exiting program"), "\n")
    cat("Could not find ", vars[i],  " files - exiting program\n") 
    q()
    
    assign(fls_list[i], fls)
  }
  
  # Extract highest quality files and assign object name to list
  fls_best <- ExtractBestPRISM(fls, forecast_data, keep_leap)[start_doy:end_doy]
  assign(fls_list[i], fls_best)
  
}

## Extract date from temperature files using regex pattern matching
dats <- unique(regmatches(tminfiles, regexpr(pattern = "[0-9]{8}", 
                                             text = tminfiles)))

# Specify sampling frequency (how many days until output maps are generated?)
# This feature may be removed in production version
if (out_option == 1) {
  sample_freq <- 30 # Monthly maps
} else if (out_option == 2) {
  sample_freq <- 14  # Biweekly maps
} else if (out_option == 3) {
  sample_freq <- 10 # Dekad maps
} else if (out_option == 4) {
  sample_freq <- 7  # Weekly maps
} else if (out_option == 5) {
  sample_freq <- 2  # Even day maps
} else if (out_option == 6) {
  sample_freq <- 1 # Daily maps
}

# * Dates and sampling ####
# Make vector of dates to use when processing results 
# The current date will be sampled if it's the current year AND if 
# the current day falls within the range of start_doy and end_doy.
# The last date of year will always be sampled.
# Using "unique" will only keep date if it doesn't already occur in vector
# This happens if the end day of year is a multiple of the sampling frequency 
# (e.g. 1 to 300, w/ a 30 day sampling frequency), or if the current date falls
# within the sampling frequency
today <- strftime(Sys.time(), format = "%Y%m%d")
onedayago <- as.character(as.numeric(today) - 1) 
twodayago <- as.character(as.numeric(today) - 2) # change when NDFD used
current_year <- strftime(Sys.time(), format = "%Y")

if (start_year == current_year & 
    yday(Sys.time()) >= start_doy &
    yday(Sys.time()) <= end_doy) {
  dats2 <- sort(as.numeric(
    unique(c(dats[seq(0, length(dats), sample_freq)], today, 
             onedayago, twodayago, last(dats))))
    )
  LogWrap("\n", paste("Sampling every", sample_freq, "days between", 
                      first(dats), "and", last(dats)), "\n")
  LogWrap("", paste0("Sampling also includes today (", today, ") one day ago (",
                    onedayago, ") and two days ago (", twodayago, ")"), "\n")
  
} else {
  dats2 <- sort(as.numeric(unique(c(dats[seq(0, length(dats), sample_freq)],
                                    last(dats)))))
  LogWrap("\n", paste("Sampling every", sample_freq, "days between", 
                      first(dats), "and", last(dats)), "\n")
}

dats2 <- as.character(dats2) # Need to be in character format for plotting
num_dats <- length(dats2) # How many sampled dates? 

# Create vector of days in the sublist that will be sampled (rasters are saved 
# for those days) in the Daily Loop, and also tack on the last day in the list. 
sample_pts <- c(sublist[seq(0, length(sublist), sample_freq)],
                last(sublist))

# Add the present day if DDRP run is being run for the current year AND if 
# the current day falls within the range of start_doy and end_doy. 
if (start_year == current_year & 
    yday(Sys.time()) >= start_doy &
    yday(Sys.time()) <= end_doy) {
  today_doy <- strftime(Sys.time(), format = "%j") # Day of year
  sample_pts <- sort(as.numeric(unique(c(sample_pts, today_doy))))
}

# Keep only unique sampling points (there may be duplicates for example
# if the present day is already in the list).
sample_pts <- unique(sample_pts)

# Log file and terminal messages
LogWrap("", paste0("Finished loading ", forecast_data, " files for ", 
    length(start_doy:end_doy), " days"), "\n")
LogWrap("", paste0("Creating template file for ", region_param), "\n")
cat("\nFinished loading ", forecast_data, " files for ", 
    length(start_doy:end_doy), " days\n\nCreating template file for ", 
    region_param, "\n", sep = "")

# * Make template ####
# This template is used for cropping the temperature (tmin, tmax) rasters
# First define with extent of the region
#### Set up regions - use switch() (works like a single use hash) ##
REGION <- switch(region_param,
                 "CONUS"        = extent(-125.0,-66.5,24.0,50.0),
                 "WEST"         = extent(-125.0, -102, 31.1892, 49.4),
                 "EAST"         = extent(-106.8, -66.5, 24.54, 49.4),
                 "MIDWEST"      = extent(-104.2,-87,30,49),
                 "NORTHWEST"    = extent(-125.1,-103.8,40.6,49.2),
                 "SOUTHWEST"    = extent(-124.6,-101.5,31.2,42.2),
                 "SOUTHCENTRAL" = extent(-83.6,-78.3,31.8,35.3),
                 "NORTHCENTRAL" = extent(-104.3,-80.2,35.7,49.7),
                 "SOUTHEAST"    = extent(-107.1,-75.0,24.1,39.6),
                 "NORTHEAST"    = extent(-84.2,-64.3,36.9,48.1),
                 "AL"           = extent(-88.5294,-84.7506,30.1186,35.1911),
                 "AR"           = extent(-94.8878,-89.5094,32.8189,36.6936),
                 "AZ"           = extent(-115, -108.98, 31.2, 37),
                 "CA"           = extent(-124.6211, -113.7428, 32.2978, 42.2931),
                 "CO"           = extent(-109.2625, -101.8625, 36.7461, 41.2214),
                 "CT"           = extent(-73.7700, -71.7870, 40.9529, 42.0355),
                 "DL"           = extent(-76.1392, -74.1761, 38.3508, 39.9919),
                 "FL"           = extent(-87.8064, -79.9003, 24.4494, 31.1214),
                 "GA"           = extent(-85.7850, -80.5917, 30.1767, 35.1594),
                 "IA"           = extent(-96.8617, -89.9697, 40.1147, 43.7353),
                 "ID"           = extent(-117.3917, -110.6167, 41.4500, 49.3583),
                 "IL"           = extent(-91.5897, -87.0461, 36.8903, 42.6375),
                 "IN"           = extent(-88.1686, -84.4686, 37.7836, 41.9794),
                 "KS"           = extent(-102.3342, -94.1756, 36.6369, 40.2836),
                 "KT"           = extent(-89.3581, -81.8425, 36.4208, 39.3347),
                 "LA"           = extent(-94.3019, -88.7758, 28.8333, 33.2994),
                 "MA"           = extent(-73.5639, -69.7961, 41.1689, 42.9525),
                 "MD"           = extent(-79.7014, -74.8833, 37.0631, 39.9075),
                 "ME"           = extent(-71.4056, -66.6667, 42.9525, 47.5228),
                 "MI"           = extent(-90.5542, -82.3047, 41.6311, 47.5739),
                 "MN"           = extent(-97.4000, -89.3786, 43.2550, 49.3506),
                 "MO"           = extent(-95.8803, -88.9883, 35.8822, 40.7058),
                 "MS"           = extent(-91.7475, -87.8522, 29.9842, 35.2631),
                 "MT"           = extent(-116.3667, -103.8250, 44.0667, 49.3917),
                 "NC"           = extent(-84.440918, -75.300293, 33.682906, 36.646092),
                 "ND"           = extent(-104.2708, -96.3075, 45.6403, 49.1817),
                 "NE"           = extent(-104.3553, -95.0464, 39.7506, 43.2022),
                 "NH"           = extent(-72.6617, -70.6142, 42.6256, 45.4700),
                 "NJ"           = extent(-75.9175, -73.1892, 38.8944, 41.5806),
                 "NM"           = extent(-109.2942, -102.6383, 31.1892, 37.2000),
                 "NV"           = extent(-120.3358, -113.6803, 34.7356, 42.2981),
                 "NY"           = extent(-80.0867, -71.7381, 40.4828, 45.1692),
                 "OH"           = extent(-85.0439, -80.2464, 38.2797, 42.0217),
                 "OK"           = extent(-103.2850, -94.1964, 33.3839, 37.2850),
                 "OR"           = extent(-124.7294, -116.2949, 41.7150, 46.4612),
                 "PA"           = extent(-80.7672, -74.5033, 39.4694, 42.5094),
                 "RI"           = extent(-71.8628, -71.1206, 41.1463, 42.0188),
                 "SC"           = extent(-83.6422, -78.3275, 31.8814, 35.3811),
                 "SD"           = extent(-104.3553, -96.0806, 42.3050, 46.2050),
                 "TN"           = extent(-90.3239, -81.5047, 34.5578, 37.1125),
                 "TX"           = extent(-107.1592, -93.2411, 25.8614, 36.7200),
                 "UT"           = extent(-114.2925, -108.7450, 36.7778, 42.2347),
                 "VA"           = extent(-83.8322, -75.6200, 36.3892, 39.7886),
                 "VT"           = extent(-73.6747, -71.4108, 42.5886, 45.1956),
                 "WA"           = extent(-124.9585, -116.8364, 45.4554, 49.1664),
                 "WI"           = extent(-93.1572, -86.6822, 42.2733, 46.9914),
                 "WV"           = extent(-82.8783, -77.5114, 37.1158, 40.7836),
                 "WY"           = extent(-111.6167, -103.7333, 40.6667, 45.4833))

# Now make a template for the region
template <- crop(raster(tminfiles[1]), REGION) # Template for cropping
template[!is.na(template)] <- 0
dataType(template) <- "INT2U"

# * Crop climate rasters ####
# If CONUS or EAST, split template into tiles (and run in parallel)
# Benefit of tiles is lost for smaller regions, so these are not split
# This part differs from DDRP_v2.R because the cropped rasters are not held
# in a list because this uses way too much memory with the addition of 3
# other climate variables (tmean, tdmean, ppt)
ncores <- detectCores()
RegCluster(round(ncores/4))

if (region_param %in% c("CONUS", "EAST")) {
  
  # Split template (2 pieces per side), resulting in 4 tiles
  tile_list <- SplitRas(template, ppside = 2, save = FALSE, plot = FALSE) 
  LogWrap("", paste("Splitting template into", length(tile_list), "tiles"), "\n")
  
  # Name the 4 tiles (tile1, tile2, tile3, tile4)
  template <- mapply(function(n, t) {
    names(n) <- paste0("tile", t)
    return(n)
  }, n = tile_list, t = 1:length(tile_list))
  
  rm(tile_list)
  
  # Crop climate files by each template tile
  LogWrap("", paste("Cropping tmax, tmin, tmean, tdmean and ppt tiles for", 
                    region_param), "\n")
  cat("\nCropping tmax, tmin, tmean, tdmean and ppt for", region_param, "\n")

  tmin_list <- foreach(tile = template, .packages = "raster") %:% 
    foreach(tmin = tminfiles, .packages = "raster", .inorder = TRUE) %dopar% { 
      m <- as.matrix(crop(raster(tmin), tile))
    }
    
  tmax_list <- foreach(tile = template, .packages = "raster") %:% 
    foreach(tmax = tmaxfiles, .packages = "raster", .inorder = TRUE) %dopar% { 
      m <- as.matrix(crop(raster(tmax), tile))
    }
  
  tmean_list <- foreach(tile = template, .packages = "raster") %:% 
    foreach(tmean = tmeanfiles, .packages = "raster", .inorder = TRUE) %dopar% { 
      m <- as.matrix(crop(raster(tmean), tile))
    }
  
  tdmean_list <- foreach(tile = template, .packages = "raster") %:% 
    foreach(tdmean = tdmeanfiles, .packages = "raster", .inorder = TRUE) %dopar% { 
      m <- as.matrix(crop(raster(tdmean), tile))
    }
  
  ppt_list <- foreach(tile = template, .packages = "raster") %:% 
    foreach(ppt = pptfiles, .packages = "raster", .inorder = TRUE) %dopar% { 
      m <- as.matrix(crop(raster(ppt), tile))
    }

# If region is not CONUS or EAST, simply crop temp files by the single template
} else {
  LogWrap("", paste("Cropping tmax, tmin, tmean, tdmean and ppt tiles for", 
                     region_param), "\n")
  cat("\nCropping tmax, tmin, tmean, tdmean and ppt tiles for", region_param, "\n")
  
  tmin_list <- foreach(tmin = tminfiles, .packages = "raster", 
                       .inorder = TRUE) %dopar% {
    m <- as.matrix(crop(raster(tmin), template))
  }

  tmax_list <- foreach(tmax = tmaxfiles, .packages = "raster", 
                       .inorder = TRUE) %dopar% {
    m <- as.matrix(crop(raster(tmax), template))
  }

  tmean_list <- foreach(tmean = tmeanfiles, .packages = "raster", 
                       .inorder = TRUE) %dopar% {
    m <- as.matrix(crop(raster(tmean), template))
                       }
  tdmean_list <- foreach(tdmean = tdmeanfiles, .packages = "raster", 
                       .inorder = TRUE) %dopar% {
    m <- as.matrix(crop(raster(tdmean), template))
  }

  ppt_list <- foreach(ppt = pptfiles, .packages = "raster", 
                       .inorder = TRUE) %dopar% {
    m <- as.matrix(crop(raster(ppt), template))
  }
}

## (7). DAILY LOOP FUNCTION -----

#### * DailyLoop function: daily time step ####
# The daily loop is the daily time step used to generate all results.
# If the study region is CONUS or EAST, then each of 4 tiles (= tile_num) are 
# run in parallel as well. The template (= template) is needed for setting up 
# data sets to populate during the model run, and to convert matrix outputs 
# into raster format (matrices are used as inputs because it runs much faster)
# Here it is written as a function because it will be applied in slightly different
# ways depending on whether "templates" is a list (CONUS/EAST) or a single template.
DailyLoop <- function(tile_num, template) {
  
  #### * Initialize rasters for all variables ####
  # Main rasters - these are built upon within daily loop
  DDaccum <- as.matrix(template)
  # Track total degree days accumulated, for a single Lifestage (larvae)
  DDtotal <- as.matrix(template)
  # Track Lifestage for each cell per day (all cells start in OW = stage1)
  Lifestage <- as.matrix(template) + 1
  # Track voltinism per cell per day, starting at 0
  NumGen <- as.matrix(template)
  
  # BOXB templates
  relhum <- as.matrix(template)
  LW <- as.matrix(template)
  DDs <- as.matrix(template)
  DDLW <- as.matrix(template)
  CUMDDs <- as.matrix(template)
  CUMDDLW <- as.matrix(template)
  
  # Additional rasters - created depending on input setting
  if (exclusions_stressunits) {
    # Create masks for each variable
    coldmask         <- as.matrix(template)  # Daily cold units mask
    coldstress       <- as.matrix(template)  # Count of daily cold units
    coldstressTHRESH  <- as.matrix(template)  # Mask for coldstrs units thres
    coldstressTHRESH  <- coldstress_threshold # Mask for coldstrs units thres
    coldunitsCUM     <- as.matrix(template)  # Cumulative cold units
    coldstressMAX1    <- as.matrix(template)  # Max cold before most die
    coldstressMAX1    <- coldstress_units_max1 # Max cold before most die
    coldstressMAX2    <- as.matrix(template)  # Max cold before all die
    coldstressMAX2    <- coldstress_units_max2 # Max cold before all die
    coldEXCL         <- as.matrix(template)  # Cold stress exclusion
    heatmask          <- as.matrix(template)  # Daily heat stress units mask
    heatstress        <- as.matrix(template)  # Count of daily heat stress units
    heatstressTHRESH  <- as.matrix(template)  # Mask for heatstress units thres
    heatstressTHRESH  <- heatstress_threshold # Mask for heatstress units thres
    heatunitsCUM      <- as.matrix(template)  # Cumulative heat stress units
    heatstressMAX1    <- as.matrix(template)  # Max heat before most die
    heatstressMAX1    <- heatstress_units_max1 # Max heat before most die
    heatstressMAX2    <- as.matrix(template)  # Max heat before all die
    heatstressMAX2    <- heatstress_units_max2 # Max heat before all die
    heatEXCL          <- as.matrix(template)  # Heat stress exclusions
    # BOXB moisture stress templates
    drymask         <- as.matrix(template)  # binary mask for daily dry units
    drystress       <- as.matrix(template)  # count of daily dry   units
    drystressTHRESH  <- as.matrix(template)  # mask for drystress units threshold
    drystressTHRESH  <- drystress_threshold # mask for drystress units threshold
    dryunitsCUM     <- as.matrix(template)  # cumulative dry units
    drystressMAX1    <- as.matrix(template)  # use for max dry before mortality??
    drystressMAX1    <- drystress_units_max1 # use for max dry before mortality??
    drystressMAX2    <- as.matrix(template)  # use for max dry before mortality??
    drystressMAX2    <- drystress_units_max2 # use for uncertainty zone max dry before mortality??
    dryEXCL         <- as.matrix(template)  # EXCL map for dry stress
    # All climate stress exclusions
    AllEXCL           <- as.matrix(template)  # Combined stress exclusions
  }
  
  if (pems) {
    
    if (mapE == 1 & eggEventDD) {  
      # Event DD for all must be specified in spp.params file for this to run
      if (PEMnumgens > 0) {
        # Egg DOYs for when cumDDs > eggEvent threshold OW generation
        PEMe0 <- as.matrix(template)  
        # Egg DOYs for when cumDDs > eggEvent threshold  1st Gen
        PEMe1 <- as.matrix(template)  
      }
      if (PEMnumgens > 1) {
        # Egg DOYs for when cumDDs > eggEvent threshold  2nd Gen
        PEMe2 <- as.matrix(template)  
      }
      if (PEMnumgens > 2) {
        # Egg DOYs for when cumDDs > eggEvent threshold  3rd Gen
        PEMe3 <- as.matrix(template)  
      }
      if (PEMnumgens > 3) {
        # Egg DOYs for when cumDDs > eggEvent threshold  4th Gen
        PEMe4 <- as.matrix(template)  
      }
    }
    
    if (mapL == 1 & larvaeEventDD) {
      if (PEMnumgens > 0) {
        # Larval DOYs for when cumDDs > larvaeEvent threshold OW generation
        PEMl0 <- as.matrix(template)  
        # Larval DOYs for when cumDDs > larvaeEvent threshold  1st Gen
        PEMl1 <- as.matrix(template)  
      }
      if (PEMnumgens > 1) {
        # Larval DOYs for when cumDDs > larvaeEvent threshold  2nd Gen
        PEMl2 <- as.matrix(template)  
      }
      if (PEMnumgens > 2) {
        # Larval DOYs for when cumDDs > larvaeEvent threshold  3rd Gen
        PEMl3 <- as.matrix(template)  
      }
      if (PEMnumgens > 3) {
        # Larval DOYs for when cumDDs > larvaeEvent threshold  4th Gen
        PEMl4 <- as.matrix(template)  
      }
    }
    
    if (mapP == 1 & pupaeEventDD) {
      if (PEMnumgens > 0) {
        # Pupal DOYs for when cumDDs > pupalEvent threshold OW generation
        PEMp0 <- as.matrix(template)  
        # Pupal DOYs for when cumDDs > pupalEvent threshold  1st Gen
        PEMp1 <- as.matrix(template)  
      }
      if (PEMnumgens > 1) {
        # Pupal DOYs for when cumDDs > pupalEvent threshold  2nd Gen
        PEMp2 <- as.matrix(template)  
      }
      if (PEMnumgens > 2) {
        # Pupal DOYs for when cumDDs > pupalEvent threshold  3rd Gen
        PEMp3 <- as.matrix(template)  
      }
      if (PEMnumgens > 3) {
        # Pupal DOYs for when cumDDs > pupalEvent threshold  4th Gen
        PEMp4 <- as.matrix(template)  
      }
    }
    
    if (mapA == 1 & adultEventDD) {
      if (PEMnumgens > 0) {
        # Adult DOYs for when cumDDs > adultEvent threshold OW generation
        PEMa0 <- as.matrix(template)  
        # Adult DOYs for when cumDDs > adultEvent threshold  1st Gen
        PEMa1 <- as.matrix(template)  
      }
      if (PEMnumgens > 1) {
        # Adult DOYs for when cumDDs > adultEvent threshold  2nd Gen
        PEMa2 <- as.matrix(template)  
      }
      if (PEMnumgens > 2) {
        # Adult DOYs for when cumDDs > adultEvent threshold  3rd Gen
        PEMa3 <- as.matrix(template)  
      }
      if (PEMnumgens > 3) {
        # Adult DOYs for when cumDDs > adultEvent threshold  4th Gen
        PEMa4 <- as.matrix(template)  
      }
    }
    
  }
  
  #### * Step through days ####
  # tryCatch(
  for (d in 1:length(sublist)) {
    #tic()
    if (region_param %in% c("CONUS", "EAST")) {
      tmax <- as.numeric(tmax_list[[tile_num]][[d]])
      tmin <- as.numeric(tmin_list[[tile_num]][[d]])
      tmean <- as.numeric(tmean_list[[tile_num]][[d]])
      tdmean <- as.numeric(tdmean_list[[tile_num]][[d]])
      ppt <- as.numeric(ppt_list[[tile_num]][[d]])
    } else {
      tmax <- as.numeric(tmax_list[[d]])
      tmin <- as.numeric(tmin_list[[d]])
      tmean <- as.numeric(tmean_list[[d]])
      tdmean <- as.numeric(tdmean_list[[d]])
      ppt <- as.numeric(ppt_list[[d]])
    }

    # Assign each cell to a Lifestage
    # These three matrices assign physiological parameters by cell
    ls_ldt <- stage_ldt[Lifestage]
    ls_udt <- stage_udt[Lifestage]
    ls_dd <- stage_dd[Lifestage]
    
    # Accumulate total DDs across year, using larvae Lifestage
    # Can not use DDaccum because this one has values reset when progress = 1
    ls_ldt_larv <- stage_ldt[which(stgorder == "L")]
    ls_udt_larv <- stage_udt[which(stgorder == "L")]
    
    # Calculate stage-specific degree-days for each cell per day
    if (calctype == "average") {
      dd_tmp <- AvgDD(tmax, tmin, ls_ldt, ls_udt)
      dd_tmp_larv <- AvgDD(tmax, tmin, ls_ldt_larv, ls_udt_larv)
    } else if (calctype == "triangle") {
      dd_tmp <- TriDD(tmax, tmin, ls_ldt, ls_udt)
      dd_tmp_larv <- TriDD(tmax, tmin, ls_ldt_larv, ls_udt_larv)
    # For BOXB
    } else if (calctype == "boxblt") {
      # This needs to be replaced. Currently Lifestage, NumGen, and StageCount
      # (combo of Lifestage and NumGen) are calculated using the "triangle"
      # DD calculation, which is relevant to daily DDs, not hourly
      dd_tmp <- TriDD(tmax, tmin, ls_ldt_larv, ls_udt_larv) # REPLACE ME
      relhum <- calcrelhum(tmean, tdmean)
      LW <- calcLW(ppt, relhum)
      DDs <- calcDDs(tmax, tmin)
      DDLW <- LW * dd_tmp/24
      CUMDDs <- CUMDDs + DDs
      CUMDDLW <- CUMDDLW + DDLW
      # Same as dd_tmp - this should be changed for plant diseases
      dd_tmp_larv <- TriDD(tmax, tmin, ls_ldt_larv, ls_udt_larv) # REPLACE ME
    }
      
    # Accumulate degree days
    DDaccum <- DDaccum + dd_tmp
    
    # Calculate total degree days for larvae
    DDtotal <- DDtotal + dd_tmp_larv 
    
    # Climate stress exclusions
    # ASSUMES: -2=severe -1=mod 0=none throughout
    if (exclusions_stressunits) {
      # Cold stress accumulation
      # Make today's cold mask and calculate today's cold stress DDs
      coldmask <- tmin < coldstressTHRESH  
      coldstress <- coldmask * abs(coldstressTHRESH - tmin) 
      coldunitsCUM <- coldunitsCUM + coldstress
      # ASSUME NEW -2=severe -1=mod 0=none throughout
      coldEXCL <- Cond(coldunitsCUM >= coldstressMAX2, -2, 
                       Cond(coldunitsCUM >= coldstressMAX1, -1, 0))
      # Heat stress accumulation
      # Make today's heat mask and calculate today's heat stress DDs
      heatmask <- tmax > heatstressTHRESH  
      heatstress <- heatmask * abs(tmax - heatstressTHRESH) 
      heatunitsCUM <- heatunitsCUM + heatstress
      heatEXCL <- Cond(heatunitsCUM >= heatstressMAX2, -2, 
                       Cond(heatunitsCUM >= heatstressMAX1, -1, 0))

      # For BOXB
      # BOXB Changes start: need to fix equations once we have relhum
      ##-- DDLW CALCULATION??
      ##-- Dry Stress Accumulation
      drymask <- relhum < drystressTHRESH  # make todays dry mask
      drystress <- drymask * ((abs(drystressTHRESH - relhum) - 500 * ppt)) # compute todays dry stress DDs
      #drystress <- drystress - 500 * ppt   # ppt is hundreds of an inch, it reduces drystess; TODO: calibrate
      drystress <- Cond(drystress >= 0, drystress, 0)
      dryunitsCUM <- dryunitsCUM + drystress
      dryEXCL <- Cond(dryunitsCUM >= drystressMAX2,-2,
                      Cond(dryunitsCUM >= drystressMAX1,-1,0))
      
      # All climate stress exclusion for BOXB 
      AllEXCL <- Cond((dryEXCL == 0) & (coldEXCL == 0) & (heatEXCL == 0),0,
              Cond((dryEXCL == -2) | (coldEXCL == -2) | (heatEXCL == -2),-2,-1))
      
    }
    
    # Calculate pest events
    # DOYs for when cumDDs > event threshold for a given generation
    if (pems) {
      
      # FIX OWEventDD
      # The OW pest event is a percentage of the OW stage DDs 
      # (e.g. OWEventP = 0.5 would be 50% of the OWlarvaeDD, or half-way 
      # through OWlarvaeDD). It is specified in the species param file.
      OWEventDD <- stage_dd[1] * OWEventP # What should this be for a plant pathogen?
      
      # Egg PEMs
      if (mapE == 1 & eggEventDD) {  
        if (PEMnumgens > 0) {
          if (owstage == "OE") {
            # If owstage = egg, then eggs finish developing after the DD of OWegg 
            PEMe0 <- Cond(PEMe0 == 0 & NumGen == 0 & (DDaccum >= OWEventDD), 
                          d * (Lifestage == which(stgorder == "OE")), PEMe0)
          }
          # Egg DOYs for when cumDDs > eggEvent threshold 1st gen
          PEMe1 <- Cond(PEMe1 == 0 & NumGen == 1 & (DDaccum >= eggEventDD), 
                        d * (Lifestage == which(stgorder == "E")), PEMe1) 
        }
        if (PEMnumgens > 1) {
          # Egg DOYs for when cumDDs > eggEvent threshold 2nd gen
          PEMe2 <- Cond(PEMe2 == 0 & NumGen == 2 & (DDaccum >= eggEventDD), 
                        d * (Lifestage == which(stgorder == "E")), PEMe2)
        }
        if (PEMnumgens > 2) {
          # Egg DOYs for when cumDDs > eggEvent threshold 3rd gen
          PEMe3 <- Cond(PEMe3 == 0 & NumGen == 3 & (DDaccum >= eggEventDD), 
                        d * (Lifestage == which(stgorder == "E")), PEMe3) 
        }
        if (PEMnumgens > 3) {
          # Egg DOYs for when cumDDs > eggEvent threshold 4th gen
          PEMe4 <- Cond(PEMe4 == 0 & NumGen == 4 & (DDaccum >= eggEventDD), 
                        d * (Lifestage == which(stgorder == "E")), PEMe4) 
        }
      }
      
      # Larvae PEMs
      if (mapL == 1 & larvaeEventDD) {
        if (PEMnumgens > 0) {
          if (owstage == "OL") {
            # If owstage = larvae, then larvae finish developing after the DD of 
            # OWlarvae 
            PEMl0 <- Cond(PEMl0 == 0 & NumGen == 0 & (DDaccum >= OWEventDD), 
                          d * (Lifestage == which(stgorder == "OL")), PEMl0) 
            # If owstage = egg, then larvae of this OW gen will have to go through 
            # full development 
          } else if (owstage %in% c("OE")) {
            PEMl0 <- Cond(PEMl0 == 0 & NumGen == 0 & (DDaccum >= larvaeEventDD), 
                          d * (Lifestage == which(stgorder == "L")), PEMl0) 
          }
          # Larvae DOYs for when cumDDs > larvaeEvent threshold 1st gen
          PEMl1 <- Cond(PEMl1 == 0 & NumGen == 1 & (DDaccum >= larvaeEventDD), 
                        d * (Lifestage == which(stgorder == "L")), PEMl1) 
        }
        if (PEMnumgens > 1) {
          # Larvae DOYs for when cumDDs > larvaeEvent threshold 2nd gen
          PEMl2 <- Cond(PEMl2 == 0 & NumGen == 2 & (DDaccum >= larvaeEventDD), 
                        d * (Lifestage == which(stgorder == "L")), PEMl2) 
        }
        if (PEMnumgens > 2) {
          # Larvae DOYs for when cumDDs > larvaeEvent threshold 3rd gen
          PEMl3 <- Cond(PEMl3 == 0 & NumGen == 3 & (DDaccum >= larvaeEventDD), 
                        d * (Lifestage == which(stgorder == "L")), PEMl3) 
        }
        if (PEMnumgens > 3) {
          # Larvae DOYs for when cumDDs > larvaeEvent threshold 4th gen
          PEMl4 <- Cond(PEMl4 == 0 & NumGen == 4 & (DDaccum >= larvaeEventDD), 
                        d * (Lifestage == which(stgorder == "L")), PEMl4) 
        }
      }
      
      # Pupae PEMs
      if (mapP == 1 & pupaeEventDD) {
        if (PEMnumgens > 0) {
          if (owstage == "OP") {
            # If owstage = pupae, then pupae finish developing after the DD of 
            # OWpupae 
            PEMp0 <- Cond(PEMp0 == 0 & NumGen == 0 & (DDaccum >= OWEventDD), 
                          d * (Lifestage == which(stgorder == "OP")), PEMp0)
          } else if (owstage %in% c("OE", "OL")) {
            # If owstage = egg or larvae, then pupae of the OW gen will have to 
            # go through full development 
            PEMp0 <- Cond(PEMp0 == 0 & NumGen == 0 & (DDaccum >= pupaeEventDD), 
                          d * (Lifestage == which(stgorder == "P")), PEMp0) 
          }  
          # Pupae DOYs for when cumDDs > pupaeEvent threshold 1st gen
          PEMp1 <- Cond(PEMp1 == 0 & NumGen == 1 & (DDaccum >= pupaeEventDD), 
                        d * (Lifestage == which(stgorder == "P")), PEMp1) 
        }
        if (PEMnumgens > 1) {
          # Pupae DOYs for when cumDDs > pupaeEvent threshold 2nd gen
          PEMp2 <- Cond(PEMp2 == 0 & NumGen == 2 & (DDaccum >= pupaeEventDD), 
                        d * (Lifestage == which(stgorder == "P")), PEMp2) 
        }
        if (PEMnumgens > 2) {
          # Pupae DOYs for when cumDDs > pupaeEvent threshold 3rd gen
          PEMp3 <- Cond(PEMp3 == 0 & NumGen == 3 & (DDaccum >= pupaeEventDD), 
                        d * (Lifestage == which(stgorder == "P")), PEMp3) 
        }
        if (PEMnumgens > 3) {
          # Pupae DOYs for when cumDDs > pupaeEvent threshold 4th gen
          PEMp4 <- Cond(PEMp4 == 0 & NumGen == 4 & (DDaccum >= pupaeEventDD), 
                        d * (Lifestage == which(stgorder == "P")), PEMp4) 
        }
      }
      
      # Adult PEMs
      if (mapA == 1 & adultEventDD) {
        if (PEMnumgens > 0) {
          if (owstage == "OA") {
            # If owstage = adult, then adults finish developing after the DD of
            # OWadult
            PEMa0 <- Cond(PEMa0 == 0 & NumGen == 0 & (DDaccum >= OWEventDD), 
                          d * (Lifestage == which(stgorder == "OA")), PEMa0) 
          } else if (owstage %in% c("OL", "OP")) {
            # If owstage = larvae or pupae, then adults of the OW gen will 
            # have to go through full development 
            PEMa0 <- Cond(PEMa0 == 0 & NumGen == 0 & (DDaccum >= adultEventDD), 
                          d * (Lifestage == which(stgorder == "A")), PEMa0) 
          }
          # Adult DOYs for when cumDDs > adultEvent threshold 1st gen
          PEMa1 <- Cond(PEMa1 == 0 & NumGen == 1 & (DDaccum >= adultEventDD), 
                        d * (Lifestage == which(stgorder == "A")), PEMa1) 
        } 
        if (PEMnumgens > 1) {
          # Adult DOYs for when cumDDs > adultEvent threshold 2nd gen
          PEMa2 <- Cond(PEMa2 == 0 & NumGen == 2 & (DDaccum >= adultEventDD), 
                        d * (Lifestage == which(stgorder == "A")), PEMa2) 
        }
        if (PEMnumgens > 2) {
          # Adult DOYs for when cumDDs > adultEvent threshold 3rd gen
          PEMa3 <- Cond(PEMa3 == 0 & NumGen == 3 & (DDaccum >= adultEventDD), 
                        d * (Lifestage == which(stgorder == "A")), PEMa3) 
        }
        if (PEMnumgens > 3) {
          # Adult DOYs for when cumDDs > adultEvent threshold 4th gen
          PEMa4 <- Cond(PEMa4 == 0 & NumGen == 4 & (DDaccum >= adultEventDD), 
                        d * (Lifestage == which(stgorder == "A")), PEMa4) 
        }
      }
      
    }
    
    # Calculate Lifestage progression: Is accumulation > Lifestage requirement 
    # (0 = FALSE, 1 = TRUE)? 
    progress <- as.integer(DDaccum >= ls_dd)
    
    # If reproductive adult stage progresses, then that cell has oviposition 
    # and the generation count increases. If species has OW adults, then need 
    # to change stage value to "adult" to allow NumGen to increase when it 
    # progresses. "OA" = 1, and "adult" = 5 for species with OW adults
    if (owstage == "OA") {
      Lifestage2 <- Lifestage
      Lifestage2[Lifestage2 == 1] <- 5
      NumGen <- NumGen + (progress == 1 & Lifestage2 == 5)
    } else {
      # Value for "adult" varies depending on OW stage
      NumGen <- NumGen + (progress == 1 & Lifestage == which(stgorder == "A")) 
    }
    
    # If progress is 1, then there is progression to the next life stage
    Lifestage <- Lifestage + progress
    
    # Reset the DDaccum cells to zero for cells that progressed to next stage
    DDaccum <- DDaccum - (progress * ls_dd)
    
    # Reassign cells that progressed past end of stgorder to first non-OW stage
    Lifestage <- Cond(Lifestage == (length(stgorder) + 1), 2, Lifestage)
    
    # Stage Count analysis and raster brick. The generation numbers is divided
    # by 10 so that the gen # is a decimal (e.g. Lifestage 4 of 2nd gen = 4.2)
    StageCount <- NumGen + Lifestage/10
  
    #### * Save data for certain days, specified by sampling frequency ####
    # Data from last sampling day of year is also saved
    if (sublist[d] %in% sample_pts) {
      # Convert Lifestage and Numgen matrices to rasters and put into a brick
      # Added rasters for BOXB (CUMDDs, CUMDDLW)
      mat_list <- list(Lifestage, NumGen, StageCount, DDtotal, CUMDDs, CUMDDLW)
      ext <- as.data.frame(as.matrix(extent(template)))
      rast_list <- lapply(mat_list, Mat_to_rast, ext = ext, template = template)
      names(rast_list) <- c("Lifestage_rast", "NumGen_rast", "StageCount_rast",
                            "DDtotal_rast", "CUMDDs_rast", "CUMDDLW_rast")
      
      # Lifestage brick
      if (!exists("Lifestage_brick")) {
        Lifestage_brick <- brick(rast_list$Lifestage_rast, crs = crs)
      } else {
        Lifestage_brick <- addLayer(Lifestage_brick, rast_list$Lifestage_rast)
      }
      
      # NumGen brick
      if (!exists("NumGen_brick")) {
        NumGen_brick <- brick(rast_list$NumGen_rast, crs = crs)
      } else {
        NumGen_brick <- addLayer(NumGen_brick, rast_list$NumGen_rast)
      }
      
      # StageCount brick
      if (!exists("StageCount_brick")) {
        StageCount_brick <- brick(rast_list$StageCount_rast, crs = crs)
      } else {
        StageCount_brick <- addLayer(StageCount_brick, rast_list$StageCount_rast)
      }
      
      # DDtotal
      if (!exists("DDtotal_brick")) {
        DDtotal_brick <- brick(rast_list$DDtotal_rast, crs = crs)
      } else {
        DDtotal_brick <- addLayer(DDtotal_brick, rast_list$DDtotal_rast)
      }
      
      # For BOXB - Plant disease DDs and CUMDDLW (inc. estimates of leaf wetness)
      if (!exists("CUMDDs_brick")) {
        CUMDDs_brick <- brick(rast_list$CUMDDs_rast, crs = crs)
      } else {
        CUMDDs_brick <- addLayer(CUMDDs_brick, rast_list$CUMDDs_rast)
      }
      
      if (!exists("CUMDDLW_brick")) {
        CUMDDLW_brick <- brick(rast_list$CUMDDLW_rast, crs = crs)
      } else {
        CUMDDLW_brick <- addLayer(CUMDDLW_brick, rast_list$CUMDDLW_rast)
      }
      
      rm(rast_list) # Free up memory
      
      #  If exclusions_stressunits = 1, save those results in raster bricks
      if (exclusions_stressunits) {
        
        # Do the same for climate stress exclusions
        # Convert matrices to rasters and put them into a raster brick
        mat_list3 <- list(coldunitsCUM, coldEXCL, heatunitsCUM, heatEXCL, 
                          dryunitsCUM, dryEXCL, AllEXCL)
        ext <- as.data.frame(as.matrix(extent(template)))
        rast_list2 <- lapply(mat_list3, Mat_to_rast, ext = ext, 
                             template = template)
        names(rast_list2) <- c("coldunitsCUM_rast", "coldEXCL_rast", 
                               "heatunitsCUM_rast", "heatEXCL_rast", 
                               "dryunitsCUM_rast", "dryEXCL_rast", "AllEXCL_rast")
          
        # Cold stress unit accumulation brick
        if (!exists("coldunitsCUM_brick")) {
          coldunitsCUM_brick <- brick(rast_list2$coldunitsCUM_rast, 
                                      crs = crs)
        } else {
          coldunitsCUM_brick <- addLayer(coldunitsCUM_brick, 
                                         rast_list2$coldunitsCUM_rast)
        }

        # Cold stress exclusion brick
        if (!exists("coldEXCL_brick")) {
          coldEXCL_brick <- brick(rast_list2$coldEXCL_rast, crs = crs)
        } else {
          coldEXCL_brick <- addLayer(coldEXCL_brick, 
                                     rast_list2$coldEXCL_rast)
        }
        
        # Heat stress unit accumulation brick
        if (!exists("heatunitsCUM_brick")) {
          heatunitsCUM_brick <- brick(rast_list2$heatunitsCUM_rast, crs = crs)
        } else {
          heatunitsCUM_brick <- addLayer(heatunitsCUM_brick, 
                                         rast_list2$heatunitsCUM_rast)
        }
          
        # Heat stress exclusion brick
        if (!exists("heatEXCL_brick")) {
          heatEXCL_brick <- brick(rast_list2$heatEXCL_rast, crs = crs)
        } else {
          heatEXCL_brick <- addLayer(heatEXCL_brick, rast_list2$heatEXCL_rast)
        }
          
        # For BOXB
        # Dry strss unit accumulation brick
        if (!exists("dryunitsCUM_brick")) {
          dryunitsCUM_brick <- brick(rast_list2$dryunitsCUM_rast, crs = crs)
        } else {
          dryunitsCUM_brick <- addLayer(dryunitsCUM_brick, 
                                         rast_list2$dryunitsCUM_rast)
        }

        # Dry stress exclusion brick
        if (!exists("dryEXCL_brick")) {
          dryEXCL_brick <- brick(rast_list2$dryEXCL_rast, crs = crs)
        } else {
          dryEXCL_brick <- addLayer(dryEXCL_brick, rast_list2$dryEXCL_rast)
        }
          
        # All stress exclusion brick (cold stress + heat stress exclusions)
        # For BOXB - also includes dry stress exclusions
        if (!exists("AllEXCL_brick")) {
          AllEXCL_brick <- brick(rast_list2$AllEXCL_rast, crs = crs)
        } else {
          AllEXCL_brick <- addLayer(AllEXCL_brick, rast_list2$AllEXCL_rast)
        }
        
        rm(rast_list2) # Free up memory
      }
    }
    #toc()
  }
  
  ### * Daily loop done - save raster bricks ###
  # "tile_num" has a value if multiple template tiles are being run in parallel
  # Save non-optional raster bricks = DDtotal, Lifestage, NumGen, and StageCount
  SaveRaster(DDtotal_brick, "DailyLoop", tile_num, "DDtotal", "INT2S",
             paste("-", "Degree day accum. for all", num_dats, "dates"))
  SaveRaster(Lifestage_brick, "DailyLoop", tile_num, "Lifestage", "INT1U",
             paste("-", "Lifestage for all", num_dats, "dates"))
  SaveRaster(NumGen_brick, "DailyLoop", tile_num, "NumGen", "INT1U",
             paste("-", "No. of gens for all", num_dats, "dates"))
  SaveRaster(StageCount_brick, "DailyLoop", tile_num, "StageCount", "FLT4S",
             paste("-", "Generation and stage for all", num_dats, "dates"))
  
  # For BOXB
  SaveRaster(CUMDDs_brick, "DailyLoop", tile_num, "CumDDs", "INT2S",
             paste("-", "Plant disease degree day unit accum. for all", num_dats, "dates"))
  SaveRaster(CUMDDLW_brick, "DailyLoop", tile_num, "Cum_Inf_Risk", "INT2S",
             paste("-", "Plant disease cum. infection risk for all", num_dats, "dates"))
  
  # If exclusions_stressunits = 1, then save stress unit and exclusions bricks
  if (exclusions_stressunits) {
    
    # Save climate stress exclusion and unit bricks
    brk_list1 <- c(coldunitsCUM_brick, coldEXCL_brick, heatunitsCUM_brick,
                   heatEXCL_brick, dryunitsCUM_brick, dryEXCL_brick, AllEXCL_brick)
    outnams1 <- c("Cold_Stress_Units", "Cold_Stress_Excl", "Heat_Stress_Units", 
                  "Heat_Stress_Excl", "Dry_Stress_Units", "Dry_Stress_Excl", "All_Stress_Excl")
    for (i in 1:length(brk_list1)) {
      log_capt <- str_to_sentence(gsub("_", " ", outnams1[i]))
      SaveRaster(brk_list1[[i]], "DailyLoop", tile_num, outnams1[i], 
                  "INT2S", paste("-", log_capt, "for all", num_dats, "dates"))
    }

    # For certain raster bricks from above, generate a brick that includes 
    # climate stress exclusion values (-2 = moderate stress, -1 = severe stress only)
    brk_list2 <- rep(list(Lifestage_brick, NumGen_brick, CUMDDLW_brick, StageCount_brick), 2)
    type_list <- c(rep("Excl1", 4), rep("Excl2", 4))
    outnams2 <- c("Lifestage_Excl1", "NumGen_Excl1", "Cum_Inf_Risk_Excl1", "StageCount_Excl1",
      "Lifestage_Excl2", "NumGen_Excl2", "Cum_Inf_Risk_Excl2", "StageCount_Excl2")
    formats <- c(c(rep("INT2S", 3), "FLT4S"), c(rep("INT2S", 3), "FLT4S"))
    for (i in 1:length(brk_list2)) {
      brk_excl <- brick(Rast_Subs_Excl(brk_list2[[i]], tile_num, type_list[i]))
      log_capt <- str_to_sentence(gsub("_", " ", outnams2[i]))
      SaveRaster(brk_excl, "DailyLoop", tile_num, outnams2[i], formats[i], 
                 paste("-", log_capt, "for all", num_dats, "dates"))
    }
  }
  
  # If Pest Event Maps are turned on, save rasters including (optional) ones with
  # climate stress exclusions
  if (pems) {
    pem_list <- mget(ls(pattern = "PEMe|PEMl|PEMp|PEMa"))
    
    # Convert each matrix in the list to a raster and save it
    for (i in 1:length(pem_list)) {
      pem_mat <- pem_list[[i]]
      pem_rast <- Mat_to_rast(pem_mat, ext, template)
      names(pem_rast) <- names(pem_list[i])
      nam <- names(pem_rast)
      
      if (region_param %in% c("CONUS", "EAST")) {
        SaveRaster(pem_rast, "DailyLoop", tile_num, nam, "INT2U",
                   paste("-", "Pest event map -", nam))
      } else {
        SaveRaster(pem_rast, "DailyLoop", NA, names(pem_list[i]), "INT2U",
                   paste("-", "Pest event map -", nam))
      }
      
      # PEMS with climate stress exclusions
      if (exclusions_stressunits) {
        pem_excl1 <- brick(Rast_Subs_Excl(pem_rast, tile_num, "Excl1"))
        SaveRaster(pem_excl1, "DailyLoop", tile_num, 
                   paste0(nam, "_Excl1"), "INT2S", 
                   paste("-", log_capt, "for all", num_dats, "dates"))
        pem_excl2 <- brick(Rast_Subs_Excl(pem_rast, tile_num, "Excl2"))
        SaveRaster(pem_excl2, "DailyLoop", tile_num, 
                   paste0(nam, "_Excl2"), "INT2S", 
                   paste("-", log_capt, "for all", num_dats, "dates"))
        }
      
    }
  }

  # Free up memory
  rm(list=ls(pattern = "_brick"))
  
  delfiles <- dir(pattern = "*xml") # Is this necessary? Only happens w/ some GDAL vs?
  suppressWarnings(file.remove(delfiles))
}

# (8). RUN THE DAILY LOOP ----
tic("Daily loop run time") # Start timing the daily loop run-time
LogWrap("\n", "DAILY TIME STEP", "\n")
LogWrap("", "Starting daily time step", "") 
 
# Run the DailyLoop function. If region is CONUS/EAST, then the 4 tiles will
# be run in parallel to significantly increase speed. 
tryCatch( { 
  
  if (region_param %in% c("CONUS", "EAST")) {
    # Why does foreach not work for running DailyLoop in parallel???
    # Have to "mclapply" instead, which does not run on Windows
    RegCluster(4) # 4 cores are needed
    # foreach(tile_num = 1:length(template), .packages = pkgs, 
    #         .inorder = FALSE) %dopar% {
    #   tile <- template[[tile_num]]
    #   DailyLoop(tile_num, tile) 
    #   }
    mclapply(1:length(template), function(tile_num) {
      tile <- template[[tile_num]]
      DailyLoop(tile_num, tile)
      }, mc.cores = 4)
    
    stopCluster(cl)

  } else {
    # "tile_num" is NA because there are no tiles
    DailyLoop(NA, template) 
  }
  
  },
    error = function(e) {
    LogWrap("", "Error in Daily Loop", "\n") 
    cat("\nError in Daily Loop\n")
})

# Remove memory eating climate raster lists
rm(tmin_list, tmax_list, tmean_list, tdmean_list, ppt_list, template)

# Document daily loop execution time
loop_exectime <- toc(quiet = TRUE)
loop_exectime <- (loop_exectime$toc - loop_exectime$tic) / 60 

LogWrap("\n\n", paste0("Daily loop done (run time = ", 
                     round(loop_exectime, digits = 2), " min)"), "\n") 
LogWrap("\n", "FINAL ANALYSES AND MAP PRODUCTION", "\n")
cat("\nDaily loop done (run time = ", round(loop_exectime, digits = 2), " min)",
    "\n\nFINAL ANALYSES AND MAP PRODUCTION\n", sep = "")

# (9). PROCESS DAILY LOOP RESULTS -----
setwd(output_dir)

tic("Data processing run time") # Start timing for data processing

# Create a directory ("Misc_output") to put secondary outfiles
dir.create("Misc_output")

#### * Merge and delete tiles (CONUS/EAST) ####
# If CONUS or EAST, merge the tiles
if (region_param %in% c("CONUS", "EAST")) {
  
  LogWrap("\n", paste("Merging tiles for", region_param), "\n\n")
  cat("\nMerging tiles for", region_param, "\n")
  
  # Get list of brick files for each type of files, and then split it up into
  # several parts so the tiles can be merged in parallel for increased speed
  brick_files <- list.files(pattern = glob2rx("*tile*tif$"), recursive = FALSE)
  type_list <- unique(str_split_fixed(brick_files, pattern = "_tile", 2)[,1])
  type_list_split <- split(type_list, ceiling(1:length(type_list)/3)) 
  
  # For each file type, merge the 4 tiles
  RegCluster(4)
  mrg_by_type <- foreach(type = type_list_split, .packages = pkgs, 
                         .inorder = FALSE) %dopar% {
    type_vec <- unname(unlist(type)) # Change to an unnamed vector
    
    for (t in type_vec) {
      CombineMaps(brick_files, t)
      LogWrap("", paste("Merged", t, "tiles"), "\n")
      # Delete tiles now that they're merged
      unlink(list.files(pattern = glob2rx(paste0("*", t, "_*tile*tif$"))))
    }
    
  }
  
  stopCluster(cl)
  rm(cl)
  LogWrap("\n", "Done merging and deleting tiles", "")
  cat("\nDone merging and deleting tiles\n")

}
    
# Some ggplot2 settings to use for summary (PNG) maps
# Map production in ggplot requires specifying plot.height and plot.width
# These need to be dynamic because regions have different aspect ratios, 
# which result in some warped looking maps
# Calculate bounding box (xmin, xmax, ymin, ymax) of REGION 
coord <- coord_quickmap(xlim = c(REGION@xmin, REGION@xmax), 
                        ylim = c(REGION@ymin, REGION@ymax), expand = FALSE)
asp <- coord$aspect(list(x.range = c(REGION@xmin, REGION@xmax), 
                         y.range = c(REGION@ymin, REGION@ymax))) # aspect ratio

# Adjust base_size for ggplot2 (font size) according to aspect ratio
if (asp >= 1.7) {
  base_size <- 10.5
  legend_units <- 1.4
} else if (asp >= 1.5 & asp < 1.7) {
  base_size <- 9.5
  legend_units <- 1.3
} else if (asp >= 1.2 & asp < 1.5) {
  base_size <- 8.5 
  legend_units <- 1.2
} else if (asp < 1.2 & asp >= 0.6) {
  base_size <- 7
  legend_units <- 1
} else if (asp < 0.6 & asp >= 0.45) {
  base_size <- 5.7
  legend_units <- 0.6
} else if (asp < 0.45) {
  base_size <- 5.2
  legend_units <- asp
}

# Theme to use for summary maps (ggplot)
mytheme <- theme(
  legend.text = element_text(size = rel(1)), 
  legend.title = element_text(size = rel(1.2), face = "bold"),
  legend.position = "right", 
  legend.justification = "left",
  legend.margin = margin(t = 0, r = 0.10, b = 0, l = 0.10, unit = "cm"),
  legend.key.width = unit(legend_units, "line"), 
  legend.key.height = unit(legend_units, "line"),
  plot.title = element_text(size = rel(1.55), face = "bold", hjust = 0.5, 
                            vjust = -3, lineheight = 1, 
                            margin = margin(t = 0, r = 0, b = 2, l = 0)), 
  plot.subtitle = element_text(size = rel(1.25), hjust = 0.5, vjust = -3, 
                              lineheight = 1, 
                              margin = margin(t = 5, r = 0, b = 15, l = 0)),
  plot.margin = margin(t = 0.05, r = 0.25, b = 0.05, l = 0.25, unit = "cm"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  panel.background = element_blank(), panel.border = element_blank(),
  axis.title.x = element_blank(), axis.title.y = element_blank(), 
  axis.ticks = element_blank(),
  axis.text.x = element_blank(), axis.text.y = element_blank())

### * Generate summary maps ####
if (exclusions_stressunits) {
  LogWrap("\n\n", paste0("### SUMMARY MAPS: DDTOTAL, CUMDDs, CUMDDLW, LIFESTAGE, 
  NUMGEN, CLIMATE STRESS EXCL., CLIMATE STRESS UNITS ###"), "")
  cat("\n", str_wrap("### SUMMARY MAPS: DDTOTAL, CUMDDs, CUMDDLW, LIFESTAGE, 
  NUMGEN, CLIMATE STRESS EXCL., CLIMATE STRESS UNITS ###", width = 80), "\n", sep = "")
} else {
  LogWrap("\n\n", "### SUMMARY MAPS: DDTOTAL, CUMDDs, CUMDDLW, LIFESTAGE, AND NUMGEN", "")
  cat("\nSUMMARY MAPS: SUMMARY MAPS: DDTOTAL, CUMDDs, CUMDDLW, LIFESTAGE, AND NUMGEN\n")
}

# Split up the dates into chunks (no chunks have single dates - this results in
# warning messages. Splitting up the dates avoids overloading the server w/ 
# running too many dates in parallel.
dats_list <- split(dats2, ceiling(seq_along(dats2)/(length(dats2)/4)))

# Create a data frame of life stages spelled out to be used for plotting
# Stage order is specific to the species, but also need a column for order
# that stages actually occur (life cycle steps)
stg_vals <- data.frame("stg_name" = gsub("O", "", stgorder), 
                       "stg_num" = c("1", "2", "3", "4", "5"), 
                       stringsAsFactors = FALSE) %>% 
  mutate(stg_name = mgsub(stg_name, c("E", "L", "P", "A"), 
                          c("eggs", "larvae", "pupae", "adults")),
         life_cycle = case_when((stg_name == "eggs") ~ 1, 
                                (stg_name == "larvae") ~ 2,
                                (stg_name == "pupae") ~ 3,
                                (stg_name == "adults") ~ 4)) 

# For each date in a date chunk, plot and save summary maps for main outputs
# (The only exceptions are PEMs and StageCount maps which are dealt with below)
RegCluster(round(ncores/3))

#for (dat in dats_list) {
summary_maps <- foreach(dat = dats_list, .packages = pkgs,
                         .inorder = TRUE) %:%
  foreach(d = unname(unlist(dat)), .packages = pkgs, .inorder = TRUE) %dopar% {

   dat_vec <- unname(unlist(dat)) # change to an unnamed date vector
  #for (d in dat_vec) {
      # Get position (layer) of date in raster brick
      lyr <- which(dats2 == d)
      # Make the plots
      # If plant disease (BOXB), use different titles than for insects
      PlotMap(subset(brick("DDtotal.tif"), lyr), d, 
            "Degree day (DD) accumulation", "Degree days", "Misc_output/DDtotal")
      PlotMap(subset(brick("CumDDs.tif"), lyr), d, 
            "Plant disease degree-day (DD) accum.", "Plant disease\ndegree days", 
            "Misc_output/Cum_DDs")
      PlotMap(subset(brick("Cum_Inf_Risk.tif"), lyr), d, 
            "Cumulative infection risk", "Cumulative\ninfection risk", "Cum_Inf_Risk")
      
      # Generate summary maps for Lifestage and Number of Generations
      PlotMap(subset(brick("Lifestage.tif"), lyr), d,
              "Lifestage", "Stage", "Misc_output/Lifestage")
      PlotMap(subset(brick("NumGen.tif"), lyr), d,
              "Number of generations", "No. of\ngenerations", "Misc_output/NumGen")
      
      # Climate stress exclusion and stress unit summary maps 
      if (exclusions_stressunits) {
        
        # For BOXB
        PlotMap(subset(brick("Cum_Inf_Risk_Excl1.tif"), lyr), d,
                "Cumulative infection risk w/ climate stress excl.", 
                "Cumulative infection\nrisk", "Cum_Inf_Risk_Excl1")
        PlotMap(subset(brick("Cum_Inf_Risk_Excl2.tif"), lyr), d,
                "Cumulative infection risk w/ climate stress excl.",
                "Cumulative infection\nrisk", "Cum_Inf_Risk_Excl2")

        # Generate summary maps for Lifestage and NumGen w/ climate stress excl.
        PlotMap(subset(brick("Lifestage_Excl1.tif"), lyr), d,
                "Lifestage w/ climate stress excl.", "Stage", 
                "Misc_output/Lifestage_Excl1")
        PlotMap(subset(brick("Lifestage_Excl2.tif"), lyr), d,
                "Lifestage w/ climate stress excl.", "Stage", 
                "Misc_output/Lifestage_Excl2")
        PlotMap(subset(brick("NumGen_Excl1.tif"), lyr), d,
                "No. of gens w/ climate stress excl.", 
                "No. of\ngenerations", "Misc_output/NumGen_Excl1")
        PlotMap(subset(brick("NumGen_Excl2.tif"), lyr), d,
                "NumGen w/ climate stress excl.", 
                "No. of\ngenerations", "Misc_output/NumGen_Excl2")
        
        # Exclusion maps (-1 = moderate; -2 = severe)
        excl_types <- c("Cold_Stress_Excl", "Heat_Stress_Excl", 
                        "Dry_Stress_Excl", "All_Stress_Excl")
        unit_types <- c("Cold_Stress_Units", "Heat_Stress_Units", "Dry_Stress_Units")
        max1 <- c(coldstress_units_max1, heatstress_units_max1, drystress_units_max1)
        max2 <- c(coldstress_units_max2, heatstress_units_max2, drystress_units_max2)
        
        for (i in 1:length(excl_types)) {
          titl <- str_to_sentence(gsub("_", " ", excl_types[i]))
          titl <- gsub("excl", "exclusion", titl)
          # Exclusion summary maps
          PlotMap(subset(brick(paste0(excl_types[i], ".tif")), lyr), d, 
                  titl, "Exclusion status", 
                  paste0("Misc_output/", excl_types[i]))
        }
        
        # Climate stress unit accumulation summary maps
        unit_types <- c("Cold_Stress_Units", "Heat_Stress_Units", "Dry_Stress_Units")
        max1 <- c(coldstress_units_max1, heatstress_units_max1, drystress_units_max1)
        max2 <- c(coldstress_units_max2, heatstress_units_max2, drystress_units_max2)
        
        for (i in 1:length(unit_types)) {
          titl <- str_to_sentence(gsub("_", " ", unit_types[i]))
          PlotMap_stress(subset(brick(paste0(unit_types[i], ".tif")), lyr), d, 
                         max1[i], max2[i], titl, titl, 
                         paste0("Misc_output/", unit_types[i]))
        }
        
    }
#}
}

stopCluster(cl)
rm(cl)

# Generation and stage maps (StageCount) require some extra processing

# Make list of Stage Count raster bricks for plotting
if (exclusions_stressunits) {
  StageCt_lst <- c("StageCount.tif", "StageCount_Excl1.tif", "StageCount_Excl2.tif")
} else {
  StageCt_lst <- c("StageCount.tif")
}

# Make summary maps
RegCluster(round(ncores/3))

for (i in 1:length(StageCt_lst)) {

  #for (dat in dats_list) {
  foreach(dat = dats_list, .packages = pkgs, .inorder = TRUE) %:%
    foreach(d = unname(unlist(dat)),
          .packages = pkgs, .inorder = TRUE) %dopar% {
    
   # dat_vec <- unname(unlist(dat)) # Change to an unnamed date vector
    
    #for (d in dat_vec) {
      
      # Get position (layer) of date in raster brick, subset layer from bricks,
      # convert this layer to a data frame, and mutate data to extract gen 
      # number and stage number. Then join to the stg_vals data frame to get
      # the stage name. These operations are done only if the data are not
      # climate stress exclusion values (i.e. are greater than 0). The stage
      # order (stg_order) column is for ordering stages correctly in map legend.
      lyr <- which(dats2 == d)
      
      StageCt_df <- ConvDF(brick(StageCt_lst[i])[[lyr]]) %>% 
        mutate(value_orig = value) %>% # Keep original value (for plotting)
        mutate(value = round(value, 1)) %>% # Round decimal (must do this)
        mutate(gen = as.numeric(ifelse(value > 0, 
                            str_split_fixed(value, "[.]", 2)[,1], value))) %>%
        mutate(stg_num = as.character(ifelse(value > 0, 
                              str_split_fixed(value, "[.]", 2)[,2], value))) %>%
        left_join(., stg_vals, by = "stg_num") %>%
        mutate(gen_stg = paste0(gen, ".", life_cycle)) # Correctly sort legend
          
      # Format data values for plotting - Gen 0 = OW gen and ordinal values for
      # other generations (1st, 2nd, 3rd, etc.). This step is pretty slow - 
      # see if it can be sped up in future version.
      OW_gen <- StageCt_df %>% 
        filter(gen == 0) %>%
        mutate(value = paste("OW gen.", stg_name))
      
      # Filter out climate stress values to tack on after formatting other data
      excl_df <- StageCt_df %>% 
        dplyr::filter(value < 0) %>%
        mutate(gen_stg = ifelse(value == -2, -2, -1))
          
      if (any(StageCt_df$gen > 0)) {
        StageCt_df2 <- StageCt_df %>% 
          filter(gen > 0) %>%
          mutate(value = 
                   map_chr(gen, function(x) paste(toOrdinal(x), "gen."))) %>% 
          unnest() %>% 
          mutate(value = paste(value, stg_name)) %>%
          rbind(OW_gen, .) 
      } else {
        StageCt_df2 <- OW_gen
      }
      
      # Add climate stress values back in
      if (exclusions_stressunits) {
        StageCt_df2 <- rbind(StageCt_df2, excl_df)
      }
      
      # Plot the results
      if (StageCt_lst[i] == "StageCount.tif") {
        PlotMap(StageCt_df2, d, "Generation and stage", "Gen. x stage", 
                "Misc_output/StageCount")
      } else if (StageCt_lst[i] == "StageCount_Excl1.tif") {
        PlotMap(StageCt_df2, d, 
                "Generation and stage w/ climate stress exclusion", 
                "Gen. x stage", "Misc_output/StageCount_Excl1")
      } else if (StageCt_lst[i] == "StageCount_Excl2.tif") {
        PlotMap(StageCt_df2, d, 
                "Generation and stage w/ climate stress exclusion", 
                "Gen. x stage", "Misc_output/StageCount_Excl2")
      }
  
    }
  #}
}

stopCluster(cl)
rm(cl)

# Log messages
if (exclusions_stressunits) {
  LogWrap("\n\n", "Done with DDtotal, CumDDs, CumDDLW, Lifesetage, NumGen,
          Generation and Stage, climate stress exclusions, and climate 
          stress unit maps", "\n")
  cat("\n", str_wrap("Done with DDtotal, CumDDs, CumDDLW, Lifestage, NumGen,
                     Generation and Stage, climate stress exclusions, and 
                     climate stress unit maps", width = 80), "\n", sep = "")
} else {
  LogWrap("\n", "Done with DDtotal, CumDDs, CumDDLW, Lifesetage, NumGen, 
          and Generation and Stage maps", "\n") 
  cat("\n", str_wrap("Done with DDtotal, CumDDs, CumDDLW, Lifestage, NumGen, and
                     Generation and Stage maps", width = 80), "\n", sep = "")
}

if (pems & !exclusions_stressunits) {
  LogWrap("\n", "### SUMMARY MAP OUTPUTS: PEST EVENT MAPS ###", "")
  cat("\nSUMMARY MAP OUTPUTS: PEST EVENT MAPS\n", sep = "")
} else if (pems & exclusions_stressunits) {
  LogWrap("\n", "### SUMMARY MAP OUTPUTS: PEST EVENT MAPS W/
          CLIMATE STRESS EXCL. ###", "\n")
  cat("\n\n", str_wrap("RASTER AND SUMMARY MAP OUTPUTS: PEST EVENT MAPS W/
                CLIMATE STRESS EXCL.\n", width = 80), sep = "")
}

#### * Pest Event Maps ####

# Process and plot the Pest Event Maps (PEMs) 
if (pems) {
  # Get all PEM files, split them by type (e.g., PEMe1, PEMe2) and stage 
  # (e.g., "egg", "adult")
  PEM_files <-  list.files(pattern = glob2rx("*PEM*.tif$")) # all PEM files
  PEM_types <- str_split_fixed(PEM_files, ".tif", 2)[,1] # split by type
  
  # Create a data frame with PEM labels - the labels will be joined to the 
  # appropriate PEM file below
  PEM_event_labels <- cbind(data.frame("pem_types" = PEM_types), 
    data.frame("gen" = substr(PEM_types, start = 5, stop = 5)), 
    data.frame("stg" = substr(PEM_types, start = 4, stop = 4)))
  PEM_event_labels <- PEM_event_labels %>%
    mutate(genLabel = ifelse(gen == 0, "Date of OW gen.", 
                      ifelse(gen == 1, "Date of 1st gen.", 
                      ifelse(gen == 2, "Date of 2nd gen.", 
                      ifelse(gen == 3, "Date of 3rd gen.", 
                      ifelse(gen == 4, "Date of 4th gen.", NA)))))) %>% 
    mutate(eventLabel = ifelse(stg == "e", eggEventLabel, 
                        ifelse(stg == "l", larvaeEventLabel, 
                        ifelse(stg == "p", pupaeEventLabel, 
                        ifelse(stg == "a", adultEventLabel, NA))))) %>%
    mutate(finalLabel = paste(genLabel, eventLabel, sep = " ")) %>% 
                      dplyr::select(pem_types, finalLabel)
  
  # Which PEM is for the OW stage?
  OW_pem <- paste0("PEM", tolower(substr(owstage, start = 2, stop = 2)), "0")
  
  # Analyze the PEMs of each type (PEMa0, PEMp0, etc.) and then export resulting 
  # rasters and generate summary maps
  RegCluster(length(PEM_types))
  
  foreach(type = PEM_types, .packages = pkgs, .inorder = FALSE) %dopar% {
  #for (type in PEM_types) {
  #print(type)
    files_by_type <- PEM_files[grep(pattern = paste0(type, ".tif"), 
                                    x = PEM_files, fixed = TRUE)]
    # Change 0 values to NA
    PEM_brk <- brick(raster::stack(files_by_type))
    PEM_brk[PEM_brk == 0] <- NA
    
    # Create event label to be used for making summary maps, and plot map
    eventLabel_df <- dplyr::filter(PEM_event_labels, PEM_types == type) %>% 
      dplyr::select(finalLabel) %>%
      mutate(., finalLabel = ifelse(type == OW_pem, paste("Date of OW gen.", 
                                            OWEventLabel), finalLabel))          
    eventLabel <- paste(eventLabel_df$finalLabel)
    
    # Plot summary maps; if PEM has climate stress excl. modify labels/outname
    if (grepl("Excl", type)) {
      titl <- paste0(eventLabel, " w/ climate stress exclusion")
    } else {
      titl <- eventLabel
    }
           
    PlotMap(PEM_brk, last(dats2), titl, paste(eventLabel, sep = " "), 
            paste0("Misc_output/", type))
    
  }
}
#}

stopCluster(cl)
rm(cl)
  
# Log file messages
if (pems == 1) {
  LogWrap("\n\n", "Done with Pest Event Maps", "")
  cat("\nDone with Pest Event Maps\n")  
} 

# (10). Wrap-up -----
processing_exectime <- toc(quiet = TRUE)
processing_exectime <- (processing_exectime$toc - processing_exectime$tic) / 60 

LogWrap("\n", "### Done w/ processing of daily loop results ###", "\n")
LogWrap("", paste0("Run time for analyses and map production = ", 
    round(processing_exectime, digits = 2), " min"), "\n")
cat("\nDone w/ final analyses and map production\n\n", 
    "Run time for analyses and mapping run time = ", 
    round(processing_exectime, digits = 2),
    " min\n\n", "Deleting, renaming, and moving miscellaneous files\n\n", 
    sep = "")

#### * Rename final files and move misc files ####

# Create list of files that will be kept in the main output folder. These 
# include outputs for the last day of the sampled time period, with the 
# exception of Stage Count outputs.
last_dat_fls <- list.files(pattern = glob2rx(paste0("*", last(dats2), "*.png$")))

#stgCnt_remove <- grep(pattern = glob2rx(paste0("*StageCount*", last(dats2), 
#                                              "*")), last_dat_fls, value = TRUE)
#last_dat_fls <- last_dat_fls[!last_dat_fls %in% stgCnt_remove]

# If current year was sampled then keep infection risk outfile for the current day 
# Then make a list of final output files to rename.
if (start_year == current_year) {
  today_fls <- list.files(pattern = glob2rx(paste0("*Cum_Inf_Risk*", today, "*.png$")))
  onedayago_fls <- list.files(pattern = glob2rx(paste0("*Cum_Inf_Risk*", onedayago, "*.png$")))
  twodayago_fls <- list.files(pattern = glob2rx(paste0("*Cum_Inf_Risk*", twodayago, "*.png$")))
  final_fls <- c(last_dat_fls, today_fls, onedayago_fls, twodayago_fls)
} else {
  final_fls <- c(last_dat_fls)
}

# Put species abbreviation in final output files (rename)
new_names <- paste0(spp, "_", final_fls)
if (length(final_fls) > 0) {
  invisible(file.rename(final_fls, new_names))  
} else {
  LogWrap("\n", "No PNG files for final outputs - check for errors", "\n")
  cat("\nNo PNG files for final outputs - check for errors\n")
}

LogWrap("", paste0("Renamed all final PNG files to include ", spp, " in file name"), "\n")
cat("Renamed all final PNG files to include ", spp, " in file name\n", sep = "")

# All other misc files (w/out spp name in file name) are moved to "/Misc_output"
misc_fls <- grep(list.files(path = output_dir), 
                 pattern = spp, invert = TRUE, value = TRUE) 
misc_fls <- misc_fls[!(misc_fls %in% 
                         c("Misc_output", "previous_run", "Logs_metadata"))]
invisible(file.copy(misc_fls, paste0(output_dir, "/Misc_output/")))
invisible(file.remove(misc_fls))

# Wrap up log file and report time for entire model run
LogWrap("\n", "MODEL RUN DONE", "\n")
cat("\nMODEL RUN DONE\n")
total_exectime <- toc(quiet = TRUE) # Execution time for entire run
total_exectime <- round((total_exectime$toc - total_exectime$tic) / 60, digits = 2)
LogWrap("", paste("Run time for entire model =", total_exectime, "min"), "")
cat("\nRun time for entire model =", total_exectime, "min\n\n")

# Clean up
rm(list = ls(all.names = TRUE)) # Clear all objects including hidden objects
gc()