#!/usr/bin/Rscript
.libPaths("/usr/local/lib64/R/library/")
options(echo=FALSE)

# Install pakcages if necessary and load them
pkgs <- c("sp", "rgdal", "raster", "lubridate","dplyr", "stringr", "ggplot2", "optparse",
          "ggthemes","mapdata","foreach", "parallel", "RColorBrewer","tictoc","tidyr")
ld_pkgs <- lapply(pkgs, library, character.only = TRUE) # load them

# Bring in states feature for summary maps (PNG files), extract the lower 48 states, and project it
# Requires these libraries: "mapdata" and "maptools"
cat("Downloading US states feature\n")
states <- map_data("state")

################################################################################
########                  Header Documentation Section             #############
#  Cohorts DDRP: Degree-Day, establishment Risk, and Pest event mapping system #
########  By Brittany Barker, Len Coop, Gericke Cook, Dan Upper, and ###########
########  Tyson Wepprich for APHIS PPQ and IPM needs ###########################
################################################################################

# 11v. This is for development of a plant disease risk version of DDRP so use prism data such as:
#   A) Add Drystress looking similar to Heat and Cold (Chill) stress, using relhum=func(tmean,tdmean)
#   B) Add DDLW unique to Boxwood blight since BoxDD is a piecewise regression formula based on BB DD lookup table used in site vers.
#"PRISM_tddif_early_4kmD2_MTD_20210401=(PRISM_tmean_early_4kmD2_MTD_20210401-PRISM_tdmean_early_4kmD2_MTD_20210401)"
#"LW_20210401 = if(PRISM_ppt_early_4kmD2_MTD_20210401 > 8,1,if(PRISM_tddef_early_4kmD2_MTD_20210401 < 5,2,0))"
# the thresholds values can be calibrated and included as input params from the spp_param files
#

# (1). PARAM HANDLING -----
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
  ncohort <- opts$ncohort
  odd_gen_map <- opts$odd_gen_map
} else {
  #### * Default values for params, if not provided in command line ####
  spp           <- "BOXB" # Default species to use
  forecast_data <- "PRISM" # Forecast data to use (PRISM or NMME)
  start_year    <- "2016" # Year to use
  start_doy     <- 1 # Start day of year          
  end_doy       <- 365 # End day of year - need 365 if voltinism map 
  keep_leap     <- 1 # Should leap day be kept?
  region_param  <- "OR" # Region [CONUS, EAST, WEST, or state (2-letter abbr.)]
  exclusions_stressunits    <- 1 # Turn on/off climate stress unit exclusions
  pems          <- 0 # Turn on/off pest event maps
  mapA          <- 1 # Make maps for adult stage
  mapE          <- 0 # Make maps for egg stage
  mapL          <- 0 # Make maps for larval stage
  mapP          <- 0 # Make maps for pupal stage
  out_dir       <- "BOXB_test" # Output dir
  out_option    <- 1 # Sampling frequency
}


# (2). FUNCTIONS AND PLOT SETTINGS -----
Cond <- function(condition, trueValue, falseValue) {
  return(condition * trueValue + (!condition)*falseValue)
}

#### DD Calculation Methods: Tmean (no upper threshold, Tmax+Tmin AVG w/horiz upper threshold, 
#   and Single Triangle w/horiz. upper threshold

# Simple Mean Temp DD Calc method: ((tmean > LDT) * (tmean - LDT))
# same result as  max((tmax + tmin)/2 - LDT,0)  
# so no need for tmean PRISM data. 
SimpDD <- function(tmax,tmin,LDT) {
  return(max((tmax + tmin)/2 - LDT,0))
}

# Averaging DD Calc method (max + min/2 - tlow) but with horizontal (substitution) upper threshold:
AvgDD <- function(tmax, tmin, LDT, UDT) {
  return(Cond(tmax < LDT, 0, Cond(tmin > UDT, 0, Cond(tmax > UDT, (UDT + tmin)/2 - LDT, Cond((tmax + tmin)/2-LDT < 0,0,(tmax + tmin)/2 - LDT)))))
}

# Single triangle with upper threshold (Sevachurian et al. 1977) - also a good substitution for 
# single sine method
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

# BOXB changes begin
#### Relative Humidity - from tmean and tdmean
# RH: =100*(EXP((17.625*TD)/(243.04+TD))/EXP((17.625*T)/(243.04+T)))
# calc RH for todays avg Temp and Dewpoint Temp.  -verified formula
# use ppt and RH to estimate leaf wetness (LW)
#Later add: r.mapcalc "LW=if(ppt > $ppt_T3,4,if(RH > $RH_T3,3,if(ppt > $ppt_T2 && RH > $RH_T2,2,if(ppt > $ppt_T1 && RH > $RH_T1,1,0))))"
#thresholds for DDs; DONE: replace with lookup table or piecewise regression approximation of infection curve
calcrelhum=function(tmean,tdmean){
  return(100*(exp((17.625*tdmean)/(243.04+tdmean))/exp((17.625*tmean)/(243.04+tmean))))
}

ppt_T1 <- 0.3
ppt_T2 <- 8
ppt_T3 <- 30
RH_T1 <- 68
RH_T2 <- 82
RH_T3 <- 88

calcLW <- function(ppt,relhum){
  return(Cond(ppt > ppt_T3,4,
              Cond(relhum > RH_T3,3,
                   Cond(ppt > ppt_T2 & relhum > RH_T2,2,
                        Cond(ppt > ppt_T1 & relhum > RH_T1,1,0)))))
}

calcDDs <- function(tmax,tmin){
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

#### Functions for creating summary maps and saving rasters and maps ####
# Save raster - simply the "writeRaster" function from raster library but prints progress in log file
SaveRaster <- function(r, d, outnam) {
  writeRaster(r, file = paste0(outnam, "_", d), format="GTiff", overwrite=TRUE)
  print(noquote(paste0("Saving raster: ", outnam,"_", d,".tif")))
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
  
  df <- ConvDF(r) # convert raster to data frame
  sp <- paste0(gsub(pattern = "_", replacement = " ", fullname),":")
  nam <- deparse(substitute(r)) # Get name of raster, capitalize first letter
  dat <- as.character(format(strptime(d,format="%Y%m%d"), format="%m/%d/%Y")) # format the date
  titl <- paste(titl,dat,sep=" ")
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
  #
  #### * DD accumulations ####
  if (grepl("ddtotal|CUMDDs|CUMDDLW", outfl)) {
    
    # Use pre-defined bins
    if (grepl("ddtotal|CUMDDs", outfl)) {
      bins <- c("0", "1-250", "251-500", "501-750", "751-1000", "1001-1500", "1501-2000",
                "2001-3500", "3501-5000", "5001-6000", ">6000")
      vals <- c(250, 500, 750, 1000, 1500, 2000, 3500, 5000, 6000)
    } else {
      bins <- c("0", "1-10", "11-25", "26-50", "51-100", "101-300", "301-600",
                "601-1000", "1001-2000", "2001-3000", ">3000") 
      vals <- c(10, 25, 50, 100, 300, 600, 1000, 2000, 3000)
    }
    bins <- factor(bins, levels = bins)
    
    # Assign bins to data and order factor levels
    df2 <- df %>% mutate( value = case_when(value == 0 ~ bins[1],
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
    df2$value <- factor(df2$value, levels = levels(bins))
    
    # Color bins
    #cols <- setNames(colorRampPalette(c("blue3", "deepskyblue", "lightyellow", "gold", "red3"),
    #                                   interpolate = "spline", space = "rgb")(11), levels(df2$value))
    cols <-  colorRampPalette(rev(brewer.pal(11, "Spectral")))(11)
    cols <- setNames(cols, levels(df2$value))
    
    # Make the plot
    p <- Base_map(df2) + 
      scale_fill_manual(values = cols, name = str_wrap(paste0(lgd), width = 15), drop = FALSE) +
      #scale_fill_manual(values = cols, name = str_wrap(paste0(lgd), width = 15)) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width), 
           subtitle = str_wrap(subtitl, width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme
    
    #### * Lifestage ####
  } else if ((grepl("Lifestage", nam))) {
    # recode OW in StageCount to -3 (currently -1) so there's no overlap w/ clim. exclusion values (-2 and -1)  
    if (!grepl("EXCL",nam)) { # if it's simply a lifestage raster, then recode -1 to -3
      df$value <- gsub("-1","-3", df$value)
    }
    df <- mutate(df, value = ifelse(value==-3,"OW",ifelse(value==-2,"excl.-severe",
                                                          ifelse(value==-1,"excl.-moderate",ifelse(value==0,"egg",
                                                                                                   ifelse(value==1,"larva",ifelse(value==2,"pupa","adult")))))))
    df$value <- factor(df$value, levels = c("excl.-severe","excl.-moderate","OW","egg",
                                            "larva", "pupa", "adult"))   
    # Make the color key for the legend
    col_key <- cbind(
      data.frame("cols" = c("gray40","gray70", "gray70","#E41A1C","#377EB8",
                            "#4DAF4A","#984EA3"), stringsAsFactors = FALSE),
      data.frame("value" = c("excl.-severe","excl.-moderate","OW","egg",
                             "larva", "pupa", "adult"))
    )
    col_df <- suppressWarnings(suppressMessages(semi_join(col_key,df,by="value"))) # join the colors to the data
    cols <- setNames(as.character(col_df$cols),col_df$value) # convert col_df to a named vector
    
    # Make the plot
    p <- Base_map(df) + 
      scale_fill_manual(values = cols, name = str_wrap(paste0(lgd), width = 15)) +
      labs(title = str_wrap(paste(sp, titl), width = titl_width), 
           subtitle = str_wrap(subtitl, width = subtitl_width)) +
      theme_map(base_size = base_size) + 
      mytheme
    
    #### * Number of generations ####
  } else if (grepl("NumGen|NumGenExcl1|NumGenExcl2", outfl)) {
    
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
    
    #### * Stage count ####
  } else if (grepl("StageCount", outfl)) {
    
    # recode OW in StageCount to -3 (currently -1) so there's no overlap w/ clim. exclusion values (-2 and -1)  
    if (!grepl("EXCL",nam)) { # if it's simply a stage count raster, then recode -1 to -3
      df$value <- gsub("-1","-3", df$value)
    }
    df$value_orig <- as.numeric(df$value) # create for ordering recoded values
    df$value <- as.character(as.numeric(df$value))  
    # create 3 data frames with needed codes (clim. excl. codes [climEXCL], overwintering [OW], and generations [Gens]) for each unique raster value, to be a key
    climEXCL <- cbind(data.frame("gen"=c("excl.","excl."), stringsAsFactors = FALSE), 
                      data.frame("value"=as.character(-2:-1), stringsAsFactors = FALSE),
                      data.frame("stg"=c("-severe","-moderate"), stringsAsFactors = FALSE))
    OW <- cbind(data.frame("gen"=c("OW","OW","OW","OW","OW"), stringsAsFactors = FALSE), # allow up to 12 generations
                data.frame("value"=c("-3","0","1","2","3"), stringsAsFactors = FALSE),
                data.frame("stg"=c("","eggs","larvae","pupae","adults"), stringsAsFactors = FALSE))
    Gens <- cbind(data.frame("gen"=paste0("G",rep(1:(48/4),each=4)), stringsAsFactors = FALSE), # allow up to 25 generations
                  data.frame("value"=as.character(4:51), stringsAsFactors = FALSE),
                  data.frame("stg"=rep_len(c("eggs","larvae","pupae","adults"), length.out=48), stringsAsFactors = FALSE))
    climEXCL_OW_Gens <- bind_rows(climEXCL,OW,Gens) # combine all 3
    # join this "key" to the real data, and then format the values
    df2 <- suppressWarnings(suppressMessages(left_join(df,climEXCL_OW_Gens,by="value")))
    df2$value <- paste(df2$gen,df2$stg,sep=" ") # create Gen x Stg factor, the new value to be plotted
    df2$value <- gsub("excl. ","excl.", df2$value) # remove a space
    df2$value <- factor(df2$value, levels = unique(df2$value[order(df2$value_orig)])) # order by original values so plots in numerical order
    # make the color ramp (excl.-severe, excl.-moderate, and OW are gray shades)
    col_key <- cbind(  
      data.frame("cols" = 
                   c("gray40","gray70","gray70", Colfunc("mediumpurple4","mediumpurple1", 4),
                     Colfunc("magenta4","magenta", 4),Colfunc("midnightblue","lightskyblue1", 4),
                     Colfunc("cyan4","cyan", 4), Colfunc("darkgreen","lightgreen", 4),
                     Colfunc("olivedrab4","olivedrab1", 4), Colfunc("darkgoldenrod3","yellow", 4),
                     Colfunc("sienna4","sienna1", 4),Colfunc("darkred","red", 4),
                     Colfunc("burlywood4","burlywood", 4),Colfunc("rosybrown4","rosybrown1", 4),
                     Colfunc("deeppink4","lightpink", 4),Colfunc("blue","deepskyblue", 4))),
      "value_orig" = as.numeric(climEXCL_OW_Gens$value)
    )
    col_key <- suppressWarnings(suppressMessages(semi_join(col_key, df2, by="value_orig"))) # keep only the number of colors needed for the data
    cols <- setNames(as.character(col_key$cols), levels(df2$value)) 
    # Make the plot
    p <- Base_map(df2) + 
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
    
    #### * Climate stress exclusion maps ####
  } else if (grepl("Heat_Stress_Excl|Cold_Stress_Excl|All_Stress_Excl", 
                   outfl)) {
    
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
    
    #### * Pest Event Maps ####
  } else if (grepl("PEM", outfl)) {
    start_year <- as.numeric(start_year)
    
    # Format the data value column
    df <- df %>% 
      dplyr::filter(!(value %in% c(0, 366))) # remove day 0 and day 366
    df$value <- round(df$value)
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
                          as.Date(paste0(start_year, "-01-01"))) 
    
    # The resulting values may have dates before Jan-1 (Dec-30, Dec-31) because 
    # they occur in the same week as Jan-1. The ceiling_date function (dplyr)
    # rounds up to the next month so that they begin on Jan-1 
    # (e.g., 2019-12-31 == 2020-01-01)
    df$value <- as.character(cut.POSIXt(strptime(df$value, format = "%Y-%m-%d"), 
                                        breaks = "1 weeks"))
    badDates <- df %>% filter(grepl(as.character(start_year - 1), value))
    badDates$value <- as.character(ceiling_date(as.Date(badDates$value, 
                                                        format = "%Y-%m-%d"), "month"))
    
    # Now replace old data with data that have fixed date
    # Then add week of month column
    df <- df %>% filter(!(grepl(as.character(start_year - 1), value))) %>%
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
    weeks_df <- data.frame("mnth" = c(rep("Jan", 5), rep("Feb", 5), 
                                      rep("Mar", 5), rep("Apr", 5), rep("May", 5), rep("Jun", 5), rep("Jul", 5),
                                      rep("Aug", 5), rep("Sep", 5), rep("Oct", 5), rep("Nov", 5), 
                                      rep("Dec", 5))) # 5 weeks per month
    weeks_df <- data.frame(weeks_df %>% group_by(mnth) %>% # Group by month
                             mutate(mnth_wk = row_number()) %>% # Assign unique row # to rep. week #
                             mutate(mnth_wk = paste(mnth, mnth_wk, sep = "_")))
    
    # Attach those data frames to make the key
    col_key <- cbind(cols_df, weeks_df)
    col_key$mnth <- as.character(col_key$mnth) # Add which month
    
    # Extract all unique weeks from data,and count no. of bins (months)
    dats <- df %>% distinct(value) 
    dats$mnth <- str_split_fixed(dats$value, pattern = "-", 2)[,1]
    dats <- dats %>% arrange(., value) %>%
      group_by(mnth) %>% # Group by month
      # Assign a unique value to each row in a month - used for joining later
      # This will also be done for other data frames below 
      
      #mutate(week = ceiling(day(mnth_day) / 7)) %>%
      left_join(dplyr::select(df, value, week), by = "value") %>%
      mutate(mnth_day = format(as.Date(value, "%b-%d"))) %>%
      mutate(mnth_wk = paste(mnth, week, sep = "_")) %>%
      distinct(., .keep_all = TRUE) %>%
      arrange(mnth_day)
    
    # Necessary for removing weeks that are not in the data in the col_key2
    mnth_ct <- data.frame(dats %>% group_by(mnth) %>% 
                            dplyr::mutate(freq = n_distinct(value))) %>%
      group_by(mnth) #%>% # Group by month
    
    # Finally filter unncessary weeks out of generic color key (some months 
    # have only 4 weeks)
    col_key2 <- dplyr::semi_join(col_key, dats, by = "mnth") %>% 
      dplyr::left_join(., mnth_ct, by = c("mnth_wk")) %>%
      na.omit %>%
      dplyr::select(cols, value)
    
    # Format the dates dataframe (dats2) for joining to col_key2 (color key)
    dats2 <- data.frame(dplyr::select(dats, value, mnth) %>% 
                          arrange(mnth, value))
    
    # Attach the colors to the value and format with needed colunms, etc.
    col_key2 <- left_join(col_key2, dplyr::select(dats2, -mnth), by = "value")
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
      df3$value <- factor(df3$value, levels = 
                            unique(df3$value[order(as.numeric(as.character(df3$value_orig)))]))
    } else {
      df3 <- df
      df3$value <- factor(df3$value, levels = 
                            unique(df3$value[order(as.numeric(as.character(df3$value_orig)))]))
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
  
  #### * Save the plots ####
  # Save the plot, or else report that there was an error and skip
  # See "rmessages.txt" for error report
  # save the plot, or else report that there was an error and skip - see "rmessages.txt" for error report
  tryCatch(
    {
      ggsave(p,file=paste0(outfl,"_",d,".png"), height = asp * 7, 
             units = c('in'), dpi = 250) # save the plot
      print(noquote(paste0("Saving summary map: ",outfl,"_",d,".png")))
    },
    error=function(e){
      print(noquote(paste0("Could not create plot for ",outfl,"_",d)))
    } )
}


# Create summary maps (PNG) of heat stress and chill stress units, with max1 (Stress limit 1) and max2 (Stress limit 2) shown as "countour" lines
PlotMap_stress <- function(r, d, max1, max2, titl, lgd, outfl) {
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
      ggsave(p,file=paste0(outfl,"_", d, ".png"), height = 7 * asp, 
             units = c('in'), dpi = 200) # save the plot
      print(noquote(paste0("Saving summary map: ",outfl,"_",d,".png")))
    },
    error=function(e){
      print(noquote(paste0("Could not create plot for ",outfl,"_",d)))
    } )
  
}

# Write raster (TIF) maps from within the daily loop
WriteMaps <- function(d){   
  setwd(output_dir)
  # always show lifestages, NumGen, StageCount
  if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays main maps
    Lifestage <- Cond(allminEXCL,Lifestage,0)
    NumGen <- Cond(allminEXCL,NumGen,0)
    # save the rasters - make lists and save w/ mapply
    rast_list <- list(allminEXCL,Lifestage,NumGen,StageCount,ddtotal)
    names(rast_list) <- c("All_Stg_Excl","Lifestage","NumGen","StageCount","ddtotal")
    invisible(mapply(SaveRaster, rast_list, d, paste(names(rast_list))))   
    # make and save the plots 
    PlotMap(allminEXCL,d,"All stage climate stress exclusion","Exclusion status")
    PlotMap(Lifestage,d,"Lifestage","Stage","Lifestage")
    PlotMap(NumGen,d,"Number of generations","No. generations","NumGen")
    PlotMap(StageCount,d,"Generation x stage count","Gen. x stage count","StageCount")
    PlotMap(ddtotal,d,"Degree day (DD) accumulation","DD accum.","ddtotal")
  }
  # BOXB changes begin
  else { # write regular old maps
    # save the rasters - make lists and save w/ mapply
    rast_list <- list(Lifestage,NumGen,StageCount,ddtotal,CUMDDs,CUMDDLW)
    names(rast_list) <- c("Lifestage","NumGen","StageCount","ddtotal","CUMDDs","CUMDDLW")
    invisible(mapply(SaveRaster, rast_list, d, paste(names(rast_list))))
    # make and save the plots 
    PlotMap(Lifestage,d,"Lifestage","Stage","Lifestage")
    PlotMap(NumGen,d,"Number of generations","No. generations","NumGen")
    PlotMap(StageCount,d,"Generation x stage count","Gen. x stage count","StageCount")
    PlotMap(ddtotal,d,"Degree day (DD) accumulation","DD accum.","ddtotal")
    PlotMap(CUMDDs,d,"Plant Disease degree day (DD) accumulation","PD-DDs","CUMDDs")
    PlotMap(CUMDDLW,d,"Plant Disease infection risk unit accumulation","DDLW_accum","CUMDDLW")
  }
  # BOXB changes end
  
  if (exclusions_stressunits == 1) { # make additional maps for "stageless" stress maps
    # exclusions included: -2=severe, -1=moderate stress
    # Lifestage, NumGen, and StageCount calculations
    LifestageEXCL1 <- Cond((AllEXCL > -2),Lifestage,-2)
    LifestageEXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,Lifestage))
    NumGenEXCL1 <- Cond((AllEXCL > -2),NumGen,-2)
    NumGenEXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,NumGen))
    StageCountEXCL1 <- Cond((AllEXCL > -2),StageCount,-2)
    StageCountEXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,StageCount))
    # BOXB changes begin
    CumInfRiskEXCL1 <- Cond((AllEXCL > -2),CUMDDLW,-2)
    CumInfRiskEXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,CUMDDLW))
    # save the rasters
    rast_list <- list(chillEXCL,heatEXCL,dryEXCL,AllEXCL,LifestageEXCL1,LifestageEXCL2,
                      NumGenEXCL1,NumGenEXCL2,StageCountEXCL1,StageCountEXCL2,chillunitsCUM,
                      heatunitsCUM,dryunitsCUM,CumInfRiskEXCL1,CumInfRiskEXCL2)
    names(rast_list) <- c("Chill_Stress_Excl","Heat_Stress_Excl","Dry_Stress_Excl",
                          "All_Stress_Excl","LifestageExcl1","LifestageExcl2","NumGenExcl1",
                          "NumGenExcl2","StageCountExcl1","StageCountExcl2","Chill_Stress_Units",
                          "Heat_Stress_Units","Dry_Stress_units","Cum_Inf_Risk_EXCL1",
                          "Cum_Inf_Risk_EXCL2")
    invisible(mapply(SaveRaster, rast_list, d, paste(names(rast_list))))
    # make and save the plots
    PlotMap(chillEXCL,d,"Chill stress exclusion","Exclusion status","Chill_Stress_Excl")
    PlotMap(heatEXCL,d,"Heat stress exclusion","Exclusion status","Heat_Stress_Excl")
    PlotMap(dryEXCL,d,"Dry stress exclusion","Exclusion status","Dry_Stress_Excl")
    PlotMap(AllEXCL,d,"All climate stress exclusion","Exclusion status","All_Stress_Excl")
    PlotMap(LifestageEXCL1,d,"Lifestage w/ climate stress excl.","Stage","LifestageExcl1")
    PlotMap(LifestageEXCL2,d,"Lifestage w/ climate stress excl.","Stage","LifestageExcl2")
    PlotMap(NumGenEXCL1,d,"No. of gens. w/ climate stress excl.", "No. generations","NumGenExcl1")
    PlotMap(NumGenEXCL2,d,"No. of gens. w/ climate stress excl.", "No. generations","NumGenExcl2")
    PlotMap(StageCountEXCL1,d,"Gen. x stage cnt. w/ climate stress excl.", "Gen. x stage count","StageCountExcl1")
    PlotMap(StageCountEXCL2,d,"Gen. x stage cnt. w/ climate stress excl.", "Gen. x stage count","StageCountExcl2")
    PlotMap(CumInfRiskEXCL1,d,"Cum. Infection Risk w/ sev. only clim. stress", "Cum. Infection Risk","CumInfRiskExcl1")
    PlotMap(CumInfRiskEXCL2,d,"Cum. Infection Risk w/ mod. & sev. climate stress", "Cum. Infection Risk","CumInfRiskExcl2")
    PlotMap_stress(chillunitsCUM,d,chillstress_units_max1,chillstress_units_max2,"Chill stress units","Chill stress units","Chill_Stress_Units")
    PlotMap_stress(heatunitsCUM,d,heatstress_units_max1,heatstress_units_max2,"Heat stress units","Heat stress units","Heat_Stress_Units")
    PlotMap_stress(dryunitsCUM,d,drystress_units_max1,drystress_units_max2,"Dry stress units","Dry stress units","Dry_Stress_Units")
  }
  # BOXB changes end
  # rast_label not in use delete all except this example:
  ##rast_label <- sprintf("%s%s%s",spp," Percent Egg Devel. on ",d)
  if (mapE == 1) { #show egg stage
    #rast_name <- sprintf("%s%s%s","EggDev_",d,".tif")
    #writeRaster(eggDev,paste("EggDev_",d,sep=""), format="GTiff", overwrite=TRUE)
    #PlotMap(d)
    if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
      eggP <- Cond(allminEXCL,eggP,0)
      SaveRaster(eggP,d, "EggP_Exclusions")
      PlotMap(eggP,d,"Prop. eggs w/ climate stress excl.","Prop. eggs","EggP_Exclusions")
    }
    else {    # no exclusions used this run
      SaveRaster(eggP,d, "EggP")
      PlotMap(eggP,d,"Proportion eggs","Prop. eggs","EggP")
    }
    if (exclusions == 1) {   # display eggmindays map
      SaveRaster(eggmindays,d, "Egg_chill_days")
      PlotMap(eggmindays,d,"Egg chill days","Chill days","Egg_chill_days")
      SaveRaster(eggminEXCL,d, "Egg_exclusion_areas")
      PlotMap(eggminEXCL,d,"Egg chill days w/ climate stress excl.","Exclusion status","Egg_exclusion_areas")
    }
  }
  if (mapL == 1) {
    if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
      larvP <- Cond(allminEXCL,larvP,0)
      SaveRaster(larvP,d, "LarvP_Exclusions")
      PlotMap(larvP,d,"Prop. larvae w/ climate stress excl.","Prop. larvae","LarvP_Exclusions")
    }
    else {    # no exclusions used this run
      SaveRaster(larvP,d, "LarvP")
      PlotMap(larvP,d,"Proportion larvae","Prop. larvae","LarvP")
    }
    if (exclusions == 1) {   # display larvaemindays map
      SaveRaster(larvaemindays,d, "Larv_chill_days")
      PlotMap(larvaemindays,d,"Larvae chill days","Larvae chill days","Larv_chill_days")
      SaveRaster(larvaeminEXCL,d, "Larv_exclusion_areas")
      PlotMap(larvaeminEXCL,d,"Larvae chill days w/ climate stress excl.","Exclusion status","Larv_exclusion_areas")
    }
  }
  # FIXME - Pupal exclusions should be fully fleshed out
  if (mapP == 1) {
    SaveRaster(pupP,d, "PupP")
    PlotMap(pupP,d,"Proportion pupae","Prop. pupae","PupP")
  }
  if (mapA == 1) {
    if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
      adultP <- Cond(allminEXCL,adultP,0)
      SaveRaster(adultP,d, "AdultP")
      PlotMap(adultP,d,"Proportion adults","Prop. adults","AdultP")
    }
    else {    # no exclusions used this run
      SaveRaster(adultP,d, "AdultP")
      PlotMap(adultP,d,"Proportion adults","Prop. adults","AdultP")
    }
    if (exclusions == 1) {   # display adultmindays map
      SaveRaster(adultmindays,d, "Adult_chill_days")
      PlotMap(adultmindays,d,"Adult chill days","Chill days","Adult_chill_days")
      SaveRaster(adultminEXCL,d, "Adult_exclusion_areas")
      PlotMap(adultminEXCL,d,"Adult chill days w/ climate stress excl.","Exclusion status","Adult_exclusion_areas")
    }
  }
  setwd(prism_dir)
}


# (3). DIRECTORY INIT ------

#### * Param inputs - species params; thresholds, weather, etc. ####
params_dir <- "/usr/local/dds/DDRP_B1/spp_params/"

#### * Weather inputs and outputs - climate data w/subdirs 4-digit year ####
# If outdir has 2 consec. numbers, assume webuser; otherwise just use base dir
#if (grepl("16", start_year, perl = TRUE)) {
#  base_dir <- "/mnt/ssd1/PRISM/"
#} else {
base_dir <- "/data/PRISM/"
#}
prism_dir <- paste0(base_dir, start_year)

cat("\nBASE DIR: ", base_dir, "\n")
cat("\nWORKING DIR: ", prism_dir, "\n")

#### * Output directory, log file, and error message file ####
# MUST remove .tif files or script will crash during processing because it will 
# try to analyze previously processed results. 

#output_dir <- paste0("/home/httpd/html/CAPS/", spp, "_cohorts")
output_dir <- paste0("/usr/local/dds/DDRP_B1/DDRP_results/", out_dir)

# If the directory already exists, then a backup directory will be created that
# contains the old run files. Old backup directories will be removed if present.
if (file.exists(output_dir)) {
  
  # Create a directory for the previous run to backup old run files to and 
  # delete old backup (if present)
  backup_dir <- paste0(output_dir, "/previous_run")
  unlink(backup_dir, recursive = TRUE)
  dir.create(backup_dir)
  
  # Make list of old directories from most recent run to be backed up
  old_subdirs <- list.dirs(path = output_dir, recursive = FALSE)
  old_subdirs <- old_subdirs[!grepl("previous_run", old_subdirs)]

  if (any(grepl("previous_run", old_subdirs))) {
    pos <- which(grepl("previous_run", old_subdirs))
    unlink(paste0(old_subdirs[pos], "/"), recursive = TRUE)
    old_subdirs <- old_subdirs[-pos] # Remove backup dir from subdir list
  }
  
  # Copy files from main dir and then subdirs to backup directory
  file.copy(list.files(output_dir, full.names = TRUE), backup_dir) # main dir
  for (i in 1:length(old_subdirs)) {
    subdir <-  old_subdirs[i] %>%
      str_split(pattern = "/") %>%
      unlist %>%
      last()
    subdir_copy <- paste(backup_dir, subdir, sep = "/")
    dir.create(subdir_copy)
    file.copy(list.files(old_subdirs[i], full.names = TRUE), subdir_copy)
    unlink(paste0(old_subdirs[i], "/"), recursive = TRUE)
  }
  
  # Remove old files now that they have been copied to backup folder
  unlink(paste0(output_dir, "/*"))  
  cat("\n", str_wrap(paste0("EXISTING OUTPUT DIR: ", output_dir, ";\n", 
                            "copied old run files to", backup_dir, "\n"), 
                     width = 80), sep = "") 
  
} else {
  dir.create(output_dir)
  cat("NEW OUTPUT DIR:", output_dir, "\n")
}

# Push out a rlogging file with all main messages in model
# Put all log, message, and metadata files in a separate folder
setwd(output_dir)
dir.create("Logs_metadata")
Model_rlogging <- sprintf("%s%s", "./", "/Logs_metadata/Model_rlogging.txt")

# Make header for logging file
cat(paste0(rep("#", 36), collapse = ""), "\n", 
    "### Log file for DDRP v2 ###\n", 
    paste0(rep("#", 36), collapse = ""), "\n\n", sep = "", 
    file = Model_rlogging)

# Record PRISM and output dir
cat("BASE DIR: ", base_dir, "\n", file = Model_rlogging, append = TRUE)
cat("WORKING DIR: ", prism_dir, "\n", file = Model_rlogging, append = TRUE)
cat(str_wrap(paste0("EXISTING OUTPUT DIR: ", output_dir, 
                    "; removing all files"), width = 80), "\n\n", sep = "", 
    file = Model_rlogging, append = TRUE)

# Push out a message file with all R error messages
msg <- file(paste0(output_dir, "/Logs_metadata/rmessages.txt"), open = "wt")
sink(msg, type = "message")

# (3). PARAMETER AND SETTINGS SETUP ----- 
cat("PARAMETER AND SETTINGS SETUP: getting parameters to use in model\n", 
    file = Model_rlogging, append = TRUE)
cat("\n\nPARAMETER AND SETTINGS SETUP: getting parameters to use in model\n")

# Read from source param files in ./spp_params/SPP.params
param_file <- sprintf("%s%s", spp, ".params")
spp <- gsub(".params", "", param_file) # Get species abbr.
species_params <- sprintf("%s%s", params_dir, param_file) # Location of file

if (file.exists(species_params)) {
  cat("Species params: ", species_params, "\n", file = Model_rlogging, 
      append = TRUE)
  source(species_params) # Read in species parameters
  cat("Reading params for species: ", spp, " Fullname: ", fullname, "\n", 
      file = Model_rlogging, append = TRUE)
  cat("\nReading params for species: ", spp, " Fullname: ", fullname, "\n")
} else {
  cat("Param file: ", species_params, "...not found; exiting program\n", 
      file = Model_rlogging, append = TRUE)
  cat("\nParam file: ", species_params, "...not found; exiting program\n")
  q()  # No reason to keep going without any params
}

# Change year to numeric if it's a specific year
# If using climate normals, there may be letters in folder name
if (!grepl("[A-z]", start_year)) {
  start_year <- as.numeric(start_year)
} else {
  start_year <- "30yr" # This needs to be changed depending on folder name
}

# Set up start and stop day of year depending on whether it's a leap year or
# not (need to modify the last day of the year if the year is a leap year 
# and user wants to include it)
# This does not apply to 30 yr climate data, which would not have a numeric
# start year
if (is.numeric(start_year)) {
  
  if (start_year %% 4 == 0 & keep_leap == 1) {
    cat(str_wrap(paste0(start_year, " is a leap year and leap day (2/29) will be 
                 included in the model"), width = 80), "\n", 
        file = Model_rlogging, append = TRUE)
    
    # Need to add an extra day onto year if all 365 days are being included
    if (end_doy == 365) {
      end_doy <- 366
    }
    
  } else if (start_year %% 4 == 0 & keep_leap == 0) {
    cat(str_wrap(paste0(start_year, " is a leap year but leap day (2/29) will 
                        not be included in the model"), width = 80), "\n", 
        file = Model_rlogging, append = TRUE)
  } else if (start_year %% 4 != 0) {
    cat(start_year, "is not a leap year - ignoring 'keep_leap' parameter\n", 
        file = Model_rlogging, append = TRUE)
  }
}

# Check for appropriate command line parameters
# Exit program if no pest event maps have been specified but users want pest
# event maps (pems = 1)
if (pems == 1 & !(1 %in% c(mapA, mapE, mapL, mapP))) {
  cat("\n", str_wrap("No pest event maps (mapA, mapE, mapL, mapP) specified; 
                     exiting program", width = 80), sep = "", 
      file = Model_rlogging, append = TRUE)
  cat("\n", str_wrap("No pest event maps (mapA, mapE, mapL, mapP) specified; 
          exiting program", width = 80), sep = "")
  q()
}

# Exit program if an incorrect sampling frequency has been specified
if (out_option %in% !c(1, 2, 3, 4, 5, 6)) {
  cat("Out_option =", out_option, "is unacceptable; exiting program\n", 
      file = Model_rlogging, append = TRUE)
  cat("Out_option =", out_option, "is unacceptable; exiting program\n")
  q() 
}

# Exit if end day of year is inappropriate
if (end_doy > 366) {
  cat("\n", str_wrap(paste("End day of year (end_doy) of", end_doy, "is 
                     unacceptable; exiting program"), width = 80), sep = "", 
      file = Model_rlogging, append = TRUE)
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

# (4). METADATA OUTPUT FILE -----

# Push out a metadata file with all inputs used in model
cat("\nMETADATA: creating metadata file for all inputs used in model\n", 
    file = Model_rlogging, append = TRUE)
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


cat("Done writing metadata file\n\n", forecast_data, " DATA PROCESSING\n", 
    sep = "", file = Model_rlogging, append = TRUE)
cat("\nDone writing metadata file\n\n", forecast_data, " DATA PROCESSING\n",
    sep = "")

# (5). WEATHER DATA LOADING AND PROCESSING -----

# Weather inputs and outputs - PRISM climate data w/subdirs 4-digit year
# New feature - choose whether to use PRISM or NMME for weather forecasts 
# (forecast_data = PRISM, or forecast_data = NMME)

# Loop through each needed variable and create a list of needed files
vars <- c("tmin", "tmax", "tmean", "tdmean", "ppt")
fls_list <- c("tminfiles", "tmaxfiles", "tmeanfiles", "tdmeanfiles", "ppt")

for (i in seq_along(vars)) {
  # Create list of files
  fls <- list.files(path = prism_dir, 
                    pattern = glob2rx(paste0("*PRISM_", vars[i], "*", start_year, "*.bil$*")), 
                    all.files = FALSE, full.names = TRUE, recursive = TRUE)
  
  # Exit program if files are missing
  if (length(fls) == 0) {
    cat("Could not find ", vars[i], " files - exiting program\n", 
        file = Model_rlogging, append = TRUE) 
    cat("Could not find ", vars[i],  " files - exiting program\n") 
    q()
  }
  
  # Extract highest quality files and assign object name to list
  ExtractBestPRISM(fls, forecast_data, keep_leap) [start_doy:end_doy]
  assign(fls, fls_list[i])
  
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

# Make vector of dates to use when processing results 
# The current date will be sampled if it's the current year AND if 
# the current day falls within the range of start_doy and end_doy.
# The last date of year will always be sampled.
# Using "unique" will only keep date if it doesn't already occur in vector
# This happens if the end day of year is a multiple of the sampling frequency 
# (e.g. 1 to 300, w/ a 30 day sampling frequency), or if the current date falls
# within the sampling frequency
today_dat <- strftime(Sys.time(), format = "%Y%m%d")
current_year <- strftime(Sys.time(), format = "%Y")

if (start_year == current_year & 
    yday(Sys.time()) >= start_doy &
    yday(Sys.time()) <= end_doy) {
  dats2 <- sort(as.numeric(unique(c(dats[seq(0, length(dats), sample_freq)], 
                                    today_dat, last(dats)))))
} else {
  dats2 <- sort(as.numeric(unique(c(dats[seq(0, length(dats), sample_freq)],
                                    last(dats)))))
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
cat("Finished loading ", forecast_data, " files for ", 
    length(start_doy:end_doy), " days\nCreating template file for ", 
    region_param, "\n", sep = "", file = Model_rlogging, append = TRUE)
cat("\nFinished loading ", forecast_data, " files for ", 
    length(start_doy:end_doy), " days\n\nCreating template file for ", 
    region_param, "\n", sep = "")

### * Create blank template from a temp file
# This template is used for cropping the temperature (tmin, tmax) rasters
# First define with extent of the region
#### Set up regions - use switch() (works like a single use hash) ####
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

#### 
# If CONUS or EAST, split template into tiles (and run in parallel)
# Benefit of tiles is lost for smaller regions, so these are not split

# Register DoParallel
# The "RegCluster" function creates a set of copies of R running in parallel and 
# communicating over sockets (parallel socket clusters). The value may be 
# specified manually; however, here the value is estimated based on the number 
# of available cores on the computer or server DDRP is being run on. 
# Specifying too many clusters will overload the computer.
ncores <- detectCores()

# The "RegCluster" function determines an appropriate # of cores depending on 
# the "region_param" and "ncohort" parameters, so the server doesn't become
# overloaded
RegCluster(round(ncores/4))

if (region_param %in% c("CONUS", "EAST")) {
  # Split template (2 pieces per side)
  tile_list <- SplitRas(template, ppside = 2, save = FALSE, plot = FALSE) 
  tile_n <- 1:length(tile_list) # How many tiles?
  cat("Splitting template into", length(tile_list), "tiles\n", 
      file = Model_rlogging, append = TRUE)
  
  # Name the 4 tiles (tile1, tile2, tile3, tile4)
  template <- mapply(function(n, t) {
    names(n) <- paste0("tile", t)
    return(n)
  }, n = tile_list, t = tile_n )
  
  rm(tile_list)
  
  # Crop temp files by each template tile
  cat("Cropping tmax, tmin, tmean, tdmean and ppt tiles for", region_param, "\n", 
      file = Model_rlogging, append = TRUE)
  cat("\nCropping tmax, tmin, tmean, tdmean and ppt tiles for", region_param, "\n")
  
  tmax_list <- foreach(tile = template, .packages = "raster", 
                       .inorder = FALSE) %:% 
    foreach(tmax = tmaxfiles, .packages = "raster", .inorder = TRUE) %dopar% { 
      m <- as.matrix(crop(raster(tmax), tile))
    }
  
  tmin_list <- foreach(tile = template, .packages = "raster", 
                       .inorder = FALSE) %:% 
    foreach(tmin = tminfiles, .packages = "raster", .inorder = TRUE) %dopar% { 
      m <- as.matrix(crop(raster(tmin), tile))
    }
  
  # If region is not CONUS or EAST, simply crop temp files by the single template
} else {
  cat("Cropping tmax and tmin tiles for", region_param, "\n", 
      file = Model_rlogging, append = TRUE)
  cat("\nCropping tmax and tmin tiles for", region_param, "\n")
  
  tmax_list <- foreach(t = tmaxfiles, .packages = "raster", 
                       .inorder = TRUE) %dopar% {
                         m <- as.matrix(crop(raster(t), template))
                       }
  
  tmin_list <- foreach(t = tminfiles, .packages = "raster", 
                       .inorder = TRUE) %dopar% {
                         m <- as.matrix(crop(raster(t), template))
                       }
  
}

stopCluster(cl)
rm(cl)

cat("Done processing ", forecast_data, " data\n\nDAILY LOOP\n", sep = "",
    file = Model_rlogging, append = TRUE)
cat("\nDone processing ", forecast_data, " data\n\nDAILY LOOP\n", sep = "")

# Some ggplot2 settings to use for summary (PNG) maps
# Map production in ggplot requires specifying plot.height and plot.width
# These need to be dynamic because regions have different aspect ratios, 
# which result
# Somes warped looking maps
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

# Theme to use for plots
mytheme <- theme(legend.text = element_text(size = rel(1)), 
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
                 axis.title.x = element_blank(), 
                 axis.title.y = element_blank(), 
                 axis.ticks = element_blank(),
                 axis.text.x = element_blank(), 
                 axis.text.y = element_blank())

tic("Daily loop run time") # Start timing the daily loop run-time
# cat("DAILY LOOP: daily loop log files show loop progress and 
# output file info\n", file = Model_rlogging, append = TRUE)
cat("Sampling every", sample_freq, "days between", first(dats), "and", 
    last(dats), "\n", file = Model_rlogging, append = TRUE) 
cat("\nSampling every", sample_freq, "days between", first(dats), "and", 
    last(dats), "\n") 


#### Initialize all tracking rasters as zero with the template ####
ddtotal <- template
DDaccum <- template

# BOXB changes begin
relhum <- template
LW <- template
DDs <- template
CUMDDs <- template
DDLW <- template
CUMDDLW <- template
# BOXB changes end

Lifestage <- template
Lifestage <- Lifestage -1  # Need for OW designated as stage -1
NumGen <- template
StageCount <- template

if (exclusions) {
  eggminEXCL        <- template  # final map with regions with mindays >= mindaysMAX
  eggmindays        <- template  # track no. days Tlow < eggLLT
  eggmindaysMAX     <- template
  eggmindaysMAX     <- eggLLDAYS # max no. days Tlow < eggLLT = exclusion
  larvaeminEXCL     <- template  # final map with regions with mindays >= mindaysMAX
  larvaemindays     <- template
  larvaemindaysMAX  <- template
  larvaemindaysMAX  <- larvaeLLDAYS # max no. days Tlow < larvaeLLT = exclusion
  #pupaeminEXCL      <- template  # final map with regions with mindays >= mindaysMAX
  #pupaemindays      <- template
  #pupaemindaysMAX   <- template
  #pupaemindaysMAX   <- pupaeLLDAYS # max no. days Tlow < pupaeLLT = exclusion
  adultminEXCL      <- template  # final map with regions with mindays >= mindaysMAX
  adultmindays      <- template
  adultmindaysMAX   <- template
  adultmindaysMAX   <- adultLLDAYS # max no. days Tlow < adultLLT = exclusion
  allminEXCL        <- template  # product of all EXCL maps
} else if (exclusions_stressunits) {
  # BOXB Changes begin
  # NEW approach dry/chill/heat stress units
  drymask         <- template  # binary mask for daily dry units
  drystress       <- template  # count of daily dry   units
  drystressTHRESH  <- template  # mask for drystress units threshold
  drystressTHRESH  <- drystress_threshold # mask for drystress units threshold
  dryunitsCUM     <- template  # cumulative dry units
  drystressMAX1    <- template  # use for max dry before mortality??
  drystressMAX1    <- drystress_units_max1 # use for max dry before mortality??
  drystressMAX2    <- template  # use for max dry before mortality??
  drystressMAX2    <- drystress_units_max2 # use for uncertainty zone max dry before mortality??
  dryEXCL         <- template  # EXCL map for dry stress
  # BOXB changes end
  chillmask         <- template  # binary mask for daily chill units
  chillstress       <- template  # count of daily chill units
  chillstressTHRESH  <- template  # mask for chillstress units threshold
  chillstressTHRESH  <- chillstress_threshold # mask for chillstress units threshold
  chillunitsCUM     <- template  # cumulative chill units
  chillstressMAX1    <- template  # use for max chill before mortality??
  chillstressMAX1    <- chillstress_units_max1 # use for max chill before mortality??
  chillstressMAX2    <- template  # use for max chill before mortality??
  chillstressMAX2    <- chillstress_units_max2 # use for uncertainty zone max chill before mortality??
  chillEXCL         <- template  # EXCL map for chilling
  heatmask          <- template  # binary mask for daily heat stress units
  heatstress        <- template  # use for heat stress units
  heatstressTHRESH  <- template  # mask for heatstress units threshold
  heatstressTHRESH  <- heatstress_threshold # mask for heatstress units threshold
  heatunitsCUM      <- template  # cumulative heat stress units
  heatstressMAX1    <- template  # use for max chill before mortality??
  heatstressMAX1    <- heatstress_units_max1 # use for max chill before mortality??
  heatstressMAX2    <- template  # use for max chill before mortality??
  heatstressMAX2    <- heatstress_units_max2 # use for max chill before mortality??
  heatEXCL          <- template  # EXCL map for heat stress             
  AllEXCL           <- template  # EXCL map for combined stresses (chill,heat,later: moisture)            
  # end NEW
}

#### If Pest Event Maps (PEMS) wanted then init PEM rasters ####
if (pems) { 
  if (eggEventDD) {  # must be specified in spp.params file
    eggEvent <- template # raster to track DDs for egg hatch PEMs
    eggEvent <- eggEvent + eggEventDD  # allow extra DDs for event to be triggered
    if (PEMnumgens > 0) {
      PEMe1 <- template  # egg DOYs for when cumDDs > eggEvent threshold  1st Gen
    } 
    if (PEMnumgens > 1) {
      PEMe2 <- template  # egg DOYs for when cumDDs > eggEvent threshold  2nd Gen
    } 
    if (PEMnumgens > 2) {
      PEMe3 <- template  # egg DOYs for when cumDDs > eggEvent threshold  3rd Gen
    } 
    if (PEMnumgens > 3) {
      PEMe4 <- template  # egg DOYs for when cumDDs > eggEvent threshold  4th Gen
    }
  }
  if (larvaeEventDD) {
    larvaeEvent <- template # raster to track DDs for larvae PEMs
    larvaeEvent <- larvaeEvent + larvaeEventDD  # mid-larval dev after ca. xx DDs 
    if (PEMnumgens > 0) {
      PEMl1 <- template  # larval DOYs for when cumDDs > larvaeEvent threshold  1st Gen
    } 
    if (PEMnumgens > 1) {
      PEMl2 <- template  # larval DOYs for when cumDDs > larvaeEvent threshold  2nd Gen
    } 
    if (PEMnumgens > 2) {
      PEMl3 <- template  # larval DOYs for when cumDDs > larvaeEvent threshold  3rd Gen
    } 
    if (PEMnumgens > 3) {
      PEMl4 <- template  # larval DOYs for when cumDDs > larvaeEvent threshold  4th Gen
    } 
  }
  if (pupaeEventDD) {
    pupaeEvent <- template # raster to track DDs for pupal PEMs
    pupaeEvent <- pupaeEvent + pupaeEventDD  # mid-pupal dev after ca. xx DDs 
    if (PEMnumgens > 0) {
      PEMp1 <- template  # pupal DOYs for when cumDDs > pupalEvent threshold  1st Gen
    } 
    if (PEMnumgens > 1) {
      PEMp2 <- template  # pupal DOYs for when cumDDs > pupalEvent threshold  2nd Gen
    } 
    if (PEMnumgens > 2) {
      PEMp3 <- template  # pupal DOYs for when cumDDs > pupalEvent threshold  3rd Gen
    } 
    if (PEMnumgens > 3) {
      PEMp4 <- template  # pupal DOYs for when cumDDs > pupalEvent threshold  4th Gen
    } 
  }
  if (adultEventDD) {
    adultEvent <- template # raster of DOY for adult emerge PEMs
    adultEvent <- adultEvent + adultEventDD # adult emerge occurs after xx DDs during adult stage
    if (PEMnumgens > 0) {
      PEMa1 <- template  # adult DOYs for when cumDDs > adultEvent threshold  1st Gen
    } 
    if (PEMnumgens > 1) {
      PEMa2 <- template  # adult DOYs for when cumDDs > adultEvent threshold  2nd Gen
    } 
    if (PEMnumgens > 2) {
      PEMa3 <- template  # adult DOYs for when cumDDs > adultEvent threshold  3rd Gen
    } 
    if (PEMnumgens > 3) {
      PEMa4 <- template  # adult DOYs for when cumDDs > adultEvent threshold  4th Gen
    } 
  }
}
#### END Pest Event Maps (PEMS) wanted then init PEM rasters ####

#### Using OVLP feature; init maps that track this #### 
if (ovlp > 0.0) {
  eggDev <- template # proportion egg has developed thus far
  eggP <- template # proportion in egg stage
  larvDev <- template # proportion larvae has developed thus far
  larvP <- template # proportion in larval stage
  pupDev <- template # proportion pupae have developed thus far
  pupP <- template # proportion in pupal stage
  adultDev <- template # proportion adult has developed thus far
  adultP <- template # proportion in adult stage
}

#### Init OWSTAGE maps #### 
if (owstage == "OE") {
  OWeggP <- template # proportion in OW egg stage
} else if (owstage == "OL") {
  OWlarvP <- template # proportion in OW larval stage
} else if (owstage == "OP") {
  OWpupP <- template # proportion in OW pup stage
} else if (owstage == "OA") {
  OWadultP <- template # proportion in OW adult stage
} 

#rm(template)
#   does this change with OW stage change ??
#Lifestage: [-1] = OW stage [0] = egg, [1] = larvae, [2] = pupae, [3] = adult

#### Accumulate degree days and reclass cells to NA with temperature exclusion ####
# NA for exclusion means that DD accum cannot occur anymore with incomplete
# generation development and no oviposition to next generation. This may
# be tested for sensitivity by allowing exlusion cells to reset to zero instead.

#### Set up start and stop dates for loop ####
doy <- start_doy
#start <- doy + 1
uniq_list <- unique(sortedlist)
sublist <- uniq_list[doy:end_doy]

#### Initialize stage specific lifestage binary rasters ####
# Limits operations to a mask for cells that are in that lifestage
# This is what allows for pixel by pixel tracking of what lifestage
# that cell is in - move updates for these to end of day

#### Init Lifestage tracking
if (owstage == "OE") {
  LSOW0 <- Lifestage == -1 
} else if (owstage == "OL") {
  LSOW1 <- Lifestage == -1 
} else if (owstage == "OP") {
  LSOW2 <- Lifestage == -1 
} else if (owstage == "OA") {
  LSOW3 <- Lifestage == -1 
}
LS0 <- Lifestage == 0
LS1 <- Lifestage == 1
LS2 <- Lifestage == 2
LS3 <- Lifestage == 3

# find todays date since we may want to use it for detn freq of output
today <- Sys.Date()
showtoday <- format(today, format="%Y%m%d")
yesterday <- format(today-1, format="%Y%m%d")
todayplus1 <- format(today+1, format="%Y%m%d")
todayplus2 <- format(today+2, format="%Y%m%d")
todayplus3 <- format(today+3, format="%Y%m%d")
todayplus7 <- format(today+7, format="%Y%m%d")
cat("today and tomorrow etc =",showtoday,todayplus1,todayplus2,todayplus7,"\n")  # *** 
########                         END initializations                               #########

########                 Begin of actual stepping through the model                #########
########                                                                           #########
cat("Begin stepping through model\n")  # *** 
tic() # start keeping track of model run time

for (d in sublist) {
  doy <- doy + 1
  #cat("d and doy =",d,doy,"\n")
  #cat("d,doy, filetype and file =",d," ",doy," ",quality[[d]]," ",maxfiles[[d]],"\n")  # *** 
  
  if ("d" == "showtoday") {
    cat("d, doy, showtoday =",d,doy,showtoday,"\n")  # *** 
  } else {
    #cat("d, doy, NOT showtoday =",d,doy,showtoday,"\n")  # *** 
  }
  
  if (.Platform$OS.type == "windows") flush.console()
  Sys.sleep(1)
  #Read in that day's PRISM raster files ## DEPRECATED PLAN TO REMOVE USE of PRISM_tmean
  if(calctype=="simple") { # only need these for simple DD calcs
    pattern = paste("(PRISM_tmean_)(.*)(",d,")(_bil.bil)$", sep="")
    temp <- list.files(pattern=pattern,all.files=FALSE, full.names=TRUE)
    #tmean <- raster(temp)
    tmean <- crop(raster(temp),REGION)
  }
  
  #### Prepare & crop to current region: tmax and tmin raster files ####
  temp <- maxfiles[[d]]
  tmax <- crop(raster(temp),REGION)
  #setwd(output_dir)
  #writeRaster(tmax,paste("Tmax_",d,sep=""), format="GTiff", overwrite=TRUE)
  #setwd(prism_dir)
  #dataType(tmax) <- "FLT4S"
  temp <- minfiles[[d]]
  tmin <- crop(raster(temp),REGION)
  # BOXB changes begin
  cat("d =",d,"\n")  # ***
  temp <- tmeanfiles[[d]]
  tmean <- crop(raster(temp),REGION) # NEW do I need this??
  #cat("tmeanfile =",tmean,"\n")  # ***
  temp <- tdmeanfiles[[d]]
  tdmean <- crop(raster(temp),REGION) # NEW do I need this??
  temp <- pptfiles[[d]]
  ppt <- crop(raster(temp),REGION) # NEW do I need this??
  #if(calctype=="simple") { # only need these for simple DD calcs
  #   tmean <- (tmax+tmin)/2
  #}
  #else {  # always need for NEW exclusion mapping
  #   tmean <- (tmax+tmin)/2
  #}
  # BOXB changes end
  
  #### Loop through Stages; order of stages now read from SPP param file ####
  #stgorder <- c("OE","L","P","A","E","F")
  #stgorder <- c("OL","P","A","E","L","F")
  
  for (i in stgorder) {  # Handle stages in the model
    ####  MAIN STEPS FOR EGG STAGE ####
    if (i == "E" | i == "OE") {   # Egg Stage
      if(calctype=="average") { #devel DDs (zero values for temps below LDT)
        dd0tmp <- AvgDD(tmax,tmin,eggLDT,eggUDT)
        #setwd(output_dir)
        #writeRaster(dd0tmp,paste("dd0tmp_",d,sep=""), format="GTiff", overwrite=TRUE)
        #rast_name <- sprintf("%s%s%s","dd0tmp_",d,".tif")
        #rast_label <- sprintf("%s%s%s",spp," Daily Degree-Days on ",d)
        #Sys.sleep(0.3)
        #PlotMap(rast_name) 
        #setwd(prism_dir)
      } else if(calctype=="triangle") {
        dd0tmp <- TriDD(tmax,tmin,eggLDT,eggUDT)
      } else { # assume (calctype=="simple") 
        dd0tmp <- SimpDD(tmean,eggLDT)
      }
      # BOXB changes begin
      relhum <- calcrelhum(tmean,tdmean)
      LW <- calcLW(ppt,relhum)
      DDs <- calcDDs(tmax,tmin)
      DDLW <- LW * DDs/24
      CUMDDs <- CUMDDs + DDs
      CUMDDLW <- CUMDDLW + DDLW
      # BOXB changes end
      if(exclusions) { #Calculate lower lethal threshold and exclusion mask
        # reminder: 0 means exclude; 1 means dont exclude
        #eggmin <- tmin > eggLLT  # old version eggmin
        # new section to accumulate days below threshold - how to increment
        eggmindays <- Cond(tmin <= eggLLT, eggmindays+1, eggmindays)
        eggminEXCL <- Cond(eggmindays >= eggLLDAYS, 0, 1) 
        # end new section to accumulate days below thresh.
        eggminEXCL[eggminEXCL==0] <- NA  # why does this need to be here??
        #writeRaster(eggmin,paste("EggMin_",d,sep=""), format="GTiff",overwrite=TRUE)
        #Calculate upper lethal threshold and exclusion mask (keep old method for now?)
        #eggmax <- tmax < eggULT
        #eggmax[eggmax==0] <- NA
        
      }
      #Apply exclusions and mask to daily DDs and limit to correct stage
      # NEED TO rethink how exclusions get used; should not be used here in fact
      #if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
      #  if (i == "OE") { dd0 <- dd0tmp * eggmin * eggmax * LSOW0 }
      #  else if (i == "E") { dd0 <- dd0tmp * eggmin * eggmax * LS0 }
      #}
      #else { # just use lifestage mask
      if (i == "OE") { dd0 <- dd0tmp * LSOW0 }
      else if (i == "E") { dd0 <- dd0tmp * LS0 }
      #}
      #Apply exclusions and lifestage mask to daily degree day surface
      #if (i == "OE") { dd0 <- dd0tmp * eggmin * eggmax * LSOW0
      #  } else if (i == "E") { dd0 <- dd0tmp * eggmin * eggmax * LS0
      #}
      #}
      #else { # just use lifestage mask
      # if (i == "OE") { dd0 <- dd0tmp * LSOW0
      ##   } else if (i == "E") { dd0 <- dd0tmp * LS0
      # }
      #}
      
      #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
      #count cells with dd>0
      dd.stat <- cellStats(dd0,stat='max',na.rm=TRUE)
      if (dd.stat > 0) {
        DDaccum <- DDaccum + dd0
        if (pems) {
          # set DOY in PEM if curr gen, not set yet, DDs > eventDDs, and Stage is eggs, else keep same value
          if (PEMnumgens > 0) {
            PEMe1 <- Cond(NumGen==1 & PEMe1 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe1) 
          }
          if (PEMnumgens > 1) {
            PEMe2 <- Cond(NumGen==2 & PEMe2 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe2) 
          }
          if (PEMnumgens > 2) {
            PEMe3 <- Cond(NumGen==3 & PEMe3 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe3) 
          }
          if (PEMnumgens > 3) {
            PEMe4 <- Cond(NumGen==4 & PEMe4 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe4) 
          }
        }
        #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
        if (i == "OE") {
          #cat("date and doy, i =",d,doy,i,"\n")
          progressOW0 <- (DDaccum * LSOW0) >= OWeggDD
          #writeRaster(progressOW0,paste("ProgressOW0_",d,sep=""), format="GTiff",overwrite=TRUE)
          Lifestage <- Cond(LSOW0 == 1 & progressOW0 == 1, 1, Lifestage)
          #Reset the DDaccum cells to zero for cells that progressed to next lifestage
          DDaccum <- Cond(progressOW0 == 1,(DDaccum - OWeggDD) * LSOW0, DDaccum)
        } else if (i == "E") {
          #cat("date and doy, i =",d,doy,i,"\n")
          progress0 <- (DDaccum * LS0) >= eggDD
          #writeRaster(progress0,paste("Progress0_",d,sep=""), format="GTiff",overwrite=TRUE)
          Lifestage <- Cond(LS0 == 1 & progress0 == 1, 1, Lifestage)
          #Reset the DDaccum cells to zero for cells that progressed to next lifestage
          DDaccum <- Cond(progress0 == 1,(DDaccum - eggDD) * LS0, DDaccum)
        }
      }
      
      ####  MAIN STEPS FOR LARVAL STAGE ####
    } else if (i == "L" | i == "OL") {  # Larval Stage
      #developmental degree days
      if(calctype=="average") {
        dd1tmp <- AvgDD(tmax,tmin,larvaeLDT,larvaeUDT)
      } else if(calctype=="triangle") {
        dd1tmp <- TriDD(tmax,tmin,larvaeLDT,larvaeUDT)
      } else { # assume (calctype=="simple") 
        dd1tmp <- SimpDD(tmean,larvaeLDT)
      }
      ddtotal <- ddtotal + dd1tmp #LEN : Accumulate total degree days for the year for larvae
      
      if(exclusions) { # Calculate lower lethal threshold and exclusion mask
        #larvaemin <- tmin > larvaeLLT
        # new section to accumulate days below threshold - how to increment
        larvaemindays <- Cond(tmin <= larvaeLLT, larvaemindays+1, larvaemindays)
        larvaeminEXCL <- Cond(larvaemindays >= larvaeLLDAYS, 0, 1) 
        # end new section to accumulate days below thresh.
        larvaeminEXCL[larvaeminEXCL==0] <- NA
        #Calculate upper lethal threshold and exclusion mask
        #larvaemax <- tmax < larvaeULT
        #larvaemax[larvaemax == 0] <- NA
        #writeRaster(larvaemax,paste("Larvaemax_",d,sep=""), format="GTiff",overwrite=TRUE)
      } 
      else if(exclusions_stressunits) {
        # NEW version exclusions - applies to all stages
        # BOXB Changes start: need to fix equations once we have relhum
        ##-- DDLW CALCULATION??--------------------
        ##-- Dry Stress Accumulation--------------------
        drymask <- relhum < drystressTHRESH  # make todays dry mask
        drystress <- drymask * ((abs(drystressTHRESH - relhum) - 500 * ppt)) # compute todays dry stress DDs
        #drystress <- drystress - 500 * ppt   # ppt is hundreds of an inch, it reduces drystess; TODO: calibrate
        drystress <- Cond(drystress >= 0, drystress, 0)
        dryunitsCUM <- dryunitsCUM + drystress
        # ASSUME NEW -2=severe -1=mod 0=none throughout
        dryEXCL <- Cond(dryunitsCUM >= drystressMAX2,-2,Cond(dryunitsCUM >= drystressMAX1,-1,0))
        # BOXB Changes end
        
        ##-- Chill Stress Accumulation--------------------
        chillmask <- tmin < chillstressTHRESH  # make todays chill mask
        chillstress <- chillmask * abs(chillstressTHRESH - tmin) # compute todays Chill stress DDs
        chillunitsCUM <- chillunitsCUM + chillstress
        # ASSUME NEW -2=severe -1=mod 0=none throughout
        chillEXCL <- Cond(chillunitsCUM >= chillstressMAX2,-2,Cond(chillunitsCUM >= chillstressMAX1,-1,0))
        
        ##-- Heat Stress Accumulation--------------------
        heatmask <- tmax > heatstressTHRESH  # make todays heat mask
        heatstress <- heatmask * abs(tmax - heatstressTHRESH) # compute todays heat stress DDs
        heatunitsCUM <- heatunitsCUM + heatstress
        heatEXCL <- Cond(heatunitsCUM >= heatstressMAX2,-2,Cond(heatunitsCUM >= heatstressMAX1,-1,0))
        # BOXB Changes begin
        # AllEXCL <- Cond((chillEXCL == 0) & (heatEXCL == 0),0,
        #           Cond((chillEXCL == -1) & (heatEXCL >= -1),-1,
        #           Cond((chillEXCL >= -1) & (heatEXCL == -1),-1,-2)))
        # BOXB Changes begin
        AllEXCL <- Cond((dryEXCL == 0) & (chillEXCL == 0) & (heatEXCL == 0),0,
                        Cond((dryEXCL == -2) | (chillEXCL == -2) | (heatEXCL == -2),-2,-1))
        # BOXB Changes end
        # end NEW exclusions
      }
      #Apply exclusions and mask to daily DDs and limit to correct stage
      #if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
      #  if (i == "OL") { dd1 <- dd1tmp * larvaemin * larvaemax * LSOW1 }
      #  else if (i == "L") { dd1 <- dd1tmp * larvaemin * larvaemax * LS1 }
      #}
      #else { # just use lifestage mask
      if (i == "OL") { dd1 <- dd1tmp * LSOW1 }
      else if (i == "L") { dd1 <- dd1tmp * LS1 }
      #}
      
      #Accumulate degree days, if dd1 > 0 otherwise exclusion masks get applied to this generation.
      dd.stat <- cellStats(dd1,stat='max',na.rm=TRUE)
      if (dd.stat > 0) {
        DDaccum <- DDaccum + dd1
        if (pems & larvaeEventDD) {
          # set DOY in PEM if curr. gen, not set yet, DDs > eventDDs, Stage is larvae, else keep same value
          if (PEMnumgens > 0) {
            PEMl1 <- Cond(NumGen==1 & PEMl1 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl1) 
          }
          if (PEMnumgens > 1) {
            PEMl2 <- Cond(NumGen==2 & PEMl2 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl2) 
          }
          if (PEMnumgens > 2) {
            PEMl3 <- Cond(NumGen==3 & PEMl3 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl3) 
          }
          if (PEMnumgens > 3) {
            PEMl4 <- Cond(NumGen==4 & PEMl4 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl4) 
          }
        }
        #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
        if (i == "OL") {
          progressOW1 <- (DDaccum * LSOW1) >= OWlarvaeDD
          #writeRaster(progressOW1,paste("ProgressOW1_",d,sep=""), format="GTiff",overwrite=TRUE)
          Lifestage <- Cond(LSOW1 == 1 & progressOW1 == 1, 2, Lifestage)
          #Reset the DDaccum cells to zero for cells that progressed to next lifestage
          DDaccum <- Cond(progressOW1 == 1,(DDaccum - OWlarvaeDD) * LSOW1, DDaccum)
        } else if (i == "L") {
          progress1 <- (DDaccum * LS1) >= larvaeDD
          #writeRaster(progress1,paste("Progress1_",d,sep=""), format="GTiff",overwrite=TRUE)
          Lifestage <- Cond(LS1 == 1 & progress1 == 1, 2, Lifestage)
          DDaccum <- Cond(progress1 == 1,(DDaccum - larvaeDD) * LS1, DDaccum)
        }
      }
      ####  MAIN STEPS FOR PUPAL STAGE ####
    } else if (i == "P" | i == "OP") {   # Pupal Stage
      #developmental degree days
      if(calctype=="average") {
        dd2tmp <- AvgDD(tmax,tmin,pupaeLDT,pupaeUDT)
      } else if(calctype=="triangle") {
        dd2tmp <- TriDD(tmax,tmin,pupaeLDT,pupaeUDT)
      } else { # assume (calctype=="simple") 
        dd2tmp <- SimpDD(tmean,pupaeLDT)
      }
      #Apply exclus (none for pupae stage) & lifestage mask to daily DD surface and limit to correct stage
      # new section to accumulate days below threshold - how to increment
      #pupaemindays <- Cond(tmin <= pupaeLLT, pupaemindays+1, pupaemindays)
      #pupaeminEXCL <- Cond(pupaemindays >= pupaeLLDAYS, 0, 1) 
      # end new section to accumulate days below thresh.
      #pupaemin[pupsemin==0] <- NA
      #Calculate upper lethal threshold and exclusion mask
      #pupaemax <- tmax < pupaeULT
      #pupaemax[pupaemax == 0] <- NA
      if (i == "OP") { dd2 <- dd2tmp * LSOW2
      } else if (i == "P") { dd2 <- dd2tmp * LS2
      }
      
      #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
      dd.stat <- cellStats(dd2,stat='max',na.rm=TRUE)
      if (dd.stat > 0) {
        DDaccum <- DDaccum + dd2
        if (pems & pupaeEventDD) {
          # set DOY in PEM if curr. gen, not set yet, DDs > eventDDs, Stage is pupa, else keep same value
          if (PEMnumgens > 0) {
            PEMp1 <- Cond(NumGen==1 & PEMp1 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp1) 
          }
          if (PEMnumgens > 1) {
            PEMp2 <- Cond(NumGen==2 & PEMp2 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp2) 
          }
          if (PEMnumgens > 2) {
            PEMp3 <- Cond(NumGen==3 & PEMp3 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp3) 
          }
          if (PEMnumgens > 3) {
            PEMp4 <- Cond(NumGen==4 & PEMp4 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp4) 
          }
        }
        #Calculate lifestage progression: for dd accum in correct lifestage, is it >= pupae DD threshold?
        if (i == "OP") {
          progressOW2 <- (DDaccum * LSOW2) >= OWpupDD
          #writeRaster(progressOW2,paste("ProgressOW2_",d,sep=""), format="GTiff",overwrite=TRUE)
          Lifestage <- Cond(LSOW2 == 1 & progressOW2 == 1, 3, Lifestage)
          #Reset the DDaccum cells to zero for cells that progressed to next lifestage
          DDaccum <- Cond(progressOW2 == 1,(DDaccum - OWpupDD) * LSOW2, DDaccum)
        } else if (i == "P") {
          progress2 <- (DDaccum * LS2) >= pupDD
          #writeRaster(progress2,paste("Progress2_",d,sep=""), format="GTiff",overwrite=TRUE)
          Lifestage <- Cond(LS2 == 1 & progress2 == 1, 3, Lifestage)
          DDaccum <- Cond(progress2 == 1,(DDaccum - pupDD) * LS2, DDaccum)
        }
      }
      ####  MAIN STEPS FOR ADULT STAGE ####
    } else if (i == "A" | i == "OA") {  # Adult stage, or time to 50% oviposition
      #developmental degree days
      if(calctype=="average") {
        dd3tmp <- AvgDD(tmax,tmin,adultLDT,adultUDT)
      } else if(calctype=="triangle") {
        dd3tmp <- TriDD(tmax,tmin,adultLDT,adultUDT)
      } else { # assume (calctype==simple) 
        dd3tmp <- SimpDD(tmean,adultLDT)
      }
      
      if(exclusions) { # Calculate lower lethal threshold and exclusion mask
        adultmin <- tmin > adultLLT
        # new section to accumulate days below threshold - how to increment
        adultmindays <- Cond(tmin <= adultLLT, adultmindays+1, adultmindays)
        adultminEXCL <- Cond(adultmindays >= adultLLDAYS, 0, 1) 
        adultminEXCL[adultminEXCL==0] <- NA
        #      writeRaster(adultmin,paste("Admin_",d,sep=""), format="GTiff",overwrite=TRUE)
        # NOW sum up all EXCL maps 
        allminEXCL <- (eggminEXCL * larvaeminEXCL * adultminEXCL)
        # end new section to accumulate days below thresh.
        #Apply exclusions and mask to daily DDs and limit to correct stage
      }
      #Apply exclusions and mask to daily DDs and limit to correct stage
      #if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
      #  if (i == "OA") { dd3 <- dd3tmp * adultmin * adultmax * LSOW3 }
      #  else if (i == "A") { dd3 <- dd3tmp * adultmin * adultmax * LS3 }
      #}
      #else { # just use lifestage mask
      if (i == "OA") { 
        print(LSOW3)
        dd3 <- dd3tmp * LSOW3 
      }
      else if (i == "A") { dd3 <- dd3tmp * LS3 }
      #}
      #  if (i == "OA") { dd3 <- dd3tmp * adultmin * LSOW3 }
      #  else if (i == "A") { dd3 <- dd3tmp * adultmin * LS3 }
      # }
      # else { # just use lifestage mask
      #  if (i == "OA") { dd3 <- dd3tmp * LSOW3 }
      #  else if (i == "A") { dd3 <- dd3tmp * LS3 }
      # }
      
      #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
      dd.stat <- cellStats(dd3,stat='max',na.rm=TRUE)
      print(dd.stat)
      if (dd.stat > 0) {
        DDaccum <- DDaccum + dd3
        if(pems & adultEventDD) { # record DOYs in PEMS for adult events
          if (PEMnumgens > 0) {
            PEMa1 <- Cond(NumGen==1 & PEMa1 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa1) 
          }
          if (PEMnumgens > 1) {
            PEMa2 <- Cond(NumGen==2 & PEMa2 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa2) 
          }
          if (PEMnumgens > 2) {
            PEMa3 <- Cond(NumGen==3 & PEMa3 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa3) 
          }
          if (PEMnumgens > 3) {
            PEMa4 <- Cond(NumGen==4 & PEMa4 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa4) 
          }
        }
        
        #Calculate lifestage progression: for dd accum in correct lifestage, is it >= adult DD threshold?
        if (i == "OA") {
          progressOW3 <- (DDaccum * LSOW3) >= OWadultDD
          #writeRaster(progressOW3,paste("ProgressOW3_",d,sep=""), format="GTiff",overwrite=TRUE)
          DDaccum <- Cond(progressOW3 == 1,(DDaccum - OWadultDD) * LSOW3, DDaccum)
          progressOW3[is.na(progressOW3)] <- template[is.na(progressOW3)]
          Lifestage <- Cond(LSOW3 == 1 & progressOW3 == 1, 0, Lifestage)
          Lifestage[is.na(Lifestage)] <- template[is.na(Lifestage)]
          #Increment NumGen + 1
          NumGen <- NumGen + progressOW3
          #writeRaster(NumGen,paste("NumGen",d,sep=""), format="GTiff",overwrite=TRUE)
        } else if (i == "A") {
          progress3 <- (DDaccum * LS3) >= adultDD
          #writeRaster(progress3,paste("Progress3_",d,sep=""), format="GTiff",overwrite=TRUE)
          #Reset the DDaccum cells to zero for cells that progressed to next lifestage
          DDaccum <- Cond(progress3 == 1,(DDaccum - adultDD) * LS3, DDaccum)
          #Remove masking effect from progress counter so it doesn't perpetuate through NumGen raster
          progress3[is.na(progress3)] <- template[is.na(progress3)]
          Lifestage <- Cond(LS3 == 1 & progress3 == 1,0, Lifestage)
          Lifestage[is.na(Lifestage)] <- template[is.na(Lifestage)]
          #Increment NumGen + 1
          NumGen <- NumGen + progress3
          #writeRaster(NumGen,paste("NumGen",d,sep=""), format="GTiff",overwrite=TRUE)
        }
        
      }
      
      ####  MAIN STEPS FOR END OF EACH DAY ####
    } else if (i == "F") { # end of the day placeholder 
      
      # update lifestage masks
      if (owstage == "OE") {
        LSOW0 <- Lifestage == -1
      } else if (owstage == "OL") {
        LSOW1 <- Lifestage == -1
      } else if (owstage == "OP") {
        LSOW2 <- Lifestage == -1
      } else if (owstage == "OA") {
        LSOW3 <- Lifestage == -1
      }
      LS0 <- Lifestage == 0
      LS1 <- Lifestage == 1
      LS2 <- Lifestage == 2
      LS3 <- Lifestage == 3
      
      # Need to do anything else (overlap and write maps)?
      if (doy %% 30 == 0 && monthlymaps == 1) { ## monthly
        cat("### Monthly maps output ### d, doy, i =",d,doy,i,"\n")
        do_ovlp <- 1     # OK to calc ovlp
        do_maps <- 1     # OK to write maps
      } else if (doy %% 14 == 0 && biweeklymaps == 1) { ## biweekly
        cat("### Biweekly maps output ### d, doy, i =",d,doy,i,"\n")
        do_ovlp <- 1     # OK to calc ovlp
        do_maps <- 1     # OK to write maps
      } else if (doy %% 10 == 0 && dekadmaps == 1) { ## every 10 days
        cat("### Dekad maps output ### d, doy, i =",d,doy,i,"\n")
        do_ovlp <- 1     # OK to calc ovlp
        do_maps <- 1     # OK to write maps
      } else if (doy %% 7 == 0 && weeklymaps == 1) { ## weekly
        cat("### Weekly maps output ### d, doy, i =",d,doy,i,"\n")
        do_ovlp <- 1     # OK to calc ovlp
        do_maps <- 1     # OK to write maps
      } else if (doy %% 2 == 0 && evendaymaps == 1) { ## every other day
        cat("### Even day maps output ### d, doy, i =",d,doy,i,"\n")
        do_ovlp <- 1     # OK to calc ovlp
        do_maps <- 1     # OK to write maps
      } else if (dailymaps == 1) { ## must be every day
        cat("### Daily maps output ### d, doy, i =",d,doy,i,"\n")
        do_ovlp <- 1     # OK to calc ovlp
        do_maps <- 1     # OK to write maps
      } else { ## no need to do ovlp or write maps
        do_ovlp <- 0     # NO to calc ovlp
        do_maps <- 0     # NO to write maps
      }
      
      if (ovlp > 0.0 && do_ovlp == 1) {  # use overlap feature which just ramps a transition
        OWCurr_Overlap <- function(d,c){
          return(Cond(d <= c,1,
                      Cond(d <= (1-c),1,
                           0.5*(1+(1-d)/c))))
        }
        Curr_Overlap <- function(d,c){
          return(Cond(d <= c,0.5*(1+d/c),
                      Cond(d <= (1-c),1,
                           0.5*(1+(1-d)/c))))
        }
        Prev_Overlap <- function(d,c){
          return(Cond(d <= c,0.5*(1-d/c),0))
        }
        Next_Overlap <- function(d,c){
          return(Cond(d <= (1-c),0,
                      0.5*(1-(1-d)/c)))
        }
        #progress1 <- (DDaccum * LS1) >= larvaeDD 
        if (owstage == "OE") {
          #cat("date and doy,owstage, i =",d,doy,owstage,i,"\n")
          OWeggP <- template # proportion in OW egg stage
          OWeggDev <- (DDaccum * LSOW0)/OWeggDD
          OWeggP <- OWeggP + LSOW0 * OWCurr_Overlap(OWeggDev,ovlp)
          # no prev stage overlap for OW stage
          larvP <- larvP + LSOW0 * Next_Overlap(OWeggDev,ovlp)
          #setwd(output_dir)
          #writeRaster(OWeggP,paste("OWEggP_",d,sep=""), format="GTiff",overwrite=TRUE)
          #writeRaster(OWeggDev,paste("OWEggDev_",d,sep=""), format="GTiff",overwrite=TRUE)
          #setwd(prism_dir)
        } else if (owstage == "OL") {
          OWlarvP <- template # proportion in OW larval stage
          OWlarvDev <- (DDaccum * LSOW1)/OWlarvaeDD
          OWlarvP <- OWlarvP + LSOW1 * OWCurr_Overlap(OWlarvDev,ovlp)
          pupP <- pupP + LSOW1 * Next_Overlap(OWlarvDev,ovlp)
        } else if (owstage == "OP") {
          OWpupP <- template # proportion in OW pup stage
          OWpupDev <- (DDaccum * LSOW2)/OWpupDD
          OWpupP <- OWpupP + LSOW2 * OWCurr_Overlap(OWpupDev,ovlp)
          adultP <- adultP + LSOW2 * Next_Overlap(OWpupDev,ovlp)
        } else if (owstage == "OA") {
          OWadultP <- template # proportion in OW adult stage
          OWadultDev <- (DDaccum * LSOW3)/OWadultDD
          OWadultP <- OWadultP + LSOW3 * OWCurr_Overlap(OWadultDev,ovlp)
          eggP <- eggP + LSOW3 * Next_Overlap(OWadultDev,ovlp)
        }
        eggP <- template
        larvP <- template
        pupP <- template
        adultP <- template
        
        eggDev <- (DDaccum * LS0)/eggDD
        eggP <- eggP + LS0 * Curr_Overlap(eggDev,ovlp)
        adultP <- adultP + LS0 * Prev_Overlap(eggDev,ovlp)
        larvP <- larvP + LS0 * Next_Overlap(eggDev,ovlp)
        
        larvDev <- (DDaccum * LS1)/larvaeDD
        larvP <- larvP + LS1 * Curr_Overlap(larvDev,ovlp)
        eggP <- eggP + LS1 * Prev_Overlap(larvDev,ovlp)
        pupP <- pupP + LS1 * Next_Overlap(larvDev,ovlp)
        
        pupDev <- (DDaccum * LS2)/pupDD
        pupP <- pupP + LS2 * Curr_Overlap(pupDev,ovlp)
        larvP <- larvP + LS2 * Prev_Overlap(pupDev,ovlp)
        adultP <-  adultP + LS2 * Next_Overlap(pupDev,ovlp)
        
        adultDev <- (DDaccum * LS3)/adultDD
        adultP <- adultP + LS3 * Curr_Overlap(adultDev,ovlp)
        pupP <- pupP + LS3 * Prev_Overlap(adultDev,ovlp)
        eggP <- eggP + LS3 * Next_Overlap(adultDev,ovlp)
        
        # add OW stageP to regular (curr and next) stageP
        if (owstage == "OE") {
          eggP <- eggP + LSOW0 * OWeggP
          larvP <- larvP + LSOW0 * (1 - OWeggP)
        } else if (owstage == "OL") {
          larvP <- larvP + LSOW1 * OWlarvP
          pupP <- pupP + LSOW1 * (1 - OWlarvP)
        } else if (owstage == "OP") {
          pupP <- pupP + LSOW2 * OWpupP
          adultP <- adultP + LSOW2 * (1 - OWpupP)
        } else if (owstage == "OA") {
          adultP <- adultP + LSOW3 * OWadultP
          eggP <- eggP + LSOW3 * (1 - OWadultP)
        }
        
        # new Nov 2015 tally stages plus gens
        StageCount <- Lifestage + (NumGen * 4)
        
        #cat("today and showtoday,d,doy,i =",today,showtoday,d,doy,i,"\n")
        
        if (do_maps == 1) { ## today
          cat("### its today - write maps today ### d, doy, i =",d,doy,i,"\n")
          WriteMaps(d)
        } 
      }  # ovlp > 0 && do_ovlp == 1
      else if (do_maps == 1) { # no overlap but maybe write maps today
        cat("### its today - write maps today ### d, doy, i =",d,doy,i,"\n")
        WriteMaps(d)
      }
    }  # lifestage F or 5 (end of day calcs)
  }  # lifestage for loop
}  # daily loop

exectime <- toc()
exectime <- (exectime$toc - exectime$tic) / 60 
cat("Run time for daily loop =",exectime," minutes\n")  # *** 

#### Model Done - final map production ####
setwd(output_dir)
cat("Model done - Final map production\n")  # ***

SaveRaster(ddtotal,d,paste(spp,"_DD",sep=""))
Sys.sleep(0.3)
PlotMap(ddtotal,d,"Degree day (DD) accumulation","DD accum.",paste0(spp,"_ddtotal"))

# BOXB changes begin
SaveRaster(CUMDDs,d,paste(spp,"_PD-DDs",sep=""))
Sys.sleep(0.3)
PlotMap(CUMDDs,d,"Plant Disease DD accumulation","PD-DD accum.",paste0(spp,"_CUMDDs"))
SaveRaster(CUMDDLW,d,paste(spp,"INF_RISK_UNITS",sep=""))
Sys.sleep(0.3)
PlotMap(CUMDDLW,d,"Cum. infection risk","Cum. Infection Risk",paste0(spp,"_CUMDDLW"))
# BOXB changes end

SaveRaster(NumGen,d,paste(spp,"_NumGen",sep=""))
Sys.sleep(0.3)
PlotMap(NumGen,d,"Number of generations","No. generations",paste0(spp,"_NumGen"))

if (exclusions_stressunits) {
  # Exclusions included: -2=severe stress, -1=moderate stress
  NumGenEXCL1 <- Cond((AllEXCL > -2),NumGen,-2)
  NumGenEXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,NumGen))
  # BOXB changes begin
  CumInfRiskEXCL1 <- Cond((AllEXCL > -2),CUMDDLW,-2)
  CumInfRiskEXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,CUMDDLW))
  # Final stress unit and exclusion rasters
  rast_list <- list(chillunitsCUM,chillEXCL,heatunitsCUM,heatEXCL,AllEXCL,NumGenEXCL1,NumGenEXCL2)
  names(rast_list) <- c("Chill_Stress_Units","Chill_Stress_Excl","Heat_Stress_Units","Heat_Strs_Excl","All_Strs_Excl","NumGenExcl1","NumGenExcl2")
  invisible(mapply(SaveRaster, rast_list, d, paste(spp,names(rast_list),sep="_")))
  # Final stress unit and exclusion maps
  PlotMap_stress(chillunitsCUM,d,chillstress_units_max1,chillstress_units_max2,"Chill stress units","Chill Stress Units",paste(spp,"Chill_Stress_Units",sep="_"))
  PlotMap(chillEXCL,d,"Chill stress exclusion","Exclusion status",paste(spp,"Chill_Stress_Excl",sep="_"))
  PlotMap_stress(heatunitsCUM,d,heatstress_units_max1,heatstress_units_max2,"Heat stress units","Heat Stress Units",paste(spp,"Heat_Stress_Units",sep="_"))
  PlotMap(heatEXCL,d,"Heat stress exclusion","Exclusion status",paste(spp,"Heat_Stress_Excl",sep="_"))
  PlotMap(dryEXCL,d,"Dry stress exclusion","Exclusion status",paste(spp,"Dry_Stress_Excl",sep="_"))
  PlotMap(AllEXCL,d,"All climate stress exclusion","Exclusion status",paste(spp,"All_Stress_Excl",sep="_"))
  PlotMap(NumGenEXCL1,d,"No. of gens. w/ climate stress excl.", "No. generations",paste(spp,"NumGenExcl1",sep="_"))
  PlotMap(NumGenEXCL2,d,"No. of gens. w/ climate stress excl.", "No. generations",paste(spp,"NumGenExcl2",sep="_"))
  PlotMap(CumInfRiskEXCL1,d,"Cum. Infection Risk w/ sev. only clim. stress", "Cum. Infection Risk","CumInfRiskExcl1")
  PlotMap(CumInfRiskEXCL2,d,"Cum. Infection Risk w/ mod. & sev. climate stress", "Cum. Infection Risk","CumInfRiskExcl2")
  # BOXB changes end
}


# PEM rasters and maps
if(pems) {  # should PEM maps only at end of year; may not be correct if not run full year
  # Clean up labels for plotting
  #eggEventLabel <- tolower(gsub("_"," ",eggEventLabel)) # Label for PEM egg stage, lower case
  eggEventLabel <- eggEventLabel %>% tolower(.) %>% gsub("_", " ", .) %>% gsub("beginning","beg.",.) # Label for PEM egg stage, lower case
  #larvaeEventLabel <- tolower(gsub("_"," ",larvaeEventLabel)) # Label for PEM larval stage, lower case
  larvaeEventLabel <- larvaeEventLabel %>% tolower(.) %>% gsub("_", " ", .) %>% gsub("development","dev.",.) # Label for PEM larvae stage, lower case
  #pupaeEventLabel <- tolower(gsub("_"," ",pupaeEventLabel)) # Label for PEM pupae stage, lower case
  pupaeEventLabel <- pupaeEventLabel %>% tolower(.) %>%  gsub("_", " ", .) %>% gsub("development","dev.",.) # Label for PEM pupal stage, lower case
  #adultEventLabel <- tolower(gsub("_"," ",adultEventLabel)) # Label for PEM adult stage, lower case
  adultEventLabel <- adultEventLabel %>% tolower(.) %>%  gsub("_", " ", .) %>% gsub(" by females","",.) # Label for PEM adult stage, lower case
  setwd(output_dir)
  if (mapE == 1) {
    if (PEMnumgens > 0) {
      SaveRaster(PEMe1,d,paste(spp,"_PEMe1",sep=""))
      PlotMap(PEMe1,d,paste("Date of 1st gen.",eggEventLabel,sep=" "),paste("Date of 1st gen.",eggEventLabel,sep=" "),paste(spp,"PEMe1",sep="_"))
      if (exclusions_stressunits) { # 1st Gen eggstage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMe1EXCL1 <- Cond((AllEXCL > -2),PEMe1,-2)
        SaveRaster(PEMe1EXCL1,d,paste(spp,"_PEMe1Excl1",sep=""))
        PlotMap(PEMe1EXCL1,d,paste("Date of 1st gen.",eggEventLabel,"w/ sev. clim. excl.", sep=" "),paste("Date 1st gen.",eggEventLabel,sep=" "),paste(spp,"PEMe1Excl1",sep="_"))
        PEMe1EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMe1))
        SaveRaster(PEMe1EXCL2,d,paste(spp,"_PEMe1Excl2",sep=""))
        PlotMap(PEMe1EXCL2,d,paste("Date of 1st gen.",eggEventLabel,"w/ all clim. excl.", sep=" "),paste("Date 1st gen.",eggEventLabel,sep=" "),paste(spp,"PEMe1Excl2",sep="_"))
      }
    }
    if (PEMnumgens > 1) {
      SaveRaster(PEMe2,d,paste(spp,"_PEMe2",sep=""))
      PlotMap(PEMe2,d,paste("Date of 2nd gen.", eggEventLabel, sep=" "),paste("Date of 2nd gen.",eggEventLabel,sep=" "),paste(spp,"PEMe2",sep="_"))
      if (exclusions_stressunits) { # 2nd Gen eggstage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMe2EXCL1 <- Cond((AllEXCL > -2),PEMe2,-2)
        SaveRaster(PEMe2EXCL1,d,paste(spp,"_PEMe2Excl1",sep=""))
        PlotMap(PEMe2EXCL1,d,paste("Date of 2nd gen.",eggEventLabel,"w/ sev. clim. excl.", sep=" "),paste("Date 2nd gen.",eggEventLabel,sep=" "),paste(spp,"PEMe2Excl1",sep="_"))
        PEMe2EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMe2))
        SaveRaster(PEMe2EXCL2,d,paste(spp,"_PEMe2Excl2",sep=""))
        PlotMap(PEMe2EXCL2,d,paste("Date of 2nd gen.",eggEventLabel,"w/ all clim. excl.", sep=" "),paste("Date 2nd gen.",eggEventLabel,sep=" "),paste(spp,"PEMe2Excl2",sep="_"))
      }
    }
    if (PEMnumgens > 2) {
      SaveRaster(PEMe3,d,paste(spp,"_PEMe3",sep=""))
      PlotMap(PEMe3,d,paste("Date of 3rd gen.",eggEventLabel, sep=" "),paste("Date of 3rd gen.",eggEventLabel,sep=" "),paste(spp,"PEMe3",sep="_"))
      if (exclusions_stressunits) { # 3rd Gen eggstage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMe3EXCL1 <- Cond((AllEXCL > -2),PEMe3,-2)
        SaveRaster(PEMe3EXCL1,d,paste(spp,"_PEMe3Excl1",sep=""))
        PlotMap(PEMe3EXCL1,d,paste("Date of 3rd gen.",eggEventLabel, "w/ sev. clim. excl.", sep=" "),paste("Date 3rd gen.",eggEventLabel,sep=" "),paste(spp,"PEMe3Excl1",sep="_"))
        PEMe3EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMe3))
        SaveRaster(PEMe3EXCL2,d,paste(spp,"_PEMe3Excl2",sep=""))
        PlotMap(PEMe3EXCL2,d,paste("Date of 3rd gen.",eggEventLabel,"w/ all clim. excl.", sep=" "),paste("Date 3rd gen.",eggEventLabel,sep=" "),paste(spp,"PEMe3Excl2",sep="_"))
      }
    }
    if (PEMnumgens > 3) {
      SaveRaster(PEMe4,d,paste(spp,"_PEMe4",sep=""))
      PlotMap(PEMe4,d,paste("Date of 4th gen.", eggEventLabel, sep=" "),paste("Date of 4th gen.",eggEventLabel,sep=" "),paste(spp,"PEMe4",sep="_"))
      if (exclusions_stressunits) { # 4th Gen eggstage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMe4EXCL1 <- Cond((AllEXCL > -2),PEMe4,-2)
        SaveRaster(PEMe4EXCL1,d,paste(spp,"_PEMe4Excl1",sep=""))
        PlotMap(PEMe4EXCL1,d,paste("Date of 4th gen.",eggEventLabel,"w/ sev. clim. excl.", sep=" "),paste("Date 4th gen.",eggEventLabel,sep=" "),paste(spp,"PEMe4Excl1",sep="_"))
        PEMe4EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMe4))
        SaveRaster(PEMe4EXCL2,d,paste(spp,"_PEMe4Excl2",sep=""))
        PlotMap(PEMe4EXCL2,d,paste("Date of 4th gen.",eggEventLabel,"w/ all clim. excl.", sep=" "),paste("Date 4th gen.",eggEventLabel,sep=" "),paste(spp,"PEMe4Excl2",sep="_"))
      }
    }
  }
  if (mapL == 1) {
    if (PEMnumgens > 0) {
      SaveRaster(PEMl1,d,paste(spp,"_PEMl1",sep=""))
      PlotMap(PEMl1,d,paste("Date of 1st gen. ",larvaeEventLabel,sep=" "),paste("Date 1st gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl1",sep="_"))
      if (exclusions_stressunits) { # 1st Gen larval stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMl1EXCL1 <- Cond((AllEXCL > -2),PEMl1,-2)
        SaveRaster(PEMl1EXCL1,d,paste(spp,"_PEMl1Excl1",sep=""))
        PlotMap(PEMl1EXCL1,d,paste("Date of 1st gen.",larvaeEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 1st gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl1Excl1",sep="_"))
        PEMl1EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMl1))
        SaveRaster(PEMl1EXCL2,d,paste(spp,"_PEMl1Excl2",sep=""))
        PlotMap(PEMl1EXCL2,d,paste("Date of 1st gen.",larvaeEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 1st gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl1Excl2",sep="_"))
      }
    }
    if (PEMnumgens > 1) {
      SaveRaster(PEMl2,d,paste(spp,"_PEMl2",sep=""))
      PlotMap(PEMl2,d,paste("Date of 2nd gen.", larvaeEventLabel,sep=" "),paste("Date 2nd gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl2",sep="_"))
      if (exclusions_stressunits) { # 2nd Gen larval stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMl2EXCL1 <- Cond((AllEXCL > -2),PEMl2,-2)
        SaveRaster(PEMl2EXCL1,d,paste(spp,"_PEMl2Excl1",sep=""))
        PlotMap(PEMl2EXCL1,d,paste("Date of 2nd gen.",larvaeEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 2nd gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl2Excl1",sep="_"))
        PEMl2EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMl2))
        SaveRaster(PEMl2EXCL2,d,paste(spp,"_PEMl2Excl2",sep=""))
        PlotMap(PEMl2EXCL2,d,paste("Date of 2nd gen.",larvaeEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 2nd gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl2Excl2",sep="_"))
      }
    }
    if (PEMnumgens > 2) {
      SaveRaster(PEMl3,d,paste(spp,"_PEMl3",sep=""))
      PlotMap(PEMl3,d,paste("Date of 3rd gen.", larvaeEventLabel,sep=" "),paste("Date 3rd gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl3",sep="_"))
      if (exclusions_stressunits) { # 3rd Gen larval stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMl3EXCL1 <- Cond((AllEXCL > -2),PEMl3,-2)
        SaveRaster(PEMl3EXCL1,d,paste(spp,"_PEMl3Excl1",sep=""))
        PlotMap(PEMl3EXCL1,d,paste("Date of 3rd gen.",larvaeEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 3rd gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl3Excl1",sep="_"))
        PEMl3EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMl3))
        SaveRaster(PEMl3EXCL2,d,paste(spp,"_PEMl3Excl2",sep=""))
        PlotMap(PEMl3EXCL2,d,paste("Date of 3rd gen.",larvaeEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 3rd gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl3Excl2",sep="_"))
      }
    }
    if (PEMnumgens > 3) {
      SaveRaster(PEMl4,d,paste(spp,"_PEMl4",sep=""))
      PlotMap(PEMl4,d,paste("Date of 4th gen.", larvaeEventLabel,sep=" "),paste("Date 4th gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl4",sep="_"))
      if (exclusions_stressunits) { # 4th Gen larval stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMl4EXCL1 <- Cond((AllEXCL > -2),PEMl4,-2)
        SaveRaster(PEMl4EXCL1,d,paste(spp,"_PEMl4Excl1",sep=""))
        PlotMap(PEMl4EXCL1,d,paste("Date of 4th gen.",larvaeEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 4th gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl4Excl1",sep="_"))
        PEMl4EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMl4))
        SaveRaster(PEMl4EXCL2,d,paste(spp,"_PEMl4Excl2",sep=""))
        PlotMap(PEMl4EXCL2,d,paste("Date of 4th gen.",larvaeEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 4th gen.",larvaeEventLabel,sep=" "),paste(spp,"PEMl4Excl2",sep="_"))
      }
    }
  }
  if (mapP == 1) {
    if (PEMnumgens > 0) {
      SaveRaster(PEMp1,d,paste(spp,"_PEMp1",sep=""))
      PlotMap(PEMp1,d,paste("Date of 1st gen.",pupaeEventLabel,sep=" "),paste("Date 1st gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp1",sep="_"))
      if (exclusions_stressunits) { # 1st Gen pupal stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMp1EXCL1 <- Cond((AllEXCL > -2),PEMp1,-2)
        SaveRaster(PEMp1EXCL1,d,paste(spp,"_PEMp1Excl1",sep=""))
        PlotMap(PEMp1EXCL1,d,paste("Date of 1st gen.",pupaeEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 1st gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp1Excl1",sep="_"))
        PEMp1EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMp1))
        SaveRaster(PEMp1EXCL2,d,paste(spp,"_PEMp1Excl2",sep=""))
        PlotMap(PEMp1EXCL2,d,paste("Date of 1st gen.",pupaeEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 1st gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp1Excl2",sep="_"))
      }
    }
    if (PEMnumgens > 1) {
      SaveRaster(PEMp2,d,paste(spp,"_PEMp2",sep=""))
      PlotMap(PEMp2,d,paste("Date of 2nd gen.",pupaeEventLabel,sep=" "),paste("Date 2nd gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp2",sep="_"))
      if (exclusions_stressunits) { # 2nd Gen pupal stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMp2EXCL1 <- Cond((AllEXCL > -2),PEMp2,-2)
        SaveRaster(PEMp2EXCL1,d,paste(spp,"_PEMp2Excl1",sep=""))
        PlotMap(PEMp2EXCL1,d,paste("Date of 2nd gen.",pupaeEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 2nd gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp2Excl1",sep="_"))
        PEMp2EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMp2))
        SaveRaster(PEMp2EXCL2,d,paste(spp,"_PEMp2Excl2",sep=""))
        PlotMap(PEMp2EXCL2,d,paste("Date of 2nd gen.",pupaeEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 2nd gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp2Excl2",sep="_")) 
      }
    }
    if (PEMnumgens > 2) {
      SaveRaster(PEMp3,d,paste(spp,"_PEMp3",sep=""))
      PlotMap(PEMp3,d,paste("Date of 3rd gen. ",pupaeEventLabel,sep=" "),paste("Date 3rd gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp3",sep="_"))
      if (exclusions_stressunits) { # 3rd Gen pupal stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMp3EXCL1 <- Cond((AllEXCL > -2),PEMp3,-2)
        SaveRaster(PEMp3EXCL1,d,paste(spp,"_PEMp3Excl1",sep=""))
        PlotMap(PEMp3EXCL1,d,paste("Date of 3rd gen.",pupaeEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 3rd gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp3Excl1",sep="_"))
        PEMp3EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMp3))
        SaveRaster(PEMp3EXCL2,d,paste(spp,"_PEMp3Excl2",sep=""))
        PlotMap(PEMp3EXCL2,d,paste("Date of 3rd gen.",pupaeEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 3rd gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp3Excl2",sep="_"))		
      }
    }
    if (PEMnumgens > 3) {
      SaveRaster(PEMp4,d,paste(spp,"_PEMp4",sep=""))	    
      PlotMap(PEMp4,d,paste("Date of 4th gen.",pupaeEventLabel,sep=" "),paste("Date 4th gen.",pupaeEventLabel),paste(spp,"PEMp4",sep="_"))
      if (exclusions_stressunits) { # 4th Gen pupal stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMp4EXCL1 <- Cond((AllEXCL > -2),PEMp4,-2)
        SaveRaster(PEMp4EXCL1,d,paste(spp,"_PEMp4Excl1",sep=""))
        PlotMap(PEMp4EXCL1,d,paste("Date of 4th gen.",pupaeEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 4th gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp4Excl1",sep="_"))
        PEMp4EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMp4))
        SaveRaster(PEMp4EXCL2,d,paste(spp,"_PEMp4Excl2",sep=""))
        PlotMap(PEMp4EXCL2,d,paste("Date of 4th gen.",pupaeEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 4th gen.",pupaeEventLabel,sep=" "),paste(spp,"PEMp4Excl2",sep="_"))	
      }
    }
  }
  if (mapA == 1) {
    if (PEMnumgens > 0) {
      SaveRaster(PEMa1,d,paste(spp,"_PEMa1",sep=""))	
      PlotMap(PEMa1,d,paste("Date of 1st gen.",adultEventLabel,sep=" "),paste("Date 1st gen.",adultEventLabel),paste(spp,"PEMa1",sep="_"))
      if (exclusions_stressunits) { # 1st Gen adult stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMa1EXCL1 <- Cond((AllEXCL > -2),PEMa1,-2)
        SaveRaster(PEMa1EXCL1,d,paste(spp,"_PEMa1Excl1",sep=""))
        PlotMap(PEMa1EXCL1,d,paste("Date of 1st gen.",adultEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 1st gen.",adultEventLabel,sep=" "),paste(spp,"PEMa1Excl1",sep="_"))
        PEMa1EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMa1))
        SaveRaster(PEMa1EXCL2,d,paste(spp,"_PEMa1Excl2",sep=""))
        PlotMap(PEMa1EXCL2,d,paste("Date of 1st gen.",adultEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 1st gen.",adultEventLabel,sep=" "),paste(spp,"PEMa1Excl2",sep="_"))
      }
    }
    if (PEMnumgens > 1) {
      SaveRaster(PEMa2,d,paste(spp,"_PEMa2",sep=""))	
      PlotMap(PEMa2,d,paste("Date of 2nd gen.",adultEventLabel,sep=" "),paste("Date 2nd gen.",adultEventLabel),paste(spp,"PEMa2",sep="_"))
      if (exclusions_stressunits) { # 2nd Gen adult stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMa2EXCL1 <- Cond((AllEXCL > -2),PEMa2,-2)
        SaveRaster(PEMa2EXCL1,d,paste(spp,"_PEMa2Excl1",sep=""))
        PlotMap(PEMa2EXCL1,d,paste("Date of 2nd gen.",adultEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 2nd gen.",adultEventLabel,sep=" "),paste(spp,"PEMa2Excl1",sep="_"))
        PEMa2EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMa2))
        SaveRaster(PEMa2EXCL2,d,paste(spp,"_PEMa2Excl2",sep=""))
        PlotMap(PEMa2EXCL2,d,paste("Date of 2nd gen.",adultEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 2nd gen.",adultEventLabel,sep=" "),paste(spp,"PEMa2Excl2",sep="_"))
      }
    }
    if (PEMnumgens > 2) {
      SaveRaster(PEMa3,d,paste(spp,"_PEMa3",sep=""))	
      PlotMap(PEMa3,d,paste("Date of 3rd gen.",adultEventLabel,sep=" "),paste("Date 3rd gen.",adultEventLabel),paste(spp,"PEMa3",sep="_"))
      if (exclusions_stressunits) { # 3rd Gen adult stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMa3EXCL1 <- Cond((AllEXCL > -2),PEMa3,-2)
        SaveRaster(PEMa3EXCL1,d,paste(spp,"PEMa3Excl1",sep=""))
        PlotMap(PEMa3EXCL1,d,paste("Date of 3rd gen.",adultEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 3rd gen.",adultEventLabel,sep=" "),paste(spp,"PEMa3Excl1",sep="_"))
        PEMa3EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMa3))
        SaveRaster(PEMa3EXCL2,d,paste(spp,"PEMa3Excl2",sep=""))
        PlotMap(PEMa3EXCL2,d,paste("Date of 3rd gen.",adultEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 3rd gen.",adultEventLabel,sep=" "),paste(spp,"PEMa3Excl2",sep="_"))
      }
    }
    if (PEMnumgens > 3) {
      SaveRaster(PEMa4,d,paste(spp,"_PEMa4",sep=""))	
      PlotMap(PEMa4,d,paste("Date of 4th gen.",adultEventLabel,sep=" "),paste("Date 4th gen.",adultEventLabel),paste(spp,"PEMa4",sep="_"))
      if (exclusions_stressunits) { # 4th Gen adult stage includ exclusions
        # Exclusions included: -2=sev., -1=moderate stress
        PEMa4EXCL1 <- Cond((AllEXCL > -2),PEMa4,-2)
        SaveRaster(PEMa4EXCL1,d,paste(spp,"PEMa4Excl1",sep=""))
        PlotMap(PEMa4EXCL1,d,paste("Date of 4th gen.",adultEventLabel,"w/ sev. clim. excl.",sep=" "),paste("Date 4th gen.",adultEventLabel,sep=" "),paste(spp,"PEMa4Excl1",sep="_"))
        PEMa4EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMa4))
        SaveRaster(PEMa4EXCL2,d,paste(spp,"PEMa4Excl2",sep=""))
        PlotMap(PEMa4EXCL2,d,paste("Date of 4th gen.",adultEventLabel,"w/ all clim. excl.",sep=" "),paste("Date 4th gen.",adultEventLabel,sep=" "),paste(spp,"PEMa4Excl2",sep="_"))
      }
    }
  }
} #### if PEMS
setwd(prism_dir)

warnings()
#Possibility to calculate any number of outputs. This example was for 2014
#data only, but will want to look at multi-year calculations and how we
#can express uncertainty (annual variability) for more static risk maps.

