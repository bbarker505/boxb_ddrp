#!/usr/bin/Rscript
.libPaths("/usr/local/lib64/R/library/")
options(echo=FALSE)
## phen_model.R
# this version for ENTO adapt for CGI
args <- commandArgs(trailingOnly = TRUE)

# Install pakcages if necessary and load them
pkgs <- c("sp", "rgdal", "raster", "lubridate","plyr","dplyr", "stringr", "ggplot2",
          "ggthemes","mapdata","RColorBrewer","tictoc","tidyr")
# new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])] # check to see if any are missing
# if(length(new_pkgs)) install.packages(new_pkgs, lib="~/R/x86_64-pc-linux-gnu-library/3.5",repos = "https://cloud.r-project.org") # install missing ones
ld_pkgs <- lapply(pkgs, library, character.only = TRUE) # load them

# Bring in states feature for summary maps (PNG files), extract the lower 48 states, and project it
# Requires these libraries: "mapdata" and "maptools"
cat("Downloading US states feature\n")
states <- map_data("state")

##################################################################################################
########                  BEGINNING of Header Documentation Section             ##################
########            DDRP: Degree-Day, establishment Risk, and Pest event mapping system ##########
########  By Len Coop, Gericke Cook, Dan Upper, and Brittany Barker for APHIS PPQ and IPM needs ##
##################################################################################################
 # alternative names to consider:
 #	  PestEventAndRiskMappingSystem PEARMS                                                
 #	  RiskDDAndPestEventMappingSystem RADDAPEMS                                                
 #	  DDRiskPOPS: 
 #	  LifeStageMaps
 #   LBPERMS R-based Lifestage-Based Phenology Events and Risks Mapping System
 #
 # Features in program provided by Gericke/PPQ (208 lines code):
 # 1. DD model plus exclusions (climatic limits) for combined risk maps predicting both potential
 #    number of generations and failure to establish
 # 2. Allow for multi-stage developmental and exclusion params (e.g. separate thresholds and limits
 #    for eggs larvae pupae adults)
 # 3. Step through daily PRISM 4km climate maps (Tmax,Tmin,Tmean,Precip) for calculations
 #
 # Features added/pending by Len/Dan:
 # 1. WORKING:v5: Additional DD calc formulas allowing better precision and use of upper thresholds (UDT)
 # 2. WORKING:v7: Addition of pest event map (PEM) creation
 # 3. WORKING v8: Option to switch off exclusion and PEM processing
 # 4. WORKING:v15: Addition of stage overlap factor e.g. eggs do not all hatch to larvae at same time; with
 #    use of an overlap parameter and ramp functions
 # 5. WORKING:v16: Define a subregion, set up numerous subregions for various needs
 # 6. WORKING:v18: Add OW stage initialization and flexible stage ordering e.g.1: OW adult->e,l,p,a,repeat
 #    e.g.2: OW egg->l,p,a,e,repeat
 # 7. WORKING: v20: 1) Simple parameter passing suitable for command line and cgi-bin scripts (658 lines code):
 #                  2) Stage overlap param (0.0 to 0.50) now working 
 #      Params passed in: 
 #               SPP (3letterabbrev) 
 #               YEAR (4 digit) 
 #               START_DOY (dayofyear) 
 #               REGION (CONUS,OR,MIDWEST currently)
 #               STAGE_OVERLAP (0-0.5) 
 #               EXCLUSIONS (Y/N)
 #               PEMS (Y/N)
 #      ex: ./DDRiskPEMSv21.R FCM 2010 39 OR 0.35 N N
 #   showing correct DD calcs, NumGens, PEMs for 3 stages and 2 generations as examples
 # 8. WORKING v22: 1) Read pest/insect params from external library/database file (736 lines code)
 #                 2) Add data output dir different from data input dir
 #            ex: ./DDRiskPEMSv22.R FCM 2010 39 OR 0.35 N N results3
 # 9. WORKING: v23 file handling (only need tmax, tmin, ppt): 
 #              1) a) past years - use approp year daily 4km PRISM data (stable) - yes 
 #                    (800m daily PRISM is available to us on word of PRISM group but our servers probably not
 #                    ready to handle such high resolution, hopefully new server will be)
 #					   b) for current year (2015) use PRISM in order of avail: 
 #                   b1) stable b2) provisional b3) early  
 #                c) forecast (NDFD 2.5km is avail to try), need to work on naming conventions and resolution
 #		            d) extended forecast starting w/10-year averages (800m avail), need naming conventions and resolution
 #                e)   by 2016 should have 90+ day forecasts via ARDP grant work
 #                f) FUTURE PRISM data: possibly Dewpoint, RHmean/VPD or other moisture params
 #             3) Control over outputs:
 #                a) PEMS should be automatically on if PEMS is on
 #                b) More control over stages and times of PEMS w/species param files - working
 #                c) need params passed in for freq of outputs: daily, weekly, 10 days, 14 days, monthly, just the
 #                    final maps.
 #                    normally would set a last date and only output a map on that date
 #                    also a way to specify stage to focus on?
 #                    try starting with categories:
 #                    0=no extra output (defaults in spp file?)
 #                    1=outputs last day of month
 #                d) Control over which stages to map
 #                e) Consider degrees F on input (harder?) as well as on output (easier)
 #             4) Web/CGI-BIN interface         
 #                 a) all command line params available as form elements
 #                 b) have a download link for all files
 #                 c) make .png files for viewing that includes legend, title, statelines
 #                 d) use PID or random number (web) vs user specified output directory
 #                 e) full set of CONUS states and regions
 #                 f) have full metadata files available
 # 10. IN DEVEL: v24 expand exclusions to allow:
 #                 a) accumulate chill units 
 #                 b) separate tracking for end-of-year masking in addition to realtime exclusions as in previous versions
 # 11v. This is for development of a plant disease risk version of DDRP so use prism data such as:
 #   A) Add Drystress looking similar to Heat and Cold (Chill) stress, using relhum=func(tmean,tdmean)
 #   B) Add DDLW unique to Boxwood blight since BoxDD is a piecewise regression formula based on BB DD lookup table used in site vers.
 #"PRISM_tddif_early_4kmD2_MTD_20210401=(PRISM_tmean_early_4kmD2_MTD_20210401-PRISM_tdmean_early_4kmD2_MTD_20210401)"
 #"LW_20210401 = if(PRISM_ppt_early_4kmD2_MTD_20210401 > 8,1,if(PRISM_tddef_early_4kmD2_MTD_20210401 < 5,2,0))"
 # the thresholds values can be calibrated and included as input params from the spp_param files
 #
########                END of Header Documentation Section                        #########

########                BEGINNING of Param Handling Section                        #########
#### Default values for command line params ####
spp           <- "STB"       # default species to use
start_year    <- "2018"        # only one 4 digit year currently supported 
start_doy     <- 1          # day of year to use for biofix             
end_doy       <- 365         # day of year to use for stop date - need 365 if voltinism map 
region_param  <- "SOUTHEAST"   # default REGION to use
ovlp          <- 0.40        # proportion overlap stage to stage; must be < 0.5
exclusions    <- 3          # ability to turn on/off climate limits/exclusions 1=on 2=also apply exclusions in realtime 3 = stress unit exclusions
pems          <- "Y"        # turn on/off pest event maps
out_dir       <- "ALB888_013019"   # output directory (currently appended to /data/PRISM)
out_option    <- 0           # output option category (**to be defined**)
mapA          <- 1           # make maps for adult stage
mapE          <- 0           # make maps for egg stage
mapL          <- 0           # make maps for larval stage
mapP          <- 0           # make maps for pupal stage
do_ovlp       <- 0           # added 6/1/16 to only calc ovlp if we make maps on that day
do_maps       <- 0           # added 6/1/16 to only make maps on days selected
                             # could add: ovlp differ by stage, incr each generation
#### Process command line args ####
if (length(args) < 12) {
   cat("  ","\n")
   cat("DDRP: Degree-Day, establishment Risk, and Pest Event Mapping System","\n")
   cat("Required params: SPP (3letterabbrev)       ","\n") 
   cat("                 YEAR (4 digit)            ","\n") 
   cat("                 START_DOY (dayofyear)     ","\n") 
   cat("                 END_DOY (dayofyear)       ","\n") 
   cat("                 REGION (CONUS,regions,48 states currently)","\n")
   cat("                 OVERLAP (0-0.5)           ","\n") 
   cat("                 EXCLUSIONS (0=off,1=on,2=apply in realtime,3=stressunit-based exclusions)","\n")
   cat("                 Pest Event Maps (Y/N)     ","\n")
   cat("                 Name of output dir (results)","\n")
   cat("                 Output option 1-8 (1=default)","\n")
   cat("                 Make Maps for Adults (1=yes)","\n")
   cat("                 Make Maps for Eggs   (1=yes)","\n")
   cat("                 Make Maps for Larvae (1=yes)","\n")
   cat("                 Make Maps for Pupae  (1=yes)","\n")
   cat("Run example:     DDRPv24.R FCM 2010 39 365 MIDWEST 0.25 0 N results 1 1 0 0 0","\n")
   cat("  ","\n")
   q()
} else {
	cat("PARAM LIST:  ")
   print(args)
}

#### Read in command line args ####
spp           <- args[1] 
start_year    <- as.character(args[2])
start_doy     <- as.integer(args[3])
end_doy       <- as.integer(args[4])
region_param  <- args[5]
ovlp          <- as.double(args[6])
exclusions    <- args[7]
pems          <- args[8]
out_dir       <- args[9]
out_option    <- args[10]
mapA          <- args[11]
mapE          <- args[12]
mapL          <- args[13]
mapP          <- args[14]

#### Use gsub to make secure cgi-bin args : equiv. of perl spp =~ s/[^A-Za-z0-9]//g; ####
spp           <- gsub('[^A-Za-z0-9]',"",spp,perl=TRUE)
start_year    <- gsub('[^A-Za-z0-9._]',"",start_year,perl=TRUE)
start_doy     <- gsub('[^0-9.]',"",start_doy,perl=TRUE)
end_doy       <- gsub('[^0-9.]',"",end_doy,perl=TRUE)
region_param  <- gsub('[^A-Za-z0-9."]',"",region_param,perl=TRUE)
ovlp          <- gsub('[^0-9.]',"",ovlp,perl=TRUE)
#exclusions    <- gsub('[^A-Z]',"",exclusions,perl=TRUE)
exclusions    <- gsub('[^0-9]',"",exclusions,perl=TRUE)
pems          <- gsub('[^A-Z]',"",pems,perl=TRUE)
out_dir       <- gsub('[^A-Za-z0-9./_]',"",out_dir,perl=TRUE)
out_option    <- gsub('[^0-9.]',"",out_option,perl=TRUE)
mapA          <- gsub('[^0-9.]',"",mapA,perl=TRUE)
mapE          <- gsub('[^0-9.]',"",mapE,perl=TRUE)
mapL          <- gsub('[^0-9.]',"",mapL,perl=TRUE)
mapP          <- gsub('[^0-9.]',"",mapP,perl=TRUE)

#### Convert var types for some vars ####
#start_year    <- as.integer(start_year)
start_doy     <- as.integer(start_doy)
end_doy       <- as.integer(end_doy)
ovlp          <- as.double(ovlp)
out_option    <- as.integer(out_option)
mapA          <- as.integer(mapA)
mapE          <- as.integer(mapE)
mapL          <- as.integer(mapL)
mapP          <- as.integer(mapP)

#### Set up regions - use switch() (works like a single use hash) ####
#"MIDWEST"      = extent(-92,-90,30,49),
REGION <- switch(region_param,
"CONUS"        = extent(-125.0,-66.5,24.0,50.0),
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
	
#if (exclusions == "Y") { 
#	exclusions <- 1 
#} else {
#	exclusions <- 0 
#}

if (exclusions == "2") { # OK both use climate limiting exclusions and apply realtime masks affecting all output 
	exclusions <- 1 
	exclusionsrealtime <- 1 
   exclusions_stressunits <- 0
} else if (exclusions == "1") { # use exclusions but dont affect phenology/other maps until we apply them at end of model run 
	exclusions <- 1 
	exclusionsrealtime <- 0 
   exclusions_stressunits <- 0
} else if (exclusions == "3") {  #  == "3" NEW type stress unit based exclusions
	exclusions <- 0 
	exclusionsrealtime <- 0 
   exclusions_stressunits <- 1
} else {  #  == "0" dont use any exclusions
	exclusions <- 0 
	exclusionsrealtime <- 0 
   exclusions_stressunits <- 0
}


if (pems == "Y") { 
	pems <- 1 
} else { 
	pems <- 0 
}

#### Init. mapping output freq. options ####
finalmaps    <- 0     # minimal output only when last day of model run is done
monthlymaps  <- 0
biweeklymaps <- 0
dekadmaps    <- 0    # a dekad is every 10 days FYI e.g. Aug 1, 11, 21, 31
weeklymaps   <- 0
evendaymaps  <- 0
dailymaps    <- 0

if (out_option == "0") { 
  finalmaps = 1
} else if (out_option == "1") { 
  monthlymaps  <- 1
  finalmaps    <- 1
} else if (out_option == "2") { 
  biweeklymaps <- 1
  finalmaps    <- 1
} else if (out_option == "3") { 
  dekadmaps    <- 1
  finalmaps    <- 1
} else if (out_option == "4") { 
  weeklymaps   <- 1
  finalmaps    <- 1
} else if (out_option == "5") { 
  evendaymaps  <- 1
  finalmaps    <- 1
} else if (out_option == "6") { 
  dailymaps    <- 1
  finalmaps    <- 1
}
cat("CLEAN PARAMS: spp startyr startdoy enddoy region_param ovlp exclusions pems out_dir output_option mapA mapE mapL mapP: \n              ",spp,start_year,start_doy,end_doy,region_param,ovlp,exclusions,pems,out_dir,out_option,mapA,mapE,mapL,mapP,"\n")

##########################################################################################
########                  END of Param Handling Section                         ##########
########                BEGINNING of Directory Init Section                     ##########
########           PARAM INPUTS - location for species params; thresholds etc.  ##########
##########################################################################################
params_dir <- "/usr/local/dds/DDRP/spp_params/"

#### Weather inputs and outputs - PRISM climate data w/subdirs 4-digit year ####
if (grepl("16",start_year,perl=TRUE)) {  # if outdir has 2 consec. numbers, assume webuser
  base_dir <- "/mnt/ssd1/PRISM/"
} else {  # otherwise just use base dir
  base_dir <- "/data/PRISM/"
}
cat("BASE DIR: ",base_dir,"\n")
prism_dir <- sprintf("%s%s",base_dir,start_year)
cat("WORKING DIR: ",prism_dir,"\n")

#output_dir <- sprintf("%s%s","/usr/local/dds/DDRP_B1/DDRP_results/",out_dir)
output_dir <- sprintf("%s%s","/home/httpd/html/CAPS/",out_dir)

# Map outputs
if (grepl("[0-9]{3,5}",out_dir,perl=TRUE)) {  # if outdir has 3 to 5 consec. numbers, assume webuser
   #output_dir <- sprintf("%s%s","/usr/local/dds/DDRP_B1/DDRP_results/",out_dir)
   output_dir <- sprintf("%s%s","/home/httpd/html/CAPS/",out_dir)
} else {  # otherwise just use base dir
   output_dir <- sprintf("%s%s",base_dir,out_dir)
}

if (file.exists(output_dir)) {
  cat("EXISTING OUTPUT DIR: ",output_dir,"\n")
} else {
  dir.create(output_dir)
  cat("NEW OUTPUT DIR: ",output_dir,"\n")
}

# Log file - placed in out_dir
sink(paste0("/home/httpd/html/CAPS/",out_dir,"/rlogging.txt"))
#sink(paste0("/usr/local/dds/DDRP_B1/DDRP_results/",out_dir,"/rlogging.txt"))

# Capture messages and errors to a file
msg <- file(paste0("/home/httpd/html/CAPS/",out_dir,"/rmessages.txt"), open="wt")
#msg <- file(paste0("/usr/local/dds/DDRP_B1/DDRP_results/",out_dir,"/rmessages.txt"), open="wt")
sink(msg, type="message")

setwd(prism_dir)

##########################################################################################
########                  END of Directory Init Section                         ##########
########                BEGINNING of function definitions                       ##########
##########################################################################################
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
calcLW=function(ppt,relhum){
return(Cond(ppt > ppt_T3,4,Cond(relhum > RH_T3,3,Cond(ppt > ppt_T2 & relhum > RH_T2,2,Cond(ppt > ppt_T1 & relhum > RH_T1,1,0)))))
}

calcDDs=function(tmax,tmin){
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
# BOXB changes end

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
Cut_bins <- function(f) {
  f$value_orig <- f$value # keep old value so can sort factors against it
  f$value <- ifelse(f$value < 0.1, round_any(f$value, 0.01),  # round to hundredth dec
             ifelse(f$value < 1 & f$value > 0.1, round_any(f$value, .1), # round to tenth
             round_any(f$value, 1)))
  #f$value <- ifelse(f$value < 1, round_any(f$value, .1), round_any(f$value, 1)) # round decimal places appropriately
  f2 <- f %>% mutate(value = cut(value, breaks = 10, dig.lab=4)) # cut values into values
  f2$value <- gsub(pattern = "[()]|\\[|\\]", replacement = "", f2$value) # remove brackets and parentheses from output value
  f2$value <- gsub(pattern = "[,]", replacement = "-", f2$value) # replace comma with a dash
  f2$value <- factor(f2$value, levels = unique(f2$value[order(f2$value_orig)])) # Order by values so plots in numerical order
  return(f2)
}

# Plot theme to be used for all summary (PNG) maps
mytheme <- theme(plot.title=element_text(size=9, hjust=0.5, vjust=-3),  # Specify title attributes
            plot.subtitle = element_text(size=6.5, hjust=0.5, vjust=-3, lineheight = 1.25),
            plot.margin = unit(c(0,0.5,0,0.5),"line"),
            legend.text=element_text(size = 6), # Specify legend attributes
            legend.position="right", 
            legend.justification="left", 
            legend.margin=margin(0,0,0,0),
            legend.key.width = unit(0.75,"line"), legend.key.height = unit(0.4,"line"))

# Base features used for all summary (PNG) maps in "PlotMap" function
Base_map <- function(df) {
  p <- ggplot(states, aes(x=long, y=lat)) + 
        geom_raster(data=df, aes(x=x, y=y, fill=value)) + # plot raster (raster values are now in a data frame)
        geom_path(aes(group=group), color="gray20", lwd = 0.25) +
        theme_map(base_size = 7)
}

# Creates a gradient color ramp consisting of n colors
Colfunc <- function(x,y,n) {
  colorRampPalette(c(x,y))(n)
}

# Deal with 0 vs non-zero values when plotting chill and heat stress unit rasters
Stress_Val_Conv <- function(x) {
  if(all(x$value == 0)){
    x2 <- x
  } else {
    x2 <- Cut_bins(x)
  }
  x2$value <- factor(x2$value)
  return(x2)
}

# Create and write summary maps (PNG) of all rasters except for stress units
PlotMap <- function(r,d,titl,lgd,outfl) {
  df <- ConvDF(r) # convert raster to data frame
  sp <- paste0(gsub(pattern = "_", replacement = " ", fullname),":")
  nam <- deparse(substitute(r)) # Get name of raster, capitalize first letter
  dat <- as.character(format(strptime(d,format="%Y%m%d"), format="%m/%d/%Y")) # format the date
  titl <- paste(titl,dat,sep=" ")
  subtitl <- paste("Maps and modeling",format(Sys.Date(),"%m/%d/%Y"),"by Oregon State University OIPMC USPEST.ORG \nand USDA-APHIS-PPQ; climate data from OSU PRISM Climate Group")


  # lifestage and lifestage w/ climate exclusion
  if ((grepl("Lifestage", nam))) {
    # recode OW in StageCount to -3 (currently -1) so there's no overlap w/ clim. exclusion values (-2 and -1)  
    if (!grepl("EXCL",nam)) { # if it's simply a lifestage raster, then recode -1 to -3
      df$value <- gsub("-1","-3", df$value)
    }
    df <- mutate(df, value = ifelse(value==-3,"OW",ifelse(value==-2,"excl.-severe",ifelse(value==-1,"excl.-moderate",ifelse(value==0,"egg",ifelse(value==1,"larva",ifelse(value==2,"pupa","adult")))))))
    df$value <- factor(df$value, levels = c("excl.-severe","excl.-moderate","OW","egg","larva", "pupa", "adult"))   
    col_key <- cbind(data.frame("cols" = c("gray40","gray70", "gray70","#E41A1C","#377EB8","#4DAF4A","#984EA3"), stringsAsFactors = FALSE),
                      data.frame("value" = c("excl.-severe","excl.-moderate","OW","egg","larva", "pupa", "adult")))
    col_df <- suppressWarnings(suppressMessages(semi_join(col_key,df,by="value"))) # join the colors to the data
    cols <- setNames(as.character(col_df$cols),col_df$value) # convert col_df to a named vector
    # make the plot
    p <- Base_map(df) +
          scale_fill_manual(values=cols, name = paste0(lgd)) +
          labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
          mytheme
  # number of generations
  } else if (nam == "NumGen") {
      df$value <- factor(df$value) 
      # make the color ramp - enough for 15 generations
      col_key <- cbind(data.frame("cols" = c(brewer.pal(8,"Set1"),brewer.pal(8,"Pastel2")), stringsAsFactors = FALSE),
                       data.frame("value" = c(as.character(0:15)), stringsAsFactors = FALSE))
      col_df <- suppressWarnings(suppressMessages(semi_join(col_key,df,by="value"))) # join the colors to the data
      cols <- setNames(as.character(col_df$cols),col_df$value) # convert col_df to a named vector
      # make the plot
      p <- Base_map(df) +       
            scale_fill_manual(values=cols, name = paste0(lgd)) +
            labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
            mytheme
  # number of generations with climate stress exclusion
  } else if (grepl("NumGenEXCL", nam)) {
      df$value_orig <- df$value # create for ordering recoded values
      df <- mutate(df, value = factor(ifelse(value==-2,"excl.-severe",ifelse(value==-1,"excl.-moderate",value))))    
      df$value <- factor(df$value, levels = unique(df$value[order(df$value_orig)])) # order by original values so plots in numerical order
      # make the color ramp (excl.-severe, excl.-moderate, are gray shades) - enough for 15 generations
      col_key <- cbind(data.frame("cols" = c("gray40","gray70",brewer.pal(8,"Set1"),brewer.pal(8,"Pastel2")), stringsAsFactors = FALSE),
                       data.frame("value" = c("excl.-severe","excl.-moderate",as.character(0:15)), stringsAsFactors = FALSE))
      col_df <- suppressWarnings(suppressMessages(semi_join(col_key,df,by="value"))) # keep only the number of colors needed for the data
      cols <- setNames(as.character(col_df$cols),col_df$value) # convert the two columns to a named vector
      # make the plot
      p <- Base_map(df) + 
        scale_fill_manual(values=cols, name = paste0(lgd)) +
        labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
        mytheme
  # stage count and stage count with climate stress exclusion
  } else if (grepl("StageCount", nam)) {
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
      col_key = cbind("cols" = data.frame("cols" = c("gray40","gray70","gray70",Colfunc("mediumpurple4","mediumpurple1",4),Colfunc("magenta4","magenta",4),Colfunc("midnightblue","lightskyblue1",4),Colfunc("cyan4","cyan",4),
                  Colfunc("darkgreen","lightgreen",4),Colfunc("olivedrab4","olivedrab1",4),Colfunc("darkgoldenrod3","yellow",4),Colfunc("sienna4","sienna1",4),Colfunc("darkred","red",4),
                  Colfunc("burlywood4","burlywood",4),Colfunc("rosybrown4","rosybrown1",4),Colfunc("deeppink4","lightpink",4),Colfunc("blue","deepskyblue",4))),
                  "value_orig" = as.numeric(climEXCL_OW_Gens$value))
      col_df <- suppressWarnings(suppressMessages(semi_join(col_key,df2,by="value_orig"))) # keep only the number of colors needed for the data
      cols <- setNames(as.character(col_df$cols),levels(df2$value)) # convert color ramp data frame into a named vector
      # make the plot
      p <- Base_map(df2) + 
            scale_fill_manual(values=cols, name = paste0(lgd)) +
            labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +  
            mytheme
  # Pest Event Map (PEM)
  } else if (grepl("PEM", nam)) {
      if (sum(df$value) <= 0){ # using all(df$value == 0) does not work here correctly
        print(noquote(paste0("Warning: All PEM values are 0 for ",outfl,"_",d)))
        df$value_orig <- df$value # create for ordering recoded values
        df <- mutate(df, value = ifelse(df$value_orig==-2,"exc.-sev",ifelse(df$value_orig==-1,"exc.-mod",df$value)))
        df$value <- factor(df$value,levels= unique(df$value[order(as.numeric(as.character(df$value_orig)))])) # order according to orig. vals (DOY)
        cols <- c("blue")
        if (any(df$value == "exc.-mod")){
          cols <- c("gray70",cols) # add grayshade for excl.-moderate
        }
        if (any(df$value == "exc.-sev")){
          cols <- c("gray40",cols) # add grayshade for excl.-severe
        }
        cols <- setNames(cols,levels(df$value)) # assign a color to each data value
        # make the plot
        p <- Base_map(df) + 
          scale_fill_manual(values=cols, name = paste0(lgd)) +
          labs(title = str_wrap(paste(sp,titl), width = 75), subtitle = str_wrap(subtitl, width = 75)) +
          mytheme
      } else {
        df <- df %>% filter(!(value %in% c(0,366,367))) # remove day 0, 366, 367
        df$value_orig <- df$value # create for ordering recoded values
        # cut the dates into 1 week bins and order them
        df$value <- cut.POSIXt(strptime(df$value_orig, format="%j"), breaks="1 weeks") # climate exc. values (-1 and -2 will become NA)
        df$value <- format(strptime(df$value, format="%Y-%m-%d"), format="%b-%d")
        df <- mutate(df, value = ifelse(df$value_orig==-2,"exc.-sev",ifelse(df$value_orig==-1,"exc.-mod",df$value)))
        df$value <- factor(df$value,levels= unique(df$value[order(as.numeric(as.character(df$value_orig)))])) # order according to orig. vals (DOY)
        # extract dates and months for creating color scale
        dats <- data.frame(value = as.factor(unique(df$value))) %>% filter(!value %in% c("exc.-sev","exc.-mod"))
        dats$mnth <- str_split_fixed(dats$value, pattern = "-", 2)[,1]
        levels(dats$mnth) <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
        dats <- dats[order(dats$value,dats$mnth),]
        mnth_ct <- plyr::count(dats$mnth) # count number of bins
        # make the color ramp (excl.-severe and excl.-moderate are gray shades) - enough colors for 52 weeks
        # in current version, there will be no more than 5 sampling points in a month
        col_key <- cbind(data.frame("cols" = c(Colfunc("rosybrown4","rosybrown1",5),Colfunc("royalblue1","royalblue4",5),Colfunc("mediumpurple1","mediumpurple4",5),Colfunc("magenta","magenta4",5),Colfunc("cyan","cyan4",5),
                                               Colfunc("lightgreen","darkgreen",5),Colfunc("yellow","darkgoldenrod3",5),Colfunc("sienna4","sienna1",5),Colfunc("red","darkred",5),
                                               Colfunc("lightpink","deeppink4",5),Colfunc("olivedrab1","olivedrab4",5),Colfunc("burlywood","burlywood4",5))),
                         data.frame("mnth" = c(rep("Jan",5),rep("Feb",5),rep("Mar",5),rep("Apr",5),rep("May",5),rep("Jun",5),rep("Jul",5),rep("Aug",5),rep("Sep",5),rep("Oct",5),rep("Nov",5),rep("Dec",5))))
        # retain only enough cols needed for data
        # colors are joined to data, so each month has a distinct color with up to 5 shades
        col_key2 <- suppressWarnings(suppressMessages(semi_join(col_key,dats,by="mnth") %>% 
          left_join(.,mnth_ct,by=c("mnth" = "x")) %>%
          group_by(mnth) %>% tidyr::nest() %>%
          mutate(samp = purrr::map(data, ~top_n(data.frame(cols=.$cols), unique(.$freq)))) %>% # top_n takes the first n rows; keep colors consistent across dates
          tidyr::unnest(samp)))
          col_df <- cbind(dplyr::select(dats,value),dplyr::select(col_key2,cols)) # tack colors onto data frame key 
        # add gray shades if climate stress exclusions
        if (any(df$value == "exc.-mod")){
          col_df <- rbind(data.frame("cols" = "gray70","value" = "exc.-mod"), col_df)
        }
        if (any(df$value == "exc.-sev")){
          col_df <- rbind(data.frame("cols" = "gray40","value" = "exc.-sev"), col_df)
        }
        cols <- setNames(as.character(col_df$cols),col_df$value) # convert col_df to a named vector
        # make the plot        
        p <- Base_map(df) + 
                scale_fill_manual(values=cols, name = str_wrap(paste0(lgd), width = 15)) +
                labs(title = str_wrap(paste(sp,titl), width = 75), subtitle = str_wrap(subtitl, width = 75)) +
                mytheme
        }
  # proportion of each life stage and DD accumulation
  } else if (grepl("adultP|eggP|larvP|pupP|ddtotal|mindays", nam, ignore.case = TRUE)) {
    # create plot separately for rasters where all values are 0  
    if (all(df$value == 0)){
        df$value <- factor(df$value)
        p <- Base_map(df) +       
              scale_fill_brewer(palette="Spectral", direction = -1, name = paste0(lgd)) +
              labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
              mytheme
        # else first divide values into bins and then plot
      } else {
        df <- Cut_bins(df)  
        df$value <- factor(df$value)
        p <- Base_map(df) +       
              scale_fill_brewer(palette="Spectral", direction = -1, name = paste0(lgd)) +
              labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
              mytheme
      }
# BOXB changes begin
    # other remaining climate stress exclusion data (heat, chill, dry, all, and min)
  } else if (grepl("heatEXCL|chillEXCL|dryEXCL|AllEXCL|minEXCL", nam)) {
# BOXB changes end
      df <- mutate(df, value = factor(ifelse(value==-2,"excl.-severe",ifelse(value==-1,"excl.-moderate",ifelse(value==0,"not excluded",NA)))))    
      df$value <- factor(df$value, levels = c("excl.-severe","excl.-moderate","not excluded")) # order by values so plots in numerical order
      p <- Base_map(df) +        
            scale_fill_manual(values = c("excl.-severe" = "gray40","excl.-moderate" = "gray70", "not excluded" = "green2"), name=paste0(lgd)) +
            labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
            mytheme
  }
  # save the plot, or else report that there was an error and skip - see "rmessages.txt" for error report
  tryCatch(
    {
      ggsave(p,file=paste0(outfl,"_",d,".png"), width = 6, height = 3, units = c('in'), dpi = 200) # save the plot
        print(noquote(paste0("Saving summary map: ",outfl,"_",d,".png")))
    },
      error=function(e){
        print(noquote(paste0("Error: could not create plot for ",outfl,"_",d)))
    } )
    #ggsave(p,file=paste0(outfl,"_",d,".png"), width = 6, height = 3, units = c('in'), dpi = 200) # save the plot
    #print(noquote(paste0("Saving summary map: ",outfl,"_",d,".png"))) # print progress in log file
}

# Create summary maps (PNG) of heat stress and chill stress units, with max1 (Stress limit 1) and max2 (Stress limit 2) shown as "countour" lines
PlotMap_stress <- function(r,d,max1,max2,titl,lgd,outfl){
  newproj <- "+proj=longlat +datum=WGS84"
  r <- projectRaster(r,crs=newproj)
  nam <- deparse(substitute(r))
  sp <- paste0(gsub(pattern = "_", replacement = " ", fullname),":")
  dat <- as.character(format(strptime(d,format="%Y%m%d"), format="%m/%d/%Y"))
  titl <- paste(titl,dat,sep=" ")
  subtitl <- paste("Maps and modeling",format(Sys.Date(),"%m/%d/%Y"),"by Oregon State University IPPC USPEST.ORG and USDA-APHIS-PPQ; climate data from OSU PRISM Climate Group")
  df <- ConvDF(r)
  df2 <- Stress_Val_Conv(df) # properly formats values
  # create contours from raster values greater than limit 1 and limit 2, if raster values are greater than max1 and/or max2
  max1_c <- tryCatch(
    {
      max1_c <- rasterToContour(r > max1) # object of class "SpatialLinesDataFrame"
    },
    error=function(e){
      max1_c <- 0 # if "rasterToContour" throws an error, then code max1_c as 0 (is.numeric)
    } )
  max2_c <- tryCatch(
    {
      max2_c <- rasterToContour(r > max2) # object of class "SpatialLinesDataFrame"
    },
    error=function(e){
      max2_c <- 0 # if "rasterToContour" throws an error, then code max2_c as 0 (is.numeric)
    } )
  # if all values are 0, then don't include contours (must include this code or it will throw an error)
  if(all(df$value == 0)){
    p <- ggplot(states, aes(x=long, y=lat)) + 
      geom_raster(data=df2, aes(x=x, y=y, fill=value)) +
      geom_path(aes(group=group), color="gray20", lwd = 0.2) +
      theme_map(base_size = 7) +
      scale_fill_manual(values=c("#5E4FA2"), name = paste0(lgd)) +
      labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
      mytheme
    # if any values are greater than 0, but not greater than max1, then don't plot either countour (must include or it will thrown an error)
  } else if (any(df$value > 0) & is.numeric(max1_c)) { # max1_c is class "numeric" (max1_c = 0)
    p <- ggplot(states, aes(x=long, y=lat)) + 
      geom_raster(data=df2, aes(x=x, y=y, fill=value)) + 
      geom_path(aes(group=group), color="gray20", lwd = 0.2) +
      theme_map(base_size = 7) +
      scale_fill_brewer(palette="Spectral", direction = -1, name = paste0(lgd)) +
      labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
      mytheme
    # if any values are greater than limit 1 but less than limit 2, then plot countour line for just limit 1
  } else if (any(df$value > 0) & !is.numeric(max1_c) & is.numeric(max2_c)) { # max1_c is class "SpatialLinesDataFrame" but max2_c is class "numeric" (max2_c = 0)
    p <- ggplot(states, aes(x=long, y=lat)) + 
      geom_raster(data=df2, aes(x=x, y=y, fill=value)) + 
      geom_path(aes(group=group), color="gray20", lwd = 0.2) +
      geom_path(data=max1_c, aes(x=long, y=lat, group=group, color="Stress limit 1"), lwd = 0.15) +
      theme_map(base_size = 7) +
      scale_fill_brewer(palette="Spectral", direction = -1, name = paste0(lgd)) +
      labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
      mytheme
    # if any values are greater than limit1 and limit2, then plot countour lines for both limit1 and limit2
  } else if (!is.numeric(max1_c) & !is.numeric(max2_c)) { # both max1_c and max2c are class "SpatialLinesDataFrame"
    p <- ggplot(states, aes(x=long, y=lat)) + 
      geom_raster(data=df2, aes(x=x, y=y, fill=value)) +
      geom_path(aes(group=group), color="gray20", lwd = 0.2) +
      geom_path(data=max1_c, aes(x=long, y=lat, group=group, color="Stress limit 1"), lwd = 0.15) +
      geom_path(data=max2_c, aes(x=long, y=lat, group=group, color="Stress limit 2"), lwd = 0.15) +
      theme_map(base_size = 7) +
      labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
      scale_fill_brewer(palette="Spectral", direction = -1, name = paste0(lgd)) +
      scale_color_manual(name = "Stress Limits", values = c("Stress limit 1" = "magenta", "Stress limit 2" = "black")) +
      mytheme +
      guides (color = guide_legend(order = 1)) # puts the "Stress limits" legend item on top
  }
  # save the plot, or else report that there was an error and skip - see "rmessages.txt" for error report
  tryCatch(
    {
      ggsave(p,file=paste0(outfl,"_",d,".png"), width = 6, height = 3, units = c('in'), dpi = 200) # save the plot
      print(noquote(paste0("Saving summary map: ",outfl,"_",d,".png")))
    },
    error=function(e){
      print(noquote(paste0("Could not create plot for ",outfl,"_",d)))
    } )
  #ggsave(p,file=paste0(outfl,"_",d,".png"), width = 6, height = 3, units = c('in'), dpi = 200) # save the plot
  #print(noquote(paste0("Saving summary map: ",outfl,"_",d,".png")))
}

cat("6WORKING DIR: ",prism_dir,"\n")

# Write raster (TIF) maps
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
    rast_list <- list(chillEXCL,heatEXCL,dryEXCL,AllEXCL,LifestageEXCL1,LifestageEXCL2,NumGenEXCL1,NumGenEXCL2,StageCountEXCL1,StageCountEXCL2,chillunitsCUM,heatunitsCUM,dryunitsCUM,CumInfRiskEXCL1,CumInfRiskEXCL2)
    names(rast_list) <- c("Chill_Stress_Excl","Heat_Stress_Excl","Dry_Stress_Excl","All_Stress_Excl","LifestageExcl1","LifestageExcl2","NumGenExcl1","NumGenExcl2","StageCountExcl1","StageCountExcl2","Chill_Stress_Units","Heat_Stress_Units","Dry_Stress_units","Cum_Inf_Risk_EXCL1","Cum_Inf_Risk_EXCL2")
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

##########################################################################################
########                  END of function definitions                           ##########
########                BEGINNING of Initialization Section                     ##########
##########################################################################################
#### Abbrevs used in model params:
  #LDT = lower development threshold, temp at which growth = 0 (using PRISM tmean)
  #DD = degree days, number of cumulative heat units to complete that lifestage
  #OW = overwintering (stage, DD reqs., etc)

#### Pest Specific, Multiple Life Stage Phenology Model Parameters: ####
## now read from source param files in ./spp_params/SPP.params
  param_file <- sprintf("%s%s",spp,".params")
  species_params <- sprintf("%s%s",params_dir,param_file)
  if (file.exists(species_params)) {
    cat("SPP PARAMS: ",species_params,"\n")
    source(species_params)
    cat("Reading Params for SPP: ",spp," Fullname: ",fullname,"\n")
  } else {
    cat("PARAM FILE: ",species_params,"...Not found; exiting Program...\n")
    q()  # no reason to keep going without any params
  }

#### Push out a metadata file with all inputs used in model   ####
    cat("output dir: ",output_dir,"\n")
cat("made it here line xxx \n")
setwd(output_dir)
metadata <- sprintf("%s%s","./","metadata.txt")
if (file.exists(metadata)) {
  cat("Metadata for DDRP - Degree-Day, Risk, and Pest event Maps \n",file=metadata)
} else {
  file.create(metadata)
  cat("Metadata for DDRP - Degree-Day, Risk, and Pest event Maps \n",file=metadata)
}
########                      END of Param Handling Section                        #########


########                      BEGIN Start Metadata Output File                     #########

#cat("CLEAN PARAMS: spp startyr startdoy enddoy region_param ovlp exclusions pems out_dir output_option mapA mapE mapL mapP: \n            ",spp,start_year,start_doy,end_doy,region_param,ovlp,exclusions,pems,out_dir,out_option,mapA,mapE,mapL,mapP,"\n")

cat("\n Model Species Params:\n Species Abbrev:",spp,"\n Full Name:",fullname,"\n Pest of:",pestof,file=metadata,append=TRUE)
cat("\n Overwintering Stage:",owstage,"\n Egg Lower Devel Threshold:",eggLDT,"\n Egg Upper Devel Threshold:",eggUDT,file=metadata,append=TRUE)
cat("\n Larvae Lower Devel Threshold:",larvaeLDT,"\n Larvae Upper Devel Threshold:",larvaeUDT,"\n Pupae Lower Devel Threshold:",pupaeLDT,file=metadata,append=TRUE)
cat("\n Pupae Upper Devel Threshold:",pupaeUDT,"\n Adult Lower Devel Threshold:",adultLDT,"\n Adult Upper Devel Threshold:",adultUDT,file=metadata,append=TRUE)

# BOXB changes begin
if (exclusions_stressunits) {
   if (chillstress_threshold) {
   cat("\n Lower Chill Threshold:",chillstress_threshold,"\n Upper Heat Threshold:",heatstress_threshold,"\n Dry Stress Threshold:",drystress_threshold,file=metadata,append=TRUE)
   cat("\n Max Chill Units (lower bound):",chillstress_units_max1,"\n Max Chill Units (upper bound):",chillstress_units_max2,file=metadata,append=TRUE)
   cat("\n Max Heat Stress Units (lower bound):",heatstress_units_max1,"\n Max Heat Stress Units (upper bound):",heatstress_units_max2,file=metadata,append=TRUE)
   cat("\n Max Dry Stress Units (lower bound):",drystress_units_max1,"\n Max Dry Stress Units (upper bound):",drystress_units_max2,file=metadata,append=TRUE)
   }
}
# BOXB changes end

if (exclusions) {
   if (eggLLT) {
   cat("\n Egg Lower Lethal Threshold:",eggLLT,"\n Egg Upper Lethal Threshold:",eggULT,file=metadata,append=TRUE)
   }
   if (eggLLDAYS)    {
   cat("\n Egg Lower Lethal Days Needed:",eggLLDAYS,file=metadata,append=TRUE)
   }
   cat("\n Larvae Lower Lethal Threshold:",larvaeLLT,"\n Larvae Upper Lethal Threshold:",larvaeULT,file=metadata,append=TRUE)
   if (larvaeLLDAYS) {
   cat("\n Larvae Lower Lethal Days Needed:",larvaeLLDAYS,file=metadata,append=TRUE)
   }
   cat("\n Adult Lower Lethal Threshold:",adultLLT,file=metadata,append=TRUE)
   if (adultLLDAYS) {
   cat("\n Adult Lower Lethal Days Needed:",adultLLDAYS,file=metadata,append=TRUE)
   }
}

if (pems){
  cat("\n Number of generations to make Pest Event Maps (PEMs):",PEMnumgens,file=metadata,append=TRUE)
  cat("\n Egg Event DDs and Label:",eggEventDD," ",eggEventLabel,file=metadata,append=TRUE)
  cat("\n Larvae Event DDs and Label:",larvaeEventDD," ",larvaeEventLabel,file=metadata,append=TRUE)
  cat("\n Pupae Event DDs and Label:",pupaeEventDD," ",pupaeEventLabel,file=metadata,append=TRUE)
  cat("\n Adult Event DDs and Label:",adultEventDD," ",adultEventLabel,file=metadata,append=TRUE)
}

cat("\n\n Model Input Params:\n Start Year:",start_year,"\n Start Day-of-year:",start_doy,file=metadata,append=TRUE)
cat("\n End Day-of-Year:",end_doy,"\n Region:",region_param,"\n Stage overlap:",ovlp,file=metadata,append=TRUE)
cat("\n Exclusion Maps:",exclusions,"\n PestEvent Maps:",pems,"\n Output_Dir:",out_dir,file=metadata,append=TRUE)
cat("\n Output Option:",out_option,"\n Map Adults:",mapA,"\n Map Eggs:",mapE,file=metadata,append=TRUE)
cat("\n Map Larvae:",mapL,"\n Map Pupae:",mapP,file=metadata,append=TRUE)

#spp startyr startdoy enddoy region_param ovlp exclusions pems out_dir output_option mapA mapE mapL mapP: \n            ",spp,start_year,start_doy,end_doy,region_param,ovlp,exclusions,pems,out_dir,out_option,mapA,mapE,mapL,mapP,"\n")

#cat(species_params,file=metadata,append=TRUE)

setwd(prism_dir)
########                        END Start Metadata Output File                     #########


#### Example params spp=FCM in case you dont have any param files: ####
  #fullname   <- "False Codling Moth"
  #stgorder   <- c("OA","E","L","P","A","F")
  #owstage    <- "OA"
  #eggLDT     <- 11.93
  #eggUDT     <- 40  # Daiber: 40C as upper dev. threshold
  #larvaeLDT  <- 11.6
  #larvaeUDT  <- 34 #upper dev. threshold-need to verify 
  #pupaeLDT   <- 11.9
  #pupaeUDT   <- 40
  #adultLDT   <- 12.2 #for oviposition
  #adultUDT   <- 40
  #eggDD      <- 69  # round from 69.3
  #larvaeDD   <- 156
  #pupDD      <- 174 #females
  #OWpupDD    <- 86  # text OW stage dev 39 DD "post diapause"
  #adultDD    <- 79 # round from 79.2 time to 50% eggs laid
  #OWadultDD  <- 86  # text OW stage dev 39 DD "post diapause"
  #calctype   <-"average"
  #LLT = lower lethal temperature (PRISM tmin), ULT = upper lethal temperature (PRISM tmax)
  #eggLLT     <- -3
  #eggLLDAYS  <- 7  # NEW v24 add # days to accumulate for local extinction
  #eggULT     <- 41
  #larvaeLLT  <- -12
  #larvaeLLDAYS <- 2
  #larvaeULT  <- 40
  #adultLLT   <- 0.5
  #adultLLDAYS <- 5
  # Pest Event Maps (PEMs) must be turned on for these to get used:
#if (pems) {  # init as spp_param files may not have these (but should if PEMS turned on)
#  PEMnumgens <- 2  # create PEMS for up to this many generations (max is 4)
#  eggEventDD <- 5 # PEMs for egg stage is 5 DDs into egg stage
#  eggEventLabel <- "Beginning of", eggEventLabel # Label for PEM egg stage
#  larvaeEventDD <- 10 # PEMs for larvae stage is 78 DDs into larval stage
#  larvaeEventLabel <- "Early larval development" # Label for PEM larval stage
#  pupaeEventDD <- 10 # PEMs for pupal stage is 58 DDs into pupal stage
#  pupaeEventLabel <- "Early pupal development" # Label for PEM pupal stage
#  adultEventDD <- 30 # PEMs for adult stage (1st oviposition) is 22 DDs into adult stage
#  adultEventLabel <- "Approx. 1st egglaying by females" # Label for PEM adult stage
#}
#### END Pest Specific, Multiple Life Stage Phenology Model Parameters ####

#### Search pattern for PRISM/other daily temperature grids. Load them for processing. ####
pattern <- glob2rx(paste0("*PRISM_tmax_*",start_year,"*bil$")) # glob2rx searches for multiple conditions - all must be met
files <- list.files(pattern=pattern, all.files=FALSE, full.names=TRUE)

#Check that there are enough files for the year
#length(files)
numlist <- vector()

#### New code (V23) to add to DDRPx to pick most useful/recent file types needed for current year data ####
# stable (4)  > provisional (3) > early (2)  > 7day forecast (NA) > 90day fc (NA) > 10yr Avg data (new!) (1) > none (0)
# just do for tmax, assume tmin is same (add error flagging somehow??)
filelist <- vector()
# set up R version of hash tables (new.env)
maxfiles <- new.env(hash=TRUE, parent=emptyenv(), size=100L)
minfiles <- new.env(hash=TRUE, parent=emptyenv(), size=100L)
quality  <- new.env(hash=TRUE, parent=emptyenv(), size=100L)

# init file type hash keys
for (mon in 1:12) {
   for (day in 1:31) {
      skey = sprintf("%04s%02d%02d",start_year,mon,day)
      quality[[skey]] <- 0
   }
}

fkey <- 0
########                      END of Initialization Section                        #########


########                 Begin of actual stepping through the model                #########
#### Go through all files, store best set in a hash using these 0 to 4 priority numbers ####
for (file in sort(files)) {
    status = strsplit(file,split="_")[[1]][3]  # type of data: early provisional stable etc
    dat <- file %>% str_match_all("[0-9]+") %>% unlist %>% tail(.,n=1) # extract the number of the file name, which is the date
    fkey <- dat
    if ((status == "stable") && (quality[[fkey]] < 5)) {
      filelist <-c(filelist,file)
      numlist <-c(numlist,dat)
      quality[[fkey]] <- 5
      maxfiles[[fkey]] <- file  # tmax stored now check tmin
      minfile  <- gsub("tmax","tmin",file,perl=TRUE)
      if (file.exists(minfile)) {
        minfiles[[fkey]] <- minfile  # tmin exists assume its good
		}
		else {
        cat("no tmin file found equiv to avail tmax file: ",minfile,"\n")
		  }
    # FIXME no smarts to use lower grade tmin

      # BOXB changes begin
      tmeanfile  <- gsub("tmax","tmean",file,perl=TRUE)
      if (file.exists(tmeanfile)) {
        tmeanfiles[[fkey]] <- tmeanfile  # tmean exists assume its good
      } else { cat("no tmean file found equiv to avail tmax file: ",tmeanfile,"\n")  # FIXME no smarts to use lower grade tmean
      }
      tdmeanfile  <- gsub("tmax","tdmean",file,perl=TRUE)
      if (file.exists(tdmeanfile)) {
        tdmeanfiles[[fkey]] <- tdmeanfile  # tdmean exists assume its good
      } else { cat("no tdmean file found equiv to avail tmax file: ",tdmeanfile,"\n")  # FIXME no smarts to use lower grade tdmean
      }
      pptfile  <- gsub("tmax","ppt",file,perl=TRUE)
      if (file.exists(pptfile)) {
        pptfiles[[fkey]] <- pptfile  # ppt exists assume its good
      } else { cat("no ppt file found equiv to avail tmax file: ",pptfile,"\n")  # FIXME no smarts to use lower grade ppt
      }
      # BOXB changes end

    } else if ((status == "provisional") && (quality[[fkey]] < 4)) {
      filelist <-c(filelist,file)
      numlist <-c(numlist,dat)
      quality[[fkey]] <- 4
      maxfiles[[fkey]] <- file
      minfile  <- gsub("tmax","tmin",file,perl=TRUE)
      if (file.exists(minfile)) {
        minfiles[[fkey]] <- minfile
		}
		else {
        cat("no tmin file found equiv to avail tmax file: ",minfile,"\n")
		}
    } else if ((status == "early") && (quality[[fkey]] < 3)) {
      filelist <-c(filelist,file)
      numlist <-c(numlist,dat)
      quality[[fkey]] <- 3
      maxfiles[[fkey]] <- file
      minfile  <- gsub("tmax","tmin",file,perl=TRUE)
      if (file.exists(minfile)) {
        minfiles[[fkey]] <- minfile
		}
		else {
        cat("no tmin file found equiv to avail tmax file: ",minfile,"\n")
		}
	#
	# future short and extended forecast options/filetypes get inserted here
	#
    #} If analysis is based on other types of data; e.g., the 10 yr avg 2006-2015 PRISM data
    } else if (!status %in% c("stable","provisional","early")) {
      filelist <-c(filelist,file)
      numlist <-c(numlist,dat)
      quality[[fkey]] <- 2
      maxfiles[[fkey]] <- file
      minfile  <- gsub("tmax","tmin",file,perl=TRUE)
      if (file.exists(minfile)) {
        cat("found other type of minfile: ",minfile,"\n")
        minfiles[[fkey]] <- minfile
		}
		else {
        cat("no tmin file found equiv to avail tmax file: ",minfile,"\n")
		}
   #} else if ((status == "10yr0716") && (quality[[fkey]] < 1)) {
   #   filelist <-c(filelist,file)
   #   numlist <-c(numlist,num)
   #   quality[[fkey]] <- 1
   #   maxfiles[[fkey]] <- file
   #   minfile  <- gsub("tmax","tmin",file,perl=TRUE)
   #   if (file.exists(minfile)) {
   #     #cat("found nmme minfile: ",minfile,"\n")
   #     minfiles[[fkey]] <- minfile
   #   }
   #   else {
   #     cat("no tmin file found equiv to avail tmax file: ",minfile,"\n")
   #   }

    } else {
        cat("no match for this file: ",file," ",dat,"\n")
    }
}
# cat("Done1. Full numlist: ",sort(numlist),"\n")
sortedlist <- sort(numlist)
#### END code to distinguish file type/quality  ####

#### Read in first raster as a template, crop to REGIONS setting and intialize tracking rasters ####
# SET Region Here:
  #  cat("LINE 313:",sortedlist,"\n")
template <- crop(raster(files[1]),REGION)
#template <- crop(raster(maxfiles[[1]]),REGION)
template[!is.na(template)] <- 0

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
cat("made it here line endxxxx \n")
  #  cat("LINE 323: \n")

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
			  if (i == "OA") { dd3 <- dd3tmp * LSOW3 }
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
