# Code to transform the modelevaluation.org data into a binned dataset
# where each variable is split into 50 bins. This produces a dataset of observed
# and modeled fluxes assigned to 8-dimensional cells

rm(list=ls())

# User-defined variables
DaytimeFlags <- c(TRUE, FALSE)
PhysicalFlags <- c(TRUE, FALSE)
WindyFlags <- c(TRUE, FALSE)
# The below data file is produced using modelevaluation.org
TempSaveFile = "data/PLUMBER2AnalysisScript_RawData.RData"

# Set variables expected for the MEorg output
sites = c('AR-SLu', 'AT-Neu', 'AU-ASM', 'AU-Cow', 'AU-Cpr', 'AU-Ctr', 'AU-Cum', 'AU-DaP', 'AU-DaS', 'AU-Dry', 'AU-Emr', 'AU-Gin', 'AU-GWW', 'AU-How', 'AU-Lit', 'AU-Otw', 'AU-Sam', 'AU-Stp', 'AU-TTE', 'AU-Tum', 'AU-Wrr', 'BE-Bra', 'BE-Lon', 'BE-Vie', 'BR-Sa3', 'BW-Ma1', 'CA-Qcu', 'CA-Qfo', 'CH-Cha', 'CH-Dav', 'CH-Fru', 'CH-Oe1', 'CN-Cha', 'CN-Cng', 'CN-Dan', 'CN-Din', 'CN-Du2', 'CN-HaM', 'CN-Qia', 'CZ-wet', 'DE-Bay', 'DE-Geb', 'DE-Gri', 'DE-Hai', 'DE-Kli', 'DE-Meh', 'DE-Obe', 'DE-Seh', 'DE-SfN', 'DE-Tha', 'DE-Wet', 'DK-Fou', 'DK-Lva', 'DK-Ris', 'DK-Sor', 'DK-ZaH', 'ES-ES1', 'ES-ES2', 'ES-LgS', 'ES-LMa', 'ES-VDA', 'FI-Hyy', 'FI-Kaa', 'FI-Lom', 'FI-Sod', 'FR-Fon', 'FR-Gri', 'FR-Hes', 'FR-LBr', 'FR-Lq1', 'FR-Lq2', 'FR-Pue', 'GF-Guy', 'HU-Bug', 'ID-Pag', 'IE-Ca1', 'IE-Dri', 'IT-Amp', 'IT-BCi', 'IT-CA1', 'IT-CA2', 'IT-CA3', 'IT-Col', 'IT-Cpz', 'IT-Isp', 'IT-Lav', 'IT-LMa', 'IT-Mal', 'IT-MBo', 'IT-Noe', 'IT-Non', 'IT-PT1', 'IT-Ren', 'IT-Ro1', 'IT-Ro2', 'IT-SR2', 'IT-SRo', 'JP-SMF', 'NL-Ca1', 'NL-Hor', 'NL-Loo', 'PL-wet', 'PT-Esp', 'PT-Mi1', 'PT-Mi2', 'RU-Fyo', 'SD-Dem', 'SE-Deg', 'UK-Gri', 'UK-Ham', 'US-AR1', 'US-AR2', 'US-ARM', 'US-Aud', 'US-Bar', 'US-Bkg', 'US-Blo', 'US-Bo1', 'US-Cop', 'US-FPe', 'US-GLE', 'US-Goo', 'US-Ha1', 'US-Ho1', 'US-KS2', 'US-Los', 'US-Me2', 'US-Me4', 'US-Me6', 'US-MMS', 'US-MOz', 'US-Myb', 'US-Ne1', 'US-Ne2', 'US-Ne3', 'US-NR1', 'US-PFa', 'US-Prr', 'US-SP2', 'US-SP3', 'US-SRG', 'US-SRM', 'US-Syv', 'US-Ton', 'US-Tw4', 'US-Twt', 'US-UMB', 'US-Var', 'US-WCr', 'US-Whs', 'US-Wkg', 'ZA-Kru', 'ZM-Mon')
models = c('1lin_raw','3km27_raw','6km729lag_raw','RF_raw','LSTM_raw','CABLE','CABLE-POP', 'CHTESSEL_1', 'CLM5', 'GFDL', 'JULES_GL9', 'JULES_GL9_LAI', 'MATSIRO','NoahMP','ORCHIDEE2','ORCHIDEE3') # LSMs and benchmarks
benchmarks = c('1lin_raw','3km27_raw','LSTM_raw') # user-nominated benchmarks, as though in modelevaluation.org


# Load packages and 'pals' 
library(tidyverse)
library(parallel)

# The pals data is sourced from modelevaluation.org
FunctionSources = list.files(path = "pals", full.names = TRUE, pattern="*.R")
sapply(FunctionSources,source,.GlobalEnv)

# Define function for extracting model errors
DataExtraction <- function(siteIDX, breaks, IDXlist, AllObs, AllObsSubsetV, AllMods, modsiteidx, sitesfailed, models, sites, Daytime = FALSE, Physical = FALSE, Windy = FALSE){
      
      message(paste0("        site ", siteIDX, " / ", length(sites)))
      
      site.data <- data.frame("time" = as.POSIXlt(paste0(AllObsSubsetV[[siteIDX]][[IDXlist$RainfIDX]]$timing$syear,"-1-1 00:00"), tz = "UTC") + seq(0, by = 1800, length.out = AllObsSubsetV[[siteIDX]][[IDXlist$RainfIDX]]$timing$tsteps))
      
      # First check there's good data at this site:
      if(all(is.na(AllObs[[siteIDX]][[IDXlist$QleIDX]]$data)) | all(is.na(AllObsSubsetV[[siteIDX]][[IDXlist$RainfIDX]]$data))){
            sitesfailed[1] = sitesfailed[1] + 1
      }else{
            # Extract fluxes
            site.data$obs.nee <- as.vector(AllObs[[siteIDX]][[IDXlist$NEEIDX]]$data)
            site.data$obs.qh <- as.vector(AllObs[[siteIDX]][[IDXlist$QhIDX]]$data)
            site.data$obs.qle <- as.vector(AllObs[[siteIDX]][[IDXlist$QleIDX]]$data)
            
            # Extract met vars
            site.data$qair_raw <- as.vector(AllObsSubsetV[[siteIDX]][[IDXlist$QairIDX]]$data)
            site.data$swdown_raw <- as.vector(AllObsSubsetV[[siteIDX]][[IDXlist$SWdownIDX]]$data)
            site.data$tair_raw <- as.vector(AllObsSubsetV[[siteIDX]][[IDXlist$TairIDX]]$data) - 273.15 # change K to C
            site.data$wind_raw <- as.vector(AllObsSubsetV[[siteIDX]][[IDXlist$WindIDX]]$data)
            site.data$lwdown_raw <- as.vector(AllObsSubsetV[[siteIDX]][[IDXlist$LWdownIDX]]$data)
            site.data$vpd_raw <- as.vector(AllObsSubsetV[[siteIDX]][[IDXlist$VPDIDX]]$data) / 1000 # change Pa to kPa
            site.data$lai_raw <- as.vector(AllObsSubsetV[[siteIDX]][[IDXlist$LAIIDX]]$data)
            site.data$psurf_raw <- as.vector(AllObsSubsetV[[siteIDX]][[IDXlist$PSurfIDX]]$data) / 1000 # change Pa to kPa
            site.data$precip_raw <- as.vector(AllObsSubsetV[[siteIDX]][[IDXlist$RainfIDX]]$data)
            
            
            # Prep for model predictions and errors
            site.data[, paste0(models, ".qle_prediction")] <- NA
            site.data[, paste0(models, ".qh_prediction")] <- NA
            site.data[, paste0(models, ".nee_prediction")] <- NA
            site.data[, paste0(models, ".qle_error")] <- NA
            site.data[, paste0(models, ".qh_error")] <- NA
            site.data[, paste0(models, ".nee_error")] <- NA
            
            # Check model data is good
            for(m in 1:length(models)){
                  if(all(is.na(AllMods[[ modsiteidx[siteIDX,m] ]][[IDXlist$QleIDX]]$data))){
                        sitesfailed[m+1] = sitesfailed[m+1] + 1
                  } else if (!(AllMods[[ modsiteidx[siteIDX,m] ]][[IDXlist$QleIDX]]$timing$interval %in% c('subdaily'))){
                        # Fail if the model data isn't daily
                        sitesfailed[m+1] = sitesfailed[m+1] + 1
                  }else if (AllMods [[ modsiteidx[siteIDX,m] ]][[IDXlist$QleIDX]]$metricerr == TRUE){
                        # Fail if something is squiffy
                        sitesfailed[m+1] = sitesfailed[m+1] + 1
                  } else {
                    # If good, extract model prediction
                        site.data[, paste0(models[m], ".qle_prediction")] <- as.vector(AllMods[[ modsiteidx[siteIDX,m] ]][[IDXlist$QleIDX]]$data)
                        site.data[, paste0(models[m], ".qh_prediction")] <- as.vector(AllMods[[ modsiteidx[siteIDX,m] ]][[IDXlist$QhIDX]]$data)
                        site.data[, paste0(models[m], ".nee_prediction")] <- as.vector(AllMods[[ modsiteidx[siteIDX,m] ]][[IDXlist$NEEIDX]]$data)
                    # If good, extract model error
                        site.data[, paste0(models[m], ".qle_error")] <- as.vector(AllMods[[ modsiteidx[siteIDX,m] ]][[IDXlist$QleIDX]]$data) - site.data$obs.qle
                        site.data[, paste0(models[m], ".qh_error")] <- as.vector(AllMods[[ modsiteidx[siteIDX,m] ]][[IDXlist$QhIDX]]$data) - site.data$obs.qh
                        site.data[, paste0(models[m], ".nee_error")] <- as.vector(AllMods[[ modsiteidx[siteIDX,m] ]][[IDXlist$NEEIDX]]$data) - site.data$obs.nee
                  }
            }
      }

      # Do we want only daytime data?
      if(Daytime == TRUE){
        site.data <- site.data[site.data$swdown_raw > 10, ]
      }
      
      # Do we want only high wind data?
      if(Windy == TRUE){
        site.data <- site.data[site.data$wind_raw > 2, ]
      }
      
      # Do we want tofilter by physical consistency?
      if(Physical == TRUE){
        # First (and so far only) check is for non-physical specific humidity/VPD at the given temperature and pressure
        # We cannot have VPD greater than the SVP (since VPD = SVP-AVP)
        # Calculate saturated vapour pressure
        SVP <- 610.94 * exp((17.625 * site.data$tair_raw) / (site.data$tair_raw + 243.04)) # Pa from Tair in C
        # Specific Humidity at saturation
        qair_sat <- 0.622*SVP/(site.data$psurf_raw*1000 - SVP) # Note psurf_raw is in kPa
        site.data <- site.data[site.data$qair_raw <= qair_sat & site.data$vpd_raw*1000 <= SVP, ] # Note vpd_raw is in kPa
      }
      
      # Bin the met values
      site.data$qair <- cut(site.data$qair_raw, breaks$qair, include.lowest = TRUE, dig.lab = 5)
      site.data$swdown <- cut(site.data$swdown_raw, breaks$swdown, include.lowest = TRUE, dig.lab = 5)
      site.data$tair <- cut(site.data$tair_raw, breaks$tair, include.lowest = TRUE, dig.lab = 5)
      site.data$wind <- cut(site.data$wind_raw, breaks$wind, include.lowest = TRUE, dig.lab = 5)
      site.data$lwdown <- cut(site.data$lwdown_raw, breaks$lwdown, include.lowest = TRUE, dig.lab = 5)
      site.data$vpd <- cut(site.data$vpd_raw, breaks$vpd, include.lowest = TRUE, dig.lab = 5)
      site.data$lai <- cut(site.data$lai_raw, breaks$lai, include.lowest = TRUE, dig.lab = 5)
      site.data$psurf <- cut(site.data$psurf_raw, breaks$psurf, include.lowest = TRUE, dig.lab = 5)
      site.data$precip <- cut(site.data$precip_raw, breaks$precip, include.lowest = TRUE, dig.lab = 5)
      
      site.data$site <- sites[siteIDX]
      
      return(site.data)
      
}


# Load in the MEorg output
load(TempSaveFile)

# Note the variable index in the list of loaded variables:
foundQle = FALSE
foundQh = FALSE
foundNEE = FALSE
foundRainf = FALSE
foundTair = FALSE
foundSWdown = FALSE
foundWind = FALSE
foundQair = FALSE
foundLWdown = FALSE
foundLAI = FALSE
foundPSurf = FALSE
foundVPD = FALSE

IDXlist <- list()
for(sv in 1:length(subsetvars)){
      if(subsetvars[[sv]]$Name[1]=='Rainf'){
            foundRainf = TRUE
            IDXlist$RainfIDX = sv
      } else if(subsetvars[[sv]]$Name[1]=='Tair'){
            foundTair = TRUE
            IDXlist$TairIDX = sv
      } else if(subsetvars[[sv]]$Name[1]=='SWdown'){
            foundSWdown = TRUE
            IDXlist$SWdownIDX = sv
      } else if(subsetvars[[sv]]$Name[1]=='Wind'){
            foundWind = TRUE
            IDXlist$WindIDX = sv
      } else if(subsetvars[[sv]]$Name[1]=='Qair'){
            foundQair = TRUE
            IDXlist$QairIDX = sv
      } else if(subsetvars[[sv]]$Name[1]=='LWdown'){
            foundLWdown = TRUE
            IDXlist$LWdownIDX = sv
      } else if(subsetvars[[sv]]$Name[1]=='LAI'){
        foundLAI = TRUE
        IDXlist$LAIIDX = sv
      } else if(subsetvars[[sv]]$Name[1]=='PSurf'){
        foundPSurf = TRUE
        IDXlist$PSurfIDX = sv
      } else if(subsetvars[[sv]]$Name[1]=='VPD'){
        foundVPD = TRUE
        IDXlist$VPDIDX = sv
      }
}

for(v in 1:length(vars)){
      if((vars[[v]]$Name[1]=='Qle') | (vars[[v]]$Name[1]=='Qle_cor')){
            foundQle = TRUE
            IDXlist$QleIDX = v
      }
      if((vars[[v]]$Name[1]=='Qh') | (vars[[v]]$Name[1]=='Qh_cor')){
            foundQh = TRUE
            IDXlist$QhIDX = v
      }
      if((vars[[v]]$Name[1]=='NEE') | (vars[[v]]$Name[1]=='NEE_cor')){
            foundNEE = TRUE
            IDXlist$NEEIDX = v
      }
}

# Prepare vector for failure flag
sitesfailed = c(1:(length(models)+1)) * 0 # initialise sites failed for obs + all models at 0

# Define the breaks
# Aim for 50 bins
nbins <- 50
breaks <- list("lwdown" = seq(0, 700, length.out = nbins+1),
               "qair" = seq(0, 0.06, length.out = nbins+1),
               "swdown" = seq(0, 1350, length.out = nbins+1),
               "tair" = seq(220-273.15, 325-273.15, length.out = nbins+1), # change K to C
               "wind" = seq(0, 40, length.out = nbins+1),
               "vpd" = seq(-100 / 1000, 8200 / 1000, length.out = nbins+1), # change Pa to kPa
               "lai" = seq(0, 6.5, length.out = nbins+1),
               "psurf" = seq(59000 / 1000, 145000 / 1000, length.out = nbins+1), # change Pa to kPa
               "precip" = seq(0, 0.05, length.out = nbins+1)) 


# Extract the data for all combos of filter flags
for (DaytimeFlag in DaytimeFlags){
  for (PhysicalFlag in PhysicalFlags){
    for (WindyFlag in WindyFlags){
      
      flags <- paste0(ifelse(DaytimeFlag,"_Daytime",""),
                      ifelse(PhysicalFlag,"_Physical",""),
                      ifelse(WindyFlag,"_Windy",""))
      
      if (!file.exists(paste0("data/ExtractedBinnedData",
                              flags,
                              ".rds"))){
        
        message(paste0("Calculating ExtractedBinnedData",flags))
        
        data <- lapply(1:length(sites), function(x) DataExtraction(x, 
                                                                   breaks,
                                                                   IDXlist,
                                                                   AllObs,
                                                                   AllObsSubsetV,
                                                                   AllMods,
                                                                   modsiteidx,
                                                                   sitesfailed,
                                                                   models, 
                                                                   sites,
                                                                   Daytime = DaytimeFlag,
                                                                   Physical = PhysicalFlag,
                                                                   Windy = WindyFlag))
        # Bind the results together
        data <- data %>% bind_rows()
        
        # Save it for later
        saveRDS(data, 
        file = paste0("data/ExtractedBinnedData",
                      flags,
                      ".rds"))

        rm(data)
        gc()
      } else {
        message(paste0("ExtractedBinnedData",
                       flags,
                       ".rds exists!"))
      }
    }
  }
}


