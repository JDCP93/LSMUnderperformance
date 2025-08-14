# Code to calculate the seven metrics at each site for each flux and model

# Tidy up
rm(list = ls())

# Limit to models of interest and sites that are high quality
fluxes <- c("nee", "qh", "qle")
models <- c('1lin_raw','3km27_raw','6km729lag_raw','RF_raw','LSTM_raw','CABLE', 'CABLE-POP', 'CHTESSEL_1', 'CLM5', 'GFDL', 'JULES_GL9', 'JULES_GL9_LAI', 'MATSIRO','NoahMP','ORCHIDEE2','ORCHIDEE3') # these are all the locally forced LSMs in PLUMBER2
sites <- c('AR-SLu', 'AT-Neu', 'AU-ASM', 'AU-Cow', 'AU-Cpr', 'AU-Ctr', 'AU-Cum', 'AU-DaP', 'AU-DaS', 'AU-Dry', 'AU-Emr', 'AU-Gin', 'AU-GWW', 'AU-How', 'AU-Lit', 'AU-Otw', 'AU-Sam', 'AU-Stp', 'AU-TTE', 'AU-Tum', 'AU-Wrr', 'BE-Bra', 'BE-Lon', 'BE-Vie', 'BR-Sa3', 'BW-Ma1', 'CA-Qcu', 'CA-Qfo', 'CH-Cha', 'CH-Dav', 'CH-Fru', 'CH-Oe1', 'CN-Cha', 'CN-Cng', 'CN-Dan', 'CN-Din', 'CN-Du2', 'CN-HaM', 'CN-Qia', 'CZ-wet', 'DE-Bay', 'DE-Geb', 'DE-Gri', 'DE-Hai', 'DE-Kli', 'DE-Meh', 'DE-Obe', 'DE-Seh', 'DE-SfN', 'DE-Tha', 'DE-Wet', 'DK-Fou', 'DK-Lva', 'DK-Ris', 'DK-Sor', 'DK-ZaH', 'ES-ES1', 'ES-ES2', 'ES-LgS', 'ES-LMa', 'ES-VDA', 'FI-Hyy', 'FI-Kaa', 'FI-Lom', 'FI-Sod', 'FR-Fon', 'FR-Gri', 'FR-Hes', 'FR-LBr', 'FR-Lq1', 'FR-Lq2', 'FR-Pue', 'GF-Guy', 'HU-Bug', 'ID-Pag', 'IE-Ca1', 'IE-Dri', 'IT-Amp', 'IT-BCi', 'IT-CA1', 'IT-CA2', 'IT-CA3', 'IT-Col', 'IT-Cpz', 'IT-Isp', 'IT-Lav', 'IT-LMa', 'IT-Mal', 'IT-MBo', 'IT-Noe', 'IT-Non', 'IT-PT1', 'IT-Ren', 'IT-Ro1', 'IT-Ro2', 'IT-SR2', 'IT-SRo', 'JP-SMF', 'NL-Ca1', 'NL-Hor', 'NL-Loo', 'PL-wet', 'PT-Esp', 'PT-Mi1', 'PT-Mi2', 'RU-Fyo', 'SD-Dem', 'SE-Deg', 'UK-Gri', 'UK-Ham', 'US-AR1', 'US-AR2', 'US-ARM', 'US-Aud', 'US-Bar', 'US-Bkg', 'US-Blo', 'US-Bo1', 'US-Cop', 'US-FPe', 'US-GLE', 'US-Goo', 'US-Ha1', 'US-Ho1', 'US-KS2', 'US-Los', 'US-Me2', 'US-Me4', 'US-Me6', 'US-MMS', 'US-MOz', 'US-Myb', 'US-Ne1', 'US-Ne2', 'US-Ne3', 'US-NR1', 'US-PFa', 'US-Prr', 'US-SP2', 'US-SP3', 'US-SRG', 'US-SRM', 'US-Syv', 'US-Ton', 'US-Tw4', 'US-Twt', 'US-UMB', 'US-Var', 'US-WCr', 'US-Whs', 'US-Wkg', 'ZA-Kru', 'ZM-Mon')

# Set the loss limits (use "None" for no loss filtering)
LossLimitValues <- c("None", 5, 20, 50, 80, 95)

# Load libraries
library(tidyverse)
library(parallel)
library(data.table)

Metrics <- function(Data, Flux, ModelName, SiteName, DaytimeFlag, PhysicalFlag, WindyFlag, LossLimit){
  # Function to calculate metrics and save them in a list 
  # Stolen almost entirely from meorg code
  
  # Create masks to make sure metrics are calculated over the same time steps for obs and model
  obsNotNAmask <- !is.na(as.vector(Data$obs))
  NotNAmask <- ! ( is.na(as.vector(Data$obs)) | is.na(as.vector(Data[, ModelName])) )
  sum(NotNAmask)/nrow(Data)
  
  # Calculate the PLUMBER2 metrics
  Correlation <- cor(as.vector(Data[, ModelName]), as.vector(Data$obs), use='complete.obs')
  
  Q5Diff <- abs(quantile(as.vector(Data[, ModelName])[NotNAmask],0.05) - quantile(as.vector(Data$obs)[NotNAmask],0.05))
  Q95Diff <- abs(quantile(as.vector(Data[, ModelName])[NotNAmask],0.95) - quantile(as.vector(Data$obs)[NotNAmask],0.95))
  
  MBE <- mean(as.vector(Data[, ModelName])-as.vector(Data$obs), na.rm=TRUE)
  
  NME <- sum(abs(as.vector(Data[, ModelName])-as.vector(Data$obs)),na.rm=TRUE)/sum(abs(mean(Data$obs,na.rm=TRUE)-Data$obs), na.rm=TRUE)
  
  SDdiff <- abs(1 - (sd(as.vector(Data[, ModelName])[NotNAmask])/sd(as.vector(Data$obs)[NotNAmask])))
  
  nbins <- 2048 # Number of bins for density plots, taken from pals/Constants.R
  minsample <- 1000 # Number of data we want before trying to calculate a density, taken from pals/Constants.R
  if(sum(NotNAmask) < minsample){ # threshold for PDFoverlap to mean something
    PDFoverlap <- NA
  } else {
    xmax <- max(as.vector(Data[, ModelName])[NotNAmask], as.vector(Data$obs)[NotNAmask])
    xmin <- min(as.vector(Data[, ModelName])[NotNAmask], as.vector(Data$obs)[NotNAmask])
    modden <- density(as.vector(Data[, ModelName])[NotNAmask], from=xmin, to=xmax, n=nbins)
    obsden <- density(as.vector(Data$obs)[NotNAmask], from=xmin, to=xmax, n=nbins)
    # Calculate overlap with observational pdf:
    intersection <- c()
    for(b in 1:nbins){
      intersection[b] <- min(modden[[2]][b], obsden[[2]][b])
    }
    PDFoverlap <- min(signif(sum(intersection) * (xmax - xmin) / nbins * 100, 2), 100) # as a percentage
  }
  
  # Nice dataframe for output
  metrics <- data.frame(Flux = Flux,
                        Model = ModelName,
                        Site = SiteName,
                        LossFilter = LossLimit,
                        Daytime = DaytimeFlag,
                        Physical = PhysicalFlag,
                        Windy = WindyFlag,
                        Correlation = Correlation,
                        Q5Diff = Q5Diff,
                        Q95Diff = Q95Diff,
                        MBE = MBE,
                        NME = NME,
                        SDdiff = SDdiff,
                        PDFOverlap = PDFoverlap)
  
  return(metrics)
}


MetricCalculation <- function(Flux, ModelName, SiteName, DaytimeFlag, PhysicalFlag, WindyFlag, LossLimit, ExtractedData, LossData){
  
  # Define the run name since we have so many variables
  run_name <- paste0(Flux,
                     "_",
                     ModelName,
                     "_",
                     SiteName,
                     "_",
                     ifelse(LossLimit == "None","NoLRF",paste0(LossLimit,"pcLRF")),
                     ifelse(DaytimeFlag, "_Daytime",""),
                     ifelse(PhysicalFlag,"_Physical",""),
                     ifelse(WindyFlag,"_Windy",""))

  # In many circumstances, we need to add blank metrics as the data doesn't exist
  blank_df <- data.frame(Flux = Flux,
                         Model = ModelName,
                         Site = SiteName,
                         LossFilter = LossLimit,
                         Daytime = DaytimeFlag,
                         Physical = PhysicalFlag,
                         Windy = WindyFlag,
                         Correlation = NA,
                         Q5Diff = NA,
                         Q95Diff = NA,
                         MBE = NA,
                         NME = NA,
                         SDdiff = NA,
                         PDFOverlap = NA)
    
    message("    ",SiteName," - ",ModelName)
    
    # Define the variable names and a useful flux selector variable
    VarNames <- c("lai","lwdown","precip","qair","swdown","tair","vpd","wind")
    FluxSelector <- c("nee", "qh", "qle")[c(!(Flux == "nee"), !(Flux == "qh"), !(Flux == "qle"))]
    
    # Make the general data more usable
    ExtractedData_modified <- ExtractedData %>%
      filter(site == SiteName) %>% # Filter for the site
      select(-contains("_eb"), # Ignore energy corrected models
             -contains(FluxSelector)) %>% # Ignore other fluxes
      rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
      select(-contains("_error")) %>% # No longer need the error, only the predictions
      rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
      select(time, obs, all_of(contains("_raw")), any_of(c(VarNames, ModelName))) # Select what we need
    
    
    if(LossLimit == "None"){
      # The original data
      MetricData <- ExtractedData_modified
    } else {
      setDT(LossData)
      setkey(LossData, Model)
      # Find the data we want to kick out
      LossData_filtered <- LossData[.(ModelName)][Loss >= as.numeric(LossLimit)]
      setDT(ExtractedData_modified)
      # Data begone!
      MetricData <- ExtractedData_modified[!LossData_filtered, on = VarNames]
      MetricData <- as.data.frame(MetricData)
    }
    
    if (!(ModelName %in% colnames(MetricData))){
      # Some models don't have certain data
      message("      - Data does not exist!")
      Metrics <- blank_df
    } else if (sum(!is.na(MetricData[, ModelName])) == 0){
      # Some models are all NAs, rather than totally missing data
      message("      - Data is all NA!")
      Metrics <- blank_df
    } else if (nrow(MetricData) == 0){
      # If every time step is a loss, then we don't have any timesteps left to calculate metrics for!
      message("      - Every timestep loses!")
      Metrics <- blank_df
    } else {
      # Calculate the metrics and save them!
      Metrics <- Metrics(MetricData,
                         Flux,
                         ModelName,
                         SiteName,
                         DaytimeFlag,
                         PhysicalFlag,
                         WindyFlag,
                         LossLimit)
    }
  return(Metrics)
}




# Combine the different characteristics we have
Combos <- expand.grid("Model" = models,
                      "Site" = sites,
                      stringsAsFactors = FALSE)

# Apply across our 8 flags
for(DaytimeFlag in c(FALSE, TRUE)){
  for(WindyFlag in c(FALSE, TRUE)){
    for(PhysicalFlag in c(FALSE, TRUE)){

      flags <- paste0(ifelse(DaytimeFlag,"_Daytime",""),
                      ifelse(PhysicalFlag,"_Physical",""),
                      ifelse(WindyFlag,"_Windy",""))
      
      # Load the general data
      ExtractedData_Location <- paste0("data/ExtractedBinnedData",
                                       flags,
                                       ".rds")
      ExtractedDataframe <- readRDS(ExtractedData_Location)
      
      for (Flux in fluxes){
        # Load the data where we calculated the LSM losses
        # Different file for each flux
        LossData_Location <- paste0("data/LossByModelAndVars_",
                                    Flux, 
                                    flags,
                                    ".rds")
        LossDataframe <- readRDS(LossData_Location)
        
        for(LossLimitValue in LossLimitValues){
          
          file.name <- paste0(Flux,
                              "_", 
                              ifelse(LossLimitValue == "None",
                                               "NoLRF",
                                               paste0(LossLimitValue,"pcLRF")),
                              flags)
          
          if(file.exists(paste0("data/Metrics/SiteModelMetrics_",
                                file.name,
                                ".rds"))){
            message(paste0(file.name,
                           " data exists!"))
          } else {
            message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    " - ",
                    file.name)
            # Calculate the metrics for each site, model, and loss limit
            # (Obviously we already splitting by daytime/physical/windy and fluxes)
            result <- mclapply(1:nrow(Combos), 
                             function(x) MetricCalculation(Flux,
                                                           Combos[x, ]$Model,
                                                           Combos[x, ]$Site,
                                                           DaytimeFlag, 
                                                           PhysicalFlag,
                                                           WindyFlag,
                                                           LossLimit = LossLimitValue,
                                                           ExtractedDataframe, 
                                                           LossDataframe),
                             mc.cores = 8)
            # Combine the results and save them
            result <- result %>% bind_rows()
            rownames(result) <- NULL
            saveRDS(result, file = paste0("data/Metrics/SiteModelMetrics_",
                                          file.name,
                                          ".rds"))
          }
        }
      }
    }
  }
}



