# A script to check the distribution of fluxes, both from original data and also
# after various filters are applied

# The code is not pretty - much copy and paste

rm(list = ls())

# Load packages
library(tidyverse)
library(data.table)
library(dtplyr)
library(RColorBrewer)

# Define different things to look at
fluxes <- c("qh", "qle", "nee")
VarNames <- c("lai","lwdown","precip","qair","swdown","tair","vpd","wind")
ModelNames <- c('CABLE', 'CABLE-POP', 'CHTESSEL_1', 'CLM5', 'GFDL', 'JULES_GL9', 'JULES_GL9_LAI', 'MATSIRO','NoahMP','ORCHIDEE2','ORCHIDEE3') # these are all the locally forced LSMs in PLUMBER2

FullData <- readRDS("data/ExtractedBinnedData.rds") 
setDT(FullData)

PercentRemoved <- data.table()
RemovedQuantiles <- data.table()

for (Flux in fluxes){
  message(Flux)
  FluxSelector <- setdiff(fluxes, Flux)
  
  FullDataQuantiles <- FullData %>% 
    select(-contains("_eb")) %>% # Ignore energy corrected models
    select(-contains(FluxSelector)) %>% # Ignore other fluxes
    rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
    select(-contains("_error")) %>% # No longer need the error, only the predictions
    rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
    select(obs, any_of(c(VarNames))) %>% # Select what we need
    reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
            Q = seq(0,1,by=0.001)) %>%
    mutate(Flux = Flux,
           LRF = "None")
           
    FullDataQuantiles <- crossing(Model = c("All",ModelNames), FullDataQuantiles)
  

  WindyData <- readRDS("data/ExtractedBinnedData_Windy.rds")
  setDT(WindyData)
  NotWindyData <- FullData[!WindyData, on=c("time","site")]
  PercentRemoved <- rbind(PercentRemoved,
                          data.table(Model = "All",
                                     Flux = Flux,
                                     LRF = "Windy",
                                     PercentRemoved = nrow(NotWindyData)/nrow(FullData)))
  NotWindyQuantiles <- NotWindyData %>% 
    select(-contains("_eb")) %>% # Ignore energy corrected models
    select(-contains(FluxSelector)) %>% # Ignore other fluxes
    rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
    select(-contains("_error")) %>% # No longer need the error, only the predictions
    rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
    select(obs, any_of(c(VarNames))) %>% # Select what we need
    reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
              Q = seq(0,1,by=0.001)) %>%
    mutate(Flux = Flux,
           LRF = "Windy",
           Model = "All")
  rm(WindyData, NotWindyData)
  
  DaytimeData <- readRDS("data/ExtractedBinnedData_Daytime.rds")
  setDT(DaytimeData)
  NotDaytimeData <- FullData[!DaytimeData, on=c("time","site")]
  PercentRemoved <- rbind(PercentRemoved,
                          data.table(Model = "All",
                                     Flux = Flux,
                                     LRF = "Daytime",
                                     PercentRemoved = nrow(NotDaytimeData)/nrow(FullData)))
  NotDaytimeQuantiles <- NotDaytimeData %>% 
    select(-contains("_eb")) %>% # Ignore energy corrected models
    select(-contains(FluxSelector)) %>% # Ignore other fluxes
    rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
    select(-contains("_error")) %>% # No longer need the error, only the predictions
    rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
    select(obs, any_of(c(VarNames))) %>% # Select what we need
    reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
              Q = seq(0,1,by=0.001)) %>%
    mutate(Flux = Flux,
           LRF = "Daytime",
           Model = "All")
  rm(DaytimeData, NotDaytimeData)
  
  
  PhysicalData <- readRDS("data/ExtractedBinnedData_Physical.rds")
  setDT(PhysicalData)
  NotPhysicalData <- FullData[!PhysicalData, on=c("time","site")]
  PercentRemoved <- rbind(PercentRemoved,
                          data.table(Model = "All",
                                     Flux = Flux,
                                     LRF = "Physical",
                                     PercentRemoved = nrow(NotPhysicalData)/nrow(FullData)))
  
  NotPhysicalQuantiles <- NotPhysicalData %>% 
    select(-contains("_eb")) %>% # Ignore energy corrected models
    select(-contains(FluxSelector)) %>% # Ignore other fluxes
    rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
    select(-contains("_error")) %>% # No longer need the error, only the predictions
    rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
    select(obs, any_of(c(VarNames))) %>% # Select what we need
    reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
              Q = seq(0,1,by=0.001)) %>%
    mutate(Flux = Flux,
           LRF = "Physical",
           Model = "All")
  rm(PhysicalData, NotPhysicalData)
  gc()
  
  LossData <- readRDS(paste0("data/LossByModelAndVars_",
                     Flux,".rds"))
  setDT(LossData)
  LossData_50_RemovedQuantiles <- data.table()
  LossData_95_RemovedQuantiles <- data.table()
  for (ModelName in ModelNames){
    message(" - ",ModelName)
    LossData_50 <- LossData[Loss >= 50 & Model == ModelName]
    LossData_50_Removed <- FullData[LossData_50, on=VarNames]
    PercentRemoved <- rbind(PercentRemoved,
                             data.table(Model = ModelName,
                                        Flux = Flux,
                                        LRF = "50",
                                        PercentRemoved = nrow(LossData_50_Removed)/nrow(FullData)))
    
    LossData_50_RemovedQuantiles <- rbind(LossData_50_RemovedQuantiles,
                                          LossData_50_Removed %>% 
      select(-contains("_eb")) %>% # Ignore energy corrected models
      select(-contains(FluxSelector)) %>% # Ignore other fluxes
      rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
      select(-contains("_error")) %>% # No longer need the error, only the predictions
      rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
      select(obs, any_of(c(VarNames))) %>% # Select what we need
      reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
              Q = seq(0,1,by=0.001)) %>%
      mutate(Flux = Flux,
             LRF = "50",
             Model = ModelName))
    
    LossData_95 <- LossData[Loss >= 95 & Model == ModelName]
    LossData_95_Removed <- FullData[LossData_95, on=VarNames]
    PercentRemoved <- rbind(PercentRemoved,
                            data.table(Model = ModelName,
                                       Flux = Flux,
                                       LRF = "95",
                                       PercentRemoved = nrow(LossData_95_Removed)/nrow(FullData)))
    
    LossData_95_RemovedQuantiles <- rbind(LossData_95_RemovedQuantiles,
                                          LossData_95_Removed %>% 
                                            select(-contains("_eb")) %>% # Ignore energy corrected models
                                            select(-contains(FluxSelector)) %>% # Ignore other fluxes
                                            rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
                                            select(-contains("_error")) %>% # No longer need the error, only the predictions
                                            rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
                                            select(obs, any_of(c(VarNames))) %>% # Select what we need
                                            reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
                                                    Q = seq(0,1,by=0.001)) %>%
                                            mutate(Flux = Flux,
                                                   LRF = "95",
                                                   Model = ModelName))
  }
  
  
  RemovedQuantiles <- rbind(RemovedQuantiles,
                     FullDataQuantiles,
                     NotPhysicalQuantiles,
                     NotDaytimeQuantiles,
                     NotWindyQuantiles,
                     LossData_95_RemovedQuantiles,
                     LossData_50_RemovedQuantiles)
}

RemovedQuantiles <- RemovedQuantiles %>%
  mutate(LRF = factor(LRF, levels = c("None", "50", "95", "Physical", "Windy", "Daytime")))


RemovedBoxPlots <- RemovedQuantiles %>% mutate(Value = ifelse(Model %in% c("CLM5", "JULES_GL9", "JULES_GL9_LAI", "MATSIRO") & Flux == "nee", NA, Value)) %>%
  filter(LRF %in% c("None", "Physical", "Daytime", "Windy", "95", "50")) %>%
  pivot_wider(names_from = Q, values_from = Value, names_prefix = "Q") %>%
  ggplot(aes(x = LRF)) +
  geom_point(aes(y = Q0)) +
  geom_point(aes(y = Q1)) +
  geom_errorbar(aes(ymin = Q0.01, ymax = Q0.99)) +
  geom_crossbar(aes(ymin = Q0.1,y = Q0.5, ymax = Q0.9, fill = LRF), width = 0.5) +
  geom_crossbar(aes(ymin = Q0.25,y = Q0.5, ymax = Q0.75, fill = LRF)) +
  scale_y_continuous(transform = "pseudo_log") +
  facet_grid(Flux ~ Model, scales = "free", space = "free_x") +
  scale_fill_manual(values = c("grey50",brewer.pal(7,"RdYlBu")[c(1,2,5,6,7)])) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Removed Fluxes",
       subtitle = "Distribution of flux values removed by filters")



RemovedCDFs <- RemovedQuantiles %>% mutate(Value = ifelse(Model %in% c("CLM5", "JULES_GL9", "JULES_GL9_LAI", "MATSIRO") & Flux == "nee", NA, Value)) %>%
  ggplot() +
  geom_line(aes(x = Value, y = Q, color = LRF)) +
  facet_grid(Model ~ Flux, scales = "free") +
  scale_color_manual(values = c("grey50",brewer.pal(7,"RdYlBu")[c(1,2,5,6,7)])) +
  scale_x_continuous(transform = "pseudo_log") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  labs(title = "Removed Fluxes",
       subtitle = "Distribution of flux values removed by filters")


PercentRemaining<- data.table()
RemainingQuantiles <- data.table()

for (Flux in fluxes){
  message(Flux)
  FluxSelector <- setdiff(fluxes, Flux)
  
  FullDataQuantiles <- FullData %>% 
    select(-contains("_eb")) %>% # Ignore energy corrected models
    select(-contains(FluxSelector)) %>% # Ignore other fluxes
    rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
    select(-contains("_error")) %>% # No longer need the error, only the predictions
    rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
    select(obs, any_of(c(VarNames))) %>% # Select what we need
    reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
            Q = seq(0,1,by=0.001)) %>%
    mutate(Flux = Flux,
           LRF = "None")
  
  FullDataQuantiles <- crossing(Model = c("All",ModelNames), FullDataQuantiles)
  
  
  WindyData <- readRDS("data/ExtractedBinnedData_Windy.rds")
  setDT(WindyData)
  WindyQuantiles <- WindyData %>% 
    select(-contains("_eb")) %>% # Ignore energy corrected models
    select(-contains(FluxSelector)) %>% # Ignore other fluxes
    rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
    select(-contains("_error")) %>% # No longer need the error, only the predictions
    rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
    select(obs, any_of(c(VarNames))) %>% # Select what we need
    reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
            Q = seq(0,1,by=0.001)) %>%
    mutate(Flux = Flux,
           LRF = "Windy",
           Model = "All")
  PercentRemaining <- rbind(PercentRemaining,
                            data.table(Model = "All",
                                       Flux = Flux,
                                       LRF = "Windy",
                                       PercentRemaining = nrow(WindyData)/nrow(FullData)))
  rm(WindyData)

  DaytimeData <- readRDS("data/ExtractedBinnedData_Daytime.rds")
  setDT(DaytimeData)
  DaytimeQuantiles <- DaytimeData %>% 
    select(-contains("_eb")) %>% # Ignore energy corrected models
    select(-contains(FluxSelector)) %>% # Ignore other fluxes
    rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
    select(-contains("_error")) %>% # No longer need the error, only the predictions
    rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
    select(obs, any_of(c(VarNames))) %>% # Select what we need
    reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
            Q = seq(0,1,by=0.001)) %>%
    mutate(Flux = Flux,
           LRF = "Daytime",
           Model = "All")
  PercentRemaining <- rbind(PercentRemaining,
                            data.table(Model = "All",
                                       Flux = Flux,
                                       LRF = "Daytime",
                                       PercentRemaining = nrow(DaytimeData)/nrow(FullData)))
  rm(DaytimeData)

  
  PhysicalData <- readRDS("data/ExtractedBinnedData_Physical.rds")
  setDT(PhysicalData)
  PhysicalQuantiles <- PhysicalData %>% 
    select(-contains("_eb")) %>% # Ignore energy corrected models
    select(-contains(FluxSelector)) %>% # Ignore other fluxes
    rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
    select(-contains("_error")) %>% # No longer need the error, only the predictions
    rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
    select(obs, any_of(c(VarNames))) %>% # Select what we need
    reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
            Q = seq(0,1,by=0.001)) %>%
    mutate(Flux = Flux,
           LRF = "Physical",
           Model = "All")
  PercentRemaining <- rbind(PercentRemaining,
                            data.table(Model = "All",
                                       Flux = Flux,
                                       LRF = "Physical",
                                       PercentRemaining = nrow(PhysicalData)/nrow(FullData)))
  rm(PhysicalData)
  gc()
  
  LossData <- readRDS(paste0("data/LossByModelAndVars_",
                             Flux, 
                             ".rds"))
  setDT(LossData)
  LossData_50_RemainingQuantiles <- data.table()
  LossData_95_RemainingQuantiles <- data.table()
  for (ModelName in ModelNames){
    message(" - ",ModelName)
    LossData_50 <- LossData[Loss >= 50 & Model == ModelName]
    LossData_50_Remaining <- FullData[!LossData_50, on=VarNames]

    LossData_50_RemainingQuantiles <- rbind(LossData_50_RemainingQuantiles,
                                          LossData_50_Remaining %>% 
                                            select(-contains("_eb")) %>% # Ignore energy corrected models
                                            select(-contains(FluxSelector)) %>% # Ignore other fluxes
                                            rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
                                            select(-contains("_error")) %>% # No longer need the error, only the predictions
                                            rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
                                            select(obs, any_of(c(VarNames))) %>% # Select what we need
                                            reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
                                                    Q = seq(0,1,by=0.001)) %>%
                                            mutate(Flux = Flux,
                                                   LRF = "50",
                                                   Model = ModelName))
    PercentRemaining <- rbind(PercentRemaining,
                              data.table(Model = ModelName,
                                         Flux = Flux,
                                         LRF = "50",
                                         PercentRemaining = nrow(LossData_50_Remaining)/nrow(FullData)))
    
    LossData_95 <- LossData[Loss >= 95 & Model == ModelName]
    LossData_95_Remaining <- FullData[!LossData_95, on=VarNames]
    LossData_95_RemainingQuantiles <- rbind(LossData_95_RemainingQuantiles,
                                          LossData_95_Remaining %>% 
                                            select(-contains("_eb")) %>% # Ignore energy corrected models
                                            select(-contains(FluxSelector)) %>% # Ignore other fluxes
                                            rename_with(~ sub(paste0(".",Flux), "", .x), everything()) %>% # Rename for ease
                                            select(-contains("_error")) %>% # No longer need the error, only the predictions
                                            rename_with(~ sub("_prediction", "", .x), everything()) %>% # Rename predictions 
                                            select(obs, any_of(c(VarNames))) %>% # Select what we need
                                            reframe(Value = quantile(obs, seq(0,1,by=0.001), na.rm = TRUE, names = F),
                                                    Q = seq(0,1,by=0.001)) %>%
                                            mutate(Flux = Flux,
                                                   LRF = "95",
                                                   Model = ModelName))
    PercentRemaining <- rbind(PercentRemaining,
                            data.table(Model = ModelName,
                                       Flux = Flux,
                                       LRF = "95",
                                       PercentRemaining = nrow(LossData_95_Remaining)/nrow(FullData)))
  }
  
  
  RemainingQuantiles <- rbind(RemainingQuantiles,
                               FullDataQuantiles,
                               PhysicalQuantiles,
                               DaytimeQuantiles,
                               WindyQuantiles,
                               LossData_95_RemainingQuantiles,
                               LossData_50_RemainingQuantiles)
}

RemainingQuantiles <- RemainingQuantiles %>%
  mutate(LRF = factor(LRF, levels = c("None", "50", "95", "Physical", "Windy", "Daytime")))


RemainingBoxPlots <- RemainingQuantiles %>% mutate(Value = ifelse(Model %in% c("CLM5", "JULES_GL9", "JULES_GL9_LAI", "MATSIRO") & Flux == "nee", NA, Value)) %>%
  filter(LRF %in% c("None", "Physical", "Daytime", "Windy", "95", "50")) %>%
  pivot_wider(names_from = Q, values_from = Value, names_prefix = "Q") %>%
  ggplot(aes(x = LRF)) +
  geom_point(aes(y = Q0)) +
  geom_point(aes(y = Q1)) +
  geom_errorbar(aes(ymin = Q0.01, ymax = Q0.99)) +
  geom_crossbar(aes(ymin = Q0.1,y = Q0.5, ymax = Q0.9, fill = LRF), width = 0.5) +
  geom_crossbar(aes(ymin = Q0.25,y = Q0.5, ymax = Q0.75, fill = LRF)) +
  scale_y_continuous(transform = "pseudo_log") +
  facet_grid(Flux ~ Model, scales = "free", space = "free_x") +
  scale_fill_manual(values = c("grey50",brewer.pal(7,"RdYlBu")[c(1,2,5,6,7)])) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Remaining Fluxes",
       subtitle = "Distribution of flux values remaining after filters")



RemainingCDFs <-RemainingQuantiles %>% mutate(Value = ifelse(Model %in% c("CLM5", "JULES_GL9", "JULES_GL9_LAI", "MATSIRO") & Flux == "nee", NA, Value)) %>%
  ggplot() +
  geom_line(aes(x = Value, y = Q, color = LRF)) +
  facet_grid(Model ~ Flux, scales = "free") +
  scale_color_manual(values = c("grey50",brewer.pal(7,"RdYlBu")[c(1,2,5,6,7)])) +
  scale_x_continuous(transform = "pseudo_log") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  labs(title = "Remaining Fluxes",
       subtitle = "Distribution of flux values remaining after filters")

ggsave(filename = "plots/RemainingFluxDistribution.png",
       RemainingBoxPlots, width = 45, height = 30, create.dir = T, bg = "white", units = "cm", limitsize = FALSE)

ggsave(filename = "plots/RemovedFluxDistribution.png",
       RemovedBoxPlots, width = 45, height = 30, create.dir = T, bg = "white", units = "cm", limitsize = FALSE)
