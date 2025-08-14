# Code to group the data produced from 1_BinnedDataExtraction.R and plot it as 
# 2-dimensional fingerprints
rm(list = ls())

# Load libraries
library(ggplot2)
library(tidyverse)
library(reshape2)

# Load the loss data
LossFiles <- list.files(path = "data/", pattern = "LossByModelAndVars")


for (File in LossFiles){
  message(File)
  # Attributes are saved in the file name
  attr <- unlist(str_split(File, pattern ="[[:punct:]]"))
  LossData_raw <- readRDS(paste0("data/",File))
  
  LossData_part <- LossData_raw %>% 
    mutate(Flux = attr[2],
           Daytime = ifelse("Daytime" %in% attr, TRUE, FALSE),
           Physical = ifelse("Physical" %in% attr, TRUE, FALSE),
           Windy = ifelse("Windy" %in% attr, TRUE, FALSE),
           LossToNearestInteger = floor(Loss)) %>% 
    group_by(Model, Flux, Daytime, Physical, Windy, LossToNearestInteger) %>%
    summarise(Count = sum(Count),
              .groups = "drop_last") %>%
    arrange(desc(LossToNearestInteger)) %>%
    mutate(CumCount = cumsum(Count),
           Percent = CumCount / max(CumCount) * 100) %>%
    # Loss data has some NaNs
    mutate(LossToNearestInteger = ifelse(is.na(LossToNearestInteger), NA, as.character(LossToNearestInteger))) %>%
    drop_na(LossToNearestInteger) %>%
    mutate(LossToNearestInteger = factor(LossToNearestInteger, levels = c(seq(0,100, by = 1), "None"))) 
  
  # Add the No filter data
  LossData_part <- rbind(LossData_part,
                    data.frame(Model = unique(LossData_part$Model),
                               Flux = attr[2],
                               Daytime = ifelse("Daytime" %in% attr, TRUE, FALSE),
                               Physical = ifelse("Physical" %in% attr, TRUE, FALSE),
                               Windy = ifelse("Windy" %in% attr, TRUE, FALSE),
                               LossToNearestInteger = factor("None", levels = c(seq(0,100, by = 1), "None")),
                               Count = 0,
                               CumCount = 0,
                               Percent = 0))
  
  # Bind if this isn't the first time, otherwise just start the df
  if (exists("LossData")){
    LossData <- rbind(LossData, as.data.frame(LossData_part, row.names = NULL))
  } else {
    LossData <- as.data.frame(LossData_part, row.names = NULL)
  }
}
# Make the flux names pretty and plot in the right order
LossData <- LossData %>%
  mutate(Flux = case_when(Flux == "nee" ~ "NEE",
                          Flux == "qh" ~ "Qh",
                          Flux == "qle" ~ "Qle"),
         Flux = factor(Flux, levels = c("Qh", "Qle", "NEE")))
# We'll need this again
saveRDS(LossData, file = "data/PercentRemovedPerLossFilter.rds")