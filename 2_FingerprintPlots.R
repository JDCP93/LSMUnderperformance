# Code to group the data produced from 1_BinnedDataExtraction.R and plot it as 
# 2-dimensional fingerprints

rm(list = ls())

# User-defined variables
DaytimeFlags = c(FALSE, TRUE)
PhysicalFlags = c(FALSE, TRUE)
WindyFlags = c(FALSE, TRUE)

# Load the libraries (these may not all be needed)
library(tidyverse)
library(colorspace)
library(scales)
library(parallel)
library(ggh4x)
library(ggnewscale)
library(ggpubr)
library(patchwork)
library(ggtext)
library(data.table)

# Define fluxes, models, etc. etc.
fluxes <- c("qh", "qle", "nee")
variables <- c("lwdown","swdown","wind","vpd","precip", "lai","tair")
models = c('1lin_raw','3km27_raw','6km729lag_raw','RF_raw','LSTM_raw','CABLE','CABLE-POP', 'CHTESSEL_1', 'CLM5', 'GFDL', 'JULES_GL9', 'JULES_GL9_LAI', 'MATSIRO','NoahMP','ORCHIDEE2','ORCHIDEE3') # complete list of models in Experiment
benchmarks = c('1lin_raw','3km27_raw','LSTM_raw') # user-nominated benchmarks, as though in modelevaluation.org
EFMs = c('1lin_raw','3km27_raw',"6km729lag_raw", "RF_raw","LSTM_raw")
bestEFMs = c("6km729lag", "RF","LSTM")
nBestEFMs <- length(bestEFMs)
nEFMs <- length(EFMs)
nLSMs <- length(models) - nEFMs

# Run it through for each flag combo and met variable you wish to plot on x-axis
# This is in a for loop since it was ran locally on a computer with limited resources
# Feel free to parallelise if resources are available.
for (DaytimeFlag in DaytimeFlags){
  for (PhysicalFlag in PhysicalFlags){
    for (WindyFlag in WindyFlags){
      for (var1 in variables){
        for (Flux in fluxes){
          
          flags <- paste0(ifelse(DaytimeFlag,"_Daytime",""),
                          ifelse(PhysicalFlag,"_Physical",""),
                          ifelse(WindyFlag,"_Windy",""))
          gc()
          message(Sys.time(),
                  ": BEGIN - ",
                  Flux,
                  flags,
                  "_",
                  var1)
        
          # We check save files don't already exist
          lossdataname <- paste0("data/LossByModelAndVars_",
                                 Flux,
                                 flags,
                                 ".rds")
          
          plotname <- paste0("plots/FingerprintPlots/",
                             Flux,
                             flags,
                             "_",
                             var1,
                             ".png")
          
          if (file.exists(lossdataname) & file.exists(plotname)){
            
            message(" - Done already ")
            
          } else {
            
            # Get that mostly-pre-processed data
            Data <- readRDS(paste0("data/ExtractedBinnedData",
                                   flags,".rds"))
            
            # Calculate the number of timesteps nice and easy
            Timesteps <- nrow(Data)
            
            # Make sure we have the slightly-more-pre-processed data available
            FluxSelector <- setdiff(fluxes, Flux)
            FluxTitle <- c("Qh","Qle","NEE")[Flux==fluxes]
            
            # Calculate the EFM prediction errors
            Data_modified <- Data %>% select(-contains("_eb"),
                                        -contains(FluxSelector),
                                        -ends_with("_raw"),
                                        -time,
                                        -site) %>%
              rename_with(~ sub(paste0("\\.", Flux), "", .x)) %>%
              rename_with(~ sub("_raw", "", .x)) %>%
              mutate(
                MeanEFM_prediction = rowMeans(select(., `6km729lag_prediction`, RF_prediction, LSTM_prediction)),
                MeanEFM_error = MeanEFM_prediction - obs
                ) %>% 
              select(-contains("_prediction")) %>%
              rename_with(~ sub("_error", "", .x))
            
            # For each LSM, check if its error is bigger than all of the good benchmarks
            for (LSM in setdiff(models,EFMs)){
              LSMLoss_name <- paste0("Loss_", LSM)
              Data_modified[[LSMLoss_name]] <- abs(Data_modified[[LSM]]) > do.call(pmax, abs(Data_modified[, bestEFMs]))
            }
            
            setDT(Data_modified)       
            
            if (!file.exists(lossdataname)){
              # We will need this data later 
              AllVarGroupedData_DT <- data.table::melt(Data_modified,
                                           measure.vars = patterns("^Loss_"),
                                           variable.name = "Model",
                                           value.name = "ModelLoss",
                                           variable.factor = FALSE
                                          )[, Model := sub("^Loss_","",Model)
                                            ][, .(
                                              Loss = mean(ModelLoss, na.rm = TRUE) * 100,
                                              Count = sum(!is.na(ModelLoss))
                                            ),
                                            by = c("Model", variables)
                                            ]
    
              setDF(AllVarGroupedData_DT)
              
              saveRDS(AllVarGroupedData_DT, 
                       file = lossdataname)
              
              rm(AllVarGroupedData_DT)
            }
          
            # Now we can plot
            # Let's make the outputs a bit nicer 
            xLabels <- list(Var = c("<br>Downwelling Longwave (W/m<sup>2</sup>)",
                                    "<br>Specific Humidity (kg/kg)", 
                                    "<br>Downwelling Shortwave (W/m<sup>2</sup>)",
                                    "<br>Air Temperature (&deg;C)",
                                    "<br>Wind Speed (m/s)",
                                    "<br>LAI (m<sup>2</sup>/m<sup>2</sup>)", 
                                    "<br>VPD (kPa)",
                                    "<br>Rainfall (mm/s)"))
            
            yLabels <- list(Var = c("Downwelling<br>Longwave<br>(W/m<sup>2</sup>)",
                                    "Specific<br>Humidity<br>(kg/kg)", 
                                    "Downwelling <br>Shortwave<br>(W/m<sup>2</sup>)",
                                    "Air<br>Temperature<br>(&deg;C)",
                                    "Wind<br>Speed<br>(m/s)",
                                    "LAI<br>(m<sup>2</sup>/m<sup>2</sup>)", 
                                    "VPD<br>(kPa)",
                                    "Rainfall<br>(mm/s)"))
    
            # Start plot counter
            i = 0
            for(var2 in setdiff(variables, var1)){
              
              message(paste0("    ",var1,"-",var2," plot"))
              i = i + 1
                
                # Nice labels and names
                varName1 = xLabels$Var[which(c("lwdown", "qair", "swdown", "tair", "wind", "lai", "vpd", "precip") == var1)]
                varName2 = yLabels$Var[which(c("lwdown", "qair", "swdown", "tair", "wind", "lai", "vpd", "precip") == var2)]
                metric_names = c("Loss" = "LSM Loss Ratio",
                                 "BinNo" = "Timesteps in Bin")
                
                # Arrange loss data
                setnames(Data_modified, var1, "var1")
                setnames(Data_modified, var2, "var2")
                plotdata <- data.table::melt(Data_modified,
                                             measure.vars = patterns("^Loss_"),
                                             variable.name = "Model",
                                             value.name = "ModelLoss",
                                             variable.factor = FALSE
                )[, Model := sub("^Loss_","",Model)
                ][, .(
                  Loss = mean(ModelLoss, na.rm = TRUE) * 100
                ),
                by = c("Model", "var1", "var2")
                ]
                
                # Include the root mean square of LSM Losses as an aggregate
                RMSmodeldata <- plotdata[,
                                 .(Loss = sqrt(mean((Loss^2), na.rm = TRUE))),
                                 by = .(var1, var2)
                                 ][,Model := "bold(atop('All LSMs','(RMS)'))"]
                
                # Also useful if we can see where the majority of timesteps are located
                timestepdata <- Data_modified[,
                                              .(Loss = .N * 100 / Timesteps),
                                              by = .(Model = rep("bold(Timesteps)", Timesteps), var1, var2)]
                
                
                setnames(Data_modified, "var1", var1)
                setnames(Data_modified, "var2", var2)
                
                # Whack it all together
                plotdata <- rbind(plotdata, timestepdata, RMSmodeldata)
                
                plot <- ggplot(plotdata) +
                  geom_tile(aes(x = var1, y = var2, fill = Loss)) +
                  scale_fill_binned_divergingx(name = paste0("LSM Loss Ratio for ",FluxTitle," in Cell (%)"),
                                               palette = 'RdYlBu', 
                                               breaks = c(0, 5, 20, 30, 50, 80, 95, 100), 
                                               limits = c(0,100),
                                               mid = 25, 
                                               rev = TRUE) +
                  ggnewscale::new_scale_fill() +
                  geom_tile(data = timestepdata,
                            aes(x = var1,
                                y = var2,
                                fill = Loss)) +
                  scale_fill_continuous_sequential(palette = 'YlGnBu', 
                                                   guide = "none", 
                                                   limits = c(0,8), 
                                                   transform = "sqrt") + 
                  scale_x_discrete(varName1,
                                   breaks = levels(Data[[var1]])[c(T, rep(F, 9))],
                                   drop = FALSE) +
                  scale_y_discrete(varName2,
                                   breaks = levels(Data[[var2]])[c(T, rep(F, 9))],
                                   drop = FALSE) +
                  theme_bw() +
                  theme(text = element_text(size = 36),
                        axis.text.x = element_text(angle = 270, vjust = 0.2, hjust = 0),
                        legend.key.width = unit(10,"cm"),
                        legend.position = "bottom",
                        legend.title.position = "top",
                        legend.title = element_text(hjust = 0.5),
                        plot.title = element_text(hjust = 0.5),
                        axis.title.y=element_markdown(angle=0, vjust = 0.5),
                        axis.title.x = element_markdown())
                
                # If this is the first plot, we want the facet labels along the top
                if(i == 1){
                  plot <- plot +
                    facet_wrap(Model ~ ., nrow = 1, strip.position = "top", labeller = label_parsed)
                  
                } else {
                  # otherwise we want to get rid of them 
                  plot <- plot +
                    theme(strip.background = element_blank(),
                          strip.text.x = element_blank()) +
                    facet_wrap(Model ~ ., nrow = 1)
                }
                # Assign the plot a name ready for combining
                assign(paste0("plot",i), plot)
              }
            
              # Use get_legend so that we can have a single consistent legend
              LSM_Loss_legend <- get_legend(plot1)
              # We also need a timesteps legend too
              BinNo_legend <- get_legend(ggplot(data.frame(x = c(0,1),
                                                           y = c(0,1),
                                                           z = c(0,100))) +
                                           geom_tile(aes(x = x, y = y, fill = z)) +
                                           scale_fill_continuous_sequential(name = "Timesteps in Cell (%)",
                                                                            palette = 'YlGnBu',
                                                                            limits = c(0,8),
                                                                            transform = "sqrt") +
                                           theme(text = element_text(size = 36), 
                                                 legend.key.width = unit(10,"cm"),
                                                 legend.position = "bottom",
                                                 legend.title.position = "bottom",
                                                 legend.title = element_text(hjust = 0.5)))
              
              # Arrange the plots
              # First we add the actual plots together using wrap_plots
              # This ensures the plot areas align despite plot1 having the facet names
              arrangedPlot_nolegend <- wrap_plots(list(plot1 + rremove("legend") + rremove("x.title") + rremove("x.text"),
                                                       plot2 + rremove("legend") + rremove("x.title") + rremove("x.text"),
                                                       plot3 + rremove("legend") + rremove("x.title") + rremove("x.text"), 
                                                       plot4 + rremove("legend") + rremove("x.title") + rremove("x.text"),
                                                       plot5 + rremove("legend") + rremove("x.title") + rremove("x.text"),
                                                       plot6 + rremove("legend")), ncol = 1, align = "hv")
              # Having removed the legends, we add one of each back in!
              arrangedPlot <- ggarrange(LSM_Loss_legend, arrangedPlot_nolegend, BinNo_legend, nrow = 3, heights = c(1, 12, 1))
              
              # Save our glorious creation
              ggsave(filename = plotname,
                     arrangedPlot, width = 116, height = 72, create.dir = T, bg = "white", units = "cm", limitsize = FALSE)
          }
        }
      }
    }
  }
}
