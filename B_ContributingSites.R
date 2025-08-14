# Code to plot how many sites contribute to each cell

rm(list = ls())

# User-defined variables
DaytimeFlags = c(FALSE) #, TRUE)
PhysicalFlags = c(FALSE) #, TRUE)
WindyFlags = c(FALSE) #,TRUE)

# Load libraries
library(tidyverse)
library(colorspace)
library(scales)
library(parallel)
library(ggh4x)
library(ggnewscale)
library(ggpubr)
library(patchwork)
library(ggtext)

# Define different things to look at
fluxes <- c("qh", "qle", "nee")
variables <- c("lai","lwdown","swdown","tair","vpd","wind","precip")
models = c('1lin_raw','3km27_raw','6km729lag_raw','RF_raw','LSTM_raw','CABLE','CABLE-POP', 'CHTESSEL_1', 'CLM5', 'GFDL', 'JULES_GL9', 'JULES_GL9_LAI', 'MATSIRO','NoahMP','ORCHIDEE2','ORCHIDEE3') # complete list of models in Experiment
benchmarks = c('1lin_raw','3km27_raw','LSTM_raw') # user-nominated benchmarks, as though in modelevaluation.org
EFMs = c('1lin_raw','3km27_raw',"6km729lag_raw", "RF_raw","LSTM_raw")
bestEFMs = c("6km729lag_raw", "RF_raw","LSTM_raw")
nBestEFMs <- length(bestEFMs)
nEFMs <- length(EFMs)
nLSMs <- length(models) - nEFMs

# All fluxes are identical
Flux <- "qle"

for (DaytimeFlag in DaytimeFlags){
  for (PhysicalFlag in PhysicalFlags){
    for (WindyFlag in WindyFlags){
      for (var2 in variables){

        flags <- paste0(ifelse(DaytimeFlag,"_Daytime",""),
                        ifelse(PhysicalFlag,"_Physical",""),
                        ifelse(WindyFlag,"_Windy",""))
        
        message(Sys.time(),
                ": BEGIN - ",
                var2,
                " - ",
                flags)
        
        plotname <- paste0("plots/FingerprintPlots/",
                           "ContributingSites",
                           flags,
                           "_",
                           var2,
                           ".png")
        
        if (file.exists(plotname)){
          message(" - Done already ")
        } else {
          # Get that mostly-pre-processed data
          Data <- readRDS(paste0("data/ExtractedBinnedData",
                                 flags,
                                 ".rds"))
          # Calculate the number of timesteps nice and easy
          Timesteps <- nrow(Data)
          
          # Make sure we have the slightly-more-pre-processed data available
          FluxSelector <- setdiff(fluxes, Flux)

          Data_modified <- Data %>% select(-contains("_eb"),
                                           -contains(FluxSelector),
                                           -ends_with("_raw"),
                                           -time) %>%
            rename_with(~ sub(paste0("\\.", Flux), "", .x)) %>%
            rename_with(~ sub("_raw", "", .x)) %>%
            select(-contains("_prediction")) %>%
            rename_with(~ sub("_error", "", .x))
          
          # Now we can plot
          # 
          # Let's make the outputs a bit nicer 
          xLabels <- list(Var = c("<br>Downwelling<br>Longwave<br>(W/m<sup>2</sup>)",
                                  "<br>Specific<br>Humidity<br>(kg/kg)", 
                                  "<br>Downwelling <br>Shortwave<br>(W/m<sup>2</sup>)",
                                  "<br>Air<br>Temperature<br>(&deg;C)",
                                  "<br>Wind<br>Speed<br>(m/s)",
                                  "<br>LAI<br>(m<sup>2</sup>/m<sup>2</sup>)", 
                                  "<br>VPD<br>(kPa)",
                                  "<br>Rainfall<br>(mm/s)"))
          
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
          for(var1 in variables[variables != var2]){
            i = i + 1
            
            # Nice labels and names
            varName1 = xLabels$Var[which(c("lwdown", "qair", "swdown", "tair", "wind", "lai", "vpd", "precip") == var1)]
            varName2 = yLabels$Var[which(c("lwdown", "qair", "swdown", "tair", "wind", "lai", "vpd", "precip") == var2)]
            metric_names = c("Loss" = "LSM Loss Ratio",
                             "BinNo" = "Timesteps in Bin")
    
            # Calculate the percentage of timesteps in each bin ready for the plot
            sitedata <- Data_modified %>% rename(var1 = {{var1}},
                                                     var2 = {{var2}}) %>%
              mutate(Model = "bold(Contributing~Sites)") %>%
              group_by(Model, var1, var2) %>%
              summarise(Sites = n_distinct(site, na.rm = TRUE),
                        .groups = "drop")
            
            # Whack it all together
            plotdata <- sitedata
            
            plot <- ggplot(plotdata) +
              geom_tile(aes(x = var1, y = var2, fill = Sites)) +
              scale_fill_binned_sequential(name = "Contributing Sites",
                                           palette = 'YlGnBu', 
                                           rev = TRUE,
                                           breaks = c(1,2,5,10,20,50,100,153),
                                           limits = c(1,153),
                                           transform = "log") +
              scale_x_discrete(varName1,
                               breaks = levels(Data[[var1]])[c(T, rep(F, 9))],
                               drop = FALSE) +
              scale_y_discrete(varName2,
                               breaks = levels(Data[[var2]])[c(T, rep(F, 9))],
                               drop = FALSE) +
              theme_bw() +
              theme(text = element_text(size = 28),
                    axis.text.x = element_text(angle = 270, vjust = 0.2, hjust = 0),
                    legend.key.width = unit(4,"cm"),
                    legend.position = "bottom",
                    legend.title.position = "top",
                    legend.title = element_text(hjust = 0.5),
                    plot.title = element_text(hjust = 0.5),
                    axis.title.y=element_markdown(angle=90, vjust = 0.5),
                    axis.title.x = element_markdown()) +
              theme(strip.background = element_blank(),
                    strip.text.x = element_blank()) +
              facet_wrap(. ~ Model)
            
            # Assign the plot a name ready for combining
            assign(paste0("plot",i), plot)
          }
          
          # Use get_legend so that we can have a single consistent legend
          Site_legend <- get_legend(plot1)
          
          # Arrange the plots
          # First we add the actual plots together using wrap_plots
          # This ensures the plot areas align despite plot1 having the facet names
          arrangedPlot_nolegend <- wrap_plots(list(plot1 + rremove("legend"),
                                                   plot2 + rremove("legend") + rremove("ylab") + rremove("y.text"),
                                                   plot3 + rremove("legend") + rremove("ylab") + rremove("y.text"), 
                                                   plot4 + rremove("legend") + rremove("ylab") + rremove("y.text"),
                                                   plot5 + rremove("legend") + rremove("ylab") + rremove("y.text"),
                                                   plot6 + rremove("legend") + rremove("ylab") + rremove("y.text")), nrow = 1, align = "hv")
          # Having removed the legends, we add one of each back in!
          arrangedPlot <- ggarrange(Site_legend, arrangedPlot_nolegend, nrow = 2, heights = c(1, 7))
          
          # Save our glorious creation
          ggsave(filename = plotname,
                 arrangedPlot, width = 60, height = 30, create.dir = T, bg = "white", units = "cm", limitsize = FALSE)
          
          rm(plotdata,arrangedPlot, arrangedPlot_nolegend,
             plot, plot1, plot2, plot3, plot4, plot5, plot6)
        }
      }
    }
  }
}


