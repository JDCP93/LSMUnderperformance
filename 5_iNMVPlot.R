# Code to calculate the iNMV and plot it

rm(list = ls())

# Load the libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(gtable)
library(grid)
library(cowplot)

# If we change the number of filters we want to plot, we also need to change the 
# manual color scales in the plots below.
FiltersToPlot = c("50", "95", "Physical", "Windy", "Daytime", "None")

MetricFiles <- list.files(path = "data/Metrics", pattern = "SiteModelMetrics")

for (File in MetricFiles){
  attr <- unlist(str_split(File, pattern ="[[:punct:]]"))
  metrics <- readRDS(paste0("data/Metrics/",File))
  if (exists("Metrics")){
    Metrics <- rbind(Metrics, as.data.frame(metrics, row.names = NULL))
  } else {
    Metrics <- as.data.frame(metrics, row.names = NULL)
  }
}

# Make the metrics dataframe nice
Metrics <- Metrics %>% distinct()
Metrics$LossFilter <- factor(Metrics$LossFilter, levels = c("5", "20", "50", "80", "95", "None"))
Metrics <- Metrics %>%  mutate(Flux = case_when(Flux == "nee" ~ "NEE",
                                                Flux == "qh" ~ "Qh",
                                                Flux == "qle" ~ "Qle"),
                               Flux = factor(Flux, levels = c("Qh", "Qle", "NEE")))

Metric_Names <-  c("Correlation", 
                   "Q5Diff",
                   "Q95Diff", 
                   "MBE",
                   "NME", 
                   "SDdiff",
                   'PDFOverlap')

# Define the models to use for calculating min and max
Benchmarks <- c('1lin_raw', '3km27_raw', 'LSTM_raw')

# We only want to be comparing the LSMs to the original, raw benchmark model metrics

# Add benchmark values for different LRF values
for (Benchmark in Benchmarks){
  
  # First, remove any Physical=TRUE, Daytime=TRUE or Windy=TRUE for the benchmarks
  Metrics <- Metrics %>% filter(!(Model == Benchmark & (Daytime == TRUE | Physical == TRUE | Windy == TRUE)))
  
  Metrics[Metrics$Model == Benchmark & Metrics$LossFilter == "5", Metric_Names] <- Metrics[Metrics$Model == Benchmark & Metrics$LossFilter == "None", Metric_Names]
  Metrics[Metrics$Model == Benchmark & Metrics$LossFilter == "20", Metric_Names] <- Metrics[Metrics$Model == Benchmark & Metrics$LossFilter == "None", Metric_Names]
  Metrics[Metrics$Model == Benchmark & Metrics$LossFilter == "50", Metric_Names] <- Metrics[Metrics$Model == Benchmark & Metrics$LossFilter == "None", Metric_Names]
  Metrics[Metrics$Model == Benchmark & Metrics$LossFilter == "80", Metric_Names] <- Metrics[Metrics$Model == Benchmark & Metrics$LossFilter == "None", Metric_Names]
  Metrics[Metrics$Model == Benchmark & Metrics$LossFilter == "95", Metric_Names] <- Metrics[Metrics$Model == Benchmark & Metrics$LossFilter == "None", Metric_Names]
  
  # Add back in for Daytime=TRUE, Physical=TRUE and Windy=TRUE
  df_Daytime <- Metrics[Metrics$Model == Benchmark, ] %>% mutate(Daytime = TRUE)
  df_Physical <- Metrics[Metrics$Model == Benchmark, ] %>% mutate(Physical = TRUE)
  df_Windy <- Metrics[Metrics$Model == Benchmark, ] %>% mutate(Windy = TRUE)
  df_Daytime_Physical <- Metrics[Metrics$Model == Benchmark, ] %>% mutate(Daytime = TRUE,
                                                                              Physical = TRUE)
  df_Daytime_Windy <- Metrics[Metrics$Model == Benchmark, ] %>% mutate(Daytime = TRUE,
                                                                              Windy = TRUE)
  df_Windy_Physical <- Metrics[Metrics$Model == Benchmark, ] %>% mutate(Windy = TRUE,
                                                                              Physical = TRUE)
  df_Daytime_Windy_Physical <- Metrics[Metrics$Model == Benchmark, ] %>% mutate(Daytime = TRUE,
                                                                                Windy = TRUE,
                                                                        Physical = TRUE)
  
  Metrics <- rbind(Metrics,
                   df_Daytime, df_Windy, df_Physical,
                   df_Daytime_Physical, df_Daytime_Windy, df_Windy_Physical,
                   df_Daytime_Windy_Physical)
}




# Functions to normalise a metric
normalise_metric_lowerBetter <- function(metric, min_val, max_val) {
  result <- (metric - min_val) / (max_val - min_val)
  ifelse(abs(result) == Inf, NA, result)
}
# Functions to normalise a metric
normalise_metric_higherBetter <- function(metric, min_val, max_val) {
  result <- (max_val - metric) / (max_val - min_val)
  ifelse(abs(result) == Inf, NA, result)
}

# Apply normalisation
Metrics_normalised <- Metrics %>%
  group_by(Flux, Site, LossFilter, Daytime, Physical, Windy) %>%
  mutate(across(
    all_of(c('Q5Diff', 'Q95Diff', 'MBE', 'NME', 'SDdiff')),
    ~ normalise_metric_lowerBetter(
      .,
      min(.[Model %in% Benchmarks], na.rm = TRUE),
      max(.[Model %in% Benchmarks], na.rm = TRUE)
    ),
    .names = "{.col}"
  ),
  across(
    all_of(c('Correlation', 'PDFOverlap')),
    ~ normalise_metric_higherBetter(
      .,
      min(.[Model %in% Benchmarks], na.rm = TRUE),
      max(.[Model %in% Benchmarks], na.rm = TRUE)
    ),
    .names = "{.col}"
  ),
  LossFilter = factor(LossFilter)) %>%
  ungroup()

# Calculate the iNMV, remove combinations of filters that we don't care about,
# and make it pretty for easier plotting
plot_df <- Metrics_normalised %>% 
  distinct() %>%
  filter(!((!LossFilter == "None") & Physical == TRUE)) %>%
  filter(!((!LossFilter == "None") & Daytime == TRUE)) %>%
  filter(!((!LossFilter == "None") & Windy == TRUE)) %>%
  filter(! (Physical == TRUE & Daytime == TRUE)) %>% 
  filter(! (Physical == TRUE & Windy == TRUE)) %>% 
  filter(! (Windy == TRUE & Daytime == TRUE)) %>% 
  pivot_longer(cols = Correlation:PDFOverlap, names_to = "Metric", values_to = "Value") %>%
  group_by(Flux, Model, LossFilter, Daytime, Physical, Windy) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  pivot_wider(names_from = Model, values_from = Value) %>%
  pivot_longer(cols = c(`6km729lag_raw`, "CABLE", "CABLE-POP", "CHTESSEL_1", "CLM5", "GFDL", "JULES_GL9", "JULES_GL9_LAI", "MATSIRO", "NoahMP", "ORCHIDEE2", "ORCHIDEE3", "RF_raw"), names_to = "Model") %>%
  filter(!(Model %in% c("6km729lag_raw", "RF_raw"))) %>%
  mutate(LossFilter = case_when(Physical == FALSE ~ LossFilter,
                                Physical == TRUE ~ "Physical"),
         LossFilter = case_when(Windy == FALSE ~ LossFilter,
                                Windy == TRUE ~ "Windy"),
         LossFilter = case_when(Daytime == FALSE ~ LossFilter,
                                Daytime == TRUE ~ "Daytime"),
         LossFilter = factor(LossFilter, levels = c("5", "20", "50", "80", "95","Physical","Windy","Daytime","None")))

# We now bring in the timesteps removed per loss filter data and merge it with the plot data 
# This means we can incorporate it into the plot

# Read data file from "local_PercentRemoved.R"
LossFilterPercentRemoved <- readRDS("data/PercentRemovedPerLossFilter.rds") %>%
  mutate(LossFilter = as.character(LossToNearestInteger)) %>%
  filter(!((!LossFilter == "None") & Physical == TRUE)) %>%
  filter(!((!LossFilter == "None") & Daytime == TRUE)) %>%
  filter(!((!LossFilter == "None") & Windy == TRUE)) %>%
  filter(! (Physical == TRUE & Daytime == TRUE)) %>% 
  filter(! (Physical == TRUE & Windy == TRUE)) %>% 
  filter(! (Windy == TRUE & Daytime == TRUE)) %>% 
  mutate(LossFilter = case_when(Physical == FALSE ~ LossFilter,
                                Physical == TRUE ~ "Physical"),
         LossFilter = case_when(Windy == FALSE ~ LossFilter,
                                Windy == TRUE ~ "Windy"),
         LossFilter = case_when(Daytime == FALSE ~ LossFilter,
                                Daytime == TRUE ~ "Daytime"),
         LossFilter = factor(LossFilter, levels = c("5", "20", "50", "80", "95","Physical","Windy","Daytime","None")))

# We need to calculate the percentage removed for daytime and physical filters
ModelTimesteps_NoFilter <- LossFilterPercentRemoved %>%
  filter(Daytime == FALSE & Physical == FALSE & Windy == FALSE & Percent == "100")

ModelTimesteps_Physical <- readRDS("data/PercentRemovedPerLossFilter.rds") %>%
  filter(Daytime == FALSE & Physical == TRUE & Windy == FALSE & Percent == "100") %>%
  mutate(Percent = (ModelTimesteps_NoFilter$CumCount - CumCount) / ModelTimesteps_NoFilter$CumCount * 100)

ModelTimesteps_Daytime <- readRDS("data/PercentRemovedPerLossFilter.rds") %>%
  filter(Daytime == TRUE & Physical == FALSE & Windy == FALSE & Percent == "100") %>%
  mutate(Percent = (ModelTimesteps_NoFilter$CumCount - CumCount) / ModelTimesteps_NoFilter$CumCount * 100)

ModelTimesteps_Windy <- readRDS("data/PercentRemovedPerLossFilter.rds") %>%
  filter(Daytime == FALSE & Physical == FALSE & Windy == TRUE & Percent == "100") %>%
  mutate(Percent = (ModelTimesteps_NoFilter$CumCount - CumCount) / ModelTimesteps_NoFilter$CumCount * 100)

LossFilterPercentRemoved$Percent[LossFilterPercentRemoved$Daytime == TRUE & LossFilterPercentRemoved$Physical == FALSE & LossFilterPercentRemoved$Windy == FALSE] <- ModelTimesteps_Daytime$Percent
LossFilterPercentRemoved$Percent[LossFilterPercentRemoved$Daytime == FALSE & LossFilterPercentRemoved$Physical == TRUE & LossFilterPercentRemoved$Windy == FALSE] <- ModelTimesteps_Physical$Percent
LossFilterPercentRemoved$Percent[LossFilterPercentRemoved$Daytime == FALSE & LossFilterPercentRemoved$Physical == FALSE & LossFilterPercentRemoved$Windy == TRUE] <- ModelTimesteps_Windy$Percent

plot_df_withLoss <- merge(plot_df, LossFilterPercentRemoved, all.x = TRUE)

# Plot with the data proportion information
plot_withLoss <- plot_df_withLoss %>%
  filter(LossFilter %in% FiltersToPlot) %>%
  ggplot() +
  geom_crossbar(data = plot_df_withLoss %>% filter(LossFilter %in% FiltersToPlot[-length(FiltersToPlot)]), aes(x = Flux, ymin = -1, y = -1, ymax = -1 + Percent/100, group = LossFilter, fill = LossFilter, color = LossFilter), position = "dodge", alpha = 0.5, width = 0.75, linewidth = 0) +
  geom_text(data = plot_df_withLoss %>% filter(LossFilter %in% FiltersToPlot[-length(FiltersToPlot)]), aes(x = Flux, y = -1.1, group = LossFilter, label = round(Percent)), color = "black", position = position_dodge(width = 0.75), size = 3) +
  geom_line(aes(x = Flux, y = `1lin_raw`, group = "'1lin_raw'"), color = "grey70", linewidth = 1.5) +
  geom_point(aes(x = Flux, y = `1lin_raw`, group = "'1lin_raw'"), color = "grey70", fill = "grey70", size = 5, shape = 25) +
  geom_line(aes(x = Flux, y = `3km27_raw`, group = "'3km27_raw'"), color = "grey50", linewidth = 1.5) +
  geom_point(aes(x = Flux, y = `3km27_raw`, group = "'3km27_raw'"), color = "grey50", fill = "grey50", size = 5, shape = 23) +
  geom_line(aes(x = Flux, y = `LSTM_raw`, group = "'LSTM_raw'"), color = "grey30", linewidth = 1.5) +
  geom_point(aes(x = Flux, y = `LSTM_raw`, group = "'LSTM_raw'"), color = "grey30", fill = "grey30", size = 5, shape = 24) +
  geom_line(aes(x = Flux, y = value, group = interaction(Model, LossFilter), color = LossFilter), linewidth = 1.5, alpha = 0.75) +
  geom_point(aes(x = Flux, y = value, group = interaction(Model, LossFilter), fill = LossFilter), size = 6, shape = 21, color = "black") +
  scale_x_discrete(expand = c(0,0.1)) +
  scale_y_continuous("") +
  scale_color_manual("Loss Ratio Filter (%)", values = c(brewer.pal(7, "RdYlBu")[c(1,2)], "black", brewer.pal(7, "RdYlBu")[c(5,6,7)]),
                     breaks = c("50", "95", "None", "Physical","Windy","Daytime")) +
  scale_fill_manual("Loss Ratio Filter (%)", values = c(brewer.pal(7, "RdYlBu")[c(1,2,5,6,7)], "black"), guide = "none") +
  facet_wrap(Model ~ ., scale = "free_y", nrow = 3) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.byrow = TRUE)

# We want a legend for the benchmarks
benchmark_legend <- get_legend(
  ggplot(data.frame(Benchmarks = Benchmarks,
                    Value = c(1,1,1))) +
  geom_point(aes(x = 1, 
                 y = Value, 
                 color = Benchmarks, 
                 shape = Benchmarks,
                 fill = Benchmarks),
             size = 3) +
  geom_line(aes(x = 1, 
                 y = Value, 
                 color = Benchmarks)) +
  scale_color_manual(values = c("grey70", "grey50", "grey30"), labels = c("1lin", "3km27", "LSTM")) +
  scale_fill_manual(values = c("grey70", "grey50", "grey30"), labels = c("1lin", "3km27", "LSTM")) +
  scale_shape_manual(values = c(25, 23, 24), labels = c("1lin", "3km27", "LSTM")) +
  theme(legend.key = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        legend.box.background = element_rect(colour = "black", linewidth = 0.5),
        legend.key.spacing.y = unit(0.5, 'cm'))
)

# Add the legend
plot_withLoss_Grob <- ggplotGrob(plot_withLoss)
plot_withLoss_BenchmarkLegend <- gtable_add_grob(plot_withLoss_Grob, benchmark_legend, nrow(plot_withLoss_Grob)-7, ncol(plot_withLoss_Grob)-6)
grid.draw(plot_withLoss_BenchmarkLegend)

# Save it
ggsave(filename = paste0("plots/PLUMBER2PlotwithLoss.png"),
       plot_withLoss_BenchmarkLegend, width = 30, height = 24, create.dir = T, bg = "white", units = "cm", limitsize = FALSE)

# We can also calculate some relative differences in the iNMV across each filter
Improvements <- plot_df %>% filter(LossFilter %in% FiltersToPlot) %>%
  ungroup() %>%
  select(-Daytime, -Windy, -Physical) %>%
  pivot_wider(names_from = LossFilter, values_from = value) %>%
  mutate(Imp50 = (None-`50`),
         Imp95 = (None-`95`),
         ImpWindy = (None-Windy),
         ImpPhysical= (None-Physical),
         ImpDaytime = None-Daytime,
         RelImp50 = (None-`50`)/None,
         RelImp95 = (None-`95`)/None,
         RelImpWindy = (None-Windy)/None,
         RelImpPhysical= (None-Physical)/None,
         RelImpDaytime = (None-Daytime)/None)


SummaryImprovements <- Improvements %>% group_by(Flux) %>%
  summarise(MeanRaw = mean(None, na.rm = T),
            Mean50 = mean(`50`, na.rm = T),
            Mean95 = mean(`95`, na.rm = T),
            MeanWindy = mean(Windy, na.rm = T),
            MeanPhysical = mean(Physical, na.rm = T),
            MeanDaytime = mean(Daytime, na.rm = T),
            MeanImp50 = mean(Imp50, na.rm = T),
            MeanImp95 = mean(Imp95, na.rm = T),
            MeanImpWindy = mean(ImpWindy, na.rm = T),
            MeanImpPhysical = mean(ImpPhysical, na.rm = T),
            MeanImpDaytime = mean(ImpDaytime, na.rm = T),
            MeanRelImp50 = mean(RelImp50, na.rm = T),
            MeanRelImp95 = mean(RelImp95, na.rm = T),
            MeanRelImpWindy = mean(RelImpWindy, na.rm = T),
            MeanRelImpPhysical = mean(RelImpPhysical, na.rm = T),
            MeanRelImpDaytime = mean(RelImpDaytime, na.rm = T),)

write_csv(Improvements, file = "data/iNMVImprovements.csv")
write_csv(SummaryImprovements, file = "data/iNMVImprovementsSummary.csv")
