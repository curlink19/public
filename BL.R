`
## Author: Filatkina Kristina
## Data: 2025-04-01
## Description: NCA analysis seminar

#####--------------- Load packages and source functions ---------------#####

install.packages("tidyverse")
install.packages("readxl")
install.packages("cowplot")
install.packages("scales")
install.packages("PKNCA")

library(tidyverse)
library(readxl)
library(cowplot)
library(scales)
library(PKNCA)

theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())
MSDcol <- c("#1a1866", "#f2b93b", "#b73b58", "#a2d620", "#5839bb", "#9c4ec7", "#3a6eba", "#efdd3c", "#69686d")

funSum <- list(N = ~n(),
               mean   = ~mean(., na.rm = T),
               median = ~median(., na.rm = T),
               min    = ~min(., na.rm = T),
               max    = ~max(., na.rm = T),
               sd     = ~sd(., na.rm = T),
               se     = ~sd(., na.rm = T)/sqrt(n()))


#####--------------- Load datasets ---------------#####

x_iv <- read.csv("WorkDir/NCA/Drug_X_iv.csv")
x_po <- read.csv("WorkDir/NCA/Drug_X_po.csv")
y_iv <- read.csv("WorkDir/NCA/Drug_Y_iv.csv")
bioeq <- read.csv("WorkDir/NCA/BioEq_data.csv")

#####--------------- Task 1: Perform exploratory data analysis ---------------#####

# Create and characterize a plot for PK time profiles

# data_X_iv

ggplot(data_X_iv, aes(x = TIME, y = log(DV), by = as.factor(ID))) +
  geom_line() + geom_point() +
  labs(title = "Зависимость DV от Time",
       x = "TIME",
       y = "DV",
       color = "Id") +
  theme_minimal() 
# data_X_po

ggplot(data_X_po, aes(x = TIME, y = DV, by = as.factor(ID))) +
  geom_line() + geom_point() +
  labs(title = "Зависимость DV от Time",
       x = "TIME",
       y = "DV",
       color = "Id") +
  theme_minimal() 

# data_Y_iv


ggplot(data_Y_iv, aes(x = TIME, y = DV, by = as.factor(ID))) +
  geom_line() + geom_point() +
  labs(title = "Зависимость DV от Time",
       x = "TIME",
       y = "DV",
       color = "Id") +
  theme_minimal() 




# Set consistent plotting theme
theme_set(theme_bw(base_size = 12) + 
            theme(panel.grid.minor = element_blank(),
                  legend.position = "bottom"))


### ----------------------- Exploratory Analysis ----------------------- ###

# Function to create PK profile plots
plot_pk_profile <- function(data, title) {
  ggplot(data, aes(x = TIME, y = DV, group = factor(ID), color = factor(ID))) +
    geom_line(linewidth = 0.8) + 
    geom_point(size = 2) +
    labs(title = paste("PK Profile:", title),
         x = "Time (hours)",
         y = "Drug Concentration",
         color = "Subject ID") +
    scale_color_manual(values = custom_colors) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Generate and display all PK profiles
pk_plots <- list(
  plot_pk_profile(data$x_iv, "Drug X IV"),
  plot_pk_profile(data$x_po, "Drug X PO"),
  plot_pk_profile(data$y_iv, "Drug Y IV")
)

# Display plots in a grid
plot_grid(plotlist = pk_plots, ncol = 1)

### ----------------------- NCA Calculations ----------------------- ###

# Function to perform NCA using PKNCA
perform_nca <- function(data, route, dose_time = 0) {
  # Create concentration object
  conc_obj <- PKNCAconc(data, DV~TIME|ID)
  
  # Create dose object
  dose_data <- data %>% filter(TIME == dose_time)
  dose_obj <- PKNCAdose(dose_data, AMT~TIME|ID, route = route)
  
  # Define intervals for calculation
  intervals <- data.frame(
    start = min(data$TIME[!is.na(data$DV)]),
    end = max(data$TIME[!is.na(data$DV)]),
    auclast = TRUE,
    aucinf.obs = TRUE,
    aumclast = TRUE,
    aumcinf.obs = TRUE,
    mrt.last = TRUE
  )
  
  # Combine data and perform NCA
  data_obj <- PKNCAdata(conc_obj, dose_obj, intervals = intervals)
  results <- pk.nca(data_obj)
  
  return(as.data.frame(results) %>% spread(PPTESTCD, PPORRES))
}

# Calculate NCA parameters for each scenario
nca_results <- list(
  x_iv = perform_nca(data$x_iv, "intravascular"),
  x_po = perform_nca(data$x_po, "extravascular"),
  y_iv = perform_nca(data$y_iv, "intravascular")
)

### ----------------------- Manual Calculations ----------------------- ###

# Function for manual trapezoidal AUC calculation
calculate_auc <- function(time, conc) {
  auc <- 0
  for (i in 2:length(time)) {
    auc <- auc + (conc[i-1] + conc[i]) * (time[i] - time[i-1]) / 2
  }
  return(auc)
}

# Example manual calculation for Subject 1 (Drug X IV)
subject_data <- data$x_iv %>% filter(ID == 1) %>% arrange(TIME)

# Calculate AUC manually
manual_auc <- calculate_auc(subject_data$TIME, subject_data$DV)

# Terminal phase analysis for lambda_z
terminal_points <- tail(subject_data, 3)
lambda_z_model <- lm(log(DV) ~ TIME, data = terminal_points)
lambda_z <- -coef(lambda_z_model)[2]

# AUCinf calculation
auc_inf <- tail(subject_data$DV, 1)/lambda_z
total_auc <- manual_auc + auc_inf

# Compare with PKNCA results
comparison <- data.frame(
  Parameter = c("AUClast", "AUCinf"),
  Manual = c(manual_auc, total_auc),
  PKNCA = c(filter(nca_results$x_iv, ID == 1)$auclast,
            filter(nca_results$x_iv, ID == 1)$aucinf.obs)
)

### ----------------------- Bioequivalence Analysis ----------------------- ###

# Prepare bioequivalence data
bioeq_data <- data$bioeq %>%
  mutate(
    Period1_TRT = substr(Seq, 1, 1),
    Period2_TRT = substr(Seq, 2, 2)
  ) %>%
  pivot_longer(
    cols = starts_with("Per"),
    names_to = "Period",
    values_to = "AUC"
  ) %>%
  mutate(
    Treatment = ifelse(Period == "Per1", Period1_TRT, Period2_TRT)
  ) %>%
  select(ID, Seq, Period, Treatment, AUC)

# ANOVA model for bioequivalence
be_model <- aov(log(AUC) ~ Period + Seq + Treatment + ID, data = bioeq_data)

# Calculate geometric mean ratio and 90% CI
calculate_ci <- function(model, alpha = 0.1) {
  mse <- summary(model)[[1]]["Residuals", "Mean Sq"]
  n <- table(bioeq_data$Seq[bioeq_data$Period == "Per1"])
  
  trt_means <- tapply(log(bioeq_data$AUC), bioeq_data$Treatment, mean)
  diff <- trt_means["T"] - trt_means["R"]
  
  ci <- diff + c(-1, 1) * qt(1 - alpha/2, sum(n) - 2) * sqrt(mse/2 * (1/n[1] + 1/n[2]))
  
  return(list(
    GMR = exp(diff),
    CI = exp(ci),
    CV = sqrt(exp(mse) - 1) * 100
  ))
}

be_results <- calculate_ci(be_model)

### ----------------------- Results Presentation ----------------------- ###

# Print key results
cat("Drug X IV NCA Results Summary:\n")
print(summary(nca_results$x_iv[, c("ID", "auclast", "aucinf.obs", "half.life")]))

cat("\nBioequivalence Analysis Results:\n")
cat(sprintf("Geometric Mean Ratio (T/R): %.2f%%\n", be_results$GMR * 100))
cat(sprintf("90%% Confidence Interval: %.2f%% - %.2f%%\n", 
            be_results$CI[1] * 100, be_results$CI[2] * 100))
cat(sprintf("Intrasubject CV: %.1f%%\n", be_results$CV))

# Conclusion based on bioequivalence criteria
if (all(be_results$CI >= 0.8, be_results$CI <= 1.25)) {
  cat("\nConclusion: The formulations are bioequivalent.\n")
} else {
  cat("\nConclusion: The formulations are NOT bioequivalent.\n")
}