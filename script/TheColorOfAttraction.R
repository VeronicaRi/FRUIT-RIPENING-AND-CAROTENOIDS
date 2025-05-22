############ Capulines #################
###########  22 may del 2025 ########## 
#### Author: Verónica Rincón Rubio ####
############# Capulinero ##############


setwd("//Users/veronicari/Documentos/GitHub/Capulines")


# COLORS
#8bc34a, #fdd835, #ff6f00, #b71c1c, #581845



#Step 1: Data Preparation

fruit_data <- read.csv("ColorsInTrees.csv", header = T)
str(fruit_data)
library(tidyverse)

# Check unique values in the wavelength column
unique(fruit_data %>%
         pivot_longer(
           cols = starts_with("X"),
           names_to = "wavelength",
           names_prefix = "X",
           values_to = "reflectance"
         ) %>%
         pull(wavelength))

# Proper pivot_longer transformation
spectral_data <- fruit_data %>%
  pivot_longer(
    cols = starts_with("X"),    # Select columns X400nm to X700nm
    names_to = "wavelength",    # Create a column for wavelength
    names_prefix = "X",         # Remove the "X" prefix
    values_to = "reflectance"   # Values for reflectance
  ) %>%
  mutate(wavelength = parse_number(wavelength))  # Convert to numeric safely

# Verify the structure of the transformed data
str(spectral_data)

# Check for NAs in the wavelength column
sum(is.na(spectral_data$wavelength))

# Check the transformed data
glimpse(spectral_data)

##Step 2: Summary Statistics
# Summarize mean and standard deviation of reflectance by ColorCategory and wavelength
reflectance_summary <- spectral_data %>%
  group_by(RipeningStage, wavelength) %>%
  summarize(
    mean_reflectance = mean(reflectance, na.rm = TRUE),
    sd_reflectance = sd(reflectance, na.rm = TRUE),
    .groups = "drop"
  )

# View summary
head(reflectance_summary)


###Step 3: Visualization
##3.1 Line Plot for Reflectance Curves
# Reorder ColorCategory levels
reflectance_summary$RipeningStage <- factor(reflectance_summary$RipeningStage, 
                                            levels = c("G", "YG", "RY", "R", "P"))

# Define custom colors
custom_colors <- c("G" = "#8bc34a", 
                   "YG" = "#fdd835", 
                   "RY" = "#ff6f00", 
                   "R" = "#b71c1c", 
                   "P" = "#581845")

# Line plot with ribbons
ggplot(reflectance_summary, aes(x = wavelength, y = mean_reflectance, color = RipeningStage)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_reflectance - sd_reflectance,
                  ymax = mean_reflectance + sd_reflectance,
                  fill = RipeningStage), alpha = 0.2) +
  # Apply custom colors
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  # Add labels and theme
  labs(
    title = "Reflectance Curves by Ripening Stage",
    x = "Wavelength (nm)",
    y = "Mean Reflectance"
  ) +
  theme_minimal()



###Step 4: PCA for Spectral Data
# Extract spectral columns for PCA
spectral_matrix <- fruit_data %>%
  select(starts_with("X")) %>%
  as.matrix()

# Perform PCA
pca_result <- prcomp(spectral_matrix, center = TRUE, scale. = TRUE)

# Print PCA summary
summary(pca_result)

# Add PCA scores back to the dataset
pca_scores <- as_tibble(pca_result$x) %>%
  bind_cols(fruit_data %>% select(RipeningStage))

# Visualize PCA
# Reorder ColorCategory levels
pca_scores$RipeningStage  <- factor(pca_scores$RipeningStage, 
                                    levels = c("G", "YG", "RY", "R", "P"))

# Define custom colors
custom_colors <- c("G" = "#8bc34a", 
                   "YG" = "#fdd835", 
                   "RY" = "#ff6f00", 
                   "R" = "#b71c1c", 
                   "P" = "#581845")

# PCA plot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = RipeningStage)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = custom_colors) +  # Use custom colors
  labs(
    title = "PCA of Spectral Data",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal()

# View PCA loadings
loadings <- pca_result$rotation
print(loadings)

###Step 5: Statistical Testing
# Load mgcv package
library(mgcv)

# Fit a GAM model
gam_model1 <- gam(reflectance ~ s(wavelength, bs = "cs") + RipeningStage, data = spectral_data)
gam_model2 <- gam(reflectance ~ s(wavelength, bs = "cs") + Stage, data = spectral_data)


# Summary of the model
summary(gam_model1)
summary(gam_model2)


#################################################
################### NUTRIENTS ###################
#################################################

setwd("//Users/veronicari/Documentos/GitHub/Capulines")

# Load required libraries
library(emmeans) # For post-hoc comparisons
library(ggplot2) # For visualizations
library(patchwork)


# Ensure ColorCategory is a factor
nutrient_data <- read.csv("Nutrients.csv", header = T)
nutrient_data 

nutrient_data$RipeningStage <- as.factor(nutrient_data$RipeningStage)
nutrient_data$Site <- as.factor(nutrient_data$Site)
nutrient_data$Tree <- as.factor(nutrient_data$Tree)
nutrient_data$Carbohydrates <- as.numeric(nutrient_data$Carbohydrates)

str(nutrient_data)

#### Summary

# Run ANOVA for each variable separately
summary(aov(Moisture ~ RipeningStage + Tree, data = nutrient_data))
summary(aov(Ashes ~ RipeningStage + Tree, data = nutrient_data))
summary(aov(Protein ~ RipeningStage + Tree, data = nutrient_data))
summary(aov(Fat ~ RipeningStage + Tree, data = nutrient_data))
summary(aov(Fiber ~ RipeningStage + Tree, data = nutrient_data))
summary(aov(Carbohydrates ~ RipeningStage + Tree, data = nutrient_data))
summary(aov(TotalPolyphenols ~ RipeningStage + Tree, data = nutrient_data))
summary(aov(PolyphenolConcentration ~ RipeningStage + Tree, data = nutrient_data))
summary(aov(AntioxidantActivity ~ RipeningStage + Tree, data = nutrient_data))


TukeyHSD(aov(Moisture ~ RipeningStage + Tree, data = nutrient_data))
TukeyHSD(aov(Protein ~ RipeningStage + Tree, data = nutrient_data))
TukeyHSD(aov(Fiber ~ RipeningStage + Tree, data = nutrient_data))
TukeyHSD(aov(TotalPolyphenols ~ RipeningStage + Tree, data = nutrient_data))

# Run ANOVA for each variable without trees
summary(aov(Moisture ~ RipeningStage, data = nutrient_data))
summary(aov(Ashes ~ RipeningStage, data = nutrient_data))
summary(aov(Protein ~ RipeningStage, data = nutrient_data))
summary(aov(Fat ~ RipeningStage, data = nutrient_data))
summary(aov(Fiber ~ RipeningStage, data = nutrient_data))
summary(aov(Carbohydrates ~ RipeningStage, data = nutrient_data))
summary(aov(TotalPolyphenols ~ RipeningStage, data = nutrient_data))
summary(aov(PolyphenolConcentration ~ RipeningStage, data = nutrient_data))
summary(aov(AntioxidantActivity ~ RipeningStage, data = nutrient_data))



##Perform ANOVA for Each Nutrient
####extended

# ANOVA for Moisture
moisture_model <- aov(Moisture ~ RipeningStage  + Tree, data = nutrient_data)
summary(moisture_model)

# Simplify the model by removing the interaction term
moisture_model_simple <- lm(Moisture ~ RipeningStage  + Tree, data = nutrient_data)
summary(moisture_model_simple)
anova(moisture_model_simple)

# Boxplot for Moisture
# Reorder the levels of the Color factor
nutrient_data$RipeningStage  <- factor(nutrient_data$RipeningStage , levels = c("YG", "RY", "R", "P"))

# Define custom colors for the categories
custom_colors <- c("YG" = "#FFC300", "RY" = "#FF5733", "R" = "#C70039", "P" = "#581845")

# Calculate summary statistics (mean and standard error)
moisture_summary <- nutrient_data %>%
  group_by(RipeningStage ) %>%
  summarize(
    mean_moisture = mean(Moisture, na.rm = TRUE),
    se_moisture = sd(Moisture, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Create the refined plot
plot_moisture<-ggplot(nutrient_data, aes(x = RipeningStage , y = Moisture, color = RipeningStage , shape = Tree)) +
  # Individual data points with jitter
  geom_jitter(width = 0, size = 3, alpha = 0.8) +
  
  # Mean markers
  geom_point(data = moisture_summary, aes(x = RipeningStage , y = mean_moisture), 
             shape = 18, size = 3, color = "black", alpha = 0.4, inherit.aes = FALSE) +
  
  # Error bars for standard error
  geom_errorbar(data = moisture_summary, aes(x = RipeningStage , 
                                             ymin = mean_moisture - se_moisture, 
                                             ymax = mean_moisture + se_moisture), 
                width = 0.2, color = "black", alpha = 0.2, size = 0.5, inherit.aes = FALSE) +
  
  # Custom colors
  scale_color_manual(values = custom_colors) +
  
  # Labels and theme
  labs(title = "Moisture Across Ripening Stages",
       x = "Ripening Stage",
       y = "Moisture (g)"
  ) +
  theme_minimal()


# ANOVA for Ashes
ashes_model <- aov(Ashes ~ RipeningStage  + Tree, data = nutrient_data)
summary(ashes_model)

# Simplify the model by removing the interaction term
ashes_model_simple <- lm(Ashes ~ RipeningStage  + Tree, data = nutrient_data)
summary(ashes_model_simple)
anova(ashes_model_simple)

# Boxplot

# Calculate summary statistics (mean and standard error)
Ashes_summary <- nutrient_data %>%
  group_by(RipeningStage ) %>%
  summarize(
    mean_Ashes = mean(Ashes, na.rm = TRUE),
    se_Ashes = sd(Ashes, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Create the refined plot
plot_ashes<-ggplot(nutrient_data, aes(x = RipeningStage , y = Ashes, color = RipeningStage , shape = Tree)) +
  # Individual data points with jitter
  geom_jitter(width = 0, size = 3, alpha = 0.8) +
  
  # Mean markers
  geom_point(data = Ashes_summary, aes(x = RipeningStage , y = mean_Ashes), 
             shape = 18, size = 3, color = "black", alpha = 0.4, inherit.aes = FALSE) +
  
  # Error bars for standard error
  geom_errorbar(data = Ashes_summary, aes(x = RipeningStage , 
                                          ymin = mean_Ashes - se_Ashes, 
                                          ymax = mean_Ashes + se_Ashes), 
                width = 0.2, color = "black", alpha = 0.2, size = 0.5, inherit.aes = FALSE) +
  
  # Custom colors
  scale_color_manual(values = custom_colors) +
  
  # Labels and theme
  labs(
    title = "Ashes Across Ripening Stages",
    x = "Ripening Stage",
    y = "Ashes (g)"
  ) +
  theme_minimal()


# ANOVA for Protein
protein_model <- aov(Protein ~ RipeningStage  + Tree, data = nutrient_data)
summary(protein_model)

# Simplify the model by removing the interaction term
protein_model_simple <- lm(Protein ~ RipeningStage  + Tree, data = nutrient_data)
summary(protein_model_simple)
anova(protein_model_simple)


# Boxplot
# Calculate summary statistics (mean and standard error)
Protein_summary <- nutrient_data %>%
  group_by(RipeningStage ) %>%
  summarize(
    mean_Protein = mean(Protein, na.rm = TRUE),
    se_Protein = sd(Protein, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Create the refined plot
plot_protein<-ggplot(nutrient_data, aes(x = RipeningStage , y = Protein, color = RipeningStage , shape = Tree)) +
  # Individual data points with jitter
  geom_jitter(width = 0, size = 3, alpha = 0.8) +
  
  # Mean markers
  geom_point(data = Protein_summary, aes(x = RipeningStage , y = mean_Protein), 
             shape = 18, size = 3, color = "black", alpha = 0.4, inherit.aes = FALSE) +
  
  # Error bars for standard error
  geom_errorbar(data = Protein_summary, aes(x = RipeningStage , 
                                            ymin = mean_Protein - se_Protein, 
                                            ymax = mean_Protein + se_Protein), 
                width = 0.2, color = "black", alpha = 0.2, size = 0.5, inherit.aes = FALSE) +
  
  # Custom colors
  scale_color_manual(values = custom_colors) +
  
  # Labels and theme
  labs(
    title = "Protein Across Ripening Stages",
    x = "Ripening Stage",
    y = "Protein(g)"
  ) +
  theme_minimal()



# ANOVA for Fat
fat_model <- aov(Fat ~ RipeningStage  + Tree, data = nutrient_data)
summary(fat_model)

# Simplify the model by removing the interaction term
fat_model_simple <- lm(Fat ~ RipeningStage  + Tree, data = nutrient_data)
summary(fat_model_simple)
anova(fat_model_simple)


# Boxplot
# Calculate summary statistics (mean and standard error)
fat_summary <- nutrient_data %>%
  group_by(RipeningStage ) %>%
  summarize(
    mean_fat = mean(Fat, na.rm = TRUE),
    se_fat = sd(Fat, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Create the refined plot
plot_fat<-ggplot(nutrient_data, aes(x = RipeningStage , y = Fat, color = RipeningStage , shape = Tree)) +
  # Individual data points with jitter
  geom_jitter(width = 0, size = 3, alpha = 0.8) +
  
  # Mean markers
  geom_point(data = fat_summary, aes(x = RipeningStage , y = mean_fat), 
             shape = 18, size = 3, color = "black", alpha = 0.4, inherit.aes = FALSE) +
  
  # Error bars for standard error
  geom_errorbar(data = fat_summary, aes(x = RipeningStage , 
                                        ymin = mean_fat - se_fat, 
                                        ymax = mean_fat + se_fat), 
                width = 0.2, color = "black", alpha = 0.2, size = 0.5, inherit.aes = FALSE) +
  
  # Custom colors
  scale_color_manual(values = custom_colors) +
  
  # Labels and theme
  labs(
    title = "Fat Across Ripening Stages",
    x = "Ripening Stage",
    y = "Fat (g)"
  ) +
  theme_minimal()


# ANOVA for Fiber
fiber_model <- aov(Fiber ~ RipeningStage  + Tree, data = nutrient_data)
summary(fiber_model)

# Simplify the model by removing the interaction term
fiber_model_simple <- lm(Fiber ~ RipeningStage  + Tree, data = nutrient_data)
summary(fiber_model_simple)
anova(fiber_model_simple)


# Boxplot

# Calculate summary statistics (mean and standard error)
fiber_summary <- nutrient_data %>%
  group_by(RipeningStage ) %>%
  summarize(
    mean_fiber = mean(Fiber, na.rm = TRUE),
    se_fiber = sd(Fiber, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Create the refined plot
plot_fiber<-ggplot(nutrient_data, aes(x = RipeningStage , y = Fiber, color = RipeningStage , shape = Tree)) +
  # Individual data points with jitter
  geom_jitter(width = 0, size = 3, alpha = 0.8) +
  
  # Mean markers
  geom_point(data = fiber_summary, aes(x = RipeningStage , y = mean_fiber), 
             shape = 18, size = 3, color = "black", alpha = 0.4, inherit.aes = FALSE) +
  
  # Error bars for standard error
  geom_errorbar(data = fiber_summary, aes(x = RipeningStage , 
                                          ymin = mean_fiber - se_fiber, 
                                          ymax = mean_fiber + se_fiber), 
                width = 0.2, color = "black", alpha = 0.2, size = 0.5, inherit.aes = FALSE) +
  
  # Custom colors
  scale_color_manual(values = custom_colors) +
  
  # Labels and theme
  labs(
    title = "Fiber Across Ripening Stages",
    x = "Ripening Stage",
    y = "Fiber (g)"
  ) +
  theme_minimal()


# ANOVA for Carbohydrates
carb_model <- aov(Carbohydrates ~ RipeningStage  + Tree, data = nutrient_data)
summary(carb_model)

# Simplify the model by removing the interaction term
carb_model_simple <- lm(Carbohydrates ~ RipeningStage  + Tree, data = nutrient_data)
summary(carb_model_simple)
anova(carb_model_simple)


# Boxplot 
# Calculate summary statistics (mean and standard error)
carb_summary <- nutrient_data %>%
  group_by(RipeningStage ) %>%
  summarize(
    mean_carb = mean(Carbohydrates, na.rm = TRUE),
    se_carb = sd(Carbohydrates, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Create the refined plot
plot_carb<-ggplot(nutrient_data, aes(x = RipeningStage , y = Carbohydrates, color = RipeningStage , shape = Tree)) +
  # Individual data points with jitter
  geom_jitter(width = 0, size = 3, alpha = 0.8) +
  
  # Mean markers
  geom_point(data = carb_summary, aes(x = RipeningStage , y = mean_carb), 
             shape = 18, size = 3, color = "black", alpha = 0.4, inherit.aes = FALSE) +
  
  # Error bars for standard error
  geom_errorbar(data = carb_summary, aes(x = RipeningStage , 
                                         ymin = mean_carb - se_carb, 
                                         ymax = mean_carb + se_carb), 
                width = 0.2, color = "black", alpha = 0.2, size = 0.5, inherit.aes = FALSE) +
  
  # Custom colors
  scale_color_manual(values = custom_colors) +
  
  # Labels and theme
  labs(
    title = "Carbohydrates Across Ripening Stages",
    x = "Ripening Stage",
    y = "Carbohydrates (g)"
  ) +
  theme_minimal()



# ANOVA for Polyphenols
polyphenols_model <- aov(TotalPolyphenols ~ RipeningStage  + Tree, data = nutrient_data)
summary(polyphenols_model)

# Simplify the model by removing the interaction term
poly_model_simple <- lm(TotalPolyphenols ~ RipeningStage  + Tree, data = nutrient_data)
summary(poly_model_simple)
anova(poly_model_simple)


# Boxplot
# Calculate summary statistics (mean and standard error)
polyphenols_summary <- nutrient_data %>%
  group_by(RipeningStage ) %>%
  summarize(
    mean_polyphenols = mean(TotalPolyphenols, na.rm = TRUE),
    se_polyphenols = sd(TotalPolyphenols, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Create the refined plot
plot_poly<-ggplot(nutrient_data, aes(x = RipeningStage , y = TotalPolyphenols, color = RipeningStage , shape = Tree)) +
  # Individual data points with jitter
  geom_jitter(width = 0, size = 3, alpha = 0.8) +
  
  # Mean markers
  geom_point(data = polyphenols_summary, aes(x = RipeningStage , y = mean_polyphenols), 
             shape = 18, size = 3, color = "black", alpha = 0.4, inherit.aes = FALSE) +
  
  # Error bars for standard error
  geom_errorbar(data = polyphenols_summary, aes(x = RipeningStage , 
                                                ymin = mean_polyphenols - se_polyphenols, 
                                                ymax = mean_polyphenols + se_polyphenols), 
                width = 0.2, color = "black", alpha = 0.2, size = 0.5, inherit.aes = FALSE) +
  
  # Custom colors
  scale_color_manual(values = custom_colors) +
  
  # Labels and theme
  labs(
    title = "Polyphenols Across Ripening Stages",
    x = "Ripening Stage",
    y = "Polyphenols (mg/100mL)"
  ) +
  theme_minimal()




# ANOVA for Antioxidant Activity
antioxidant_model <- aov(AntioxidantActivity ~ RipeningStage + Tree, data = nutrient_data)
summary(antioxidant_model)

# Simplify the model by removing the interaction term
anti_model_simple <- lm(AntioxidantActivity ~ RipeningStage  + Tree, data = nutrient_data)
summary(anti_model_simple)
anova(anti_model_simple)

# Boxplot

# Calculate summary statistics (mean and standard error)
antioxidant_summary <- nutrient_data %>%
  group_by(RipeningStage ) %>%
  summarize(
    mean_antioxidant = mean(AntioxidantActivity, na.rm = TRUE),
    se_antioxidant = sd(AntioxidantActivity, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Create the refined plot
plot_antioxidant<-ggplot(nutrient_data, aes(x = RipeningStage , y = AntioxidantActivity, color = RipeningStage , shape = Tree)) +
  # Individual data points with jitter
  geom_jitter(width = 0, size = 3, alpha = 0.8) +
  
  # Mean markers
  geom_point(data = antioxidant_summary, aes(x = RipeningStage , y = mean_antioxidant), 
             shape = 18, size = 3, color = "black", alpha = 0.4, inherit.aes = FALSE) +
  
  # Error bars for standard error
  geom_errorbar(data = antioxidant_summary, aes(x = RipeningStage , 
                                                ymin = mean_antioxidant - se_antioxidant, 
                                                ymax = mean_antioxidant + se_antioxidant), 
                width = 0.2, color = "black", alpha = 0.2, size = 0.5, inherit.aes = FALSE) +
  
  # Custom colors
  scale_color_manual(values = custom_colors) +
  
  # Labels and theme
  labs(
    title = "Antioxidant Activity Across Ripening Stages",
    x = "Ripening Stage",
    y = "Antioxidant Activity  (%)"
  ) +
  theme_minimal()


(plot_moisture | plot_ashes | plot_protein) / 
  (plot_fat | plot_fiber | plot_carb) / (plot_poly | plot_antioxidant) + 
  plot_annotation(title = "Nutrient Composition Across Ripening Stages")



# Load emmeans package
library(emmeans)

# Post-hoc comparisons for Moisture
tukey_moisture <- TukeyHSD(moisture_model)
print(tukey_moisture)
TukeyHSD(moisture_model, "RipeningStage")
TukeyHSD(moisture_model, "Tree")

# Post-hoc comparisons for Protein
tukey_protein <- TukeyHSD(protein_model)
print(tukey_protein)
TukeyHSD(protein_model, "RipeningStage")
TukeyHSD(protein_model, "Tree")

# Post-hoc comparisons for Protein
tukey_fiber <- TukeyHSD(fiber_model)
print(tukey_fiber)
TukeyHSD(fiber_model, "RipeningStage")
TukeyHSD(fiber_model, "Tree")



####Just carotenoides

carotenoid_data<-read.csv("Carotenoids.csv", header = T)
str(carotenoid_data)
summary(carotenoid_data)

library(dplyr)

# Summary by Color
carotenoid_summary <- carotenoid_data %>%
  group_by(RipeningStage ) %>%
  summarize(
    mean_carotenoids = mean(TotalCarotenoids, na.rm = TRUE),
    sd_carotenoids = sd(TotalCarotenoids, na.rm = TRUE),
    mean_absorbance = mean(Absorbance, na.rm = TRUE),
    sd_absorbance = sd(Absorbance, na.rm = TRUE)
  )

print(carotenoid_summary)

carotenoid_model <- aov(TotalCarotenoids ~ RipeningStage , data = carotenoid_data)
summary(carotenoid_model)
TukeyHSD(carotenoid_model)

absorbance_model <- lm(TotalCarotenoids ~ Absorbance, data = carotenoid_data)
summary(absorbance_model)


# Reorder Color levels
carotenoid_data$RipeningStage  <- factor(carotenoid_data$RipeningStage , 
                                         levels = c("G", "YG", "RY", "R", "P"))

ggplot(carotenoid_data, aes(x = RipeningStage , y = TotalCarotenoids, fill = RipeningStage )) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.8) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.6) +
  scale_fill_manual(values = c("G" = "#8bc34a", "YG" = "#fdd835", 
                               "RY" = "#ff6f00", "R" = "#b71c1c", 
                               "P" = "#581845")) +
  labs(title = "Total Carotenoids by Ripening Stage", x = "Ripening Stage", y = "Total Carotenoids") +
  theme_minimal()


ggplot(carotenoid_data, aes(x = Absorbance, y = TotalCarotenoids, color = RipeningStage )) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c("G" = "#8bc34a", "YG" = "#fdd835", 
                                "RY" = "#ff6f00", "R" = "#b71c1c", 
                                "P" = "#581845")) +
  labs(title = "Absorbance vs. Total Carotenoids", x = "Absorbance", y = "Total Carotenoids") +
  theme_minimal()





###Germination 
germination_data<-read.csv("Germination.csv", header = T)
# Load required libraries
library(dplyr)
library(ggplot2)

(germination_data)


# Aggregate data: count total seeds and germinated seeds per ripening stage
germination_summary <- germination_data %>%
  group_by(RipeningStage) %>%  # Group by ripening stage
  summarise(
    total_seeds = n(),  # Total seeds per ripening stage
    germinated_seeds = sum(Germination, na.rm = TRUE)  # Count germinated seeds
  )

# View the transformed data
print(germination_summary)


# Fit a binomial GLM
germination_model <- glm(
  cbind(germinated_seeds, total_seeds - germinated_seeds) ~ RipeningStage, 
  data = germination_summary, 
  family = binomial
)

# View the model results
summary(germination_model)


###WITH TREES
# Load required libraries
library(dplyr)
library(ggplot2)

# Aggregate data: count total seeds and germinated seeds per ripening stage and tree
germination_summary <- germination_data %>%
  group_by(RipeningStage, Tree) %>%  # Group by both ripening stage and tree
  summarise(
    total_seeds = n(),  # Total seeds per ripening stage and tree
    germinated_seeds = sum(Germination, na.rm = TRUE),  # Count germinated seeds
    .groups = "drop"
  )

# View the transformed data
print(germination_summary)

# Fit a binomial GLM accounting for Tree
germination_model <- glm(
  cbind(germinated_seeds, total_seeds - germinated_seeds) ~ RipeningStage + Tree, 
  data = germination_summary, 
  family = binomial
)

# View the model results
summary(germination_model)


########################
## GERMINATION TIME
########################

# Load necessary library
library(survival)
germination_data<-read.csv("Germination.csv", header = T)

# Define the observation period (Censoring at 130 days)
censoring_time <- 130
germination_data$GerminationTimeDays[is.na(germination_data$GerminationTimeDays)] <- censoring_time

# Create a survival object
surv_object <- Surv(time = germination_data$GerminationTimeDays, 
                    event = germination_data$Germination == 1)  # 1 = Germinated, 0 = Censored

# Fit Cox model with Tree as a frailty term (random effect)
cox_model <- coxph(surv_object ~ RipeningStage + frailty(Tree), data = germination_data)

# Display model summary
summary(cox_model)

# Ensure RipeningStage is a factor with the desired order
germination_summary$RipeningStage <- factor(germination_summary$RipeningStage, 
                                            levels = c("G", "YG", "RY", "R", "P"))

# Create the bar plot
# ANOVA for germination timing
# Summarize data for germination success
germination_summary <- germination_data %>%
  group_by(RipeningStage, Tree) %>%
  summarize(
    germination_rate = mean(Germination, na.rm = TRUE),
    median_germination_time = median(GerminationTimeDays, na.rm = TRUE),
    .groups = "drop"
  )

print(germination_summary)

timing_model <- aov(GerminationTimeDays ~ RipeningStage * Tree, data = germination_data)

# Summary of the model
summary(timing_model)

# Post-hoc test
TukeyHSD(timing_model)
ggplot(germination_summary, aes(x = RipeningStage, y = germination_rate, fill = RipeningStage)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  facet_wrap(~ Tree) +
  scale_fill_manual(values = c("G" = "#8bc34a", "YG" = "#fdd835", 
                               "RY" = "#ff6f00", "R" = "#b71c1c", 
                               "P" = "#581845")) +
  labs(title = "Germination Success by Ripening Stage and Tree", 
       x = "Ripening Stage", y = "Proportion of Germinated Seeds") +
  theme_minimal()
View(germination_summary)

ggplot(germination_data, aes(x = RipeningStage, y = GerminationTimeDays, color = RipeningStage)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
  facet_wrap(~ Tree) +
  scale_color_manual(values = c("G" = "#8bc34a", "YG" = "#fdd835", 
                                "RY" = "#ff6f00", "R" = "#b71c1c", 
                                "P" = "#581845")) +
  labs(title = "Germination Timing by Ripening Stage and Tree", 
       x = "Ripening Stage", y = "Days to Germination") +
  theme_minimal()




######################
###### TREE 29 ######
######################

data<-read.csv("tree_a29_full_data.csv", header = T)
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Example: Cleaning the dataset
data_clean <- data %>%
  mutate(RipeningStage = factor(RipeningStage, levels = c("G", "YG", "RY", "R", "P")))

# Melting data for easier plotting of nutrients and antioxidants
library(tidyr)
data_long <- data_clean %>%
  pivot_longer(cols = c(Protein, TotalPolyphenols, AntioxidantActivity),
               names_to = "BenefitType",
               values_to = "Value")

# Plot
ggplot(data_long, aes(x = RipeningStage, y = Value, fill = BenefitType)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~BenefitType, scales = "free_y") +
  labs(title = "Benefits of Birds and Plants Across Ripening Stages",
       x = "Ripening Stage", y = "Value") +
  theme_minimal()



# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Example dataset (replace `data` with your actual dataset)
data_clean <- data %>%
  mutate(RipeningStage = factor(RipeningStage, levels = c("G", "YG", "RY", "R", "P")))

# Pivot data for nutrients into long format
data_long <- data_clean %>%
  pivot_longer(cols = c(Protein, Fat, Carbohydrates, Fiber, Ashes),
               names_to = "Nutrient",
               values_to = "Value")

# Plot nutrients with germination overlay
ggplot(data_long, aes(x = RipeningStage, y = Value, fill = Nutrient)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_line(data = data_clean, aes(x = RipeningStage, y = Germination * 100, group = 1, color = "Germination"), 
            size = 1, inherit.aes = FALSE) +
  geom_point(data = data_clean, aes(x = RipeningStage, y = Germination * 100, color = "Germination"), 
             size = 3, inherit.aes = FALSE) +
  scale_y_continuous(
    name = "Nutrient Value",
    sec.axis = sec_axis(~./100, name = "Germination (%)")
  ) +
  labs(
    title = "Nutrients and Germination Across Ripening Stages",
    x = "Ripening Stage",
    fill = "Nutrient",
    color = "Germination"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  scale_color_manual(values = c("Germination" = "red"))