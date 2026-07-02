# 29-08-2025
# Analyses of ABM simulation outputs

# Load R packages
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(patchwork)
library(Rmisc)
library(tidyverse)
library(scales)
library(lme4)



# Create Fig. 2
# Compare all models based on RMSE and degree z-scores

model1 <- read.csv("Output_model/abserrors_zscore_randomproof.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_2") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE ) %>% 
  group_by(plot_week) %>%
  dplyr::summarise(
    mean_zscore = sum(z_range, na.rm = TRUE),
    sd_zscore = sd(z_range, na.rm = TRUE),
    plot = first(plot),
    week = first(week),
    model = first(model),
    .groups = "drop"
  )


model2 <- read.csv("Output_model/abserrors_zscore_model2.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_1") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE ) %>% 
  group_by(plot_week) %>%
  dplyr::summarise(
    mean_zscore = sum(z_range, na.rm = TRUE),
    sd_zscore = sd(z_range, na.rm = TRUE),
    plot = first(plot),
    week = first(week),
    model = first(model),
    .groups = "drop"
  )

model3 <- read.csv("Output_model/abserrors_zscore_1regular.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_3") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE ) %>% 
  group_by(plot_week) %>%
  dplyr::summarise(
    mean_zscore = sum(z_range, na.rm = TRUE),
    sd_zscore = sd(z_range, na.rm = TRUE),
    plot = first(plot),
    week = first(week),
    model = first(model),
    .groups = "drop"
  )

model4 <- read.csv("Output_model/abserrors_zscore_1XY.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_4") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE ) %>% 
  group_by(plot_week) %>%
  dplyr::summarise(
    mean_zscore = sum(z_range, na.rm = TRUE),
    sd_zscore = sd(z_range, na.rm = TRUE),
    plot = first(plot),
    week = first(week),
    model = first(model),
    .groups = "drop"
  )


models <- rbind(model1, model2, model3, model4)



# Ensure models are in the correct order
models <- models %>%
  mutate(
    model = factor(model, levels = paste0("model_", 1:4)),
    week = factor(week)
  )


ggplot(models, aes(x = model, y = mean_zscore, 
                   color = plot, shape = week, 
                   group = interaction(plot, week))) +
  geom_point(size = 2, alpha=0.8) +
  geom_line(linewidth = 0.7) +
  geom_errorbar(aes(ymin = mean_zscore - sd_zscore,
                    ymax = mean_zscore + sd_zscore),
                width = 0.2, linewidth = 0.6) +
  scale_color_manual(values = c("#264653","#2a9d8f","#8ab17d","#e9c46a","#f4a261","#e76f51")) +
  scale_shape_manual(values = c(15, 17, 19)) +  # square, triangle, circle
  labs(
    x = "",
    y = "Observed vs simulated degree (z-score)",
    color = "Plot",
    shape = "Flowering \nstage"
  ) +
  facet_wrap( ~plot, scales="free_y") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
    strip.text=element_text(face="bold", size=14),
    axis.text.x = element_text(size=10, angle = 45, hjust = 1)
  ) + scale_x_discrete(breaks=c("model_1","model_2","model_3", "model_4"),
                       labels=c("Baseline", "Species \nabundance", 
                                "Pollinator \nbehaviour", "Spatial \nconfiguration")) + 
  scale_y_log10()



# Fit the mixed model
hist(log(models$mean_zscore))
mod_lmm <- lmer(log(mean_zscore) ~ model + (1 | plot_week), data = models)
summary(mod_lmm)

# Estimated marginal means for each model level
emm <- emmeans(mod_lmm, ~ model)

# Compare adjacent levels only
contrast(emm, method = "consec", adjust = "none")



# Pollinators' r values

# Define parameters for the r distribution
shape_e <- 2
shape_g <- 2
scale_e <- 5
scale_g <- 1
n_draws <- 1000  # Sample size

# Generate prior
name_prior <- paste0("_shapes_eg_", shape_e, "_", shape_g, "_scales", scale_e, "_", scale_g)

# Generate random samples with Gamma distribution
set.seed(123)  # set seed
prior_e <- rgamma(n_draws, shape = shape_e, scale = scale_e)

prior_g <- rgamma(n_draws, shape = shape_g, scale = scale_g)

# Transform to dataframe and take a look at the first rows
prior_e_df <- data.frame(prior_e = prior_e)
prior_g_df <- data.frame(prior_g = prior_g)

head(prior_e_df)
head(prior_g_df)

hist(prior_e_df$prior_e)
hist(prior_g_df$prior_g)

prior_g_df %<>% mutate(type="prior_g") %>% dplyr::rename(r_value= prior_g)
prior_e_df %<>% mutate(type="prior_e") %>% dplyr::rename(r_value= prior_e)

# Load posteriors
temp <- list.files("Output_model/Accepted/", full.names = T)
myfiles <- lapply(temp, read.csv)

myfiles[[1]]
names.temp <- gsub("Output_model/Accepted//|_gm_e\\(2,5\\)_gm_g\\(2,1\\)_th4_4|_accepted|\\.csv", "", temp)

  
for(i in 1:length(myfiles)){
  names(myfiles)[[i]] <- names.temp[[i]]
  temp1 <- myfiles[[i]] 
  temp1 %<>% dplyr::select(r_esp, r_gen) %>% 
    pivot_longer(cols=c(r_esp, r_gen), names_to = "type", values_to = "r_value") 
  temp2 <- rbind(prior_g_df, prior_e_df, temp1)
  temp2 %<>% dplyr::mutate(plot_week=sub("_.*", "", names.temp[[i]]),
                           scenario=sub(".*_", "", names.temp[[i]]))
  myfiles[[i]] <- temp2
  
}

all.plots.weeks <- do.call("rbind", myfiles)

all.plots.weeks %<>% dplyr::filter(scenario=="random")


glimpse(all.plots.weeks)
levels(as.factor(all.plots.weeks$type))


# Create Fig. S2

ggplot(all.plots.weeks, aes(x = r_value, color = type, fill = type)) +
  geom_density(alpha = 0.4) +
  scale_color_manual(values = c("#7a9aaf","#b0d8b5", "#2e4057", "#66a182"),
                     labels = c("Occasional \nforagers (prior)", 
                                "Frequent \nforagers (prior)", 
                                "Occasional \nforagers (posterior)", 
                                "Frequent \nforagers (posterior)")) +
  scale_fill_manual(values = c("#7a9aaf","#b0d8b5", "#2e4057", "#66a182"),
                    labels = c("Occasional \nforagers (prior)", 
                               "Frequent \nforagers (prior)", 
                               "Occasional \nforagers (posterior)", 
                               "Frequent \nforagers (posterior)")) +
  facet_wrap(~plot_week, ncol = 3, scales = "free") +
  theme_bw() + 
  xlim(0, 10) +  
  xlab(expression(paste(italic("r "), "value"))) + 
  ylab("Density") + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        legend.position = "bottom",
        strip.background = element_rect(color = "white", fill = "white", size = 1.5),
        strip.text = element_text(size = 14, face = "bold")) +  # Adjust facet label text
  guides(color = guide_legend(ncol = 2, byrow = TRUE, 
                              override.aes = list(order = c(4, 2, 3, 1))),
         fill = guide_legend(ncol = 2, byrow = TRUE, 
                             override.aes = list(order = c(4, 2, 3, 1))))


# Create Fig. 3 panel A

r.pol <- ggplot(filter(all.plots.weeks, plot_week=="A2"), 
                aes(x = r_value, color = type, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = c("#7a9aaf","#b0d8b5", "#2e4057", "#66a182"),
                     labels = c("Occasional \nforagers (prior)", 
                                "Frequent \nforagers (prior)", 
                                "Occasional \nforagers (posterior)", 
                                "Frequent \nforagers (posterior)")) +
  scale_fill_manual(values = c("#7a9aaf","#b0d8b5", "#2e4057", "#66a182"),
                    labels = c("Occasional \nforagers (prior)", 
                               "Frequent \nforagers (prior)", 
                               "Occasional \nforagers (posterior)", 
                               "Frequent \nforagers (posterior)")) +
  theme_bw() + 
  xlim(0, 10) +  
  xlab(expression(paste(italic("r "), "value"))) + 
  ylab("Density") + 
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=9),
        strip.background = element_rect(color = "white", fill = "white", size = 1.5),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = c(0.98, 0.98),      
        legend.justification = c(1, 1)) +  
  guides(color = guide_legend(ncol = 2, byrow = TRUE, 
                              keyheight = unit(1, "lines"),  # Increase key height for spacing
                              override.aes = list(order = c(4, 2, 3, 1))),
         fill = guide_legend(ncol = 2, byrow = TRUE, 
                             keyheight = unit(1, "lines"),  # Increase key height for spacing
                             override.aes = list(order = c(4, 2, 3, 1))))



# Estimate pollinator activity density based on interaction data

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

temp <- list.files("Output_model/Temporal_slices_individuals/", full.names = T)
myfiles <- lapply(temp, read.csv, sep=";")

names.temp <- gsub("Output_model/Temporal_slices_individuals//network_plot_|_slice_|\\.csv", "", temp)

myfiles[[1]]

n.ind.plant <- list()
rel.abun.pol <- list()
for(i in 1:length(myfiles)){
  names(myfiles)[[i]] <- names.temp[[i]]
  temp1 <- myfiles[[i]] %<>% dplyr::select(-Plant_id)
  temp2 <- data.frame(n_plants=nrow(temp1), 
                      n_pol=ncol(temp1),
                      n_int=sum(temp1 != 0),
                      plot_week= names(myfiles)[[i]])
  temp2 %<>% mutate(plot = substr(plot_week, 1, 1),  
                    slice = as.integer(substr(plot_week, 2, 2))) 
  
  n.ind.plant[[i]] <- temp2
  
  temp3 <- temp1
  temp4 <- colSums(temp3) / sum(temp3)
  
  temp5 <- temp4*temp2$n_plants*5
  temp6 <- (temp5)
  
  temp7 <- data.frame(pollinator_sp = names(temp6),
                      rel_abundance = as.numeric(temp4),
                      total_abundance = as.numeric(temp6),
                      plot_week = names(myfiles)[[i]])
  
  
  rel.abun.pol[[i]] <- temp7
}

all.n.ind <- do.call("rbind", n.ind.plant)

table.net <- all.n.ind %>% dplyr::select(-plot_week)
table.net <- table.net[,c(4,5,1,2,3)]
table.net$n_plant_sp <- c(5,5,3,4,6,3,4,7,3,4,6,3,5,8,2,4,7,2)

mean(table.net$n_plants)
sd(table.net$n_plants)

mean(table.net$n_pol)
sd(table.net$n_pol)

mean(table.net$n_int)
sd(table.net$n_int)

mean(table.net$n_plant_sp)
sd(table.net$n_plant_sp)


all.abun.pol <- do.call("rbind", rel.abun.pol)


# P abundance histogram
all.abun.pol %>% group_by(plot_week) %>% 
  summarise(sum=sum(rel_abundance))


# Create Fig. S1

bins <- seq(0, 1, by = 0.05)

# Create bins for rel_abundance using custom breaks
df <- all.abun.pol %>%
  dplyr::mutate(rel_abundance_bin = cut(rel_abundance, breaks = bins, include.lowest = TRUE))

# Summarize the data to count the number of distinct pollinator species per bin
df_summary <- df %>%
  group_by(rel_abundance_bin, plot_week) %>%
  dplyr::summarize(count = n_distinct(pollinator_sp)) %>%
  # Convert bin labels to numeric midpoints for continuous x-axis
  mutate(rel_abundance_mid = (as.numeric(gsub("[^0-9\\.]", "", sub(",.*", "", rel_abundance_bin))) +
                                as.numeric(gsub("[^0-9\\.]", "", sub(".*,", "", rel_abundance_bin)))) / 2)
glimpse(df_summary)

ggplot(df_summary, aes(x = rel_abundance_mid, y = count)) +
  geom_col(fill = "#2e4057") +  
  labs(x = "Relative interaction frequency", 
       y = "Number of pollinator species") +
  theme_bw() +
  facet_wrap(~plot_week, ncol=3, scales="free_y") +
  theme(
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.text = element_text(size = 12),
    strip.background = element_rect(color = "white", fill = "white", size = 1.5),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1))



# Create Fig. 3 panel B

#####
temp <- list.files("Output_model/Accepted/", full.names = T)
myfiles <- lapply(temp, read.csv)

myfiles[[1]]
names.temp <- gsub("Output_model/Accepted//|_gm_e\\(2,5\\)_gm_g\\(2,1\\)_th4_4|_accepted|\\.csv", "", temp)

for(i in 1:length(myfiles)){
  names(myfiles)[[i]] <- names.temp[[i]]
  temp1 <- myfiles[[i]] 
  names(temp1)[2] <- "pollinator_sp"
  
  pol.abun <- temp1 %>% dplyr::select(pollinator_sp, Tipo) %>%
    dplyr::mutate(plot_week=sub("_.*", "", names.temp[[i]]))
  pol.abun %<>% merge(all.abun.pol, by=c("pollinator_sp", "plot_week")) %>% unique()
  abundance_sum <- aggregate(total_abundance ~ Tipo, data = pol.abun, sum)
  abundance_sum$Tipo <- c("esp", "gen")
  abundance_sum %<>% dplyr::rename(sp_gen_esp=Tipo)
  
  temp2 <- temp1 %>% dplyr::select(r_esp, r_gen) %>% 
    pivot_longer(cols=c(r_esp, r_gen), names_to = "type", values_to = "r_value") 
  #temp2 <- rbind(prior_g_df, prior_e_df, temp1)
  temp2$scenario <- sub(".*_", "", names.temp[[i]])
  temp2 %<>% dplyr::mutate(plot_week=sub("_.*", "", names.temp[[i]]),
                           scenario=sub(".*_", "", names.temp[[i]]))
  
  temp2 %<>% dplyr::mutate(sp_gen_esp = substr(type, 3, 5)) %>% 
    merge(abundance_sum, by="sp_gen_esp")
  
  myfiles[[i]] <- temp2
}


all.plots.weeks <- do.call("rbind", myfiles)

all.plots.weeks <- unique(all.plots.weeks)


all.plots.weeks %<>% dplyr::filter(scenario=="random")
test <- aggregate(r_value ~ plot_week + type + sp_gen_esp + total_abundance, data = all.plots.weeks, median)

test %<>% filter(!plot_week %in% c("E3", "F3"))
r.pol.abun <- ggplot(test, aes(x = total_abundance, y = r_value, color = type, label=plot_week)) +
  geom_point(size=3, alpha=0.5) + geom_text(hjust=-0.5, vjust=0, show_guide = F, size=2.5) +
  stat_ellipse(type = "t", linetype = 2) +
  scale_color_manual(values=c("#2e4057", "#66a182"), 
                     labels=c("Occasional \nforagers", "Frequent \nforagers")) +
  labs(x = "Interaction frequency", y = "r scale (median)", color = "type") +
  theme_bw() + guides(label="none") + theme(legend.title=element_blank(),
                                            legend.key.size = unit(0.5, "cm"),
                                            legend.text = element_text(size=9),
                                            axis.title = element_text(size=16),
                                            axis.text = element_text(size=14),
                                            legend.position = c(0.98, 0.98),      
                                            legend.justification = c(1, 1)) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE, 
                              keyheight = unit(2, "lines"),
                              fill = guide_legend(ncol = 1, byrow = TRUE, 
                                                  keyheight = unit(2, "lines"))))



# Create Fig. 3

ggarrange(
  NULL, r.pol, NULL,         
  r.pol.abun, 
  nrow = 1, widths = c(0.1, 1, 0.1, 1), vjust=1.5, hjust=0.5,
  labels = c("", "A","", "B"), font.label=list(color="black",size=20) 
) 

ggarrange(
  r.pol,       
  r.pol.abun, 
  nrow = 2, 
  labels = c( "A", "B"), font.label=list(color="black",size=20) 
) 



# Create Fig. 4

# Read acceptance rates of the Bayesian Approximate Computation (ABC)
acc.rate <- read.csv("Output_model/metrics.csv", dec = ".")
glimpse(acc.rate)

acc.rate %<>% dplyr::select(-Random,-Regular,-XY) %>% 
  pivot_longer(cols=c(Tasa_Random, Tasa_XY, Tasa_Regular), 
               names_to = "scenario", values_to = "rate") %>%
  dplyr::filter(Plot_slice!="F1")

##### test overall differences between spatial scenarios

rand_xy <- subset(acc.rate, scenario %in% c("Tasa_Random", "Tasa_XY")) %>% arrange(desc(scenario))
rand_reg <- subset(acc.rate, scenario %in% c("Tasa_Random", "Tasa_Regular"))
xy_reg <- subset(acc.rate, scenario %in% c("Tasa_Regular", "Tasa_XY")) %>% arrange(desc(scenario))

wilcox.test(rate ~ scenario, data = rand_xy, paired = TRUE)
wilcox.test(rate ~ scenario, data = rand_reg, paired = TRUE)
wilcox.test(rate ~ scenario, data = xy_reg, paired = TRUE)
rand_xy %>% friedman.test(rate ~ scenario)

friedman.test(rate ~ scenario | Plot_slice, data = acc.rate)


# Panel A

test1 <- acc.rate %>%
  dplyr::group_by(Plot_slice) %>%
  dplyr::mutate(rate = scale(rate)) %>%
  ungroup()

test1$scenario <- factor(test1$scenario, levels = c("Tasa_XY", "Tasa_Random", "Tasa_Regular"))

scaled.plot <- ggplot(test1, aes(x = scenario, y = rate, fill=scenario, color=scenario)) +
  geom_violin(trim = FALSE, alpha=0.2) +
  geom_point(size = 2, position = position_jitter(w=0.05), alpha=0.5) +
  geom_boxplot(width=.05, outlier.shape=NA, alpha=0) +
  scale_color_manual(values=c("#D26635", "#971E2E", "#561942")) +
  scale_fill_manual(values=c("#D26635", "#971E2E", "#561942")) +
  theme_bw(base_size=20) + 
  labs(x = "Scenario", y = "Acceptance rate") +
  scale_x_discrete(labels=c("Observed", "Random", "Regular")) + 
  theme(legend.title=element_blank(), axis.title.x = element_blank(),
        axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        axis.text.x = element_text(vjust = -1),
        plot.margin = margin(0, 0, 10, 0, "pt")) +
  guides(color="none", fill="none") 


# Panel B

# Plot aggregated spatial distribution
set.seed(3)
n_aggregated <- 30
x_aggregated <- c(rnorm(15, mean = 10, sd = 3), rnorm(15, mean = 20, sd = 3))
y_aggregated <- c(rnorm(15, mean = 10, sd = 3), rnorm(15, mean = 20, sd = 3))

aggregated_data <- data.frame(x = x_aggregated, y = y_aggregated)
aggregated_data$sp <- sample(rep(c("a", "b", "c"), each = 10))

p1 <- ggplot(aggregated_data, aes(x, y, shape=sp)) +
  geom_point(color = "#D26635", size = 2.5, alpha=0.6, stroke = 4) +
  xlab("Longitude") + ylab("Latitude") +
  theme_minimal(base_size = 20) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(),  # Remove axis ticks
    plot.background = element_rect(fill = "white", color = "white"),  # White background with black margins
    panel.border = element_rect(color = "black", fill=NA),
    axis.title = element_text(size=16),
    plot.margin = margin(0, 4, 0, 0, "pt")  # Black border around the plot
  )+ guides(shape=FALSE)

# Plot random spatial distribution
set.seed(2)
n_random <- 30
x_random <- runif(n_random, 1, 30)
y_random <- runif(n_random, 1, 30)

random_data <- data.frame(x = x_random, y = y_random)
random_data$sp <- sample(rep(c("a", "b", "c"), each = 10))

p2 <- ggplot(random_data, aes(x, y, shape=sp)) +
  geom_point(color = "#971E2E", size = 2.5, alpha=0.6, stroke = 4) +
  xlab("Longitude") + ylab("Latitude") +
  theme_minimal(base_size = 20) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(),  # Remove axis ticks
    plot.background = element_rect(fill = "white", color = "white"),  # White background with black margins
    panel.border = element_rect(color = "black", fill=NA),
    axis.title = element_text(size=16),
    plot.margin = margin(0, 4, 0, 0, "pt")  # Black border around the plot
  )+ guides(shape=FALSE)

# Plot regular spatial distribution (grid)
n <- 5  # Number of points along each axis
grid_data <- expand.grid(x = seq(1, n), y = seq(1, n))
grid_data$sp <- c(sample(rep(c("a", "b", "c"), each = 8)), "a")

p3 <- ggplot(grid_data, aes(x, y, shape=sp)) +
  geom_point(color = "#561942", size = 2, alpha=0.6, stroke = 4) +
  xlab("Longitude") + ylab("Latitude") +
  theme_minimal(base_size = 20) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(),   # Remove axis ticks
    plot.background = element_rect(fill = "white", color = "white"),  # White background with black margins
    panel.border = element_rect(color = "black", fill=NA),
    axis.title = element_text(size=16),
    plot.margin = margin(0, 4, 0, 0, "pt")  # Black border around the plot
  ) + guides(shape=FALSE)

# Display plots
print(p1)
print(p2)
print(p3)

ggarrange(scaled.plot,
          ggarrange(p1, p2, p3, nrow=1, labels="B"),
          nrow=2, labels="A"
)

scaled.plot / ((p1+ p2+ p3) + plot_layout(axes="collect"))


#### Temporal variation in predictive ability based on acceptance rates

temp.data <- test1 %>%
  mutate(
    plot = substr(Plot_slice, 1, 1),  
    slice = as.integer(substr(Plot_slice, 2, 2))  
  ) %>% arrange(slice) %>%
  mutate(slice=as.factor(slice)) %>% filter(scenario=="Tasa_XY")


temp.data$rate <- as.vector(temp.data$rate)
glimpse(temp.data)


summary_data <- temp.data %>%
  dplyr::group_by(slice) %>%
  dplyr::summarise(mean_rate = mean(rate),
                   min_rate = mean(rate) -sd(rate),
                   max_rate = mean(rate) +sd(rate))


plot.time <- ggplot(summary_data, aes(x = slice, color=slice)) +
  geom_pointrange(aes(y = mean_rate, ymin = min_rate, ymax = max_rate), size=1.5, fatten = 2) +
  theme_bw() +
  scale_color_manual(values=c("#35757F", "#35757F", "#35757F")) + 
  labs(x = "Scenario", y = "Acceptance rate") +
  scale_x_discrete(labels=c("Early", "Peak", "Late")) + 
  theme(legend.title=element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = -1),
        axis.text = element_text(size=15), axis.title.y = element_text(size=15)) +
  guides(color="none")


temp.data <- test1 %>%
  mutate(
    plot = substr(Plot_slice, 1, 1),  # Extract the letter (e.g., A, B, C)
    slice = as.integer(substr(Plot_slice, 2, 2))  # Extract the number (e.g., 1, 2, 3)
  ) %>% arrange(slice) %>%
  mutate(slice=as.factor(slice)) %>%
  dplyr::group_by(scenario, slice) %>%
  dplyr::summarise(mean_rate = mean(rate))

temp.data$scenario <- factor(temp.data$scenario, levels = c("Tasa_XY", "Tasa_Random", "Tasa_Regular"))


# Create Fig. 5

col1 = "#F8F1D8"
col2 = "#35757F"

plot.space.time <- ggplot(temp.data, aes(x = slice, y = scenario, fill = mean_rate)) +
  geom_tile(color = "white") +  # Adding tile border for better visibility
  scale_fill_gradient(low = col1, high = col2,
                      name = "Mean Rate") +
  labs(x = "Slice",
       y = "Scenario") +
  scale_x_discrete(labels = c("Early", "Peak", "Late")) +
  scale_y_discrete(labels = c("Observed", "Random", "Regular")) +
  theme_minimal() +
  theme(axis.text = element_text(size=15, hjust = 0.5), 
        legend.title= element_text(size=15, hjust = 0.5),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12),
        axis.title = element_blank()) +
  coord_fixed()  # Keep the aspect ratio square


plot.space.time + plot.time + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size=20))



# Example data for conceptual figure
df <- data.frame(
  degree = c(
    1.0, 3.3, 5.0, 6.2, 7.3, 8.0, 8.7, 9.3, 10.0, 10.7,
    11.3, 11.9, 12.5, 13.0, 13.5, 14.0, 14.4, 14.9, 15.4,
    15.9, 16.3, 16.8, 17.2, 17.8, 18.3, 18.8, 19.4, 20.0,
    20.8, 21.6, 22.2, 22.8, 23.5, 24.2, 25.0, 26.0, 27.0,
    28.0, 31.5
  ),
  n_species = c(
    39, 33, 29, 26, 23, 21, 20, 19, 17, 16,
    15, 13, 13, 12, 11, 9, 8, 10, 8,
    6, 7, 5, 6, 5, 2, 3, 2, 2,
    1, 1, 1, 1, 1, 1, 1, 4, 1,
    1, 1
  )
)

ggplot(df, aes(degree, n_species)) +
  geom_point(size = 3.2, colour = "grey40", alpha = 0.7) +
  labs(
    x = "Species' degree\n(i.e., number of unique interactions)",
    y = "Number of species"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 18),
    axis.line = element_line(linewidth = 0.8),
    panel.grid = element_blank()
  )

