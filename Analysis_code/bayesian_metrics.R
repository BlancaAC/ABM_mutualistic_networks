# 21-05-2025

model1 <- read.csv("Data/bayesian_metrics/abserrors_randomproof.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_1") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE )

model2 <- read.csv("Data/bayesian_metrics/abserrors_model2.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_2") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE )

model3 <- read.csv("Data/bayesian_metrics/abserrors_1random.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_3") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE )

model4 <- read.csv("Data/bayesian_metrics/abserrors_1XY.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_4") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE )

model4.reg <- read.csv("Data/bayesian_metrics/abserrors_1regular.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_4_reg") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE )

models <- rbind(model1, model2, model3, model4)

models.spat <- rbind(model3, model4, model4.reg)

ggplot(models.spat, aes(x = model, y = total_abs_error, fill=model, color=model)) +
  geom_violin(trim = FALSE, alpha=0.2) +
  geom_point(size = 2, position = position_jitter(w=0.05), alpha=0.5) +
  geom_boxplot(width=.05, outlier.shape=NA, alpha=0) +
  scale_color_manual(values=c("#D26635", "#971E2E", "#561942")) +
  scale_fill_manual(values=c("#D26635", "#971E2E", "#561942")) +
  theme_bw(base_size=20) + 
  labs(x = "Scenario", y = "Acceptance rate") +
  scale_x_discrete(labels=c("XY", "Random", "Regular")) + 
  theme(legend.title=element_blank(), axis.title.x = element_blank(),
        axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        axis.text.x = element_text(vjust = -1),
        plot.margin = margin(0, 0, 6, 0, "pt")) +
  guides(color="none", fill="none") 



# Ensure models are in the correct order
models <- models %>%
  mutate(
    model = factor(model, levels = paste0("model_", 1:4)),
    week = factor(week)
  )

ggplot(models, aes(x = model, y = total_abs_error, 
                   color = plot, shape = week, 
                   group = interaction(plot, week))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#264653","#2a9d8f","#8ab17d","#e9c46a","#f4a261","#e76f51")) +
  scale_shape_manual(values = c(15, 17, 19)) +  # square, triangle, circle
  labs(
    x = "",
    y = "Total Absolute Error (RMAE)",
    color = "Plot",
    shape = "Week"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size=14),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )


models <- models %>%
  mutate(model_num = as.numeric(factor(model, levels = paste0("model_", 1:4))))

# Random intercepts for plot and week
lmm <- lmer(total_abs_error ~ model_num + (1 | plot) + (1 | week), data = models)
summary(lmm)




# Make sure model is a factor in the correct order
models <- models %>%
  mutate(model = factor(model, levels = paste0("model_", 1:4)))

# Fit the mixed model
mod_lmm <- lmer(total_abs_error ~ model + (1 | plot_week), data = models)
summary(mod_lmm)

# Estimated marginal means for each model level
emm <- emmeans(mod_lmm, ~ model)

# Compare adjacent levels only
contrast(emm, method = "consec", adjust = "none")



# with errors within network


model1 <- read.csv("Data/bayesian_metrics/with_errors/abserrors_randomproof.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_1") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE ) %>% 
  group_by(plot_week) %>%
  #dplyr::filter(mean_error==min(mean_error)) %>%
  dplyr::summarise(
    mean_mean_error = sum(mean_error, na.rm = TRUE),
    mean_std_error = sum(std_error, na.rm = TRUE),
    plot = first(plot),
    week = first(week),
    model = first(model),
    .groups = "drop"
  )



model2 <- read.csv("Data/bayesian_metrics/with_errors/abserrors_model2.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_2") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE ) %>% 
  group_by(plot_week) %>%
  #dplyr::filter(mean_error==min(mean_error)) %>%
  dplyr::summarise(
    mean_mean_error = sum(mean_error, na.rm = TRUE),
    mean_std_error = sum(std_error, na.rm = TRUE),
    plot = first(plot),
    week = first(week),
    model = first(model),
    .groups = "drop"
  )

model3 <- read.csv("Data/bayesian_metrics/with_errors/abserrors_regular.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_3") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE ) %>% 
  group_by(plot_week) %>%
  #dplyr::filter(mean_error==min(mean_error)) %>%
  dplyr::summarise(
    mean_mean_error = sum(mean_error, na.rm = TRUE),
    mean_std_error = sum(std_error, na.rm = TRUE),
    plot = first(plot),
    week = first(week),
    model = first(model),
    .groups = "drop"
  )

model4 <- read.csv("Data/bayesian_metrics/with_errors/abserrors_XY.csv", dec = ".") %>%
  dplyr::filter(metric=="RMAE") %>%
  mutate(plot_week = str_extract(plot, "(?<=_)[^_]+(?=_)"), model= "model_4") %>%
  separate(plot_week, into = c("plot", "week"), sep = "(?<=^[A-Z])", remove = FALSE ) %>% 
  group_by(plot_week) %>%
  #dplyr::filter(mean_error==min(mean_error)) %>%
  dplyr::summarise(
    mean_mean_error = sum(mean_error, na.rm = TRUE),
    mean_std_error = sum(std_error, na.rm = TRUE),
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



ggplot(models, aes(x = model, y = mean_mean_error, 
                   color = plot, shape = week, 
                   group = interaction(plot, week))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = mean_mean_error - mean_std_error,
                    ymax = mean_mean_error + mean_std_error),
                width = 0.2, linewidth = 0.6) +
  scale_color_manual(values = c("#264653","#2a9d8f","#8ab17d","#e9c46a","#f4a261","#e76f51")) +
  scale_shape_manual(values = c(15, 17, 19)) +  # square, triangle, circle
  labs(
    x = "",
    y = "Total Absolute Error (RMAE)",
    color = "Plot",
    shape = "Week"
  ) +
  facet_grid(plot ~week) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 14),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) + scale_y_log10()


# z-scores

model1 <- read.csv("Data/bayesian_metrics/with_errors/abserrors_zscore_randomproof.csv", dec = ".") %>%
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


model2 <- read.csv("Data/bayesian_metrics/with_errors/abserrors_zscore_model2.csv", dec = ".") %>%
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

model3 <- read.csv("Data/bayesian_metrics/with_errors/abserrors_zscore_1regular.csv", dec = ".") %>%
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

model4 <- read.csv("Data/bayesian_metrics/with_errors/abserrors_zscore_1XY.csv", dec = ".") %>%
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
  geom_point(size = 1, alpha=0.8) +
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
    axis.text.x = element_text(size=10, angle=45, hjust = 1)
  ) + scale_x_discrete(breaks=c("model_1","model_2","model_3", "model_4"),
                       labels=c("Baseline", "Species \nabundance", 
                                "Pollinator \nbehavior", "Spatial \nconfiguration")) + 
  scale_y_log10()



# Fit the mixed model
hist(log(models$mean_zscore))
mod_lmm <- lmer(log(mean_zscore) ~ model + (1 | plot_week), data = models)
summary(mod_lmm)

# Estimated marginal means for each model level
emm <- emmeans(mod_lmm, ~ model)

# Compare adjacent levels only
contrast(emm, method = "consec", adjust = "none")
