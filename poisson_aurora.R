library(ggplot2)
library(dplyr)
library(tidyr)
library(isotone)
library(REBayes)

## Code for Figure showing replicates produced by data fission

double_poissons <- function(z, tau = 0.1) {
  thinning <- rbinom(length(z), z, tau)
  x <- (z - thinning) / (1 - tau)
  y <- thinning / tau
  data.frame(x = x, y = y)
}

three_component_poisson_sim <- function(n, effect_sizes = c(1, 4, 7)) {
  mus <- sample(effect_sizes, size = n, replace = TRUE)
  zs <- rpois(n, lambda = mus)
  list(true_mu = mus, z = zs)
}

bayes_three_component_poisson <- function(zs, tau = 0, effect_sizes = c(1, 4, 7)) {
  zs <- zs * (1-tau) # rescale, zs ~ Pois((1-tau)*theta_i)
  effect_sizes <- effect_sizes * (1-tau)
  num <- sapply(zs, function(z) sum(effect_sizes * dpois(z, effect_sizes)))
  denom <- sapply(zs, function(z) sum(dpois(z, effect_sizes)))
  num / (denom * (1-tau))
}

set.seed(1000)
n <- 1000
effect_sizes <- c(1, 4, 7)
tau_values <- c(0.04, 0.14)
sim_df <- three_component_poisson_sim(n, effect_sizes)

plot_data <- lapply(tau_values, function(tau) {
  doubled <- double_poissons(sim_df$z, tau = tau)
  
  oracle_pred <- bayes_three_component_poisson(doubled$x, tau = tau, effect_sizes = effect_sizes)
  
  iso_fit <- gpava(doubled$x, doubled$y, ties = "secondary")
  iso_pred <- iso_fit$x
  
  data.frame(
    x = doubled$x,  
    y = doubled$y,
    oracle = oracle_pred,
    isotonic = iso_pred,
    tau = tau
  )
})
all_data <- do.call(rbind, plot_data)

summary_data <- all_data %>%
  group_by(tau, x, y) %>%
  summarize(
    y = mean(y),
    oracle = mean(oracle),
    isotonic = mean(isotonic),
    count = n()
  ) %>%
  ungroup()

x_lim <- c(0, max(summary_data$x) * 1.1)
y_lim <- c(0, max(c(summary_data$y, summary_data$oracle, summary_data$isotonic)) * 1.2)


base_size <- 26

tau_labeller <- function(variable, value) {
  return(lapply(value, function(x) bquote(tau == .(x))))
}

identity_line_data <- summary_data %>%
  group_by(tau) %>%
  summarize(x_min = min(x), x_max = max(x), y_min = min(x), y_max = max(x))


p <- ggplot(summary_data, aes(x = x)) +
  geom_point(aes(y = y, size = count), alpha = 0.5, color = "grey40") +
  geom_line(aes(y = oracle, color = "Oracle", linetype = "Oracle", linewidth = "Oracle"), alpha = 0.75) +
  geom_line(aes(y = isotonic, color = "Isotonic", linetype = "Isotonic", linewidth = "Isotonic"),  alpha = 0.75) +
  geom_segment(data = identity_line_data, aes(x = x_min, xend = x_max, y = y_min, yend = y_max, 
                                              color = "Identity", linetype = "Identity", linewidth = "Identity"), 
               alpha = 0.75) +
  facet_wrap(~tau, labeller = tau_labeller, scales = "free", ncol = 2) +
  theme_minimal(base_size = base_size) +
  labs(
    x = expression(f[tau](X[i])),
    y = expression(g[tau](X[i])),
    color = NULL,
    linetype = NULL,
    linewidth = NULL,
    size = "Count"
  ) +
  scale_color_manual(
    values = c("Oracle" = "black", "Isotonic" = "#1B9E77", "Identity" = "#666666"),
    labels = c("Oracle" = expression(paste("", E, "[", theta[i], " | ", f[tau](X[i]), "]")),
               "Isotonic" = "Isotonic Regression",
               "Identity" = "Identity Line"),
    breaks = c("Oracle", "Isotonic", "Identity") 
  ) +
  scale_linetype_manual(
    values = c("Oracle" = "11", "Isotonic" = "solid", "Identity" = "dashed"),
    labels = c("Oracle" = expression(paste("", E, "[", theta[i], " | ", f[tau](X[i]), "]")),
               "Isotonic" = "Isotonic Regression",
               "Identity" = "Identity Line"),
    breaks = c("Oracle", "Isotonic", "Identity") 
  ) +
  scale_linewidth_manual(
    values = c("Oracle" = 1.9, "Isotonic" = 2.0, "Identity" = 1.2),
    labels = c("Oracle" = expression(paste("", E, "[", theta[i], " | ", f[tau](X[i]), "]")),
               "Isotonic Regression" = "Isotonic Regression",
               "Identity" = "Identity Line"),
    breaks = c("Oracle", "Isotonic", "Identity") 
  ) +
  scale_size_continuous(range = c(0.5, 6)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = base_size * 0.8),  
    legend.text = element_text(size = base_size * 0.8),
    legend.key.width = unit(3, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.spacing.y = unit(0.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", color = "black", size = 1.5),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = base_size * 1.4),
    axis.text = element_text(size = base_size * 0.9),
    axis.title = element_text(size = base_size * 1.3, face = "bold"),
    axis.line = element_line(size = 1, color = "black"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
  ) +
  coord_cartesian(xlim = x_lim, ylim = y_lim, expand = TRUE) +
  guides(
    color = guide_legend(order = 1, title = NULL, override.aes = list(linetype=c("11","solid","dashed"), linewidth = c(1.9, 2.0, 1.2))),
    linewidth = "none",
    linetype = "none",
    size = guide_legend(order = 2, title = "Count")  
  )

ggsave("double_poisson_plot.png", plot = p, width = 18/1.11, height = 7.5/1.1, units = "in", dpi = 600, bg = "white")


## Simulation study (for Table)


rao_blackwell_isotonic <- function(z, tau, n_iterations = 100) {
  predictions <- replicate(n_iterations, {
    doubled <- double_poissons(z, tau = tau)
    iso_fit <- gpava(doubled$x, doubled$y, ties = "secondary")
    iso_fit$x
  })
  rowMeans(predictions)
}

npmle_poisson <- function(z) {
  npmle_fit <- Pmix(z, exposure = rep(1, length(z)))
  npmle_fit$dy
}

run_simulation <- function(n = 1000, tau_values, n_iterations = 100) {
  sim_df <- three_component_poisson_sim(n)
  z <- sim_df$z
  true_mu <- sim_df$true_mu
  
  results <- list()
  
  results[["Bayes_Rule"]] <- mean((bayes_three_component_poisson(z) - true_mu)^2)
  results[["MLE"]] <- mean((z - true_mu)^2)
  results[["NPMLE"]] <- mean((npmle_poisson(z) - true_mu)^2)
  
  for (tau in tau_values) {
    doubled <- double_poissons(z, tau = tau)
    
    iso_fit <- gpava(doubled$x, doubled$y, ties = "secondary")
    results[[paste0("Isotonic_", tau)]] <- mean((iso_fit$x - true_mu)^2)
    
    rb_iso_pred <- rao_blackwell_isotonic(z, tau, n_iterations)
    results[[paste0("RB_Isotonic_", tau)]] <- mean((rb_iso_pred - true_mu)^2)
  }
  
  unlist(results)
}

set.seed(42)
n_sims <-100

results <- replicate(n_sims, run_simulation(n=1000, tau_values = tau_values), simplify = FALSE)

process_results <- function(results) {
  df <- data.frame(
    Method = names(results[[1]]),
    Mean_MSE = rowMeans(do.call(cbind, results))
  )
}

result_summary <- process_results(results)
write.table(result_summary[c(1,2,3,5,7),], file="double_poisson_MSE.csv")