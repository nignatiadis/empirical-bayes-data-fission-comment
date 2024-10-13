library(ggplot2)
library(dplyr)
library(tidyr)
library(isotone)
library(REBayes)

# Helper functions
double_gaussians <- function(z, tau = 0.1, epsilon_noise = rnorm(length(z))) {
  data.frame(x = z + epsilon_noise * tau, y = z - epsilon_noise / tau)
}

three_component_normal_sim <- function(n, sigma = 1, effect_sizes = 3 * (-1:1) + 4) {
  mus <- sample(effect_sizes, size = n, replace = TRUE)
  zs <- mus + rnorm(n, sd = sigma)
  list(true_mu = mus, z = zs)
}

bayes_three_component <- function(zs, sigma = 1, tau = 0, effect_sizes = 3 * (-1:1) + 4) {
  sigma_marginal <- sqrt(sigma^2 + tau^2)
  num <- sapply(zs, function(z) sum(effect_sizes * dnorm(z, effect_sizes, sigma_marginal)))
  denom <- sapply(zs, function(z) sum(dnorm(z, effect_sizes, sigma_marginal)))
  num / denom
}

# Main simulation
set.seed(1000)
n <- 1000
tau_values <- c(0.2, 0.4)
sim_df <- three_component_normal_sim(n)
epsilon_noise <- rnorm(n)


plot_data <- lapply(tau_values, function(tau) {
  doubled <- double_gaussians(sim_df$z, tau = tau, epsilon_noise = epsilon_noise)
  
  oracle_pred <- data.frame(
    x = doubled$x,
    y = bayes_three_component(doubled$x, tau = tau),
    method = "Oracle"
  )
  
  iso_fit <- gpava(doubled$x, doubled$y, ties="secondary")
  iso_pred <- data.frame(
    x = doubled$x,
    y = iso_fit$x,
    method = "Isotonic Regression"
  )
  
  doubled_preds <- rbind(oracle_pred, iso_pred)
  
  list(
    points = mutate(doubled, tau = tau),
    lines = mutate(doubled_preds, tau = tau)
  )
})

all_points <- do.call(rbind, lapply(plot_data, function(x) x$points))
all_lines <- do.call(rbind, lapply(plot_data, function(x) x$lines))

x_range <- range(all_points$x, all_lines$x)
y_range <- range(all_points$y, all_lines$y)

x_padding <- (max(x_range) - min(x_range)) * 0.05 
x_range <- c(min(x_range) - x_padding, max(x_range) + x_padding)

# Plotting
base_size <- 26

tau_labeller <- function(variable, value) {
  return(lapply(value, function(x) bquote(tau == .(x))))
}

identity_line <- group_by(all_lines, tau) %>% 
  summarize(
    x = range(c(x,y)),
    y = x,
    method = "Identity"
  ) %>%
  ungroup()


all_lines <- rbind(all_lines, identity_line) %>%
  mutate(method = factor(method, levels=c("Oracle", "Isotonic Regression", "Identity")))

p <- ggplot(all_points, aes(x = x, y = y)) + 
  geom_point(alpha = 0.28, color = "grey40", size = 1.5, stroke=NA) +
  geom_line(data = all_lines, aes(x = x, y = y, color = method, linetype = method, linewidth = method), alpha=0.75) +
  facet_wrap(~tau, labeller = tau_labeller, scales = "free", ncol = 2) +
  theme_minimal(base_size = base_size) +
  labs(
    x = expression(f[tau](X[i])),
    y = expression(g[tau](X[i])),
    color = NULL,
    linetype = NULL,
    linewidth = NULL
  ) +
  scale_color_manual(
    values = c("Oracle" = "black", "Isotonic Regression" = "#1B9E77" , "Identity" = "#666666"),
    labels = c("Oracle" = expression(paste("", E, "[", theta[i], " | ", f[tau](X[i]), "]")),
               "Isotonic Regression" = "Isotonic Regression",
               "Identity" = "Identity Line")
  ) +
  scale_linetype_manual(
    values = c("11", "solid", "dashed"),
    labels = c("Oracle" = expression(paste("", E, "[", theta[i], " | ", f[tau](X[i]), "]")),
               "Isotonic Regression" = "Isotonic Regression",
               "Identity" = "Identity Line")
  ) +
  scale_linewidth_manual(
    values = c(1.9, 2.0, 1.2),
    labels = c("Oracle" = expression(paste("", E, "[", theta[i], " | ", f[tau](X[i]), "]")),
               "Isotonic Regression" = "Isotonic Regression",
               "Identity" = "Identity Line")
  ) + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
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
  coord_cartesian(xlim = x_range, ylim = c(-11, 20), expand = TRUE) +
  guides(
    color = guide_legend(override.aes = list(linewidth = c(1.9, 2.0, 1.2))),
    linewidth = "none"
  )

ggsave("double_gaussian_plot.png", plot = p, width = 18/1.1, height = 7.5/1.1, units = "in", dpi = 600, bg = "white")


# Simulation study
rao_blackwell_isotonic <- function(z, tau, n_iterations = 100) {
  predictions <- replicate(n_iterations, {
    epsilon_noise <- rnorm(length(z))
    doubled <- double_gaussians(z, tau = tau, epsilon_noise = epsilon_noise)
    iso_fit <- gpava(doubled$x, doubled$y, ties="secondary")
    iso_fit$x
  })
  rowMeans(predictions)
}

npmle_gaussian <- function(z) {
  npmle_fit <- GLmix(z)
  predict(npmle_fit, z)
}

run_simulation <- function(n = 1000, tau_values, n_iterations = 100) {
  sim_df <- three_component_normal_sim(n)
  z <- sim_df$z
  true_mu <- sim_df$true_mu
  
  results <- list()
  
  results[["Bayes_Rule"]] <- mean((bayes_three_component(z) - true_mu)^2)
  results[["NPMLE"]] <- mean((npmle_gaussian(z) - true_mu)^2)
  results[["MLE"]] <- mean((z - true_mu)^2)
  
  for (tau in tau_values) {
    doubled <- double_gaussians(z, tau = tau)
    iso_fit <- gpava(doubled$x, doubled$y, ties="secondary")
    iso_pred <- iso_fit$x
    results[[paste0("Isotonic_", tau)]] <- mean((iso_pred - true_mu)^2)
    
    rb_iso_pred <- rao_blackwell_isotonic(z, tau, n_iterations)
    results[[paste0("RB_Isotonic_", tau)]] <- mean((rb_iso_pred - true_mu)^2)
  }
  
  unlist(results)
}

set.seed(42)
n_sims <- 100
tau_values <- c(0.2, 0.4)

results <- replicate(n_sims, run_simulation(tau_values = tau_values), simplify = FALSE)

process_results <- function(results) {
  df <- data.frame(
    Method = names(results[[1]]),
    Mean_MSE = rowMeans(do.call(cbind, results))
  )
  df
}

result_summary <- process_results(results)
write.table(result_summary[c(1,2,3,5,7),], file="double_gaussian_MSE.csv")