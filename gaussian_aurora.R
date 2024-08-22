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

p <- ggplot(all_points, aes(x = x, y = y)) + 
  geom_point(alpha = 0.3, color = "grey50", size = 1.8) +
  geom_line(data = all_lines, aes(x = x, y = y, color = method), size = 2, alpha=0.8) + 
  geom_abline(aes(intercept = 0, slope = 1, linetype = "Identity"), color = "black", alpha = 0.5, size = 1.8) +
  facet_wrap(~tau, labeller = tau_labeller, scales = "free", ncol = 2) +
  theme_minimal(base_size = base_size) +
  labs(
    x = expression(f[tau](X[i])),
    y = expression(g[tau](X[i])),
    color = NULL,
    linetype = NULL
  ) +
  scale_color_manual(
    values = c("Oracle" = "#0072B2", "Isotonic Regression" = "#D55E00"),
    labels = c("Oracle" = expression(paste("", E, "[", theta[i], " | ", f[tau](X[i]), "]")),
               "Isotonic Regression" = "Isotonic Regression"),
    breaks = c("Oracle", "Isotonic Regression")
  ) +
  scale_linetype_manual(values = c("Identity" = "dashed"), 
                        labels = c("Identity" = "Identity Line")) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 20),
    legend.title = element_text(size = base_size * 1.1),
    legend.text = element_text(size = base_size * 1.1),
    legend.key.size = unit(2.2, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.6),
    panel.background = element_rect(fill = "white", color = "black", size = 1.8),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = base_size * 1.4),
    plot.title = element_text(hjust = 0.5, face = "bold", size = base_size * 1.6),
    plot.subtitle = element_text(hjust = 0.5, size = base_size * 1.3),
    plot.margin = margin(35, 35, 35, 35),
    axis.text = element_text(size = base_size * 1.1, face = "bold"),
    axis.title = element_text(size = base_size * 1.3, face = "bold"),
    axis.line = element_line(size = 1.2, color = "black")
  ) +
  coord_cartesian(xlim = x_range, ylim = y_range, expand = TRUE) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 5)),
         linetype = guide_legend(order = 2))

ggsave("double_gaussian_plot.png", plot = p, width = 20, height = 8, units = "in", dpi = 300, bg = "white")

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
print(result_summary)