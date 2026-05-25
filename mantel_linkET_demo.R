# ======================================================================
# Mantel Test Visualization using linkET
# Demo script with synthetic data
# 
# Author:  [Your Name]
# GitHub:  [Your GitHub URL]
# License: MIT
#
# Description:
#   This script demonstrates how to perform pairwise Mantel tests
#   between two variable groups (e.g., physiological parameters and
#   gene expression data) and visualize the results using the linkET
#   package. Data used here are synthetic (demo only) and do not
#   represent any real experimental measurements.
#
# Requirements:
#   R >= 4.1.0
#   Packages: tidyverse, linkET, RColorBrewer, vegan
#
# Install packages (run once):
#   install.packages(c("tidyverse", "RColorBrewer", "vegan"))
#   install.packages("devtools")
#   devtools::install_github("Hy4m/linkET")
# ======================================================================


# ----------------------------------------------------------------------
# 0. Load required libraries
# ----------------------------------------------------------------------
required_pkgs <- c("tidyverse", "linkET", "RColorBrewer", "vegan")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0(
      "Package '", pkg, "' is not installed.\n",
      "Please run: install.packages('", pkg, "')\n",
      "For linkET: devtools::install_github('Hy4m/linkET')"
    ))
  }
}

library(tidyverse)
library(linkET)
library(RColorBrewer)
library(vegan)

cat("All packages loaded successfully.\n")


# ----------------------------------------------------------------------
# 1. Define experimental structure
#    - Two treatment levels: Control (Ctrl) and Treatment (Trt)
#    - Four dose levels: Ctrl, Dose_25, Dose_50, Dose_75
#    - Three biological replicates per combination
# ----------------------------------------------------------------------

treatment_levels <- c("Ctrl", "Trt")
dose_levels      <- c("Ctrl", "Dose_25", "Dose_50", "Dose_75")
n_reps           <- 3

n_total <- length(treatment_levels) * length(dose_levels) * n_reps  # 24 observations


# ----------------------------------------------------------------------
# 2. Synthetic physiological parameter data
#    (Param_1 through Param_5)
#    Values are illustrative and drawn from plausible biological ranges.
# ----------------------------------------------------------------------

set.seed(42)  # for full reproducibility

# Helper: small jitter around a mean value
jitter_vals <- function(mean_val, n = 3, sd = 2) {
  round(rnorm(n, mean = mean_val, sd = sd), 2)
}

physio_data <- data.frame(
  Treatment = c(rep("Ctrl", 12), rep("Trt", 12)),
  Dose      = rep(rep(dose_levels, n_reps), 2),
  Rep       = rep(rep(1:n_reps, each = length(dose_levels)), 2),

  # Param_1: increases with treatment and dose (e.g., antioxidant enzyme)
  Param_1 = c(
    jitter_vals(88),  jitter_vals(95),  jitter_vals(97),  jitter_vals(102), # Ctrl
    jitter_vals(140), jitter_vals(150), jitter_vals(168), jitter_vals(172)  # Trt
  ),

  # Param_2: moderate increase with treatment
  Param_2 = c(
    jitter_vals(25, sd = 1.5), jitter_vals(31, sd = 1.5),
    jitter_vals(33, sd = 1.5), jitter_vals(37, sd = 1.5),
    jitter_vals(38, sd = 1.5), jitter_vals(43, sd = 1.5),
    jitter_vals(46, sd = 1.5), jitter_vals(50, sd = 1.5)
  ),

  # Param_3: relatively stable, small treatment effect
  Param_3 = c(
    jitter_vals(155, sd = 1), jitter_vals(158, sd = 1),
    jitter_vals(159, sd = 1), jitter_vals(163, sd = 1),
    jitter_vals(172, sd = 1), jitter_vals(181, sd = 1),
    jitter_vals(184, sd = 1), jitter_vals(196, sd = 1)
  ),

  # Param_4: decreases with dose under treatment (e.g., lipid peroxidation)
  Param_4 = c(
    jitter_vals(13.5, sd = 0.5), jitter_vals(11.8, sd = 0.5),
    jitter_vals(11.3, sd = 0.5), jitter_vals(9.8, sd = 0.5),
    jitter_vals(28.5, sd = 0.5), jitter_vals(22.5, sd = 0.5),
    jitter_vals(21.0, sd = 0.5), jitter_vals(16.5, sd = 0.5)
  ),

  # Param_5: strong increase with treatment and dose
  Param_5 = c(
    jitter_vals(63, sd = 2),  jitter_vals(81, sd = 2),
    jitter_vals(83, sd = 2),  jitter_vals(97, sd = 2),
    jitter_vals(143, sd = 3), jitter_vals(187, sd = 3),
    jitter_vals(193, sd = 3), jitter_vals(235, sd = 3)
  )
)


# ----------------------------------------------------------------------
# 3. Synthetic gene expression data
#    (Gene_1, Gene_2, Gene_3)
#    Values represent relative expression (e.g., fold-change vs. reference).
# ----------------------------------------------------------------------

gene_data <- data.frame(
  Treatment = c(rep("Ctrl", 12), rep("Trt", 12)),
  Dose      = rep(rep(dose_levels, n_reps), 2),
  Rep       = rep(rep(1:n_reps, each = length(dose_levels)), 2),

  # Gene_1: downregulated with dose, upregulated under treatment
  Gene_1 = c(
    jitter_vals(1.05, sd = 0.04), jitter_vals(0.97, sd = 0.04),
    jitter_vals(0.92, sd = 0.04), jitter_vals(0.86, sd = 0.04),
    jitter_vals(1.82, sd = 0.04), jitter_vals(1.77, sd = 0.04),
    jitter_vals(1.74, sd = 0.04), jitter_vals(1.68, sd = 0.04)
  ),

  # Gene_2: modest upregulation under treatment
  Gene_2 = c(
    jitter_vals(1.04, sd = 0.03), jitter_vals(0.97, sd = 0.03),
    jitter_vals(0.93, sd = 0.03), jitter_vals(0.89, sd = 0.03),
    jitter_vals(1.45, sd = 0.03), jitter_vals(1.38, sd = 0.03),
    jitter_vals(1.33, sd = 0.03), jitter_vals(1.27, sd = 0.03)
  ),

  # Gene_3: upregulation with treatment, stable with dose
  Gene_3 = c(
    jitter_vals(1.02, sd = 0.03), jitter_vals(0.92, sd = 0.03),
    jitter_vals(0.94, sd = 0.03), jitter_vals(0.89, sd = 0.03),
    jitter_vals(1.53, sd = 0.03), jitter_vals(1.46, sd = 0.03),
    jitter_vals(1.44, sd = 0.03), jitter_vals(1.36, sd = 0.03)
  )
)


# ----------------------------------------------------------------------
# 4. Merge datasets and create unique sample IDs
# ----------------------------------------------------------------------

physio_data <- physio_data %>%
  mutate(Sample = paste(Treatment, Dose, Rep, sep = "_")) %>%
  select(-Treatment, -Dose, -Rep)

gene_data <- gene_data %>%
  mutate(Sample = paste(Treatment, Dose, Rep, sep = "_")) %>%
  select(-Treatment, -Dose, -Rep)

all_data <- inner_join(physio_data, gene_data, by = "Sample")

cat(sprintf("Merged dataset: %d samples x %d variables\n",
            nrow(all_data), ncol(all_data) - 1))  # -1 for Sample column


# ----------------------------------------------------------------------
# 5. Prepare matrices for Mantel analysis
#    env   = physiological parameters (predictor block)
#    spec  = gene expression data (response block)
# ----------------------------------------------------------------------

env  <- all_data %>% select(Param_1, Param_2, Param_3, Param_4, Param_5)
spec <- all_data %>% select(Gene_1, Gene_2, Gene_3)

# Standardize env variables (zero mean, unit variance)
# This is important when variables are on different scales
env_scaled <- scale(env)

env_names  <- colnames(env_scaled)
spec_names <- colnames(spec)


# ----------------------------------------------------------------------
# 6. Pairwise Mantel tests (all env x gene combinations)
#    Method: Pearson correlation on Euclidean distance matrices
#    Permutations: 9999 (for reliable p-values at alpha = 0.05)
# ----------------------------------------------------------------------

set.seed(123)  # reproducibility of permutation tests

mantel_res <- data.frame()

for (env_i in env_names) {
  for (spec_j in spec_names) {

    env_vec  <- env_scaled[, env_i]
    spec_vec <- spec[, spec_j]

    dist_env  <- dist(env_vec,  method = "euclidean")
    dist_spec <- dist(spec_vec, method = "euclidean")

    m_test <- mantel(dist_env, dist_spec,
                     method       = "pearson",
                     permutations = 9999)

    mantel_res <- rbind(mantel_res, data.frame(
      from    = spec_j,            # gene (left-hand side in the plot)
      to      = env_i,             # parameter (bottom axis)
      r       = round(m_test$statistic, 4),
      p.value = round(m_test$signif, 4),
      stringsAsFactors = FALSE
    ))
  }
}

# 'rd' column is used by geom_couple() for the size aesthetic
mantel_res$rd <- mantel_res$r

cat("\n========== Mantel Test Results ==========\n")
print(mantel_res)


# ----------------------------------------------------------------------
# 7. Pearson correlation matrix among physiological parameters
# ----------------------------------------------------------------------

env_corr <- correlate(as.data.frame(env_scaled))


# ----------------------------------------------------------------------
# 8. Build the linkET visualization
# ----------------------------------------------------------------------

# Diverging colour palette for Pearson correlation heatmap
col_pal <- rev(brewer.pal(11, "RdBu"))

p <- qcorrplot(env_corr, type = "lower", diag = FALSE) +

  # Heatmap tiles for pairwise parameter correlations
  geom_square() +

  # Curved lines connecting gene expression to physiological parameters
  # Colour encodes Mantel p-value; size encodes Mantel r
  geom_couple(
    data       = mantel_res,
    aes(colour = p.value, size = rd),
    curvature  = nice_curvature()
  ) +

  # Colour scale: heatmap fill (Pearson r among parameters)
  scale_fill_gradientn(
    colours = col_pal,
    name    = "Pearson's r",
    limits  = c(-1, 1)
  ) +

  # Size scale: Mantel r
  scale_size_continuous(
    range = c(0.5, 2.5),
    name  = "Mantel's r"
  ) +

  # Colour scale: Mantel p-value (reversed so low p = darker/warm colour)
  scale_colour_gradientn(
    colours = c("#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B"),
    name    = "Mantel's p",
    trans   = "reverse",
    limits  = c(1, 0)
  ) +

  guides(
    fill   = guide_colorbar(title = "Pearson's r", order = 3),
    colour = guide_colorbar(title = "Mantel's p",  order = 1),
    size   = guide_legend(title   = "Mantel's r",  order = 2)
  ) +

  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    plot.margin      = margin(10, 10, 10, 10)
  )

print(p)


# ----------------------------------------------------------------------
# 9. Save outputs
# ----------------------------------------------------------------------

output_dir <- "outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Save figure
fig_path <- file.path(output_dir, "mantel_linkET_demo_plot.png")
ggsave(
  filename = fig_path,
  plot     = p,
  width    = 9,
  height   = 6,
  dpi      = 300,
  bg       = "white"
)
cat(sprintf("\nFigure saved: %s\n", fig_path))

# Save Mantel results table
csv_path <- file.path(output_dir, "mantel_results_demo.csv")
write.csv(mantel_res, csv_path, row.names = FALSE)
cat(sprintf("Results table saved: %s\n", csv_path))

cat("\nDone.\n")
