# ==============================================================================
#  f-OFDM Spectral Shaping Filter Analysis — Figures 11–14 & Tables VII–X
# ==============================================================================
#
#  Authors:  Pramoth Kumar M, Rajasingh J, Abhishek J, Sundararaj K
#  Paper:    "Spectral Shaping Filter Architectures for Filtered-OFDM in 5G
#             and Beyond: A Comparative Analysis of NGF, Sinc², TASM, and RRC"
#  Journal:  ETRI Journal (submitted)
#
#  Repository: https://github.com/j-abhishek26/filtered-OFDM
#  License:    MIT
#
#  Description:
#    Generates Figures 11–14 and supplementary figures for the above paper.
#    Evaluates NGF, Sinc², TASM, and RRC filters for f-OFDM across BER vs.
#    Eb/N₀ performance, power spectral density (PSD) and out-of-band emission
#    (OOBE), PAPR/CCDF analysis, and computational complexity modelling.
#
#  Usage:
#    source("fofdm_analysis_FIXED.R")
#
#  Outputs:
#    ./figures/Fig11_*.pdf/png  ...  Fig14_*.pdf/png
#    ./figures/FigS1_Filter_Responses.pdf/png
#    ./figures/FigSummary_All_Metrics.pdf/png
#    ./results/Table_VII_BER_EbN0.csv  ...  Table_X_Complexity_*.csv
#    ./results/BER_all_filters.csv, PSD_all_filters.csv, PAPR_CCDF_all_filters.csv
#
#  Requirements:
#    R >= 4.3, packages: ggplot2, dplyr, tidyr, patchwork, scales, pracma,
#    latex2exp (auto-installed if missing)
# ==============================================================================

# ── 0.  Package loading  ──────────────────────────────────────────────────────
# Load pracma first so dplyr::filter() takes precedence
if (!requireNamespace("pracma",   quietly = TRUE)) install.packages("pracma")
if (!requireNamespace("ggplot2",  quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",    quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr",    quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("patchwork",quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("scales",   quietly = TRUE)) install.packages("scales")
if (!requireNamespace("latex2exp",quietly = TRUE)) install.packages("latex2exp")

suppressPackageStartupMessages({
  library(pracma)     # erfc(), linspace()  — load before dplyr
  library(ggplot2)
  library(dplyr)      # dplyr::filter() now masks pracma; that is fine
  library(tidyr)
  library(patchwork)
  library(scales)
  library(latex2exp)
})

# Namespace aliases to avoid dplyr/signal conflicts
df_filter   <- dplyr::filter          # always use this for data-frame filtering
df_select   <- dplyr::select
df_mutate   <- dplyr::mutate
df_group_by <- dplyr::group_by
df_summarise<- dplyr::summarise

# ── Output directories ────────────────────────────────────────────────────────
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ── IEEE-style ggplot2 theme ──────────────────────────────────────────────────
ieee_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey88", linewidth = 0.4),
    axis.title        = element_text(size = 10, face = "bold"),
    axis.text         = element_text(size = 9),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 9),
    legend.key.width  = unit(1.4, "cm"),
    plot.title        = element_text(size = 11, face = "bold",  hjust = 0.5),
    plot.subtitle     = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
    strip.text        = element_text(size = 9,  face = "bold"),
    strip.background  = element_rect(fill = "grey94", colour = "grey70")
  )

# ── Colour / linetype / shape palette ────────────────────────────────────────
FILTER_COLS <- c(
  "NGF"           = "#1b7837",
  "RRC"           = "#2166ac",
  "TASM"          = "#d6604d",
  "Sinc2"         = "#8073ac",
  "Ideal AWGN"    = "black",
  "Unfiltered OFDM" = "grey50"
)
FILTER_LINES <- c(
  "NGF"           = "solid",
  "RRC"           = "dashed",
  "TASM"          = "dotdash",
  "Sinc2"         = "dotted",
  "Ideal AWGN"    = "longdash",
  "Unfiltered OFDM" = "twodash"
)
FILTER_SHAPES <- c(
  "NGF"           = 16,
  "RRC"           = 17,
  "TASM"          = 15,
  "Sinc2"         = 18,
  "Ideal AWGN"    = 1,
  "Unfiltered OFDM" = 4
)

# =============================================================================
#  SECTION 1 — SIMULATION PARAMETERS
# =============================================================================

SIM <- list(
  N          = 512L,        # OFDM subcarriers
  k_sub      = 64L,         # sub-stream length for BER sim
  L          = 129L,        # filter length (taps) — must be odd
  fs         = 30.72e6,     # sampling frequency Hz (3GPP NR FR1)
  delta_f    = 15e3,        # subcarrier spacing Hz
  alpha_rrc  = 0.22,        # RRC roll-off factor (3GPP NR)
  EbN0_dB    = seq(0, 30, by = 1),
  QAM_BER    = c(16L, 64L, 256L),
  N_sym_ber  = 800L,        # Monte-Carlo blocks per SNR point (↑ for accuracy)
  N_sym_psd  = 80L,         # blocks for PSD estimation
  N_blk_papr = 3000L,       # blocks for PAPR CCDF (↑ for smoother tails)
  gamma_dB   = seq(4, 14, by = 0.25)
)

cat(sprintf("[INFO] Simulation parameters initialized (N=%d, L=%d, alpha_RRC=%.2f).\n",
            SIM$N, SIM$L, SIM$alpha_rrc))

# =============================================================================
#  SECTION 2 — FILTER COEFFICIENT GENERATION
# =============================================================================

# 2a. Normalized Gaussian Filter (NGF)
make_ngf <- function(L) {
  t <- seq(-(L - 1) / 2, (L - 1) / 2)
  h <- exp(-(pi * t / L)^2)
  h / sum(h)
}

# 2b. Sinc^2 / Bartlett window
make_sinc2 <- function(L) {
  t <- seq(-(L - 1) / 2, (L - 1) / 2) / ((L - 1) / 2)
  h <- pmax(1 - abs(t), 0)
  h / sum(h)
}

# 2c. Taylor-Approximated Sinc Main-Lobe (TASM, N_terms = 7)
make_tasm <- function(L, N_terms = 7L) {
  t <- seq(-(L - 1) / 2, (L - 1) / 2) / ((L - 1) / 2)
  h <- Reduce("+", lapply(seq(0L, N_terms), function(n) {
    ((-1)^n) * (pi * t)^(2 * n) / factorial(2 * n + 1)
  }))
  h <- pmax(h, 0)      # retain main lobe only (discard negative sidelobes)
  h / sum(h)
}

# 2d. Root-Raised-Cosine (RRC), Proakis & Salehi closed-form
make_rrc <- function(L, alpha = 0.22, Tsym = 1.0) {
  t <- seq(-(L - 1) / 2, (L - 1) / 2)
  h <- numeric(L)
  for (i in seq_along(t)) {
    ti <- t[i] / Tsym
    if (abs(ti) < 1e-9) {
      # t = 0
      h[i] <- (1 / Tsym) * (1 - alpha + 4 * alpha / pi)
    } else if (abs(abs(4 * alpha * ti) - 1) < 1e-6) {
      # t = +/- T/(4*alpha)
      h[i] <- (alpha / (Tsym * sqrt(2))) *
        ((1 + 2 / pi) * sin(pi / (4 * alpha)) +
           (1 - 2 / pi) * cos(pi / (4 * alpha)))
    } else {
      num  <- sin(pi * ti * (1 - alpha)) +
        4 * alpha * ti * cos(pi * ti * (1 + alpha))
      den  <- pi * ti * (1 - (4 * alpha * ti)^2)
      h[i] <- (1 / Tsym) * num / den
    }
  }
  h / sum(abs(h))
}

# 2e. Identity (unfiltered OFDM — delta at centre)
make_identity <- function(L) {
  h <- rep(0.0, L)
  h[ceiling(L / 2)] <- 1.0
  h
}

# ── Build filter bank ─────────────────────────────────────────────────────────
L <- SIM$L
filters <- list(
  NGF      = make_ngf(L),
  RRC      = make_rrc(L, alpha = SIM$alpha_rrc),
  TASM     = make_tasm(L),
  Sinc2    = make_sinc2(L),
  Identity = make_identity(L)
)

cat("[INFO] Filter bank generated:\n")
for (nm in names(filters))
  cat(sprintf("  %-10s  sum(h)=%.6f  energy=%.6f\n",
              nm, sum(filters[[nm]]), sum(filters[[nm]]^2)))

# =============================================================================
#  SECTION 3 — QAM MODULATION HELPERS
# =============================================================================

# Square M-QAM constellation on normalised [-1,+1] grid
qam_constellation <- function(M) {
  sM <- sqrt(M)
  stopifnot(sM == round(sM))
  lv  <- seq(-(sM - 1), sM - 1, by = 2) / (sM - 1)
  as.matrix(expand.grid(I = lv, Q = lv))   # M x 2
}

# 0-based index → complex IQ
qam_map <- function(idx, M) {
  C <- qam_constellation(M)
  complex(real = C[idx + 1, 1], imaginary = C[idx + 1, 2])
}

# Hard decision: complex IQ → 0-based index
qam_demap <- function(rx, M) {
  C  <- qam_constellation(M)
  Cc <- complex(real = C[, 1], imaginary = C[, 2])
  # Vectorised: for each received sample find nearest constellation point
  vapply(rx, function(s) which.min(Mod(s - Cc)) - 1L, integer(1))
}

# Average symbol energy of normalised constellation
qam_sym_energy <- function(M) {
  C <- qam_constellation(M)
  mean(C[, 1]^2 + C[, 2]^2)
}

# =============================================================================
#  SECTION 4 — BER vs Eb/N0
# =============================================================================

# 4a. Theoretical BER — Gray-coded square M-QAM over AWGN
ber_theory <- function(M, EbN0_dB) {
  EbN0   <- 10^(EbN0_dB / 10)
  k      <- log2(M)
  sM     <- sqrt(M)
  qf     <- function(x) 0.5 * erfc(x / sqrt(2))   # Q-function via erfc
  ber    <- (4 / k) * (1 - 1 / sM) * qf(sqrt(3 * k * EbN0 / (M - 1)))
  pmax(ber, 1e-10)
}

# 4b. Simulated BER — frequency-domain filtering model
#     TX: QAM → IFFT → freq-filter
#     RX: freq-matched-filter → IFFT → hard-decision
ber_simulated <- function(filt_h, M, EbN0_dB_vec, N_sym = 500L) {
  set.seed(42L)
  k_bits   <- log2(M)
  N_fft    <- SIM$k_sub                       # sub-stream FFT size
  Es       <- qam_sym_energy(M)               # normalised symbol energy

  # Frequency-domain filter response
  # Truncate or zero-pad filter to exactly N_fft points, then FFT
  L_h   <- length(filt_h)
  if (L_h >= N_fft) {
    h_pad <- filt_h[seq_len(N_fft)]           # truncate (rare)
  } else {
    h_pad <- c(filt_h, rep(0.0, N_fft - L_h)) # zero-pad
  }
  H     <- Mod(fft(h_pad))
  H     <- H / max(H)                         # unit-peak normalisation
  H        <- H / max(H)                      # unit-peak normalisation

  ber_vec  <- numeric(length(EbN0_dB_vec))

  for (j in seq_along(EbN0_dB_vec)) {
    EbN0   <- 10^(EbN0_dB_vec[j] / 10)
    EsN0   <- EbN0 * k_bits
    # Per-real-dimension noise std-dev
    sigma  <- sqrt(Es / (2 * EsN0))

    n_err  <- 0L
    n_tot  <- 0L

    for (blk in seq_len(N_sym)) {
      # --- Transmitter ---
      tx_idx <- sample.int(M, N_fft, replace = TRUE) - 1L
      tx_iq  <- qam_map(tx_idx, M)
      TX_f   <- fft(tx_iq)                    # frequency domain
      TX_f   <- TX_f * H                      # apply filter

      # --- Channel: AWGN ---
      noise  <- complex(
        real      = rnorm(N_fft, 0, sigma),
        imaginary = rnorm(N_fft, 0, sigma)
      )
      RX_f   <- TX_f + noise

      # --- Receiver: matched filter + IFFT ---
      RX_f   <- RX_f * H                      # matched filter (|H|^2)
      rx_iq  <- fft(RX_f, inverse = TRUE) / N_fft

      # --- Hard decision ---
      rx_idx <- qam_demap(rx_iq, M)
      n_err  <- n_err + sum(tx_idx != rx_idx)
      n_tot  <- n_tot + N_fft
    }

    ber_vec[j] <- max(n_err / (n_tot * k_bits), 1e-10)
  }
  ber_vec
}

# 4c. Run BER for all filters × 3 QAM orders
cat("\n[1/4] BER vs. Eb/N0 simulation...\n")
EbN0_vec    <- SIM$EbN0_dB
ber_results <- list()

for (M in SIM$QAM_BER) {
  cat(sprintf("       %d-QAM\n", M))
  # Theoretical bound
  ber_results[[paste0("Ideal AWGN_", M)]] <- ber_theory(M, EbN0_vec)
  # Simulated filters
  for (fn in c("NGF", "RRC", "TASM", "Sinc2")) {
    cat(sprintf("       %d-QAM | %s\n", M, fn))
    ber_results[[paste0(fn, "_", M)]] <-
      ber_simulated(filters[[fn]], M, EbN0_vec, N_sym = SIM$N_sym_ber)
  }
}

# 4d. Build tidy data frame
ber_df <- do.call(rbind, lapply(names(ber_results), function(key) {
  # Key format:  "FilterName_M"  where FilterName may contain underscores
  # Strategy: last token after split on "_" is M; rest is filter name
  parts  <- strsplit(key, "_")[[1]]
  M_val  <- as.integer(parts[length(parts)])
  fname  <- paste(parts[-length(parts)], collapse = " ")   # rejoin with space
  data.frame(Filter    = fname,
             M         = M_val,
             EbN0      = EbN0_vec,
             BER       = ber_results[[key]],
             stringsAsFactors = FALSE)
}))

ber_df$Filter    <- factor(ber_df$Filter,
                           levels = c("NGF","RRC","TASM","Sinc2","Ideal AWGN"))
ber_df$QAM_label <- factor(paste0(ber_df$M, "-QAM"),
                           levels = c("16-QAM","64-QAM","256-QAM"))

# 4e. Table VII — Eb/N0 required at BER = 1e-3
#     Compute for each filter×M independently
compute_ber_table <- function(df, target_ber = 1e-3) {
  all_filters <- c("NGF","RRC","TASM","Sinc2","Ideal AWGN")
  all_M       <- sort(unique(df$M))

  result <- do.call(rbind, lapply(all_M, function(m) {
    row <- data.frame(M = m)
    for (fn in all_filters) {
      # Extract rows for this filter and M directly — no pipe filter() needed
      sub  <- df[df$Filter == fn & df$M == m, ]
      idx  <- which(sub$BER <= target_ber)
      val  <- if (length(idx) == 0) NA_real_ else round(min(sub$EbN0[idx]), 1)
      row[[fn]] <- val
    }
    row
  }))
  result
}

ber_table_VII <- compute_ber_table(ber_df)
cat("       Table VII exported to ./results/Table_VII_BER_EbN0.csv\n")
write.csv(ber_table_VII, "results/Table_VII_BER_EbN0.csv", row.names = FALSE)

# 4f. Figure 11 — BER vs Eb/N0
#     Use df_filter() alias to avoid signal::filter masking
ber_points <- df_filter(ber_df, EbN0 %% 5 == 0, Filter != "Ideal AWGN")

fig11 <- ggplot(ber_df,
                aes(x = EbN0, y = BER,
                    colour = Filter, linetype = Filter, shape = Filter)) +
  geom_line(linewidth = 0.70) +
  geom_point(data = ber_points, size = 2.0) +
  scale_y_log10(
    breaks = 10^seq(-6, 0),
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(1e-6, 1.2)
  ) +
  scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 30)) +
  scale_colour_manual(values  = FILTER_COLS)  +
  scale_linetype_manual(values = FILTER_LINES) +
  scale_shape_manual(values   = FILTER_SHAPES) +
  facet_wrap(~ QAM_label, ncol = 3) +
  labs(
    title    = "BER vs. Eb/N0 — f-OFDM Filter Comparison",
    subtitle = "NGF | RRC (\u03b1=0.22) | TASM | Sinc\u00b2 | Ideal AWGN bound",
    x        = TeX("$E_b/N_0$ (dB)"),
    y        = "Bit Error Rate (BER)"
  ) +
  ieee_theme +
  guides(colour   = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1),
         shape    = guide_legend(nrow = 1))

ggsave("figures/Fig11_BER_vs_EbN0.pdf", fig11,
       width = 7.2, height = 3.8, device = cairo_pdf)
ggsave("figures/Fig11_BER_vs_EbN0.png", fig11,
       width = 7.2, height = 3.8, dpi = 300)
cat("       Saved.\n")

# =============================================================================
#  SECTION 5 — POWER SPECTRAL DENSITY (PSD) & OOBE
# =============================================================================

# 5a. Generate time-domain f-OFDM signal
gen_fofdm_signal <- function(filt_h, N = 512L, M = 64L,
                             n_sym = 50L, seed = 123L) {
  set.seed(seed)
  out <- vector("complex", n_sym * (N + length(filt_h) - 1L))
  pos <- 1L
  for (s in seq_len(n_sym)) {
    idx    <- sample.int(M, N, replace = TRUE) - 1L
    tx_iq  <- qam_map(idx, M)
    td     <- fft(tx_iq, inverse = TRUE) / N    # IFFT → time domain
    # Linear convolution (type="open") for correct transient handling
    sig_re <- convolve(Re(td), rev(filt_h), type = "open")
    sig_im <- convolve(Im(td), rev(filt_h), type = "open")
    seg    <- complex(real = sig_re, imaginary = sig_im)
    len_s  <- length(seg)
    out[pos:(pos + len_s - 1L)] <- seg
    pos    <- pos + len_s
  }
  out[seq_len(pos - 1L)]
}

# 5b. Welch PSD (Hann window, 50% overlap, zero-padded FFT)
welch_psd <- function(x, nfft = 4096L, win_len = 512L) {
  hann <- 0.5 * (1 - cos(2 * pi * seq(0, win_len - 1) / (win_len - 1)))
  hop  <- win_len %/% 2L          # 50% overlap
  n    <- length(x)
  starts <- seq(1L, n - win_len + 1L, by = hop)
  psd_acc <- numeric(nfft)
  cnt  <- 0L
  for (s in starts) {
    seg      <- as.complex(x[s:(s + win_len - 1L)]) * hann
    seg_pad  <- c(seg, rep(0+0i, nfft - win_len))
    psd_acc  <- psd_acc + Mod(fft(seg_pad))^2
    cnt      <- cnt + 1L
  }
  psd_acc <- psd_acc / (cnt * sum(hann^2))
  # One-sided (0 to fs/2)
  idx_os  <- seq_len(nfft %/% 2L + 1L)
  psd_os  <- psd_acc[idx_os] * 2.0
  psd_os[1]             <- psd_os[1] / 2.0   # DC: not doubled
  psd_os[length(psd_os)]<- psd_os[length(psd_os)] / 2.0   # Nyquist
  freq <- (idx_os - 1L) / nfft               # normalised 0..0.5
  list(freq = freq, psd = psd_os)
}

# 5c. Compute PSD for all filters
cat("\n[2/4] Power spectral density (PSD) and OOBE...\n")
NFFT_PSD <- 4096L
psd_list <- list()

for (fn in c("NGF","RRC","TASM","Sinc2","Identity")) {
  cat(sprintf("       %s\n", fn))
  sig   <- gen_fofdm_signal(filters[[fn]], N = SIM$N,
                            n_sym = SIM$N_sym_psd)
  wp    <- welch_psd(sig, nfft = NFFT_PSD, win_len = 512L)
  label <- if (fn == "Identity") "Unfiltered OFDM" else fn
  psd_list[[label]] <- data.frame(
    Filter = label,
    freq   = wp$freq,
    psd_lin= wp$psd,
    stringsAsFactors = FALSE
  )
}

psd_df <- do.call(rbind, psd_list)
# Normalise each filter's PSD to its own maximum (0 dBr reference)
psd_df <- psd_df %>%
  df_group_by(Filter) %>%
  df_mutate(psd_dB = 10 * log10(psd_lin / max(psd_lin))) %>%
  ungroup()

psd_df$Filter <- factor(psd_df$Filter,
  levels = c("NGF","RRC","TASM","Sinc2","Unfiltered OFDM"))

# 5d. OOBE table (Table VIII) — computed from linear PSD
cat("       Computing OOBE...\n")
oobe_rows <- lapply(c("NGF","RRC","TASM","Sinc2","Unfiltered OFDM"),
  function(fn) {
    sub <- psd_df[psd_df$Filter == fn, ]
    f   <- sub$freq
    p   <- sub$psd_lin
    P_in  <- sum(p[f <= 0.250])
    P_1st <- sum(p[f > 0.250 & f <= 0.375])
    P_2nd <- sum(p[f > 0.375 & f <= 0.500])
    # Guard against zero power
    o1 <- if (P_in > 0 & P_1st > 0) round(10*log10(P_1st/P_in),1) else NA_real_
    o2 <- if (P_in > 0 & P_2nd > 0) round(10*log10(P_2nd/P_in),1) else NA_real_
    data.frame(Filter = fn, OOBE_1st_dBc = o1, OOBE_2nd_dBc = o2,
               stringsAsFactors = FALSE)
  })
oobe_table <- do.call(rbind, oobe_rows)
cat("       Table VIII exported to ./results/Table_VIII_OOBE.csv\n")
write.csv(oobe_table, "results/Table_VIII_OOBE.csv", row.names = FALSE)

# 5e. Figure 12 — PSD
psd_plot <- df_filter(psd_df, freq <= 0.50)

fig12 <- ggplot(psd_plot,
                aes(x = freq, y = psd_dB,
                    colour = Filter, linetype = Filter)) +
  geom_line(linewidth = 0.72) +
  annotate("rect",
    xmin = 0, xmax = 0.25, ymin = -95, ymax = 5,
    fill = "#4393c3", alpha = 0.06) +
  geom_vline(xintercept = 0.25, linetype = "dashed",
             colour = "grey55", linewidth = 0.45) +
  annotate("text", x = 0.12,  y = -4,
           label = "In-band",  size = 2.8, colour = "#2c7fb8") +
  annotate("text", x = 0.375, y = -4,
           label = "1st adj.", size = 2.8, colour = "grey40") +
  scale_y_continuous(limits = c(-95, 5), breaks = seq(-90, 0, 10)) +
  scale_x_continuous(
    breaks = seq(0, 0.5, 0.1),
    labels = c("0","0.1","0.2","0.3","0.4","0.5")
  ) +
  scale_colour_manual(values  = FILTER_COLS)  +
  scale_linetype_manual(values = FILTER_LINES) +
  labs(
    title    = "Power Spectral Density — f-OFDM Filter Comparison",
    subtitle = "Welch PSD, Hann window, 4096-pt FFT, 64-QAM",
    x        = "Normalised Frequency (f / fs)",
    y        = "PSD (dBr)"
  ) +
  ieee_theme

ggsave("figures/Fig12_PSD_Comparison.pdf", fig12,
       width = 6.5, height = 3.8, device = cairo_pdf)
ggsave("figures/Fig12_PSD_Comparison.png", fig12,
       width = 6.5, height = 3.8, dpi = 300)
cat("       Saved.\n")

# =============================================================================
#  SECTION 6 — PAPR / CCDF
# =============================================================================

# 6a. Compute PAPR for many random OFDM blocks
compute_papr_vec <- function(filt_h, M = 256L, N = 512L,
                             n_blocks = 2000L, seed = 7L) {
  set.seed(seed)
  papr_vec <- numeric(n_blocks)
  for (b in seq_len(n_blocks)) {
    idx    <- sample.int(M, N, replace = TRUE) - 1L
    tx_iq  <- qam_map(idx, M)
    td     <- fft(tx_iq, inverse = TRUE) / N
    sig_re <- convolve(Re(td), rev(filt_h), type = "open")
    sig_im <- convolve(Im(td), rev(filt_h), type = "open")
    pwr    <- sig_re^2 + sig_im^2
    papr_vec[b] <- max(pwr) / mean(pwr)
  }
  papr_vec
}

# 6b. CCDF = Pr(PAPR > gamma_0)
compute_ccdf <- function(papr_vec, gamma_dB_seq) {
  glin <- 10^(gamma_dB_seq / 10)
  vapply(glin, function(g) mean(papr_vec > g), numeric(1))
}

cat("\n[3/4] PAPR / CCDF analysis...\n")
gamma_dB <- SIM$gamma_dB
papr_data <- list()

for (fn in c("NGF","RRC","TASM","Sinc2","Identity")) {
  cat(sprintf("       %s\n", fn))
  pv    <- compute_papr_vec(filters[[fn]], M = 256L,
                            N = SIM$N, n_blocks = SIM$N_blk_papr)
  label <- if (fn == "Identity") "Unfiltered OFDM" else fn
  ccdf  <- compute_ccdf(pv, gamma_dB)
  papr_data[[label]] <- data.frame(
    Filter   = label,
    gamma_dB = gamma_dB,
    CCDF     = pmax(ccdf, 1e-5),
    stringsAsFactors = FALSE
  )
}

papr_df <- do.call(rbind, papr_data)
papr_df$Filter <- factor(papr_df$Filter,
  levels = c("NGF","RRC","TASM","Sinc2","Unfiltered OFDM"))

# 6c. Table IX — PAPR at CCDF = 1e-3
papr_table_IX <- do.call(rbind, lapply(
  levels(papr_df$Filter), function(fn) {
    sub <- papr_df[papr_df$Filter == fn, ]
    idx <- which(sub$CCDF <= 1e-3)
    val <- if (length(idx) == 0) NA_real_ else round(min(sub$gamma_dB[idx]), 1)
    data.frame(Filter = fn, PAPR_at_CCDF_1e3_dB = val,
               stringsAsFactors = FALSE)
  }))

cat("       Table IX exported to ./results/Table_IX_PAPR.csv\n")
write.csv(papr_table_IX, "results/Table_IX_PAPR.csv", row.names = FALSE)

# 6d. Figure 13 — PAPR CCDF
fig13 <- ggplot(papr_df,
                aes(x = gamma_dB, y = CCDF,
                    colour = Filter, linetype = Filter)) +
  geom_line(linewidth = 0.72) +
  geom_hline(yintercept = 1e-3, linetype = "dotted",
             colour = "grey45", linewidth = 0.5) +
  annotate("text", x = 13.6, y = 1.6e-3,
           label = "CCDF = 10\u207b\u00b3", size = 2.8, colour = "grey35") +
  scale_y_log10(
    breaks = 10^seq(-5, 0),
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(1e-5, 1.2)
  ) +
  scale_x_continuous(breaks = seq(4, 14, 2), limits = c(4, 14)) +
  scale_colour_manual(values  = FILTER_COLS)  +
  scale_linetype_manual(values = FILTER_LINES) +
  labs(
    title    = "PAPR CCDF — f-OFDM Filter Comparison",
    subtitle = TeX("256-QAM, $N=512$ subcarriers"),
    x        = TeX("PAPR threshold $\\gamma_0$ (dB)"),
    y        = TeX("$\\Pr(\\mathrm{PAPR} > \\gamma_0)$")
  ) +
  ieee_theme

ggsave("figures/Fig13_PAPR_CCDF.pdf", fig13,
       width = 5.8, height = 3.8, device = cairo_pdf)
ggsave("figures/Fig13_PAPR_CCDF.png", fig13,
       width = 5.8, height = 3.8, dpi = 300)
cat("       Saved.\n")

# =============================================================================
#  SECTION 7 — COMPUTATIONAL COMPLEXITY
# =============================================================================

# Real-multiply model: N_x = 8*N*log2(N) + 4*N + delta_RRC * L
build_complexity_df <- function(N_vec, L_tap = 129L) {
  rows <- list()
  for (fn in c("Unfiltered OFDM","NGF","Sinc2","TASM","RRC")) {
    for (N in N_vec) {
      fft_m  <- 8 * N * log2(N)
      filt_m <- if (fn == "Unfiltered OFDM") 0 else 4 * N
      rrc_m  <- if (fn == "RRC") L_tap else 0
      total  <- fft_m + filt_m + rrc_m
      rows[[length(rows) + 1L]] <- data.frame(
        Filter       = fn,
        N            = N,
        FFT_mults    = fft_m,
        Filter_mults = filt_m,
        RRC_extra    = rrc_m,
        Total        = total,
        Normalised   = round(total / fft_m, 4),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

N_seq <- c(64, 128, 256, 512, 1024, 2048)
cx_df <- build_complexity_df(N_seq, L_tap = SIM$L)
cx_df$Filter <- factor(cx_df$Filter,
  levels = c("NGF","RRC","TASM","Sinc2","Unfiltered OFDM"))

# Table X at N = 512
cx_512 <- df_filter(cx_df, N == 512)
cat("\n[4/4] Computational complexity analysis...\n")
cat("       Table X exported to ./results/Table_X_Complexity_*.csv\n")
write.csv(cx_512, "results/Table_X_Complexity_N512.csv", row.names = FALSE)
write.csv(cx_df,  "results/Table_X_Complexity_All_N.csv", row.names = FALSE)

# Figure 14 — Complexity vs N
fig14 <- ggplot(cx_df,
                aes(x = N, y = Total / 1e3,
                    colour = Filter, linetype = Filter, shape = Filter)) +
  geom_line(linewidth = 0.72) +
  geom_point(size = 2.2) +
  scale_x_log10(breaks = N_seq, labels = as.character(N_seq)) +
  scale_y_continuous(
    labels = label_comma(),
    breaks = pretty(cx_df$Total / 1e3, n = 6)
  ) +
  scale_colour_manual(values  = FILTER_COLS)  +
  scale_linetype_manual(values = FILTER_LINES) +
  scale_shape_manual(values   = FILTER_SHAPES) +
  labs(
    title    = "Computational Complexity vs. Number of Subcarriers",
    subtitle = TeX("Real multiplications per symbol  (L = 129,  \\times 10^3)"),
    x        = "Number of Subcarriers N  (log scale)",
    y        = TeX("Real multiplications ($\\times 10^3$)")
  ) +
  ieee_theme

ggsave("figures/Fig14_Complexity_vs_N.pdf", fig14,
       width = 5.8, height = 3.8, device = cairo_pdf)
ggsave("figures/Fig14_Complexity_vs_N.png", fig14,
       width = 5.8, height = 3.8, dpi = 300)
cat("       Saved.\n")

# =============================================================================
#  SECTION 8 — FILTER TIME/FREQUENCY RESPONSES  (Supplementary Figure S1)
# =============================================================================

L_fr   <- 129L
t_axis <- seq(-(L_fr - 1) / 2, (L_fr - 1) / 2)
NFFT_FR<- 2048L

td_df <- do.call(rbind, lapply(c("NGF","RRC","TASM","Sinc2"), function(fn) {
  h <- switch(fn,
    NGF   = make_ngf(L_fr),
    RRC   = make_rrc(L_fr, SIM$alpha_rrc),
    TASM  = make_tasm(L_fr),
    Sinc2 = make_sinc2(L_fr))
  data.frame(Filter = fn, t = t_axis, h = h / max(abs(h)),
             stringsAsFactors = FALSE)
}))

fd_df <- do.call(rbind, lapply(c("NGF","RRC","TASM","Sinc2"), function(fn) {
  h    <- switch(fn,
    NGF   = make_ngf(L_fr),
    RRC   = make_rrc(L_fr, SIM$alpha_rrc),
    TASM  = make_tasm(L_fr),
    Sinc2 = make_sinc2(L_fr))
  H    <- Mod(fft(c(h, rep(0, NFFT_FR - L_fr))))
  H    <- H / max(H)
  # FFT-shift to centre DC
  half <- NFFT_FR %/% 2L
  H_s  <- c(H[(half + 1):NFFT_FR], H[1:half])
  freq <- seq(-0.5, 0.5 - 1/NFFT_FR, by = 1/NFFT_FR)
  data.frame(Filter = fn, freq = freq,
             H_dB   = 20 * log10(H_s + 1e-12),
             stringsAsFactors = FALSE)
}))

td_df$Filter <- factor(td_df$Filter, levels = c("NGF","RRC","TASM","Sinc2"))
fd_df$Filter <- factor(fd_df$Filter, levels = c("NGF","RRC","TASM","Sinc2"))

fig_s1a <- ggplot(td_df, aes(x = t, y = h, colour = Filter, linetype = Filter)) +
  geom_line(linewidth = 0.72) +
  scale_colour_manual(values  = FILTER_COLS[c("NGF","RRC","TASM","Sinc2")]) +
  scale_linetype_manual(values = FILTER_LINES[c("NGF","RRC","TASM","Sinc2")]) +
  labs(title = "Time-Domain Impulse Responses",
       x = "Sample index n", y = "Normalised amplitude h[n]") +
  ieee_theme

fig_s1b <- ggplot(df_filter(fd_df, abs(freq) <= 0.5),
                  aes(x = freq, y = H_dB, colour = Filter, linetype = Filter)) +
  geom_line(linewidth = 0.72) +
  scale_y_continuous(limits = c(-80, 3), breaks = seq(-80, 0, 20)) +
  scale_colour_manual(values  = FILTER_COLS[c("NGF","RRC","TASM","Sinc2")]) +
  scale_linetype_manual(values = FILTER_LINES[c("NGF","RRC","TASM","Sinc2")]) +
  labs(title = "Frequency-Domain Magnitude Responses",
       x = "Normalised frequency (f/fs)", y = "|H(f)| (dB)") +
  ieee_theme

fig_s1 <- fig_s1a / fig_s1b +
  plot_annotation(
    title    = "Filter Responses: NGF | RRC | TASM | Sinc\u00b2",
    subtitle = "L = 129 taps  |  RRC roll-off \u03b1 = 0.22",
    theme    = theme(
      plot.title   = element_text(size = 11, face = "bold",  hjust = 0.5),
      plot.subtitle= element_text(size =  9, hjust = 0.5, colour = "grey40"))
  )

ggsave("figures/FigS1_Filter_Responses.pdf", fig_s1,
       width = 6.5, height = 6.2, device = cairo_pdf)
ggsave("figures/FigS1_Filter_Responses.png", fig_s1,
       width = 6.5, height = 6.2, dpi = 300)
cat("       Supplementary figure S1 saved.\n")

# =============================================================================
#  SECTION 9 — COMBINED 4-PANEL SUMMARY FIGURE
# =============================================================================

fig_summary <-
  (fig11 + theme(legend.position = "none",
                 plot.title = element_text(size = 9),
                 plot.subtitle = element_blank())) /
  ((fig12 + theme(legend.position = "none",
                  plot.title = element_text(size = 9),
                  plot.subtitle = element_blank())) |
   (fig13 + theme(legend.position = "none",
                  plot.title = element_text(size = 9),
                  plot.subtitle = element_blank())) |
   (fig14 + theme(legend.position = "none",
                  plot.title = element_text(size = 9),
                  plot.subtitle = element_blank()))) +
  plot_layout(heights = c(1.5, 1.0)) +
  plot_annotation(
    title    = "f-OFDM Filter Comparison: BER | PSD | PAPR | Complexity",
    subtitle = "NGF  \u25cf  RRC (\u03b1=0.22)  \u25cf  TASM  \u25cf  Sinc\u00b2  \u25cf  Unfiltered OFDM",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey40"))
  )

ggsave("figures/FigSummary_All_Metrics.pdf", fig_summary,
       width = 11.0, height = 8.5, device = cairo_pdf)
ggsave("figures/FigSummary_All_Metrics.png", fig_summary,
       width = 11.0, height = 8.5, dpi = 300)
cat("       Summary figure saved.\n")

# =============================================================================
#  SECTION 10 — EXPORT ALL NUMERICAL DATA
# =============================================================================

write.csv(ber_df,   "results/BER_all_filters.csv",        row.names = FALSE)
write.csv(psd_df,   "results/PSD_all_filters.csv",        row.names = FALSE)
write.csv(papr_df,  "results/PAPR_CCDF_all_filters.csv",  row.names = FALSE)

# =============================================================================
#  DONE
# =============================================================================

cat("\n", strrep("=", 65), "\n")
cat("  SIMULATION COMPLETE \u2014 Figures 11\u201314 & Supplementary\n")
cat(strrep("=", 65), "\n\n")
cat("Output directory: ./figures/\n")
cat("  Fig11  BER vs. Eb/N0 (16 / 64 / 256-QAM)\n")
cat("  Fig12  Power spectral density comparison\n")
cat("  Fig13  PAPR CCDF (256-QAM)\n")
cat("  Fig14  Computational complexity vs. N\n")
cat("  FigS1  Filter time/frequency responses (supplementary)\n")
cat("  Summary  Combined 4-panel overview\n\n")
cat("Results directory: ./results/\n")
cat("  Table_VII_BER_EbN0.csv\n")
cat("  Table_VIII_OOBE.csv\n")
cat("  Table_IX_PAPR.csv\n")
cat("  Table_X_Complexity_N512.csv\n")
cat("  Table_X_Complexity_All_N.csv\n")
cat("  BER_all_filters.csv\n")
cat("  PSD_all_filters.csv\n")
cat("  PAPR_CCDF_all_filters.csv\n\n")
