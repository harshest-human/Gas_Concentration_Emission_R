######################################################################
# 03_lins_ccc.R   -- Tier A.2
# ---------------------------------------------------------------------
# Lin's concordance correlation coefficient (CCC) per analyzer pair,
# with bootstrap 95 % CI. CCC penalises both location shift (mean
# bias) and scale shift (slope != 1), so it complements Bland-Altman.
#
# Reviewer linkage: R2-L353, R2-L391  -- "use Lin's CCC".
#
# We use epiR::epi.ccc when available; otherwise compute CCC manually.
# A 1000-sample non-parametric bootstrap gives a percentile 95 % CI.
######################################################################

stopifnot(exists("concentration_reshaped"), exists("emission_reshaped"))

v7_lins_ccc <- function(x, y) {
        ok <- is.finite(x) & is.finite(y)
        x  <- x[ok]; y <- y[ok]
        n <- length(x)
        if (n < 5L) return(c(rho_c = NA_real_, n = n))
        mx <- mean(x); my <- mean(y)
        sx2 <- mean((x - mx)^2)
        sy2 <- mean((y - my)^2)
        sxy <- mean((x - mx) * (y - my))
        rho_c <- (2 * sxy) / (sx2 + sy2 + (mx - my)^2)
        c(rho_c = rho_c, n = n)
}

v7_ccc_boot_ci <- function(x, y, R = 1000L, seed = 42L) {
        ok <- is.finite(x) & is.finite(y)
        x  <- x[ok]; y <- y[ok]; n <- length(x)
        if (n < 5L) return(c(lo = NA_real_, hi = NA_real_))
        set.seed(seed)
        boot_vals <- replicate(R, {
                idx <- sample.int(n, n, replace = TRUE)
                v7_lins_ccc(x[idx], y[idx])["rho_c"]
        })
        unname(stats::quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE)) |>
                setNames(c("lo", "hi"))
}

v7_ccc_pair <- function(df_long, var_name, pair) {
        wide <- v7_to_wide(df_long, var_name)
        v <- v7_pair_vectors(wide, pair[1], pair[2])
        if (is.null(v)) return(NULL)
        ccc <- v7_lins_ccc(v$x, v$y)
        ci  <- v7_ccc_boot_ci(v$x, v$y, R = 1000L)
        # Pearson r decomposition
        r   <- stats::cor(v$x, v$y, use = "pairwise.complete.obs")
        Cb  <- ifelse(is.finite(r) && r != 0, ccc["rho_c"] / r, NA_real_)
        tibble::tibble(
                variable    = var_name,
                analyzer_x  = pair[1],
                analyzer_y  = pair[2],
                n           = v$n,
                rho_c       = ccc["rho_c"],
                ci_lo       = ci["lo"],
                ci_hi       = ci["hi"],
                pearson_r   = r,
                Cb_accuracy = Cb       # "bias correction factor"
        )
}

v7_ccc_var <- function(df_long, var_name, pairs = v7_pairs_all) {
        purrr::map_dfr(pairs, ~ v7_ccc_pair(df_long, var_name, .x))
}

v7_ccc_table <- dplyr::bind_rows(
        v7_ccc_var(concentration_reshaped, "CO2_mgm3"),
        v7_ccc_var(concentration_reshaped, "CH4_mgm3"),
        v7_ccc_var(concentration_reshaped, "NH3_mgm3"),
        v7_ccc_var(concentration_reshaped, "delta_CO2"),
        v7_ccc_var(concentration_reshaped, "delta_CH4"),
        v7_ccc_var(concentration_reshaped, "delta_NH3"),
        v7_ccc_var(emission_reshaped,      "Q_vent"),
        v7_ccc_var(emission_reshaped,      "e_CH4_ghLU"),
        v7_ccc_var(emission_reshaped,      "e_NH3_ghLU")
) %>%
        dplyr::mutate(
                # McBride 2005 strength-of-agreement bands (poor / moderate / substantial / almost perfect)
                strength = dplyr::case_when(
                        rho_c < 0.90               ~ "poor",
                        rho_c < 0.95               ~ "moderate",
                        rho_c < 0.99               ~ "substantial",
                        rho_c >= 0.99              ~ "almost perfect",
                        TRUE                       ~ NA_character_
                )
        ) %>%
        dplyr::mutate(across(where(is.numeric), ~ v7_round(.x, 4)))

v7_write_csv(v7_ccc_table, "03_lins_ccc.csv")

# ---- compact lattice plot of CCC by pair / variable ----------------
v7_ccc_plot_df <- v7_ccc_table %>%
        dplyr::mutate(pair = paste0(analyzer_x, " / ", analyzer_y))

p_ccc <- ggplot2::ggplot(v7_ccc_plot_df,
                         ggplot2::aes(x = pair, y = variable, fill = rho_c)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho_c)), size = 2.5) +
        ggplot2::scale_fill_gradientn(
                colours = c("darkred","red","orange2","gold2","yellow",
                            "greenyellow","green1","green3","darkgreen"),
                limits = c(0, 1), name = "Lin CCC"
        ) +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
                       axis.title  = ggplot2::element_blank(),
                       legend.position = "bottom")
v7_save_plot(p_ccc, "03_lins_ccc_lattice.png", width = 12, height = 6)

message("v7 module 03 (Lin's CCC): done.")
