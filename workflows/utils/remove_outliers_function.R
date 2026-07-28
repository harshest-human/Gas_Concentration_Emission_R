# =============================================================================
# remove_outliers_function.R
#
# 1.5 x IQR outlier removal, optionally grouped (e.g. per location). Sets
# outliers to NA so downstream code can keep row structure intact.
#
# Arguments:
#   df              data.frame / tibble
#   exclude_cols    character vector of numeric columns to leave alone
#   group_cols      character vector of grouping columns; outliers are
#                   computed per group when supplied
#   summary_path    optional CSV path; when given, the per-column removed
#                   count is written there (in addition to the console
#                   print). The CSV has columns: Variable, Removed, plus
#                   analyzer_label and lab_label when supplied.
#   analyzer_label  optional analyser identifier (e.g. "FTIR.1"); written
#                   as an `analyzer` column in the summary CSV
#   lab_label       optional lab identifier (e.g. "ATB"); written as a `lab`
#                   column in the summary CSV
#
# Returns the input data.frame / tibble with outlier cells set to NA.
# =============================================================================
remove_outliers <- function(df,
                            exclude_cols   = NULL,
                            group_cols     = NULL,
                            summary_path   = NULL,
                            analyzer_label = NULL,
                            lab_label      = NULL) {
        # Helper function for a single numeric vector
        remove_vec_outliers <- function(x) {
                qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
                H <- 1.5 * IQR(x, na.rm = TRUE)
                x[x < (qnt[1] - H) | x > (qnt[2] + H)] <- NA
                return(x)
        }

        # Columns to process (numeric and not excluded)
        cols_to_process <- names(df)[sapply(df, is.numeric)]
        if (!is.null(exclude_cols)) {
                cols_to_process <- setdiff(cols_to_process, exclude_cols)
        }

        # Before counts
        before_counts <- colSums(!is.na(df[cols_to_process]))

        # Apply outlier removal
        if (!is.null(group_cols)) {
                df <- df %>%
                        group_by(across(all_of(group_cols))) %>%
                        mutate(across(all_of(cols_to_process), remove_vec_outliers)) %>%
                        ungroup()
        } else {
                df[cols_to_process] <- lapply(df[cols_to_process], remove_vec_outliers)
        }

        # After counts
        after_counts <- colSums(!is.na(df[cols_to_process]))

        # Build removed-count summary
        removed_counts <- before_counts - after_counts
        summary_df <- data.frame(Variable = names(removed_counts),
                                 Before   = unname(before_counts),
                                 After    = unname(after_counts),
                                 Removed  = unname(removed_counts),
                                 stringsAsFactors = FALSE)
        if (!is.null(analyzer_label)) summary_df$analyzer <- analyzer_label
        if (!is.null(lab_label))      summary_df$lab      <- lab_label

        print(summary_df)

        # Optional CSV export (Stage 1 dropout log)
        if (!is.null(summary_path)) {
                write.csv(summary_df, summary_path, row.names = FALSE)
        }

        return(df)
}
