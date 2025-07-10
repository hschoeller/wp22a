# library(changepoint)
library(sandwich)
library(nlme)
library(strucchange)

#--- Function to detect change points ----------------------------------------
detect_change_points <- function(data,
                                 dates,
                                 max_cp_vec = 1:4) {
    # Initialize an empty data frame to store change point info:
    # cp_index: index in the time series
    # cp_date: the corresponding date
    # max_cp: the hierarchical level (iteration when first found)
    cp_df <- data.frame(
        cp_index = integer(),
        cp_date = as.Date(character()),
        max_cp = integer(),
        stringsAsFactors = FALSE
    )

    for (i in seq_along(max_cp_vec)) {
        max_cp <- max_cp_vec[i]
        cpt_obj <- cpt.meanvar(data,
            method = "BinSeg",
            Q = max_cp,
            class = TRUE
        )
        # cpts() returns indices of change points (including final index)
        cps <- cpts(cpt_obj)
        cps <- cps[cps < length(data)]

        # Save only newly found change points
        new_cps <- cps[!cps %in% cp_df$cp_index]

        if (length(new_cps) > 0) {
            for (cp in new_cps) {
                cp_df <- rbind(
                    cp_df,
                    data.frame(
                        cp_index = cp,
                        cp_date = dates[cp],
                        max_cp = max_cp,
                        stringsAsFactors = FALSE
                    )
                )
            }
            cat(
                "Max change points:", max_cp,
                "-> New cp detected at indices:", paste(new_cps, collapse = ", "),
                "which correspond to dates:",
                paste(as.character(dates[new_cps]), collapse = ", "),
                "\n"
            )
        } else {
            cat(
                "Max change points:", max_cp,
                "-> No new change points detected.\n"
            )
        }
    }
    return(cp_df)
}
