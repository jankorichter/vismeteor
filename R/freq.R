#' @title Convert a contingency table of meteor magnitudes into integer values
#' @description
#' The meteor magnitude contingency table of VMDB contains half meteor counts (e.g. `3.5`).
#' This function converts these values to integer values.
#' @param mt table; A two-dimensional contingency table of meteor magnitudes.
#' @param expand logical; If `TRUE` expand meteor magnitudes.
#' @details
#' The contingency table of meteor magnitudes `mt` must be two-dimensional.
#' Column names must be integer meteor magnitude values.
#' The row names refer to the magnitude observations.
#' @return
#' if `expand` is `FALSE`, contingency table of meteor magnitudes is returned.
#'
#' Otherwise, a list of individual meteor magnitudes is returned.
#' The name of elements refer to the magnitude observations.
#' @examples
#' # For example, create a contingency table of meteor magnitudes
#' data <- data.frame(
#'     obs.id = rep(c("A", "B"), each=4),
#'     m = c(2, 3, 4, 5, -1, 3, 4, 5),
#'     freq = c(2.5, 0.5, 0, 3, 1, 3, 2.5, 1.5)
#' )
#' (mt <- xtabs(freq ~ obs.id + m, data=data))
#'
#' # contingency table with integer values
#' vmtable(mt)
#'
#' # list of individual meteor magnitudes
#' vmtable(mt, expand = TRUE )
#' @export
vmtable <- function(mt, expand = FALSE) {
    if (! methods::is(mt, 'table')) {
        stop(paste0('Magnitude table is not a table!'))
    }

    if (2 != length(dim(mt))) {
        stop(paste0('Magnitude table is not two-dimensional!'))
    }

    if ((is.integer(as.vector(mt)))) {
        mt.new <- mt
    } else {
        mt.list <- apply(mt, 1, function(mt) {
            list(
                m = as.numeric(names(mt)),
                freq = as.vector(mt)
            )
        }, simplify = FALSE)

        data <- do.call(
            rbind.data.frame,
            mapply(function(id, mt) {
                m <- as.numeric(mt$m)
                freq <- mt$freq

                # random order of magnitudes (asc. or desc.?)
                set.seed(19 + sum(m) + sum(freq))
                decr <- stats::rbinom(1, 1, c(0.5, 0.5))
                o <- order(m, decreasing = as.logical(decr))

                m <- m[o]
                freq.dup <- 2 * freq[o]
                i <- seq_len(length(freq.dup))
                i <- rep(i, freq.dup)
                i <- i[seq_len(length(i)) %% 2 == 0]
                mt <- table(m[i])
                freq <- as.vector(mt)
                list(
                    id = rep(id, length(freq)),
                    m = names(mt),
                    freq = as.integer(freq)
                )
            }, rownames(mt), mt.list, SIMPLIFY = FALSE)
        )

        mt.new <- stats::xtabs(freq ~ id + m, data=data)
    }

    if (expand) {
        data <- data.frame(
            id = rep(rownames(mt.new), each=ncol(mt.new)),
            m = rep(as.numeric(colnames(mt.new)), times=nrow(mt.new)),
            freq = as.vector(t(mt.new[,1:ncol(mt.new)]))
        )
        data <- subset(data, data$freq != 0)

        mt.list <- do.call(
            list,
            lapply(split(data, data$id), function(d){
                as.integer(rep(d$m, d$freq))
            })
        )

        return(mt.list)
    }

    dimnames(mt.new) <- dimnames(mt)

    return(mt.new)
}
