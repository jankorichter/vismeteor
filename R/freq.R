#' @title Form groups of equal minimum count of frequencies
#' @description
#' This function generates quantiles with a minimum frequency.
#' These quantiles are formed from a vector `freq` of frequencies.
#' Each quantile then has the minimum total frequency `freq.min`.
#' @param freq integer; A vector of frequencies.
#' @param min integer; Minimum total frequency per quantile.
#' @details
#' The frequencies `freq` are grouped in the order in which they
#' are passed as a vector.
#' The minimum `min` must be greater than `0`.
#' @return
#' A factor of indices is returned.
#' The index references the corresponding passed frequency.
#' @examples
#' freq <- c(1,2,3,4,5,6,7,8,9)
#' cumsum(freq)
#' (f <- freq.quantile(freq, 10))
#' sapply(split(freq, f), sum)
#' @export
freq.quantile <- function(freq, min) {
    if (0 >= min) {
        stop(paste0('min must be greater than 0 instead of "', min, '"!'))
    }

    n.sum <- 0
    id.last <- 1L
    id <- integer(length(freq))
    for (i in seq_along(freq)) {
        n.i <- freq[i]
        if ((n.i + n.sum) < min) {
            id[i] <- id.last
            n.sum <- n.i + n.sum
        } else {
            id[i] <- id.last
            n.sum <- 0
            id.last <- id.last + 1L
        }
    }

    if (0 == n.sum) {
        id.last <- id.last - 1L
    }

    id.max <- id[length(id)]
    if (id.max > 1 & n.sum > 0 & sum(freq[id == id.last]) < min) {
        id[id == id.max] <- id.max - 1L
    }

    factor(id, ordered = TRUE)
}

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
