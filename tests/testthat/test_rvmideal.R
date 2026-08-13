test_that("rvmideal", {
    psi <- 4.0

    with_seed <- function(seed, code) {
        code <- substitute(code)
        orig_seed <- .Random.seed
        on.exit(.Random.seed <<- orig_seed)
        set.seed(seed)
        eval.parent(code)
    }

    # meteor magnitudes <= 6.0
    m <- vismeteor::rvmideal(1000, 6.5, psi)
    expect_type(m, "double")
    expect_length(m, 1000)
    expect_false(anyNA(m))
    expect_false(any(is.infinite(m)))
    expect_true(all(m <= 6.0))
    expect_equal(as.integer(m), m)

    # meteor magnitudes equals geometric distribution
    lm <- 6.3
    m <- with_seed(8, vismeteor::rvmideal(100000, lm, 30))
    expect_type(m, "double")
    expect_length(m, 100000)
    expect_false(anyNA(m))
    expect_false(any(is.infinite(m)))
    expect_true(all(m <= 6))
    expect_equal(as.integer(m), m)
    m_geom <- with_seed(8, rvmgeom(100000, lm, 10^(1 / 2.5)))
    expect_equal(round(mean(m_geom), 2), round(mean(m), 2))
})

test_that("rvmideal with an infinite psi", {
    lm <- 6.3

    # The geometric limit is reached exactly, so the same seed has to produce
    # the same magnitudes as the geometric model rather than merely a matching
    # mean.
    set.seed(8L)
    m <- vismeteor::rvmideal(1000, lm, Inf)
    expect_type(m, "double")
    expect_length(m, 1000)
    expect_false(anyNA(m))
    expect_false(any(is.infinite(m)))
    expect_true(all(m <= 6))
    expect_equal(as.integer(m), m)

    set.seed(8L)
    expect_equal(m, vismeteor::rvmgeom(1000, lm, 10^0.4))

    expect_error(vismeteor::rvmideal(10, lm, -Inf), "psi")
})
