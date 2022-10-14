test_that("rvmgeom", {
    r <- 1.8
    perception.const <- function(m, log = FALSE) {
        rep(ifelse(log, 0.0, 1.0), length(m))
    }

    with_seed <- function(seed, code) {
        code <- substitute(code)
        orig.seed <- .Random.seed
        on.exit(.Random.seed <<- orig.seed)
        set.seed(seed)
        eval.parent(code)
    }

    ####################################
    # Tests without limiting magnitude #
    ####################################

    # meteor magnitudes >= 0.0
    m <- vismeteor::rvmgeom(1000, r)
    expect_type(m, 'double')
    expect_length(m, 1000)
    expect_false(anyNA(m))
    expect_false(anyNA(is.infinite(m)))
    expect_true(all(m >= 0.0))
    expect_equal(as.integer(m), m)

    #################################
    # Tests with limiting magnitude #
    #################################

    # meteor magnitudes <= 6.0
    m <- vismeteor::rvmgeom(1000, r, lm = 6.5)
    expect_type(m, 'double')
    expect_length(m, 1000)
    expect_false(anyNA(m))
    expect_false(anyNA(is.infinite(m)))
    expect_true(all(m <= 6.0))
    expect_equal(as.integer(m), m)

    ####################################
    # Tests with perception.fun        #
    ####################################

    m <- with_seed(7, vismeteor::rvmgeom(1000, r, perception.fun = perception.const))
    expect_type(m, 'double')
    expect_length(m, 1000)
    expect_false(anyNA(m))
    expect_false(anyNA(is.infinite(m)))
    expect_true(all(m >= 0.0))
    expect_equal(as.integer(m), m)
    m.geom <- with_seed(7, stats::rgeom(10000, 1.0 - 1/r))
    expect_equal(round(mean(m.geom), 1), round(mean(m), 1))


})