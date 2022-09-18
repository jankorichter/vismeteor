test_that("rvmgeom", {
    r <- 1.8

    ####################################
    # Tests without limiting magnitude #
    ####################################

    # meteor magnitudes >= 0.0
    m <- rvmgeom(1000, r)
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
    m <- rvmgeom(1000, r, lm = 6.5)
    expect_type(m, 'double')
    expect_length(m, 1000)
    expect_false(anyNA(m))
    expect_false(anyNA(is.infinite(m)))
    expect_true(all(m <= 6.0))
    expect_equal(as.integer(m), m)
})