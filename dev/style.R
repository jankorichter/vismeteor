# Re-style the package and run lintr.
# Convenience wrapper around the same calls used in the Makefile.
styler::style_pkg(
    indent_by = 4L,
    scope = "tokens",
    strict = FALSE
)
styler::style_dir(
    "vignettes",
    filetype = c("R", "Rmd"),
    indent_by = 4L
)
lintr::lint_package()
