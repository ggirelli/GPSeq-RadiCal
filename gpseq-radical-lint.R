#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# License: MIT - Copyright (c) 2020 Gabriele Girelli
# ------------------------------------------------------------------------------

depends_on = c("lintr", "logging")
tmp = lapply(depends_on, function(x) {
    expr = x %in% names(installed.packages()[, "Package"])
    if (!expr) {
        cat(sprintf(" Missing '%s' package.\n", x),
            sprintf("Install it with 'install.packages(\"%s\")'.\n", x))
    }
    stopifnot(expr)
})
logging::basicConfig()
logging::loginfo("Linting...")
lintr::lint("gpseq-radical.R", linters=lintr::with_defaults(
    assignment_linter=NULL,
    infix_spaces_linter=NULL,
    lintr::line_length_linter(80),
    object_usage_linter=NULL))
logging::loginfo("Done.")

################################################################################
