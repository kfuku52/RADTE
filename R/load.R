radte_module_files <- c(
    "R/version.R",
    "R/cli.R",
    "R/tree_navigation.R",
    "R/input_parsers.R",
    "R/calibration.R",
    "R/chronos_backend.R",
    "R/output.R",
    "R/mcmctree_backend.R",
    "R/generax.R",
    "R/preparation.R",
    "R/dating.R",
    "R/run_output.R",
    "R/main.R"
)

radte_source_modules <- function(root = ".", envir = parent.frame()) {
    for (module_file in radte_module_files) {
        sys.source(file.path(root, module_file), envir = envir)
    }
    invisible(envir)
}
