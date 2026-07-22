## Regression check for sdreport_native parity on a non-random model
library(TMB)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)[1]
if (is.na(file_arg) || !nzchar(file_arg)) {
    file_arg <- grep("check_sdreport_native\\.R$", args, value = TRUE)[1]
}
this_file <- if (!is.na(file_arg) && nzchar(file_arg)) sub("^--file=", "", file_arg) else ""
repo_root <- if (nzchar(this_file)) {
    normalizePath(file.path(dirname(this_file), "..", ".."), winslash = "/", mustWork = TRUE)
} else if (file.exists(file.path(getwd(), "TMB"))) {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
}

# Source local helpers so this script works even when installed TMB is older.
source(file.path(repo_root, "TMB", "R", "sdreport.R"), local = .GlobalEnv)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(file.path(repo_root, "tmb_examples"))

examples_env <- Sys.getenv("examples", unset = "")
if (nzchar(examples_env)) {
    examples <- trimws(strsplit(examples_env, ",", fixed = TRUE)[[1]])
    examples <- examples[nzchar(examples)]
} else {
    examples <- trimws(strsplit(Sys.getenv("example", unset = "linreg,linreg_parallel,longlinreg"), ",", fixed = TRUE)[[1]])
    examples <- examples[nzchar(examples)]
}
if (length(examples) == 0) {
    stop("No examples specified.")
}
cat("Running sdreport_native regression for examples:", paste(examples, collapse = ", "), "\n")

# Build/load the local native kernel when not provided by installed TMB.
if (!is.loaded("tmb_sdreport_delta_fixed")) {
    src <- file.path(repo_root, "TMB", "src", "sdreport_native.cpp")
    include_tmb <- file.path(repo_root, "TMB", "inst", "include")
    include_eigen <- system.file("include", package = "RcppEigen")
    cmd <- sprintf(
        "PKG_CPPFLAGS='-I%s -I%s' R CMD SHLIB %s",
        shQuote(include_tmb), shQuote(include_eigen), shQuote(src)
    )
    status <- system(cmd)
    if (status != 0) {
        stop("Failed to compile local sdreport_native.cpp")
    }
    dynlib <- sub("\\.cpp$", .Platform$dynlib.ext, src)
    dyn.load(dynlib)
}

if (!exists("compare_sdreport_native", mode = "function", inherits = TRUE)) {
    stop("compare_sdreport_native not found after sourcing local sdreport.R")
}

results <- vector("list", length(examples))
names(results) <- examples

for (i in seq_along(examples)) {
    ex <- examples[[i]]
    cat("\n--- Example:", ex, "---\n")
    entry <- list(pass = FALSE, error = NULL, cmp = NULL)
    tryCatch({
        expr <- sprintf(
            "library(TMB); source('%s/TMB/R/sdreport.R'); runExample('%s', thisR=TRUE, exfolder='%s'); cmp <- compare_sdreport_native(obj, compare.cov=TRUE); print(cmp); if(!isTRUE(cmp$pass$all)) stop('max_abs differences exceeded tolerance')",
            repo_root, ex, file.path(repo_root, "tmb_examples")
        )
        status <- system2("Rscript", c("-e", shQuote(expr)), stdout = "", stderr = "")
        entry$pass <- identical(status, 0L)
        if (!entry$pass) {
            entry$error <- sprintf("subprocess exited with status %d", status)
        }
    }, error = function(e) {
        entry$pass <<- FALSE
        entry$error <<- conditionMessage(e)
    })
    results[[i]] <- entry
}

cat("\nSummary:\n")
for (i in seq_along(results)) {
    ex <- names(results)[[i]]
    status <- if (isTRUE(results[[i]]$pass)) "PASS" else "FAIL"
    cat(sprintf("  [%s] %s\n", status, ex))
    if (!isTRUE(results[[i]]$pass) && !is.null(results[[i]]$error)) {
        cat("    ", results[[i]]$error, "\n", sep = "")
    }
}

failed <- names(results)[!vapply(results, function(x) isTRUE(x$pass), logical(1))]
if (length(failed) > 0) {
    stop("sdreport_native regression failed for: ", paste(failed, collapse = ", "))
}

cat("sdreport_native regression passed for all examples.\n")
