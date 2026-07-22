## Regression check for sdreport_native parity on supported models
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
fallback_example <- Sys.getenv("fallback_example", unset = "")
joint_example <- Sys.getenv("joint_example", unset = "")
bias_example <- Sys.getenv("bias_example", unset = "")

# Build/load the local native kernel when not provided by installed TMB.
if (!is.loaded("tmb_sdreport_delta_fixed") || !is.loaded("tmb_sdreport_delta_random")) {
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
    tryCatch(
        {
            expr <- sprintf(
                "library(TMB); source('%s/TMB/R/sdreport.R'); updateCholesky <- get('updateCholesky', envir=asNamespace('TMB')); solveSubset <- get('solveSubset', envir=asNamespace('TMB')); runSymbolicAnalysis <- get('runSymbolicAnalysis', envir=asNamespace('TMB')); runExample('%s', thisR=TRUE, exfolder='%s'); cmp <- compare_sdreport_native(obj, compare.cov=TRUE); print(cmp); if(!isTRUE(cmp$pass$all)) stop('max_abs differences exceeded tolerance')",
                repo_root, ex, file.path(repo_root, "tmb_examples")
            )
            status <- system2("Rscript", c("-e", shQuote(expr)), stdout = "", stderr = "")
            entry$pass <- identical(status, 0L)
            if (!entry$pass) {
                entry$error <- sprintf("subprocess exited with status %d", status)
            }
        },
        error = function(e) {
            entry$pass <<- FALSE
            entry$error <<- conditionMessage(e)
        }
    )
    results[[i]] <- entry
}

if (nzchar(fallback_example)) {
    ex <- fallback_example
    key <- paste0("fallback:", ex)
    cat("\n--- Fallback Example:", ex, "---\n")
    entry <- list(pass = FALSE, error = NULL, cmp = NULL)
    tryCatch(
        {
            expr <- sprintf(
                paste(
                    "library(TMB)",
                    "source('%s/TMB/R/sdreport.R')",
                    "updateCholesky <- get('updateCholesky', envir=asNamespace('TMB'))",
                    "solveSubset <- get('solveSubset', envir=asNamespace('TMB'))",
                    "runSymbolicAnalysis <- get('runSymbolicAnalysis', envir=asNamespace('TMB'))",
                    "runExample('%s', thisR=TRUE, exfolder='%s')",
                    "ok <- FALSE",
                    "tryCatch({ sdreport_native(obj, strict=TRUE, bias.correct=TRUE); stop('strict=TRUE unexpectedly succeeded') }, error=function(e) { ok <<- TRUE })",
                    "if (!ok) stop('strict fallback check failed')",
                    "ref <- suppressWarnings(sdreport(obj, bias.correct=TRUE))",
                    "nat <- suppressWarnings(sdreport_native(obj, strict=FALSE, bias.correct=TRUE))",
                    "m <- function(x,y){ d <- abs(as.numeric(x)-as.numeric(y)); if(length(d)==0 || all(is.na(d))) return(NA_real_); max(d, na.rm=TRUE) }",
                    "dv <- m(ref$value, nat$value)",
                    "ds <- m(ref$sd, nat$sd)",
                    "if (!is.na(dv) && dv > 1e-7) stop(sprintf('fallback value mismatch: %%g', dv))",
                    "if (!is.na(ds) && ds > 1e-7) stop(sprintf('fallback sd mismatch: %%g', ds))",
                    "cat('Fallback check passed for %s\\n')",
                    sep = "; "
                ),
                repo_root, ex, file.path(repo_root, "tmb_examples"), ex
            )
            status <- system2("Rscript", c("-e", shQuote(expr)), stdout = "", stderr = "")
            entry$pass <- identical(status, 0L)
            if (!entry$pass) {
                entry$error <- sprintf("subprocess exited with status %d", status)
            }
        },
        error = function(e) {
            entry$pass <<- FALSE
            entry$error <<- conditionMessage(e)
        }
    )
    results[[key]] <- entry
}

if (nzchar(bias_example)) {
    ex <- bias_example
    key <- paste0("bias:", ex)
    cat("\n--- Bias Correction Example:", ex, "---\n")
    entry <- list(pass = FALSE, error = NULL, cmp = NULL)
    tryCatch({
        expr <- sprintf(
            paste(
                "library(TMB)",
                "source('%s/TMB/R/sdreport.R')",
                "updateCholesky <- get('updateCholesky', envir=asNamespace('TMB'))",
                "solveSubset <- get('solveSubset', envir=asNamespace('TMB'))",
                "runSymbolicAnalysis <- get('runSymbolicAnalysis', envir=asNamespace('TMB'))",
                "runExample('%s', thisR=TRUE, exfolder='%s')",
                "ref <- suppressWarnings(TMB::sdreport(obj, bias.correct=TRUE))",
                "nat <- suppressWarnings(sdreport_native(obj, strict=TRUE, bias.correct=TRUE))",
                "if (is.null(ref$unbiased) || is.null(nat$unbiased)) stop('missing unbiased payload')",
                "m <- function(x,y){ d <- abs(as.numeric(x)-as.numeric(y)); if(length(d)==0 || all(is.na(d))) return(NA_real_); max(d, na.rm=TRUE) }",
                "du <- m(ref$unbiased$value, nat$unbiased$value)",
                "if (!is.na(du) && du > 1e-7) stop(sprintf('bias-corrected value mismatch: %%g', du))",
                "if (!is.null(ref$unbiased$sd) && !is.null(nat$unbiased$sd)) { ds <- m(ref$unbiased$sd, nat$unbiased$sd); if (!is.na(ds) && ds > 1e-7) stop(sprintf('bias-corrected sd mismatch: %%g', ds)) }",
                "cat('Bias correction check passed for %s\\n')",
                sep = "; "
            ),
            repo_root, ex, file.path(repo_root, "tmb_examples"), ex
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
    results[[key]] <- entry
}

if (nzchar(joint_example)) {
    ex <- joint_example
    key <- paste0("joint:", ex)
    cat("\n--- Joint Precision Example:", ex, "---\n")
    entry <- list(pass = FALSE, error = NULL, cmp = NULL)
    tryCatch(
        {
            expr <- sprintf(
                paste(
                    "library(TMB)",
                    "source('%s/TMB/R/sdreport.R')",
                    "updateCholesky <- get('updateCholesky', envir=asNamespace('TMB'))",
                    "solveSubset <- get('solveSubset', envir=asNamespace('TMB'))",
                    "runExample('%s', thisR=TRUE, exfolder='%s')",
                    "ref <- TMB::sdreport(obj, getJointPrecision=TRUE)",
                    "nat <- sdreport_native(obj, strict=TRUE, getJointPrecision=TRUE)",
                    "if (is.null(ref$jointPrecision) || is.null(nat$jointPrecision)) stop('jointPrecision missing')",
                    "if (!(inherits(ref$jointPrecision, 'Matrix') || is.matrix(ref$jointPrecision))) stop('reference jointPrecision type unsupported')",
                    "if (!(inherits(nat$jointPrecision, 'Matrix') || is.matrix(nat$jointPrecision))) stop('native jointPrecision type unsupported')",
                    "if (!identical(dim(ref$jointPrecision), dim(nat$jointPrecision))) stop('jointPrecision dimension mismatch')",
                    "dj <- max(abs(as.numeric(ref$jointPrecision) - as.numeric(nat$jointPrecision)), na.rm=TRUE)",
                    "if (!is.finite(dj) || dj > 1e-7) stop(sprintf('jointPrecision mismatch: %%g', dj))",
                    "cat('Joint precision check passed for %s\\n')",
                    sep = "; "
                ),
                repo_root, ex, file.path(repo_root, "tmb_examples"), ex
            )
            status <- system2("Rscript", c("-e", shQuote(expr)), stdout = "", stderr = "")
            entry$pass <- identical(status, 0L)
            if (!entry$pass) {
                entry$error <- sprintf("subprocess exited with status %d", status)
            }
        },
        error = function(e) {
            entry$pass <<- FALSE
            entry$error <<- conditionMessage(e)
        }
    )
    results[[key]] <- entry
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
