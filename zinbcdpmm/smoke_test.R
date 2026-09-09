# Smoke test for zinbcdpmm MVP
# Log written to smoke_test_log.txt next to this script's package root.

pkg <- "d:/Users/lenovo/Desktop/有序多分类-ZINB模型随机效应为CDPMM/zinbcdpmm"
logf <- file.path(pkg, "smoke_test_log.txt")
zz <- file(logf, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(zz)
}, add = TRUE)

cat("=== zinbcdpmm smoke test ===\n")
cat(R.version.string, "\n")
cat("Time:", as.character(Sys.time()), "\n\n")

need <- c("BayesLogit", "mvtnorm", "MCMCpack", "truncnorm", "loo", "coda")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("Installing", p, "...\n")
    install.packages(p, repos = "https://cloud.r-project.org")
  } else {
    cat("OK:", p, "\n")
  }
}

# Prefer load_all / document if available; else source install
has_devtools <- requireNamespace("devtools", quietly = TRUE)
has_pkgload <- requireNamespace("pkgload", quietly = TRUE)

ok <- FALSE
err <- NULL
tryCatch({
  if (has_devtools) {
    cat("\nTrying devtools::document()...\n")
    try(devtools::document(pkg), silent = TRUE)
  }
  if (has_pkgload) {
    cat("Loading with pkgload::load_all()...\n")
    pkgload::load_all(pkg, export_all = FALSE, quiet = FALSE)
  } else if (has_devtools) {
    cat("Loading with devtools::load_all()...\n")
    devtools::load_all(pkg, export_all = FALSE, quiet = FALSE)
  } else {
    cat("Installing from source (repos=NULL)...\n")
    install.packages(pkg, repos = NULL, type = "source")
    library(zinbcdpmm)
  }

  cat("\n--- simulate ---\n")
  set.seed(1)
  dat <- simulate_zinb_ordinal(n = 20, nis = 4, scenario = 2, re_dist = "normal")
  cat("n=", dat$n, " N=", dat$N, " class=", paste(class(dat), collapse = ","), "\n", sep = "")

  cat("\n--- fit (chain=50, burn=20, thin=5, G=4) ---\n")
  fit <- fit_zinb_ordinal(dat, chain = 50, burn = 20, thin = 5, G = 4)
  cat("\nprint(fit):\n")
  print(fit)
  cat("\nsummary(fit):\n")
  print(summary(fit))
  ok <- TRUE
}, error = function(e) {
  err <<- conditionMessage(e)
  cat("\nERROR:", err, "\n")
  if (!is.null(e$call)) cat("Call:", paste(deparse(e$call), collapse = " "), "\n")
})

cat("\n=== RESULT:", if (ok) "PASSED" else "FAILED", "===\n")
if (!ok && !is.null(err)) cat("Error message:", err, "\n")
