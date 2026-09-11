## install-deps.R 
## Installs the packages this study needs. Run once:
##
##   Rscript scripts/install-deps.R

pkgs <- c("mice", "smcfcs", "ranger", "ggplot2", "knitr", "rmarkdown")

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) == 0) {
  message("All packages already installed.")
} else {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  message(sprintf("%-10s %s", p, if (ok) as.character(utils::packageVersion(p)) else "NOT INSTALLED"))
}
