## tools/cran_prepare_and_check.R
## À lancer depuis la racine du package samOptiPro

message("=== samOptiPro : préparation CRAN (infra + Rd) ===")

if (!file.exists("DESCRIPTION")) {
  stop("Ce script doit être lancé dans le répertoire racine du package (où se trouve DESCRIPTION).")
}

pkgs <- c("usethis", "devtools", "tools")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}
library(usethis)
library(devtools)
library(tools)

## ------------------------------------------------------------------
## 1) LICENCE MIT + file LICENSE, URL GitHub, VignetteBuilder, Suggests
## ------------------------------------------------------------------
message("[1] Licence / DESCRIPTION / vignettes")

usethis::use_mit_license("Romuald Hounyeme")

desc <- read.dcf("DESCRIPTION")

## URL -> GitHub (évite 404 pkgdown)
desc[1, "URL"] <- "https://github.com/RomualdEcoStats/samOptiPro-core"

## BugReports si manquant
if (!"BugReports" %in% colnames(desc)) {
  desc <- cbind(desc,
                BugReports = "https://github.com/RomualdEcoStats/samOptiPro-core/issues")
}

## VignetteBuilder : knitr
if (!"VignetteBuilder" %in% colnames(desc)) {
  desc <- cbind(desc, VignetteBuilder = "knitr")
} else {
  desc[1, "VignetteBuilder"] <- "knitr"
}

## knitr / rmarkdown en Suggests
suggests <- if ("Suggests" %in% colnames(desc)) desc[1, "Suggests"] else ""
sv <- unique(strsplit(suggests, "\\s*,\\s*")[[1]])
sv <- sv[sv != ""]
add_suggest <- function(x) {
  if (!x %in% sv) sv <<- c(sv, x)
}
add_suggest("knitr")
add_suggest("rmarkdown")
desc[1, "Suggests"] <- paste(sv, collapse = ", ")

write.dcf(desc, "DESCRIPTION")

## ------------------------------------------------------------------
## 2) .Rbuildignore : ignorer Tutorial_M3b.Rmd
## ------------------------------------------------------------------
message("[2] .Rbuildignore")

rb <- if (file.exists(".Rbuildignore")) readLines(".Rbuildignore") else character(0)
if (!any(grepl("^\\^Tutorial_M3b\\\\.Rmd\\$", rb))) {
  rb <- c(rb, "^Tutorial_M3b\\.Rmd$")
  writeLines(rb, ".Rbuildignore")
}

## ------------------------------------------------------------------
## 3) Dépendances imports / suggests (Matrix, posterior, stats, utils)
## ------------------------------------------------------------------
message("[3] Dépendances imports / suggests")

usethis::use_package("Matrix",   type = "Imports")
usethis::use_package("posterior", type = "Suggests")

helper_file <- "R/zzz_monitors_helpers.R"
if (!file.exists(helper_file)) {
  writeLines(c("# Helpers for samOptiPro", ""), helper_file)
}
helper_lines <- readLines(helper_file)

add_roxy_import <- function(pattern, line) {
  if (!any(grepl(pattern, helper_lines, fixed = TRUE))) {
    helper_lines <<- c(sprintf("#' %s", line), helper_lines)
  }
}
add_roxy_import("importFrom(stats, aggregate",   "importFrom(stats, aggregate)")
add_roxy_import("importFrom(stats, reorder",     "importFrom(stats, reorder)")
add_roxy_import("importFrom(utils, capture.output", "importFrom(utils, capture.output)")
add_roxy_import("importFrom(utils, sessionInfo",    "importFrom(utils, sessionInfo)")

writeLines(helper_lines, helper_file)

## ------------------------------------------------------------------
## 4) Nettoyage Rd : enlever le @export parasite dans la Rd
## ------------------------------------------------------------------
message("[4] Nettoyage Rd (ligne '@export' parasite dans les exemples)")

rd_file <- "man/plot_strategies_from_test_result_fast.Rd"
if (file.exists(rd_file)) {
  rd <- readLines(rd_file)
  keep <- !grepl("^\\s*@export\\s*$", rd)
  if (any(!keep)) {
    rd <- rd[keep]
    writeLines(rd, rd_file)
    message("  -> '@export' retiré de man/plot_strategies_from_test_result_fast.Rd")
  } else {
    message("  -> Pas de ligne '@export' nue détectée dans ce Rd.")
  }
} else {
  message("  -> Rd man/plot_strategies_from_test_result_fast.Rd introuvable (à vérifier).")
}

## ------------------------------------------------------------------
## 5) Corriger \name dans grapes-or-or-grapes.Rd (si déjà généré)
## ------------------------------------------------------------------
message("[5] Corriger \\name dans man/grapes-or-or-grapes.Rd")

grapes_rd <- "man/grapes-or-or-grapes.Rd"
if (file.exists(grapes_rd)) {
  rd <- readLines(grapes_rd)
  idx_name <- grep("^\\\\name\\{", rd)
  if (length(idx_name) == 1L) {
    rd[idx_name] <- "\\name{grapes-or-or-grapes}"
    writeLines(rd, grapes_rd)
    message("  -> \\name{grapes-or-or-grapes} appliqué.")
  } else {
    message("  -> Pas de ligne \\name unique détectée, à inspecter à la main.")
  }
} else {
  message("  -> man/grapes-or-or-grapes.Rd introuvable (roxygen ne l'a peut-être pas reconstruit).")
}

## ------------------------------------------------------------------
## 6) Scan non-ASCII (info seulement)
## ------------------------------------------------------------------
message("[6] Scan non-ASCII (info, à corriger manuellement si tu veux zéro WARNING)")

r_files <- c("R/diagnostics.R", "R/plots.R", "R/zzz_monitors_helpers.R")
for (f in r_files) {
  if (file.exists(f)) {
    cat("\n--- Non-ASCII scan:", f, "---\n")
    print(tools::showNonASCIIfile(f))
  }
}

message(">>> Pour CRAN ultra-strict : remplacer les caractères accentués/non-ASCII dans le CODE par des \\uXXXX ou ASCII.")

## ------------------------------------------------------------------
## 7) Build du tar.gz et R CMD check --as-cran
## ------------------------------------------------------------------
message("[7] Build + R CMD check --as-cran (sans re-document)")

# Construction du tarball
tarball      <- devtools::build()
tarball_path <- normalizePath(tarball)

message("  -> Tarball construit : ", tarball_path)

# Lancer R CMD check directement sur le chemin absolu
cmd <- paste("R CMD check", shQuote(tarball_path), "--as-cran")
message("  -> Lancement : ", cmd)
system(cmd)
