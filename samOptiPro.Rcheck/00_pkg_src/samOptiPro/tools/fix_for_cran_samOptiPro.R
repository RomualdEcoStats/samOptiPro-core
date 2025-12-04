## fix_for_cran_samOptiPro.R
## À lancer depuis la racine du package samOptiPro

message("=== Fix CRAN infrastructure issues for samOptiPro ===")

if (!file.exists("DESCRIPTION")) {
  stop("Ce script doit être lancé dans le répertoire racine du package (où se trouve DESCRIPTION).")
}

## --- 0. Charge les outils nécessaires (sans planter si absent) ----
pkg_needed <- c("usethis", "devtools", "tools")
for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(usethis)
library(devtools)
library(tools)

## --- 1. LICENSE : MIT + file LICENSE correctement configuré --------
message("[1] Fix licence MIT + file LICENSE")

# Crée / écrase proprement LICENSE + LICENSE.md et synchronise DESCRIPTION
usethis::use_mit_license("Romuald Hounyeme")

# S'assure que DESCRIPTION pointe bien vers 'MIT + file LICENSE'
desc <- read.dcf("DESCRIPTION")
desc[1, "License"] <- "MIT + file LICENSE"
write.dcf(desc, "DESCRIPTION")

## --- 2. DESCRIPTION : URL, VignetteBuilder, Suggests ---------------
message("[2] Ajuste DESCRIPTION (URL, VignetteBuilder, Suggests)")

desc <- read.dcf("DESCRIPTION")

## 2.1 URL : enlever le site pkgdown 404 et garder l’URL GitHub
url_new <- "https://github.com/RomualdEcoStats/samOptiPro-core"
desc[1, "URL"] <- url_new

## 2.2 BugReports si pas encore mis
if (!"BugReports" %in% colnames(desc)) {
  desc <- cbind(desc, BugReports = paste0(url_new, "/issues"))
}

## 2.3 VignetteBuilder: knitr
if (!"VignetteBuilder" %in% colnames(desc)) {
  desc <- cbind(desc, VignetteBuilder = "knitr")
} else {
  desc[1, "VignetteBuilder"] <- "knitr"
}

## 2.4 Ajoute knitr / rmarkdown dans Suggests
suggests <- if ("Suggests" %in% colnames(desc)) desc[1, "Suggests"] else ""
suggests_vec <- unique(strsplit(suggests, "\\s*,\\s*")[[1]])
suggests_vec <- suggests_vec[suggests_vec != ""]

add_suggest <- function(x) {
  if (!x %in% suggests_vec) suggests_vec <<- c(suggests_vec, x)
}

add_suggest("knitr")
add_suggest("rmarkdown")

desc[1, "Suggests"] <- paste(suggests_vec, collapse = ", ")

write.dcf(desc, "DESCRIPTION")

## --- 3. Rbuildignore : ignorer Tutorial_M3b.Rmd au top-level -------
message("[3] Ajoute Tutorial_M3b.Rmd à .Rbuildignore (top-level)")

rb <- if (file.exists(".Rbuildignore")) readLines(".Rbuildignore") else character(0)
if (!"^Tutorial_M3b\\.Rmd$" %in% rb) {
  rb <- c(rb, "^Tutorial_M3b\\.Rmd$")
  writeLines(rb, ".Rbuildignore")
}

## --- 4. Dépendances : Matrix, posterior, stats, utils --------------
message("[4] Déclare les dépendances utilisées dans le code")

## 4.1 Matrix utilisé par loadNamespace/requireNamespace
usethis::use_package("Matrix", type = "Imports")

## 4.2 posterior : tu peux décider Imports vs Suggests. Ici Suggests (diags avancés).
usethis::use_package("posterior", type = "Suggests")

## 4.3 stats::aggregate, stats::reorder, utils::capture.output, utils::sessionInfo
## Ajoute les importFrom via roxygen, en les plaçant dans zzz_monitors_helpers.R (ou un fichier central)
helper_file <- "R/zzz_monitors_helpers.R"
if (file.exists(helper_file)) {
  helper_lines <- readLines(helper_file)
} else {
  helper_lines <- c("# Helpers for samOptiPro", "")
}

# On ajoute des directives roxygen d'import si non présentes
add_roxy_import <- function(pattern, line) {
  if (!any(grepl(pattern, helper_lines, fixed = TRUE))) {
    helper_lines <<- c(sprintf("#' %s", line), helper_lines)
  }
}

add_roxy_import("importFrom(stats,aggregate", "importFrom(stats, aggregate)")
add_roxy_import("importFrom(stats,reorder",   "importFrom(stats, reorder)")
add_roxy_import("importFrom(utils,capture.output", "importFrom(utils, capture.output)")
add_roxy_import("importFrom(utils,sessionInfo",    "importFrom(utils, sessionInfo)")

writeLines(helper_lines, helper_file)

## --- 5. @export parasite dans la Rd d’exemple ----------------------
message("[5] Corrige l'@export parasite dans man/plot_strategies_from_test_result_fast.Rd")

rd_file <- "man/plot_strategies_from_test_result_fast.Rd"
if (file.exists(rd_file)) {
  rd <- readLines(rd_file)
  # supprime toute ligne qui ne contient qu'un @export (avec ou sans espaces)
  keep <- !grepl("^\\s*@export\\s*$", rd)
  if (any(!keep)) {
    rd <- rd[keep]
    writeLines(rd, rd_file)
  }
} else {
  warning("Fichier Rd man/plot_strategies_from_test_result_fast.Rd introuvable – à vérifier.")
}

## --- 6. \name invalide dans grapes-or-or-grapes.Rd -----------------
message("[6] Corrige le \\name dans man/grapes-or-or-grapes.Rd")

grapes_rd <- "man/grapes-or-or-grapes.Rd"
if (file.exists(grapes_rd)) {
  rd <- readLines(grapes_rd)
  # Remplace la ligne \name{%||%} (ou similaire) par un nom sans |, avec alias qui doit déjà exister
  idx_name <- grep("^\\\\name\\{", rd)
  if (length(idx_name) == 1L) {
    rd[idx_name] <- "\\name{grapes-or-or-grapes}"
    writeLines(rd, grapes_rd)
  }
} else {
  warning("Fichier man/grapes-or-or-grapes.Rd introuvable – à vérifier.")
}

## --- 7. Vignettes : s’assurer que knitr est bien pris en compte -----
message("[7] Vérifie/initialise la configuration de vignettes (knitr/rmarkdown)")

if (dir.exists("vignettes")) {
  # Ajoute l’infrastructure pkgdown/knitr si besoin ; ne réécrit pas tes Rmd
  try(usethis::use_vignette("samOptiPro_overview_dummy"), silent = TRUE)
  # On supprime aussitôt le dummy si créé
  dummy <- "vignettes/samOptiPro_overview_dummy.Rmd"
  if (file.exists(dummy)) file.remove(dummy)
}

## --- 8. Non-ASCII dans le code : pointage pour correction manuelle -
message("[8] Scan des fichiers R pour caractères non-ASCII (à corriger manuellement)")

r_files_to_scan <- c(
  "R/diagnostics.R",
  "R/plots.R",
  "R/zzz_monitors_helpers.R"
)
for (f in r_files_to_scan) {
  if (file.exists(f)) {
    cat("\n--- Non-ASCII scan:", f, "---\n")
    print(tools::showNonASCIIfile(f))
  }
}

message(">>> Remplace les caractères non-ASCII par des \\uXXXX ou des équivalents ASCII (au moins dans le code, CRAN tolère en commentaires).")

## --- 9. Rappel sur la doc / Rd / export ----------------------------
message("
[9] Reste à faire à la main (pour CRAN clean) :

  - Documenter ou rendre internes (non exportés) les fonctions :
      enrich_hmc_diag_tbl_for_plots, identify_bottlenecks, sampler_df
    => soit tu ajoutes des blocs roxygen + \\link, soit tu retires @export.

  - Corriger les Roxygen/Rd pour les arguments signalés comme non documentés :
      build_fn, model, mcmc_conf, ignore_patterns, strict_sampler_only, etc.
    => regarde les WARNING 'Rd \\usage sections' dans devtools::check().

  - Corriger les vignettes si besoin (contenu), maintenant que VignetteBuilder est en place.

")

## --- 10. Regénère la doc & relance un check ------------------------
message("[10] Regénère la doc et relance un check rapide")

devtools::document()
devtools::check(cran = TRUE)
