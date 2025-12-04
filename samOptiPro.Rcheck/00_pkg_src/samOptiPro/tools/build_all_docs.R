#!/usr/bin/env Rscript
# ======================================================
# build_all_docs.R
# Génère automatiquement :
#   A. vignettes individuelles (1 par fonction)
#   B. overview global "samOptiPro-overview.Rmd"
#   C. roxygen manquant dans R/
#   D. site pkgdown complet
# ======================================================

suppressMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(fs)
})

cat("\n==============================\n")
cat("samOptiPro — Documentation Builder\n")
cat("==============================\n")


# ------------------------------------------------------
# 0 — Localisation du package
# ------------------------------------------------------
pkg <- normalizePath(".", mustWork = TRUE)
r_dir <- file.path(pkg, "R")
vig_dir <- file.path(pkg, "vignettes")
extdata <- file.path(pkg, "inst", "extdata")

if (!dir_exists(vig_dir)) dir_create(vig_dir)


# ------------------------------------------------------
# 1 — Charger l’inventaire propre
# ------------------------------------------------------
inv_file <- file.path(extdata, "functions_inventory_reordered_clean.csv")
if (!file.exists(inv_file)) stop("Inventaire manquant : ", inv_file)

df <- read_csv(inv_file, show_col_types = FALSE)

cat("[INFO] Inventaire chargé : ", nrow(df), " fonctions détectées\n")


# ======================================================
# A — Génération automatique des vignettes individuelles
# ======================================================

cat("\n[A] Génération des vignettes individuelles…\n")

df %>%
  rowwise() %>%
  do({
    fun <- .$fun_name
    desc <- .$description
    args <- .$arguments

    vignette_file <- file.path(vig_dir, paste0("fun_", fun, ".Rmd"))

    txt <- c(
      "---",
      paste0("title: \"Function: ", fun, "\""),
      "output: rmarkdown::html_vignette",
      "vignette: >",
      paste0("  %\\VignetteIndexEntry{Function: ", fun, "}"),
      "  %\\VignetteEngine{knitr::rmarkdown}",
      "  %\\VignetteEncoding{UTF-8}",
      "---",
      "",
      paste0("# `", fun, "`"),
      "",
      ifelse(desc == "", "_No description available — to be completed._", desc),
      "",
      "## Arguments",
      "```text",
      args %||% "",
      "```",
      "",
      "## Example",
      "```r",
      paste0("# ", fun, "()"),
      "```"
    )

    writeLines(txt, vignette_file)

    tibble(done = TRUE)
  })

cat("[OK] Vignettes individuelles créées.\n")


# ======================================================
# B — Création du panorama global overview
# ======================================================

cat("\n[B] Création du panorama ‘samOptiPro-overview.Rmd’…\n")

overview_file <- file.path(vig_dir, "samOptiPro-overview.Rmd")

cat("[INFO] Catégorisation automatique…\n")

get_cat <- function(file) {
  case_when(
    grepl("diagnostics", file) ~ "Diagnostics",
    grepl("plots", file)       ~ "Plots",
    grepl("runners", file)     ~ "Runners",
    grepl("utils_sampler", file) ~ "Samplers Utilities",
    grepl("check_inits", file) ~ "Inits / Robustness",
    grepl("constants", file)   ~ "Constants",
    grepl("hmc_wrappers", file)~ "HMC Wrappers",
    TRUE ~ "Misc"
  )
}

df$category <- get_cat(df$file)

body_sections <- df %>%
  group_by(category) %>%
  arrange(fun_name) %>%
  summarise(
    txt = paste0(
      "## ", first(category), "\n\n",
      paste(
        sprintf("### `%s`\n%s\n", fun_name,
                ifelse(description == "", "_No description yet._", description)),
        collapse = "\n"
      ),
      "\n"
    )
  ) %>%
  pull(txt) %>%
  paste(collapse = "\n")

overview <- c(
  "---",
  "title: \"samOptiPro — Package Overview\"",
  "output: rmarkdown::html_vignette",
  "vignette: >",
  "  %\\VignetteIndexEntry{samOptiPro — Overview}",
  "  %\\VignetteEngine{knitr::rmarkdown}",
  "  %\\VignetteEncoding{UTF-8}",
  "---",
  "",
  "# Overview",
  "",
  "This document provides a structured overview of all exported and internal functions of **samOptiPro**.",
  "",
  body_sections
)

writeLines(overview, overview_file)

cat("[OK] Overview généré.\n")


# ======================================================
# C — Auto-génération du roxygen manquant
# ======================================================

cat("\n[C] Mise à jour automatique des descriptions roxygen…\n")

auto_desc <- "#' @description FIX-ME: Auto-generated description. Please revise."

empty_desc <- df %>% filter(description == "")

for (fun in empty_desc$fun_name) {
  file <- empty_desc$df.file[empty_desc$fun_name == fun]
  path <- file.path(r_dir, file)

  if (!file.exists(path)) next

  txt <- readLines(path)

  # repérer le début de la fonction
  idx <- grep(paste0("^", fun, " *<-"), txt)

  if (length(idx) == 0) next

  # repérer si déjà un @description avant
  has_desc <- any(grepl("@description", txt[1:(idx[1]-1)], fixed = TRUE))

  if (!has_desc) {
    insert_at <- idx[1] - 1
    txt <- append(txt, auto_desc, after = insert_at)
    writeLines(txt, path)
    cat("[ADDED] @description -> ", fun, "\n")
  }
}

cat("[OK] Mise à jour des roxygen terminée.\n")


# ======================================================
# D — Construction du site pkgdown
# ======================================================

cat("\n[D] Construction du site pkgdown…\n")

pkgdown_yml <- file.path(pkg, "_pkgdown.yml")

if (!file.exists(pkgdown_yml)) {
  writeLines(c(
    "template:",
    "  bootstrap: 5",
    "",
    "reference:",
    "  - title: \"Reference\"",
    "    contents:",
    "      - \"*\""
  ), pkgdown_yml)

  cat("[INFO] _pkgdown.yml créé automatiquement.\n")
}

suppressWarnings({
  if (!requireNamespace("pkgdown", quietly = TRUE)) {
    install.packages("pkgdown")
  }
})

pkgdown::build_site(pkg = pkg, preview = FALSE)

cat("\n[FINI] Documentation complète générée.\n")
# ======================================================
# Patch pkgdown : forcer les fonctions problématiques en interne
# ======================================================

cat("\n[PATCH] Vérification des fonctions manquantes pour pkgdown…\n")

problem_funs <- c("configure_hmc_safely_bis", "test_strategy_family_fast")

for (fun in problem_funs) {

  rows <- df %>% filter(fun_name == fun)
  if (nrow(rows) == 0) next

  file <- rows$file[1]
  path <- file.path(r_dir, file)

  if (!file.exists(path)) next

  txt <- readLines(path)

  # Chercher un bloc roxygen juste avant la fonction
  idx <- grep(paste0("^", fun, " *<-"), txt)
  if (length(idx) == 0) next

  # Vérifier si @keywords internal existe déjà
  has_kw <- any(grepl("@keywords internal", txt[1:(idx[1]-1)], fixed = TRUE))

  if (!has_kw) {
    insert_at <- idx[1] - 1
    txt <- append(txt, "#' @keywords internal", after = insert_at)
    writeLines(txt, path)
    cat("  -> Ajout @keywords internal à :", fun, "\n")
  }
}

cat("[PATCH] Correction pkgdown appliquée.\n")
