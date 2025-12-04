## =============================================================
## INVENTAIRE DES FONCTIONS DU PACKAGE
## samOptiPro (ou autre package actif)
##
## Objectif :
##  - Lister toutes les fonctions du package
##  - Extraire :
##       * leur nom
##       * leurs arguments (signature)
##       * la première ligne de roxygen (description courte)
##  - Exporter sous forme de data.frame clair (+ CSV optionnel)
## =============================================================

library(stringr)

## -------- PARAMÈTRE : nom du package --------------------------
pkg <- "samOptiPro"

cat("[INFO] Scanning package:", pkg, "\n")

## -------- Trouver le chemin du dossier R/ du package ----------
pkg_path <- devtools::as.package(".")$path
rdir <- file.path(pkg_path, "R")

if (!dir.exists(rdir)) stop("Cannot find R/ directory in package.")

## -------- Lister les fichiers R -------------------------------
files <- list.files(rdir, pattern = "\\.R$", full.names = TRUE)
cat("[INFO] Found", length(files), "R files.\n")

results <- list()

## -------- Fonction utilitaire : extraire arguments ------------
extract_args <- function(fx) {
  if (!is.function(fx)) return(NA_character_)
  args <- names(formals(fx))
  paste(args, collapse = ", ")
}

## -------- Environnement du package ----------------------------
ns <- asNamespace(pkg)

## -------- Extraction pour chaque fichier ----------------------
for (f in files) {
  txt <- readLines(f, warn = FALSE)

  ## indices des lignes roxygen (#')
  roxy_idx <- grep("^#'", txt)

  ## indices des définitions de fonctions : "name <- function("
  functions_idx <- grep("^\\s*[^#].*<-\\s*function\\s*\\(", txt)

  if (length(functions_idx) == 0) next

  for (idx in functions_idx) {
    line <- txt[idx]

    ## nom de la fonction
    name <- str_match(line, "^\\s*([A-Za-z0-9_\\.]+)\\s*<-\\s*function")[, 2]
    if (is.na(name) || !nzchar(name)) next

    ## récupérer le bloc roxygen immédiatement au-dessus
    roxy_block <- character(0)
    if (length(roxy_idx)) {
      # on remonte depuis idx-1 tant que la ligne commence par #'
      start <- idx - 1L
      while (start %in% roxy_idx && start >= 1L) {
        start <- start - 1L
      }
      if (start < (idx - 1L)) {
        roxy_block <- txt[(start + 1L):(idx - 1L)]
        roxy_block <- gsub("^#' ?", "", roxy_block)
      }
    }

    ## description courte = première ligne non vide du bloc roxygen
    desc <- ""
    if (length(roxy_block)) {
      nonempty <- roxy_block[roxy_block != ""]
      if (length(nonempty)) desc <- nonempty[1]
    }

    ## extraire la fonction dans le namespace du package (si elle existe)
    fx <- tryCatch(get(name, envir = ns), error = function(e) NULL)
    args <- extract_args(fx)

    results[[length(results) + 1L]] <- list(
      file        = basename(f),
      fun_name    = name,
      arguments   = args,
      description = desc
    )
  }
}

## -------- Résultats en data.frame -----------------------------
if (length(results)) {
  df <- do.call(
    rbind,
    lapply(results, as.data.frame, stringsAsFactors = FALSE)
  )
} else {
  df <- data.frame(
    file        = character(0),
    fun_name    = character(0),
    arguments   = character(0),
    description = character(0),
    stringsAsFactors = FALSE
  )
}

cat("\n[INFO] Found", nrow(df), "functions in package.\n\n")
print(df)

## -------- Sauvegarde pour revue / vignette --------------------
out_dir  <- file.path(pkg_path, "inst", "extdata")
out_file <- file.path(out_dir, "functions_inventory.csv")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(df, out_file, row.names = FALSE)

cat("\n[INFO] Inventory saved to:\n  ", out_file, "\n")
