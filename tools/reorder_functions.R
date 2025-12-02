## =============================================================
## Réorganisation automatique des fonctions du package
## - Regroupe et réordonne les fonctions par fichier
## - Fonctionne sur le dossier R/ du package courant
## - Produit des fichiers réordonnés dans R_reordered/
##
## Logique :
##  1) Détecte tous les blocs "nom <- function(...) { ... }"
##     + roxygen (#') immédiatement au-dessus.
##  2) Lit le NAMESPACE pour savoir quelles fonctions sont exportées.
##  3) Pour chaque fichier :
##       - écrit d’abord les fonctions exportées (alphabétique),
##       - puis les fonctions internes (alphabétique),
##       - en conservant les blocs de code tels quels,
##       - en préservant aussi les lignes hors fonctions.
##
## Sécurité :
##  - Ne modifie PAS le dossier R/
##  - Écrit dans R_reordered/ que tu peux inspecter avant de remplacer.
## =============================================================

reorder_package_functions <- function(pkg_path = ".") {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required.")
  }

  pkg <- devtools::as.package(pkg_path)
  pkg_path <- pkg$path
  rdir <- file.path(pkg_path, "R")
  if (!dir.exists(rdir)) stop("Cannot find R/ directory in package.")

  out_dir <- file.path(pkg_path, "R_reordered")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("[INFO] Package path :", pkg_path, "\n")
  cat("[INFO] R directory  :", rdir, "\n")
  cat("[INFO] Output dir   :", out_dir, "\n")

  ## --- Récupérer les fonctions exportées depuis le NAMESPACE ---
  ns_path <- file.path(pkg_path, "NAMESPACE")
  exports <- character(0)
  if (file.exists(ns_path)) {
    ns_text <- readLines(ns_path, warn = FALSE)
    # lignes du type: export(foo) ou export("foo")
    export_lines <- grep("^\\s*export\\(", ns_text, value = TRUE)
    exports <- gsub("^\\s*export\\((\"|')?([^\"')]+)(\"|')?\\)\\s*$", "\\2", export_lines)
  } else {
    # fallback: essayer de charger le namespace si installé
    pkg_name <- pkg$package
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      exports <- getNamespaceExports(pkg_name)
    }
  }
  exports <- sort(unique(exports))
  cat("[INFO] Found", length(exports), "exported symbols.\n")

  ## --- utilitaire : trouver fin de fonction à partir de la ligne idx_fun ---
  find_function_end <- function(txt, start_line) {
    n <- length(txt)
    # trouver la première { à partir de start_line
    brace_open_pos <- regexpr("\\{", txt[start_line])
    if (brace_open_pos[1] == -1) {
      # si { n'est pas sur la même ligne, chercher dans les suivantes
      i <- start_line
      while (i <= n && brace_open_pos[1] == -1) {
        i <- i + 1L
        if (i > n) break
        brace_open_pos <- regexpr("\\{", txt[i])
      }
      if (i > n || brace_open_pos[1] == -1) {
        warning("No '{' found for function starting at line ", start_line)
        return(start_line)
      }
      start <- i
    } else {
      start <- start_line
    }

    # compter les accolades à partir de 'start'
    depth <- 0L
    for (i in start:n) {
      line <- txt[i]
      # retirer les chaînes ? (on reste simple : brut)
      opens  <- lengths(regmatches(line, gregexpr("\\{", line)))
      closes <- lengths(regmatches(line, gregexpr("\\}", line)))
      depth <- depth + opens - closes
      if (depth == 0L) {
        return(i)
      }
    }
    # si on arrive ici, parenthésage cassé
    warning("Could not find matching '}' for function starting at line ", start_line)
    return(n)
  }

  ## --- parcours des fichiers R ---
  r_files <- list.files(rdir, pattern = "\\.R$", full.names = TRUE)
  cat("[INFO] Found", length(r_files), "R files.\n")

  for (f in r_files) {
    rel_name <- basename(f)
    cat("\n[FILE] Processing", rel_name, "...\n")
    txt <- readLines(f, warn = FALSE)
    n   <- length(txt)

    # indices des définitions de fonctions : "name <- function("
    fun_idx <- grep("^\\s*[^#].*<-\\s*function\\s*\\(", txt)
    if (!length(fun_idx)) {
      cat("  -> No functions detected, copying file as-is.\n")
      file.copy(f, file.path(out_dir, rel_name), overwrite = TRUE)
      next
    }

    fun_blocks <- list()
    fun_order  <- character(0)

    used_lines <- rep(FALSE, n)  # marquera les lignes appartenant à des fonctions

    for (idx in fun_idx) {
      line <- txt[idx]
      # nom de fonction
      name_match <- sub("^\\s*([A-Za-z0-9_\\.]+)\\s*<-\\s*function.*$", "\\1", line)
      if (!nzchar(name_match) || grepl("<-", name_match, fixed = TRUE)) {
        next
      }
      fun_name <- name_match

      # roxygen juste au-dessus
      start <- idx
      i <- idx - 1L
      while (i >= 1L && grepl("^#'", txt[i])) {
        start <- i
        i <- i - 1L
      }

      end <- find_function_end(txt, idx)

      block_lines <- txt[start:end]
      used_lines[start:end] <- TRUE

      if (!is.null(fun_blocks[[fun_name]])) {
        warning("Function '", fun_name, "' appears multiple times in file ", rel_name,
                ". Keeping first occurence, ignoring later one(s).")
      } else {
        fun_blocks[[fun_name]] <- block_lines
        fun_order <- c(fun_order, fun_name)
      }
    }

    if (!length(fun_blocks)) {
      cat("  -> No valid function blocks found, copying file as-is.\n")
      file.copy(f, file.path(out_dir, rel_name), overwrite = TRUE)
      next
    }

    # Lignes hors fonctions (commentaires, code global, etc.)
    other_lines <- txt[!used_lines]

    # Séparer exportées / internes
    fun_names <- names(fun_blocks)
    exported_here <- sort(intersect(fun_names, exports))
    internal_here <- sort(setdiff(fun_names, exports))

    cat("  -> Functions detected:", length(fun_names),
        " (exported:", length(exported_here),
        " internal:", length(internal_here), ")\n")

    # Construction du nouveau contenu
    new_lines <- character(0)

    # 1) Lignes hors fonctions (on les met en tête)
    if (length(other_lines)) {
      new_lines <- c(new_lines, other_lines, "", "")
    }

    # 2) Fonctions exportées
    if (length(exported_here)) {
      new_lines <- c(new_lines,
                     sprintf("## ---- Exported functions (%d) ----", length(exported_here)),
                     "")
      for (nm in exported_here) {
        new_lines <- c(new_lines, fun_blocks[[nm]], "", "")
      }
    }

    # 3) Fonctions internes
    if (length(internal_here)) {
      new_lines <- c(new_lines,
                     sprintf("## ---- Internal functions (%d) ----", length(internal_here)),
                     "")
      for (nm in internal_here) {
        new_lines <- c(new_lines, fun_blocks[[nm]], "", "")
      }
    }

    # Écriture dans R_reordered/
    out_file <- file.path(out_dir, rel_name)
    writeLines(new_lines, con = out_file, useBytes = TRUE)
    cat("  -> Wrote reordered file to:", out_file, "\n")
  }

  cat("\n[INFO] Done. Reordered files are in:", out_dir, "\n")
  invisible(TRUE)
}

## Pour lancer depuis la racine du package :
## reorder_package_functions(".")
