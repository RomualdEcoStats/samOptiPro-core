## =============================================================
## Inventaire propre des fonctions du package (version R_reordered)
##  - Une ligne par fonction
##  - Arguments
##  - Description roxygen courte (sans @title, @export, etc.)
## =============================================================

library(stringr)

pkg <- "samOptiPro"
cat("[INFO] Scanning package (reordered):", pkg, "\n")

# chemin du package (point de départ = DESCRIPTION)
pkg_path <- devtools::as.package(".")$path

rdir_reordered <- file.path(pkg_path, "R_reordered")
if (!dir.exists(rdir_reordered)) {
  stop("Cannot find R_reordered/ directory in package.")
}

files <- list.files(rdir_reordered, pattern = "\\.R$", full.names = TRUE)
cat("[INFO] Found", length(files), "R files in R_reordered/.\n")

# namespace pour récupérer les formals() des fonctions déjà exportées
ns <- asNamespace(pkg)

extract_args <- function(fname) {
  fx <- tryCatch(get(fname, envir = ns), error = function(e) NULL)
  if (!is.function(fx)) return(NA_character_)
  args <- names(formals(fx))
  paste(args, collapse = ", ")
}

results <- list()

for (f in files) {
  txt <- readLines(f, warn = FALSE)

  is_roxy <- grepl("^#'", txt)
  is_def  <- grepl("^\\s*([A-Za-z0-9_.]+)\\s*<-\\s*function", txt)

  def_idx <- which(is_def)
  if (!length(def_idx)) next

  for (idx in def_idx) {
    line <- txt[idx]
    name <- str_match(line, "^\\s*([A-Za-z0-9_.]+)\\s*<-\\s*function")[, 2]
    if (is.na(name) || name == "") next

    ## bloc roxygen directement au-dessus
    j <- idx - 1L
    if (j < 1L || !is_roxy[j]) {
      roxy_block <- character(0)
    } else {
      while (j >= 1L && is_roxy[j]) j <- j - 1L
      roxy_block <- txt[(j + 1L):(idx - 1L)]
    }

    # nettoyage du bloc roxygen
    roxy_block <- gsub("^#' ?", "", roxy_block)

    # première ligne non vide qui ne commence pas par @
    desc <- ""
    if (length(roxy_block)) {
      rb <- roxy_block[roxy_block != ""]
      rb <- rb[!str_starts(rb, "@")]
      if (length(rb)) desc <- rb[1L]
    }

    args <- extract_args(name)

    results[[length(results) + 1L]] <- list(
      file        = basename(f),
      fun_name    = name,
      arguments   = args,
      description = desc
    )
  }
}

df2 <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))

## dédoublonnage : une ligne par fonction
df2_clean <- df2 %>%
  group_by(fun_name) %>%
  arrange(file) %>%
  summarise(
    file        = first(file),
    arguments   = first(arguments),
    description = {
      d <- description[description != ""]
      if (length(d) == 0) "" else d[1]
    },
    .groups = "drop"
  )

out_dir  <- file.path(pkg_path, "inst", "extdata")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_file   <- file.path(out_dir, "functions_inventory_reordered_raw.csv")
clean_file <- file.path(out_dir, "functions_inventory_reordered_clean.csv")

write.csv(df2,      raw_file,   row.names = FALSE)
write.csv(df2_clean,clean_file, row.names = FALSE)

cat("\n[INFO] Reordered inventories saved to:\n  ",
    raw_file,   "\n  ",
    clean_file, "\n")
