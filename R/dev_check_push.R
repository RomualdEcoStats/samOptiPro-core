dev_check_push <- function(commit_message = NULL,
                           run_lintr = FALSE,
                           cran_like = FALSE) {
  
  # Active renv si présent
  if (file.exists("renv/activate.R")) source("renv/activate.R")
  
  # Dépendances utilitaires
  if (!requireNamespace("devtools", quietly = TRUE)) stop("Install 'devtools'.")
  if (!requireNamespace("gert",     quietly = TRUE)) stop("Install 'gert'.")
  
  message("Documenting with roxygen2 …")
  devtools::document(quiet = TRUE)
  
  if (isTRUE(run_lintr) && requireNamespace("lintr", quietly = TRUE)) {
    message("Linting …")
    lintr::lint_package()
  }
  
  message("Checking …")
  devtools::check(
    document = FALSE,
    cran     = isTRUE(cran_like),
    quiet    = TRUE,
    error_on = "error"
  )
  message("Checks passed ✔")
  
  # Git: add/commit si nécessaire
  st <- gert::git_status()
  if (nrow(st)) {
    if (is.null(commit_message) &&
        requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      commit_message <- rstudioapi::showPrompt(
        "Commit message",
        "Enter a concise commit message:",
        default = "chore: update docs & NAMESPACE"
      )
    }
    if (is.null(commit_message) || !nzchar(commit_message)) {
      commit_message <- "chore: update package"
    }
    gert::git_add(".")
    gert::git_commit(commit_message)
  } else {
    message("Nothing to commit — working tree clean.")
  }
  
  # Push
  br <- gert::git_branch(); current <- br$name[br$head]
  message(sprintf("Pushing to origin/%s …", current))
  gert::git_push(remote = "origin", set_upstream = TRUE)
  message("Done ✔")
}
