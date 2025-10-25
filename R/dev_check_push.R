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
################################################################
# scripts/push_to_github.R
push_to_github <- function(commit_message = NULL,
                           add_paths = c("R", "man", "NAMESPACE", "DESCRIPTION",
                                         "vignettes", "inst", "tests",
                                         "README.md", "NEWS.md", ".Rbuildignore"),
                           check_large = TRUE,
                           large_threshold_mb = 45,
                           pull_before_push = TRUE) {
  
  if (!requireNamespace("gert", quietly = TRUE)) {
    stop("Please install {gert}: install.packages('gert')")
  }
  
  repo_root <- tryCatch(gert::git_find(), error = function(e) NULL)
  if (is.null(repo_root)) stop("Not inside a Git repository.")
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(repo_root)
  
  # 1) Warn on large files (common cause of push failures on GitHub)
  if (isTRUE(check_large)) {
    all_files <- list.files(".", recursive = TRUE, all.files = FALSE, full.names = TRUE)
    info <- file.info(all_files)
    big <- rownames(info)[!info$isdir & info$size > (large_threshold_mb * 1024^2)]
    if (length(big)) {
      message("⚠️ Large files detected (> ", large_threshold_mb, " MB):")
      for (b in big) message("  - ", b, " (", round(file.info(b)$size/1024^2, 1), " MB)")
      message("Consider adding them to .gitignore or using Git LFS if needed.")
    }
  }
  
  # 2) Stage files (only existing ones)
  add_paths <- unique(add_paths)
  add_paths <- add_paths[file.exists(add_paths)]
  if (length(add_paths)) {
    gert::git_add(add_paths)
  }
  
  # Stage deletions/renames and any other modifications
  gert::git_add(".")
  
  # 3) If nothing to commit -> exit early
  st <- gert::git_status()
  if (!nrow(st)) {
    message("✅ Nothing to commit — working tree is clean.")
  } else {
    # 4) Commit message (ask in RStudio if not provided)
    if (is.null(commit_message) || !nzchar(commit_message)) {
      if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
        commit_message <- rstudioapi::showPrompt(
          "Commit message",
          "Enter a concise commit message:",
          default = "chore: update package files"
        )
      }
      if (is.null(commit_message) || !nzchar(commit_message)) {
        commit_message <- paste0("chore: update (", format(Sys.time(), "%Y-%m-%d %H:%M"), ")")
      }
    }
    gert::git_commit(commit_message)
    message("📝 Commit created: ", commit_message)
  }
  
  # 5) Ensure we know current branch & remote
  br <- gert::git_branch()
  current <- br$name[br$head]
  if (is.na(current) || !nzchar(current)) current <- "main"
  
  has_origin <- any(vapply(gert::git_remote_list()$name, identical, logical(1), "origin"))
  if (!has_origin) stop("No 'origin' remote configured. Set it with: gert::git_remote_add('origin', '<url>')")
  
  # 6) Pull (rebase) before pushing, to avoid non fast-forward
  if (isTRUE(pull_before_push)) {
    message("🔄 Pulling from origin/", current, " (rebase = TRUE)…")
    try(gert::git_pull(repo = repo_root, rebase = TRUE), silent = TRUE)
  }
  
  # 7) Push (set upstream if first push)
  message("⤴️  Pushing to origin/", current, " …")
  ok <- try(gert::git_push(remote = "origin", refspec = current, set_upstream = TRUE), silent = TRUE)
  if (inherits(ok, "try-error")) {
    message("Push failed. Attempting a pull --rebase and re-push…")
    try(gert::git_pull(repo = repo_root, rebase = TRUE), silent = TRUE)
    gert::git_push(remote = "origin", refspec = current, set_upstream = TRUE)
  }
  
  message("✅ Done. Pushed branch: ", current)
  invisible(TRUE)
}
