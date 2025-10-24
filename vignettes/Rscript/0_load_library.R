suppressMessages(suppressWarnings(
  library(nimble,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(
  library(parallel,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(
  library(glue,warn.conflicts = F, quietly = T)))  # for character pasting
suppressMessages(suppressWarnings(
  library(stringr,warn.conflicts = F, quietly = T)))  # for character manipulations
suppressMessages(suppressWarnings(
  library(cli,warn.conflicts = F, quietly = T)))  # for nice communication
suppressMessages(suppressWarnings(
  library(dplyr,warn.conflicts = F, quietly = T)))

suppressMessages(suppressWarnings(
  library(rmarkdown,warn.conflicts = F, quietly = T))) # for html markdown
suppressMessages(suppressWarnings(
  library(ggplot2,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(
  library(coda,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(
  library(microbenchmark,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(
  library(rstan,warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(
  library(nimbleHMC,warn.conflicts = F, quietly = T)))