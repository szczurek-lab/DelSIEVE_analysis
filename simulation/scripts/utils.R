# Title     : Utils
# Objective : Provide utensil functions for R
# Created by: senbaikang
# Created on: 07.03.21

here::i_am("scripts/utils.R")
library(here)
here()

if(!"scales" %in% installed.packages()){
  install.packages("scales", repos = repository)
}
library(scales)

if (!"ggpubr" %in% installed.packages()) {
  install.packages("ggpubr", repos = repository)
}
library(ggpubr)

source(file = here("scripts", "magic_color.R"))

load.data <- function(
  file,
  sep = "\t",
  additional_labels = NULL,
  ...
) {
  if (!file.exists(file)) {
    warning(paste(file, "does not exist!"))
    return(NULL)
  }

  data <- read.table(file = file, header = TRUE, sep = sep, ...)

  if (!is.null(additional_labels)) {
    stopifnot(length(additional_labels[[1]]) == length(additional_labels[[2]]))

    for (i in seq_len(length(additional_labels[[1]]))) {
      data[[additional_labels[[1]][i]]] <- additional_labels[[2]][i]
    }
  }

  data
}

load.rds <- function(
  file,
  additional_labels = NULL
) {
  if (!file.exists(file)) {
    warning(paste(file, "does not exist!"))
    return(NULL)
  }

  data <- readRDS(file = file)

  if (!is.null(additional_labels)) {
    stopifnot(length(additional_labels[[1]]) == length(additional_labels[[2]]))

    for (i in seq_len(length(additional_labels[[1]]))) {
      data[[additional_labels[[1]][i]]] <- additional_labels[[2]][i]
    }
  }

  return(data)
}

scientific <- function(vals){
  sapply(
    vals,
    function(x) {
      if (x == 0)
        return("0")

      x <- ifelse(
        grepl("[eE]", x),
        x,
        label_scientific()(x)
      )

      if (grepl("^1[eE]", x)) {
        x <- gsub("^1[eE]", "10^", x)
      } else {
        x <- gsub("[eE]", "%*%10^", x)
      }

      parse(text = gsub(
        "[+]",
        "",
        x
      ))
    }
  )
}

define.legend <- function(
  use.common.legend,
  common.legend,
  original.legend
) {
  if (use.common.legend) {
    common.legend
  } else {
    original.legend
  }
}

adjust_facet_length <- function(gg.obj, dir = "h") {
  gg.obj

  # if (dir == "h") {
  #   .dir <- "widths"
  #   panel.scales <- "panel_scales_x"
  # }
  # else if (dir == "v") {
  #   .dir <- "heights"
  #   panel.scales <- "panel_scales_y"
  # }
  #
  # # convert ggplot object to grob object
  # gp <- ggplotGrob(gg.obj)
  #
  # # get gtable columns corresponding to the facets
  # facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
  #
  # # get the number of unique x-axis values per facet
  # var <- sapply(ggplot_build(gg.obj)$layout[[panel.scales]],
  #                 function(l) length(l$range$range))
  #
  # # change the relative widths of the facet columns based on
  # # how many unique x-axis values are in each facet
  # gp[[.dir]][facet.columns] <- gp[[.dir]][facet.columns] * var
  #
  # as_ggplot(gp)
}
