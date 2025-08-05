# Title     : Gather results from simulations
# Created by: senbaikang
# Created on: 07.03.21

here::i_am("scripts/gather_results.R")
library(here)
here()

repository <- "https://stat.ethz.ch/CRAN/"

if(!"dplyr" %in% installed.packages()){
  install.packages("dplyr", repos = repository)
}
library(dplyr)

if (!"ggplot2" %in% installed.packages()) {
  install.packages("ggplot2", repos = repository)
}
library(ggplot2)

if (!"ggforce" %in% installed.packages()) {
  install.packages("ggforce", repos = repository)
}
library(ggforce)

if (!"lemon" %in% installed.packages()) {
  install.packages("lemon", repos = repository)
}
library(lemon)

# library(ggpattern)


source(file = here("scripts", "utils.R"))


# Constants --------------------------------

ado.modes <- c("sa" = "single_ado")
true.ado.rate <- c("single_ado" = 0.3)
# if (length(ado.modes) == 2) {
  output.prefix <- "joint_plots"
# } else if (ado.modes == "single_ado") {
#   output.prefix <- "single_ado_plots"
# } else if (ado.modes == "locus_dropout") {
#   output.prefix <- "locus_dropout_plots"
# }

box.shape <- "box" # "box" or "violin"
if (box.shape == "box") {
  box.func <- geom_boxplot
  jitter.func <- geom_point
} else if (box.shape == "violin") {
  box.func <- geom_violin
  jitter.func <- geom_sina
}

true.eff.seq.err.rate <- 1.999E-3
true.wildtype.overdispersion <- 100.0
true.alt.overdispersion <- 2.5

isa.mu.rate.threshold <- 2.0E-6
prop.threshold <- 0.001

box.linewidth <- 1.5 / .pt
box.width <- 0.6
jitter.alpha <- 0.6
jitter.width <- 0.2
dodge.width <- 0.7
dot.size <- 0.6 / .pt

bg.color <- "gray60"

use.common.legend <- TRUE

text <- element_text(size = 7)
axis.text.x <- element_text(size = 5)
axis.text.y <- element_text(size = 6)
axis.title <- element_text(size = 7)
legend.text <- element_text(size = 7)
plot.title <- element_text(size = 8)
legend.title <- element_text(size = 7)
strip.text.x <- element_text(size = 7) #, margin = margin(0, -1, 0, -1))
strip.text.y <- element_text(size = 7)

new.col.names <- c(
  "simulated_data",
  "mutation_rate",
  "rela_deletion_rate",
  "simulated_dropout_type",
  "cnv_rate",
  "doublet_rate",
  "cnv_doublet_status"
)

# For facets
cna.doublet.status <- c(
  "CNA\nNo doublet",
  "No CNA\nDoublet",
  "CNA\nDoublet"
)
names(cna.doublet.status) <- c("c", "d", "cd")

doublet.rate.colors <- wesanderson::wes_palette("Royal2")[c(3, 4, 5)]
names(doublet.rate.colors) <- c("0", "0.02", "0.1")

cna.rate.colors <- wesanderson::wes_palette("AsteroidCity3")[c(1, 3, 4)]
names(cna.rate.colors) <- c("0", "1/3", "2/3")

cna.sites.prop.labels <- c("No CNAs", "CNAs in 1/3 of all sites", "CNAs in 2/3 of all sites")
names(cna.sites.prop.labels) <- c("0", "1/3", "2/3")

cna.status.labels <- c("None", "Excluded", "Included")
names(cna.status.labels) <- c("no_cnv", "exc_cnv", "inc_cnv")

efficiency.labels <- c("Scenario 1", "Scenario 2")
names(efficiency.labels) <- c(1, 2)

efficiency.colors <- c("#b95000", "#e64b35", "#3c5488", "#149af2")
names(efficiency.colors) <- c("DelSIEVE_stage_1", "DelSIEVE_stage_2", "SIEVE_stage_1", "SIEVE_stage_2")

# Import simulation results -------------------------

setups <- list(
  "105_cna_1" = c("105_cna_1", 8.0E-6, 0.25, "single ADO"),
  "105_cna_2" = c("105_cna_2", 8.0E-6, 0.25, "single ADO"),
  "105_doublet_1" = c("105_doublet_1", 8.0E-6, 0.25, "single ADO"),
  "105_doublet_2" = c("105_doublet_2", 8.0E-6, 0.25, "single ADO")
)

add.setups <- list(
  c("105", 8.0E-6, 0.25, "single ADO", "0", 0, "c"),
  c("105", 8.0E-6, 0.25, "single ADO", "0", 0, "d")
)

tree.info <- do.call(
  rbind,
  lapply(
    names(ado.modes),
    function(x) {
      do.call(
        rbind,
        lapply(
          names(setups),
          function(y) {
            load.data(
              file = here(ado.modes[[x]], paste0("simulated_data_", y), "trees_info_updated.tsv"),
              additional_labels = list(
                new.col.names[1:4],
                setups[[y]]
              )
            )
          }
        )
      ) %>%
        filter(
          is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
        )
    }
  )
) %>%
  mutate(
    cnv_doublet_status = case_when(
      doublet_rate == 0 & cnv_rate != "0" ~ "c",
      doublet_rate != 0 & cnv_rate == "0" ~ "d",
      doublet_rate != 0 & cnv_rate != "0" ~ "cd"
    )
  )

tree.info <- rbind(
  do.call(
    rbind,
    lapply(
      names(ado.modes),
      function(x) {
        do.call(
          rbind,
          lapply(
            add.setups,
            function(y) {
              load.data(
                file = here(ado.modes[[x]], paste0("simulated_data_", y[1]), "trees_info_updated.tsv"),
                additional_labels = list(
                  new.col.names,
                  y
                )
              )
            }
          )
        ) %>%
          filter(
            is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
          )
      }
    )
  )[colnames(tree.info)],
  tree.info
)

site.info <- do.call(
  rbind,
  lapply(
    names(ado.modes),
    function(x) {
      do.call(
        rbind,
        lapply(
          names(setups),
          function(y) {
            load.data(
              file = here(ado.modes[[x]], paste0("simulated_data_", y), "sites_info.tsv"),
              additional_labels = list(
                new.col.names[1:4],
                setups[[y]]
              ),
              check.names = FALSE
            )
          }
        )
      ) %>%
        filter(
          is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
        )
    }
  )
) %>%
  mutate(
    cnv_doublet_status = case_when(
      doublet_rate == 0 & cnv_rate != "0" ~ "c",
      doublet_rate != 0 & cnv_rate == "0" ~ "d",
      doublet_rate != 0 & cnv_rate != "0" ~ "cd"
    )
  )

site.info <- rbind(
  do.call(
    rbind,
    lapply(
      names(ado.modes),
      function(x) {
        do.call(
          rbind,
          lapply(
            add.setups,
            function(y) {
              load.data(
                file = here(ado.modes[[x]], paste0("simulated_data_", y[1]), "sites_info.tsv"),
                additional_labels = list(
                  new.col.names,
                  y
                ),
                check.names = FALSE
              )
            }
          )
        ) %>%
          filter(
            is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
          )
      }
    )
  )[colnames(site.info)],
  site.info
)

param.info <- do.call(
  rbind,
  lapply(
    names(ado.modes),
    function(x) {
      do.call(
        rbind,
        lapply(
          names(setups),
          function(y) {
            load.data(
              file = here(ado.modes[[x]], paste0("simulated_data_", y), "params_info.tsv"),
              additional_labels = list(
                new.col.names[1:4],
                setups[[y]]
              )
            )
          }
        )
      ) %>%
        filter(
          is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
        )
    }
  )
) %>%
  mutate(
    cnv_doublet_status = case_when(
      doublet_rate == 0 & cnv_rate != "0" ~ "c",
      doublet_rate != 0 & cnv_rate == "0" ~ "d",
      doublet_rate != 0 & cnv_rate != "0" ~ "cd"
    )
  )

param.info <- rbind(
  do.call(
    rbind,
    lapply(
      names(ado.modes),
      function(x) {
        do.call(
          rbind,
          lapply(
            add.setups,
            function(y) {
              load.data(
                file = here(ado.modes[[x]], paste0("simulated_data_", y[1]), "params_info.tsv"),
                additional_labels = list(
                  new.col.names,
                  y
                )
              )
            }
          )
        ) %>%
          filter(
            is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
          )
      }
    )
  )[colnames(param.info)],
  param.info
)

var.info <- do.call(
  rbind,
  lapply(
    names(ado.modes),
    function(x) {
      do.call(
        rbind,
        lapply(
          names(setups),
          function(y) {
            load.data(
              file = here(ado.modes[[x]], paste0("simulated_data_", y), "variants_info.tsv"),
              additional_labels = list(
                new.col.names[1:4],
                setups[[y]]
              )
            )
          }
        )
      ) %>%
        filter(
          is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
        )
    }
  )
) %>%
  mutate(
    cnv_doublet_status = case_when(
      doublet_rate == 0 & cnv_rate != "0" ~ "c",
      doublet_rate != 0 & cnv_rate == "0" ~ "d",
      doublet_rate != 0 & cnv_rate != "0" ~ "cd"
    )
  )

var.info <- rbind(
  do.call(
    rbind,
    lapply(
      names(ado.modes),
      function(x) {
        do.call(
          rbind,
          lapply(
            add.setups,
            function(y) {
              load.data(
                file = here(ado.modes[[x]], paste0("simulated_data_", y[1]), "variants_info.tsv"),
                additional_labels = list(
                  new.col.names,
                  y
                )
              )
            }
          )
        ) %>%
          filter(
            is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
          ) %>%
          mutate(
            doublet_rate = 0,
            cnv_rate = 0
          )
      }
    )
  )[colnames(var.info)],
  var.info
) %>%
  mutate(
    true_hetero_mu = true_positive_hetero_mu + false_negative_hetero_mu,
    hetero_mu = true_positive_hetero_mu + false_negative_hetero_mu + true_negative_hetero_mu + false_positive_hetero_mu,
    homo_mu = true_positive_homo_mu + false_negative_homo_mu + true_negative_homo_mu + false_positive_homo_mu,
    true_del_alt_left = true_positive_del_alt_left + false_negative_del_alt_left,
    del_alt_left = true_positive_del_alt_left + false_negative_del_alt_left + true_negative_del_alt_left + false_positive_del_alt_left,
    true_del_ref_left = true_positive_del_ref_left + false_negative_del_ref_left,
    del_ref_left = true_positive_del_ref_left + false_negative_del_ref_left + true_negative_del_ref_left + false_positive_del_ref_left,
    true_all_del = true_positive_all_del + false_negative_all_del,
    all_del = true_positive_all_del + false_negative_all_del + true_negative_all_del + false_positive_all_del
  ) %>%
  mutate(
    prop_true_hetero_mu = true_hetero_mu / hetero_mu,
    prop_true_homo_mu = true_homo_mu / homo_mu,
    prop_true_del_alt_left = true_del_alt_left / del_alt_left,
    prop_true_del_ref_left = true_del_ref_left / del_ref_left,
    prop_true_all_del = true_all_del / all_del
  )

ado.info <- do.call(
  rbind,
  lapply(
    names(ado.modes),
    function(x) {
      do.call(
        rbind,
        lapply(
          names(setups),
          function(y) {
            load.data(
              file = here(ado.modes[[x]], paste0("simulated_data_", y), "ado_info.tsv"),
              additional_labels = list(
                new.col.names[1:4],
                setups[[y]]
              )
            )
          }
        )
      ) %>%
        filter(
          is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
        )
    }
  )
) %>%
  mutate(
    cnv_doublet_status = case_when(
      doublet_rate == 0 & cnv_rate != "0" ~ "c",
      doublet_rate != 0 & cnv_rate == "0" ~ "d",
      doublet_rate != 0 & cnv_rate != "0" ~ "cd"
    )
  )

ado.info <- rbind(
  do.call(
    rbind,
    lapply(
      names(ado.modes),
      function(x) {
        do.call(
          rbind,
          lapply(
            add.setups,
            function(y) {
              load.data(
                file = here(ado.modes[[x]], paste0("simulated_data_", y[1]), "ado_info.tsv"),
                additional_labels = list(
                  new.col.names,
                  y
                )
              )
            }
          )
        ) %>%
          filter(
            is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
          ) %>%
          mutate(
            doublet_rate = 0,
            cnv_rate = 0
          )
      }
    )
  )[colnames(ado.info)],
  ado.info
)

# allelic.info <- do.call(
#   rbind,
#   lapply(
#     names(ado.modes),
#     function(x) {
#       do.call(
#         rbind,
#         lapply(
#           names(setups),
#           function(y) {
#             load.rds(
#               file = here(ado.modes[[x]], paste0("simulated_data_", y), "allelic_info.rds"),
#               additional_labels = list(
#                 new.col.names[1:4],
#                 c(setups[[y]], ado.modes[[x]])
#               )
#             )
#           }
#         )
#       ) %>%
#         filter(
#           is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
#         )
#     }
#   )
# )

# parallel.mut <- do.call(
#   rbind,
#   lapply(
#     names(ado.modes),
#     function(x) {
#       do.call(
#         rbind,
#         lapply(
#           names(setups),
#           function(y) {
#             load.data(
#               file = here(ado.modes[[x]], paste0("simulated_data_", y), "summary.tsv"),
#               sep = " ",
#               additional_labels = list(
#                 new.col.names[1:4],
#                 c(setups[[y]], ado.modes[[x]])
#               )
#             )
#           }
#         )
#       ) %>%
#         filter(
#           is.na(tool_setup) | grepl(paste0(".+", x, ".+"), tool_setup) | tool_setup == "true_parameters"
#         )
#     }
#   )
# )

# For efficiency benchmarking
new.eff.col.names <- c(
  "simulated_data",
  "mutation_rate",
  "rela_deletion_rate",
  "simulated_dropout_type",
  "dataset",
  "tool",
  "stage",
  "num_candidate_var_sites",
  "scenario"
)

eff.setups <- list(
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0001", "SIEVE", "stage_1", 814, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0001", "SIEVE", "stage_2", 814, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0001", "DelSIEVE", "stage_1", 814, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0001", "DelSIEVE", "stage_2", 814, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0002", "SIEVE", "stage_1", 798, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0002", "SIEVE", "stage_2", 798, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0002", "DelSIEVE", "stage_1", 798, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0002", "DelSIEVE", "stage_2", 798, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0003", "SIEVE", "stage_1", 750, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0003", "SIEVE", "stage_2", 750, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0003", "DelSIEVE", "stage_1", 750, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0003", "DelSIEVE", "stage_2", 750, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0004", "SIEVE", "stage_1", 700, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0004", "SIEVE", "stage_2", 700, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0004", "DelSIEVE", "stage_1", 700, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0004", "DelSIEVE", "stage_2", 700, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0005", "SIEVE", "stage_1", 851, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0005", "SIEVE", "stage_2", 851, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0005", "DelSIEVE", "stage_1", 851, 1),
  c("405_efficiency_1", 8.0E-6, 0.25, "single ADO", "0005", "DelSIEVE", "stage_2", 851, 1),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0001", "SIEVE", "stage_1", 1644, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0001", "SIEVE", "stage_2", 1644, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0001", "DelSIEVE", "stage_1", 1644, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0001", "DelSIEVE", "stage_2", 1644, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0002", "SIEVE", "stage_1", 1585, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0002", "SIEVE", "stage_2", 1585, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0002", "DelSIEVE", "stage_1", 1585, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0002", "DelSIEVE", "stage_2", 1585, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0003", "SIEVE", "stage_1", 1513, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0003", "SIEVE", "stage_2", 1513, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0003", "DelSIEVE", "stage_1", 1513, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0003", "DelSIEVE", "stage_2", 1513, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0004", "SIEVE", "stage_1", 1396, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0004", "SIEVE", "stage_2", 1396, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0004", "DelSIEVE", "stage_1", 1396, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0004", "DelSIEVE", "stage_2", 1396, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0005", "SIEVE", "stage_1", 1690, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0005", "SIEVE", "stage_2", 1690, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0005", "DelSIEVE", "stage_1", 1690, 2),
  c("405_efficiency_2", 8.0E-6, 0.25, "single ADO", "0005", "DelSIEVE", "stage_2", 1690, 2)
)

efficiency <- do.call(
  rbind,
  lapply(
    eff.setups,
    function(x) {
      load.data(
        file = here(
          "single_ado",
          paste0("simulated_data_", x[1]),
          "benchmark_results",
          "efficiency",
          "indelsieve",
          "mt",
          "no_cnv",
          ifelse(
            x[6] == "SIEVE",
            ifelse(
              x[7] == "stage_1",
              "smc_mu_sa_stage_1",
              "smc_mu_sa_stage_1-smc_mu_sa_stage_2"
            ),
            ifelse(
              x[7] == "stage_1",
              "smc_mudel_sa_stage_1",
              "smc_mudel_sa_stage_1-smc_mudel_sa_stage_2"
            )
          ),
          paste0(x[5], ".tsv")
        ),
        additional_labels = list(
          new.eff.col.names,
          x
        )
      ) %>%
        select(
          s, mean_load, cpu_time, all_of(new.eff.col.names)
        ) %>%
        mutate(
          mean_s = mean(s),
          mean_mean_load = mean(mean_load),
          mean_cpu_time = mean(cpu_time),
          tool_setup = paste(
            tool,
            stage,
            sep = "_"
          )
        ) %>%
        select(
          mean_s, mean_mean_load, mean_cpu_time, all_of(new.eff.col.names), tool_setup
        ) %>%
        distinct()
    }
  )
)

# Rename tools.
tree.info <- tree.info %>%
  mutate(tool = case_when(
    grepl("^sciphin$", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphin$", "SCIPhIN", tool, fixed = FALSE),
    grepl("^sciphi$", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi$", "SCIPhI", tool, fixed = FALSE),
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
    grepl("^indelsieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^indelsieve", "DelSIEVE", tool, fixed = FALSE),
    grepl("^sifit", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sifit", "SiFit", tool, fixed = FALSE)
  )) %>%
  mutate(simulated_dropout_type = case_when(
    grepl("single_ado", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "ADO",
    grepl("locus_dropout", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "LDO"
  )) %>%
  filter(grepl("^DelSIEVE$", tool, ignore.case = TRUE, fixed = FALSE))

site.info <- site.info %>%
  mutate(tool = case_when(
    grepl("^sciphin$", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphin$", "SCIPhIN", tool, fixed = FALSE),
    grepl("^sciphi$", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi$", "SCIPhI", tool, fixed = FALSE),
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
    grepl("^indelsieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^indelsieve", "DelSIEVE", tool, fixed = FALSE),
    grepl("^monovar", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^monovar", "Monovar", tool, fixed = FALSE)
  )) %>%
  mutate(simulated_dropout_type = case_when(
    grepl("single_ado", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "ADO",
    grepl("locus_dropout", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "LDO"
  )) %>%
  filter(grepl("^DelSIEVE$", tool, ignore.case = TRUE, fixed = FALSE))

param.info <- param.info %>%
  mutate(tool = case_when(
    grepl("^sciphin$", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphin$", "SCIPhIN", tool, fixed = FALSE),
    grepl("^sciphi$", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi$", "SCIPhI", tool, fixed = FALSE),
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
    grepl("^indelsieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^indelsieve", "DelSIEVE", tool, fixed = FALSE),
    grepl("^sifit", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sifit", "SiFit", tool, fixed = FALSE)
  )) %>%
  mutate(simulated_dropout_type = case_when(
    grepl("single_ado", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "ADO",
    grepl("locus_dropout", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "LDO"
  )) %>%
  filter(grepl("^DelSIEVE$", tool, ignore.case = TRUE, fixed = FALSE))

var.info <- var.info %>%
  mutate(tool = case_when(
    grepl("^sciphin$", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphin$", "SCIPhIN", tool, fixed = FALSE),
    grepl("^sciphi$", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi$", "SCIPhI", tool, fixed = FALSE),
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
    grepl("^indelsieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^indelsieve", "DelSIEVE", tool, fixed = FALSE),
    grepl("^monovar", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^monovar", "Monovar", tool, fixed = FALSE)
  )) %>%
  mutate(simulated_dropout_type = case_when(
    grepl("single_ado", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "ADO",
    grepl("locus_dropout", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "LDO"
  )) %>%
  filter(grepl("^DelSIEVE$", tool, ignore.case = TRUE, fixed = FALSE))

# fsa.var.info <- fsa.var.info %>%
#   mutate(tool = case_when(
#     grepl("^sciphin$", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphin$", "SCIPhIN", tool, fixed = FALSE),
#     grepl("^sciphi$", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi$", "SCIPhI", tool, fixed = FALSE),
#     grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
#     grepl("^indelsieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^indelsieve", "DelSIEVE", tool, fixed = FALSE),
#     grepl("^monovar", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^monovar", "Monovar", tool, fixed = FALSE)
#   )) %>%
#   mutate(simulated_dropout_type = case_when(
#     grepl("single_ado", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "Single ADO",
#     grepl("locus_dropout", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "Locus dropout"
#   ))

ado.info <- ado.info %>%
  mutate(tool = case_when(
    grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
    grepl("^indelsieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^indelsieve", "DelSIEVE", tool, fixed = FALSE),
  )) %>%
  mutate(simulated_dropout_type = case_when(
    grepl("single_ado", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "ADO",
    grepl("locus_dropout", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "LDO"
  ))

# allelic.info <- allelic.info %>%
#   mutate(tool = case_when(
#     grepl("^sieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sieve", "SIEVE", tool, fixed = FALSE),
#     grepl("^indelsieve", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^indelsieve", "DelSIEVE", tool, fixed = FALSE)
#   )) %>%
#   mutate(simulated_dropout_type = case_when(
#     grepl("single_ado", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "ADO",
#     grepl("locus_dropout", simulated_dropout_type, ignore.case = TRUE, fixed = TRUE) ~ "LDO"
#   ))

# efficiency <- efficiency %>%
#   mutate(
#     tool = case_when(
#       grepl("sieve", tool, ignore.case = TRUE) ~ paste(
#         sub("^sieve", "SIEVE", tool, fixed = FALSE),
#         if_else(grepl("-", tool_setup, fixed = TRUE), "stage 2", "stage 1")
#       ),
#       grepl("^cellphy", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^cellphy", "CellPhy", x = tool, fixed = FALSE),
#       grepl("^sciphi", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sciphi", "SCIPhI", tool, fixed = FALSE),
#       grepl("^monovar", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^monovar", "Monovar", tool, fixed = FALSE),
#       grepl("^sifit", tool, ignore.case = TRUE, fixed = FALSE) ~ sub("^sifit", "SiFit", tool, fixed = FALSE)
#     )
#   ) %>%
#   select(-tool_setup)

# process columns
tree.info$cell_num <- as.factor(tree.info$cell_num)
tree.info$coverage_mean <- as.factor(tree.info$coverage_mean)
tree.info$coverage_variance <- as.factor(tree.info$coverage_variance)
tree.info$dataset <- as.factor(tree.info$dataset)
tree.info$tool <- as.factor(tree.info$tool)
tree.info$snv_type <- as.factor(tree.info$snv_type)
tree.info$tool_setup <- as.factor(tree.info$tool_setup)
tree.info$doublet_rate <- as.factor(tree.info$doublet_rate)
tree.info$cnv_rate<- as.factor(tree.info$cnv_rate)
tree.info$data_type <- factor(tree.info$data_type, levels = names(cna.status.labels))
tree.info$max_clades <- as.numeric(tree.info$max_clades)
tree.info$RF_distance <- as.numeric(tree.info$RF_distance)
tree.info$normalized_RF_distance <- as.numeric(tree.info$normalized_RF_distance)
tree.info$weighted_RF_distance <- as.numeric(tree.info$weighted_RF_distance)
tree.info$rooted_branch_score_difference <- as.numeric(tree.info$rooted_branch_score_difference)
tree.info[[new.col.names[1]]] <- factor(tree.info[[new.col.names[1]]])
tree.info[[new.col.names[2]]] <- factor(tree.info[[new.col.names[2]]], levels = sort(as.numeric(unique(tree.info[[new.col.names[2]]]))))
tree.info[[new.col.names[3]]] <- factor(tree.info[[new.col.names[3]]])
tree.info[[new.col.names[4]]] <- factor(tree.info[[new.col.names[4]]])
tree.info[[new.col.names[7]]] <- factor(tree.info[[new.col.names[7]]], levels = names(cna.doublet.status))

site.info$cell_num <- as.factor(site.info$cell_num)
site.info$coverage_mean <- as.factor(site.info$coverage_mean)
site.info$coverage_variance <- as.factor(site.info$coverage_variance)
site.info$dataset <- as.factor(site.info$dataset)
site.info$tool <- as.factor(site.info$tool)
site.info$snv_type <- as.factor(site.info$snv_type)
site.info$tool_setup <- as.factor(site.info$tool_setup)
site.info$doublet_rate <- as.factor(site.info$doublet_rate)
site.info$cnv_rate<- as.factor(site.info$cnv_rate)
site.info$data_type <- factor(site.info$data_type, levels = names(cna.status.labels))
site.info$fine_tune_type <- as.factor(site.info$fine_tune_type)
site.info$data_type <- as.factor(site.info$data_type)
site.info$num_candidate_mutated_sites <- as.integer(site.info$num_candidate_mutated_sites)
site.info$num_background_sites <- as.integer(site.info$num_background_sites)
site.info$log10_num_background_sites <- log10(site.info$num_background_sites)
site.info$tp <- as.integer(site.info$tp)
site.info$fp <- as.integer(site.info$fp)
site.info$fn <- as.integer(site.info$fn)
site.info$recall <- as.numeric(site.info$recall)
site.info$precision <- as.numeric(site.info$precision)
site.info$f1_score <- as.numeric(site.info$f1_score)
site.info$tp_M <- as.integer(site.info$tp_M)
site.info$fn_M <- as.integer(site.info$fn_M)
site.info$recall_M <- as.numeric(site.info$recall_M)
site.info$tp_D <- as.integer(site.info$tp_D)
site.info$fn_D <- as.integer(site.info$fn_D)
site.info$recall_D <- as.numeric(site.info$recall_D)
site.info[["tp_M&D&(~I)"]] <- as.integer(site.info[["tp_M&D&(~I)"]])
site.info[["fn_M&D&(~I)"]] <- as.integer(site.info[["fn_M&D&(~I)"]])
site.info[["recall_M&D&(~I)"]] <- as.numeric(site.info[["recall_M&D&(~I)"]])
site.info[[new.col.names[1]]] <- factor(site.info[[new.col.names[1]]])
site.info[[new.col.names[2]]] <- factor(site.info[[new.col.names[2]]], levels = sort(as.numeric(unique(site.info[[new.col.names[2]]]))))
site.info[[new.col.names[3]]] <- factor(site.info[[new.col.names[3]]])
site.info[[new.col.names[4]]] <- factor(site.info[[new.col.names[4]]])
site.info[[new.col.names[7]]] <- factor(site.info[[new.col.names[7]]], levels = names(cna.doublet.status))

param.info$cell_num <- as.factor(param.info$cell_num)
param.info$coverage_mean_hline <- param.info$coverage_mean
param.info$coverage_mean <- as.factor(param.info$coverage_mean)
param.info$coverage_variance_hline <- param.info$coverage_variance
param.info$coverage_variance <- as.factor(param.info$coverage_variance)
param.info$dataset <- as.factor(param.info$dataset)
param.info$tool <- as.factor(param.info$tool)
param.info$snv_type <- as.factor(param.info$snv_type)
param.info$tool_setup <- as.factor(param.info$tool_setup)
param.info$doublet_rate <- as.factor(param.info$doublet_rate)
param.info$cnv_rate<- as.factor(param.info$cnv_rate)
param.info$data_type <- factor(param.info$data_type, levels = names(cna.status.labels))
param.info$eff_seq_err_rate <- as.numeric(param.info$eff_seq_err_rate)
param.info$allelic_seq_cov <- as.numeric(param.info$allelic_seq_cov)
param.info$allelic_seq_cov_raw_var <- as.numeric(param.info$allelic_seq_cov_raw_var)
param.info$ado_rate <- as.numeric(param.info$ado_rate)
param.info$gamma_shape <- as.numeric(param.info$gamma_shape)
param.info$wild_overdispersion <- as.numeric(param.info$wild_overdispersion)
param.info$alternative_overdispersion <- as.numeric(param.info$alternative_overdispersion)
# param.info$zygosity_rate <- as.numeric(param.info$zygosity_rate)
param.info$estimates_type <- as.factor(param.info$estimates_type)
param.info$deletion_rate <- as.numeric(param.info$deletion_rate)
param.info$insertion_rate <- as.numeric(param.info$insertion_rate)
param.info$population_size <- as.numeric(param.info$population_size)
param.info[[new.col.names[1]]] <- factor(param.info[[new.col.names[1]]])
param.info[[new.col.names[2]]] <- factor(param.info[[new.col.names[2]]], levels = sort(as.numeric(unique(param.info[[new.col.names[2]]]))))
param.info$rela_del_hline <- param.info[[new.col.names[3]]]
param.info[[new.col.names[3]]] <- factor(param.info[[new.col.names[3]]])
param.info[[new.col.names[4]]] <- factor(param.info[[new.col.names[4]]])
param.info[[new.col.names[7]]] <- factor(param.info[[new.col.names[7]]], levels = names(cna.doublet.status))

var.info$cell_num <- as.factor(var.info$cell_num)
var.info$coverage_mean <- as.factor(var.info$coverage_mean)
var.info$coverage_variance <- as.factor(var.info$coverage_variance)
var.info$dataset <- as.factor(var.info$dataset)
var.info$tool <- as.factor(var.info$tool)
var.info$snv_type <- as.factor(var.info$snv_type)
var.info$tool_setup <- as.factor(var.info$tool_setup)
var.info$doublet_rate <- as.factor(var.info$doublet_rate)
var.info$cnv_rate<- as.factor(var.info$cnv_rate)
var.info$data_type <- factor(var.info$data_type, levels = names(cna.status.labels))
var.info$true_positive <- as.numeric(var.info$true_positive)
var.info$false_positive <- as.numeric(var.info$false_positive)
var.info$true_negative <- as.numeric(var.info$true_negative)
var.info$false_negative <- as.numeric(var.info$false_negative)
var.info$recall <- as.numeric(var.info$recall)
var.info$precision <- as.numeric(var.info$precision)
var.info$fall_out <- as.numeric(var.info$fall_out)
var.info$f1_score <- as.numeric(var.info$f1_score)
var.info$true_positive_mu <- as.numeric(var.info$true_positive_mu)
var.info$false_positive_mu <- as.numeric(var.info$false_positive_mu)
var.info$true_negative_mu <- as.numeric(var.info$true_negative_mu)
var.info$false_negative_mu <- as.numeric(var.info$false_negative_mu)
var.info$recall_mu <- as.numeric(var.info$recall_mu)
var.info$precision_mu <- as.numeric(var.info$precision_mu)
var.info$fall_out_mu <- as.numeric(var.info$fall_out_mu)
var.info$f1_score_mu <- as.numeric(var.info$f1_score_mu)
var.info$true_positive_hetero_mu <- as.numeric(var.info$true_positive_hetero_mu)
var.info$false_positive_hetero_mu <- as.numeric(var.info$false_positive_hetero_mu)
var.info$true_negative_hetero_mu <- as.numeric(var.info$true_negative_hetero_mu)
var.info$false_negative_hetero_mu <- as.numeric(var.info$false_negative_hetero_mu)
var.info$recall_hetero_mu <- as.numeric(var.info$recall_hetero_mu)
var.info$precision_hetero_mu <- as.numeric(var.info$precision_hetero_mu)
var.info$fall_out_hetero_mu <- as.numeric(var.info$fall_out_hetero_mu)
var.info$f1_score_hetero_mu <- as.numeric(var.info$f1_score_hetero_mu)
var.info$true_positive_homo_mu <- as.numeric(var.info$true_positive_homo_mu)
var.info$false_positive_homo_mu <- as.numeric(var.info$false_positive_homo_mu)
var.info$true_negative_homo_mu <- as.numeric(var.info$true_negative_homo_mu)
var.info$false_negative_homo_mu <- as.numeric(var.info$false_negative_homo_mu)
var.info$recall_homo_mu <- as.numeric(var.info$recall_homo_mu)
var.info$precision_homo_mu <- as.numeric(var.info$precision_homo_mu)
var.info$fall_out_homo_mu <- as.numeric(var.info$fall_out_homo_mu)
var.info$f1_score_homo_mu <- as.numeric(var.info$f1_score_homo_mu)
var.info$true_homo_mu <- as.numeric(var.info$true_homo_mu)
var.info$true_positive_del <- as.numeric(var.info$true_positive_del)
var.info$false_positive_del <- as.numeric(var.info$false_positive_del)
var.info$true_negative_del <- as.numeric(var.info$true_negative_del)
var.info$false_negative_del <- as.numeric(var.info$false_negative_del)
var.info$recall_del <- as.numeric(var.info$recall_del)
var.info$precision_del <- as.numeric(var.info$precision_del)
var.info$fall_out_del <- as.numeric(var.info$fall_out_del)
var.info$f1_score_del <- as.numeric(var.info$f1_score_del)
var.info$true_positive_all_del <- as.numeric(var.info$true_positive_all_del)
var.info$false_positive_all_del <- as.numeric(var.info$false_positive_all_del)
var.info$true_negative_all_del <- as.numeric(var.info$true_negative_all_del)
var.info$false_negative_all_del <- as.numeric(var.info$false_negative_all_del)
var.info$recall_all_del <- as.numeric(var.info$recall_all_del)
var.info$precision_all_del <- as.numeric(var.info$precision_all_del)
var.info$fall_out_all_del <- as.numeric(var.info$fall_out_all_del)
var.info$f1_score_all_del <- as.numeric(var.info$f1_score_all_del)
var.info$true_positive_del_alt_left <- as.numeric(var.info$true_positive_del_alt_left)
var.info$false_positive_del_alt_left <- as.numeric(var.info$false_positive_del_alt_left)
var.info$true_negative_del_alt_left <- as.numeric(var.info$true_negative_del_alt_left)
var.info$false_negative_del_alt_left <- as.numeric(var.info$false_negative_del_alt_left)
var.info$recall_del_alt_left <- as.numeric(var.info$recall_del_alt_left)
var.info$precision_del_alt_left <- as.numeric(var.info$precision_del_alt_left)
var.info$fall_out_del_alt_left <- as.numeric(var.info$fall_out_del_alt_left)
var.info$f1_score_del_alt_left <- as.numeric(var.info$f1_score_del_alt_left)
var.info$true_positive_del_ref_left <- as.numeric(var.info$true_positive_del_ref_left)
var.info$false_positive_del_ref_left <- as.numeric(var.info$false_positive_del_ref_left)
var.info$true_negative_del_ref_left <- as.numeric(var.info$true_negative_del_ref_left)
var.info$false_negative_del_ref_left <- as.numeric(var.info$false_negative_del_ref_left)
var.info$recall_del_ref_left <- as.numeric(var.info$recall_del_ref_left)
var.info$precision_del_ref_left <- as.numeric(var.info$precision_del_ref_left)
var.info$fall_out_del_ref_left <- as.numeric(var.info$fall_out_del_ref_left)
var.info$f1_score_del_ref_left <- as.numeric(var.info$f1_score_del_ref_left)
var.info[[new.col.names[1]]] <- as.factor(var.info[[new.col.names[1]]])
var.info <- var.info %>%
  mutate(prop_true_homo_ref_as_hetero_mu_in_called_pos = true_homo_ref_as_hetero_mu / (false_positive_hetero_mu + true_positive_hetero_mu)) %>%
  mutate(prop_true_homo_mu_as_hetero_mu_in_called_pos = true_homo_mu_as_hetero_mu / (false_positive_hetero_mu + true_positive_hetero_mu)) %>%
  mutate(prop_true_del_ref_left_as_hetero_mu_in_called_pos = true_del_ref_left_as_hetero_mu / (false_positive_hetero_mu + true_positive_hetero_mu)) %>%
  mutate(prop_true_del_alt_left_as_hetero_mu_in_called_pos = true_del_alt_left_as_hetero_mu / (false_positive_hetero_mu + true_positive_hetero_mu))

fsa.var.info <- var.info[var.info[[new.col.names[2]]] > isa.mu.rate.threshold, ]
var.info[[new.col.names[2]]] <- factor(var.info[[new.col.names[2]]], levels = sort(as.numeric(unique(var.info[[new.col.names[2]]]))))
var.info[[new.col.names[3]]] <- factor(var.info[[new.col.names[3]]])
var.info[[new.col.names[4]]] <- factor(var.info[[new.col.names[4]]])
var.info[[new.col.names[7]]] <- factor(var.info[[new.col.names[7]]], levels = names(cna.doublet.status))

fsa.var.info[[new.col.names[2]]] <- factor(fsa.var.info[[new.col.names[2]]], levels = sort(as.numeric(unique(fsa.var.info[[new.col.names[2]]]))))
fsa.var.info[[new.col.names[3]]] <- factor(fsa.var.info[[new.col.names[3]]])
fsa.var.info[[new.col.names[4]]] <- factor(fsa.var.info[[new.col.names[4]]])
fsa.var.info[[new.col.names[7]]] <- factor(fsa.var.info[[new.col.names[7]]], levels = names(cna.doublet.status))

ado.info$cell_num <- as.factor(ado.info$cell_num)
ado.info$coverage_mean <- as.factor(ado.info$coverage_mean)
ado.info$coverage_variance <- as.factor(ado.info$coverage_variance)
ado.info$dataset <- as.factor(ado.info$dataset)
ado.info$tool <- as.factor(ado.info$tool)
ado.info$snv_type <- as.factor(ado.info$snv_type)
ado.info$tool_setup <- as.factor(ado.info$tool_setup)
ado.info$doublet_rate <- as.factor(ado.info$doublet_rate)
ado.info$cnv_rate<- as.factor(ado.info$cnv_rate)
ado.info$data_type <- factor(ado.info$data_type, levels = names(cna.status.labels))
ado.info$true_positive <- as.numeric(ado.info$true_positive)
ado.info$false_positive <- as.numeric(ado.info$false_positive)
ado.info$true_negative <- as.numeric(ado.info$true_negative)
ado.info$false_negative <- as.numeric(ado.info$false_negative)
ado.info$recall <- as.numeric(ado.info$recall)
ado.info$precision <- as.numeric(ado.info$precision)
ado.info$fall_out <- as.numeric(ado.info$fall_out)
ado.info$f1_score <- as.numeric(ado.info$f1_score)
ado.info$true_positive_single_ado <- as.numeric(ado.info$true_positive_single_ado)
ado.info$false_positive_single_ado <- as.numeric(ado.info$false_positive_single_ado)
ado.info$true_negative_single_ado <- as.numeric(ado.info$true_negative_single_ado)
ado.info$false_negative_single_ado <- as.numeric(ado.info$false_negative_single_ado)
ado.info$recall_single_ado <- as.numeric(ado.info$recall_single_ado)
ado.info$precision_single_ado <- as.numeric(ado.info$precision_single_ado)
ado.info$fall_out_single_ado <- as.numeric(ado.info$fall_out_single_ado)
ado.info$f1_score_single_ado <- as.numeric(ado.info$f1_score_single_ado)
ado.info[[new.col.names[1]]] <- factor(ado.info[[new.col.names[1]]])
ado.info[[new.col.names[2]]] <- factor(ado.info[[new.col.names[2]]], levels = sort(as.numeric(unique(ado.info[[new.col.names[2]]]))))
ado.info[[new.col.names[3]]] <- factor(ado.info[[new.col.names[3]]])
ado.info[[new.col.names[4]]] <- factor(ado.info[[new.col.names[4]]])
ado.info[[new.col.names[7]]] <- factor(ado.info[[new.col.names[7]]], levels = names(cna.doublet.status))

# allelic.info$cell_num <- factor(allelic.info$cell_num)
# allelic.info$coverage_mean <- factor(allelic.info$coverage_mean)
# allelic.info$coverage_variance <- factor(allelic.info$coverage_variance)
# allelic.info$dataset <- factor(allelic.info$dataset)
# allelic.info$tool <- factor(allelic.info$tool)
# allelic.info$tool_setup <- factor(allelic.info$tool_setup)
# allelic.info$cell_names <- factor(allelic.info$cell_names)
# allelic.info[[new.col.names[1]]] <- factor(allelic.info[[new.col.names[1]]])
# allelic.info[[new.col.names[2]]] <- factor(allelic.info[[new.col.names[2]]], levels = sort(as.numeric(unique(allelic.info[[new.col.names[2]]]))))
# allelic.info[[new.col.names[3]]] <- factor(allelic.info[[new.col.names[3]]])
# allelic.info[[new.col.names[4]]] <- factor(allelic.info[[new.col.names[4]]])

efficiency$mean_s <- as.numeric(efficiency$mean_s)
efficiency$mean_m <- efficiency$mean_s / 60
efficiency$mean_mean_load <- as.numeric(efficiency$mean_mean_load)
efficiency$mean_cpu_time <- as.numeric(efficiency$mean_cpu_time)
efficiency$dataset <- as.factor(efficiency$dataset)
efficiency$tool_setup <- as.factor(efficiency$tool_setup)
efficiency$tool <- as.factor(efficiency$tool)
efficiency$stage <- as.factor(efficiency$stage)
efficiency$num_candidate_var_sites <- as.numeric(efficiency$num_candidate_var_sites)
efficiency$scenario <- factor(efficiency$scenario, levels = c(1, 2))

# Filter out some tools.
tree.info <- tree.info %>%
  filter(tool %in% c("SCIPhIN", "SCIPhI", "SiFit") | (grepl("SIEVE|DelSIEVE", tool) & grepl(".+_stage_2$", tool_setup))) # grepl("SIEVE", tool))

param.info <- param.info %>%
  filter(tool %in% c("SCIPhIN", "SCIPhI", "SiFit") | !grepl("-", tool_setup))

var.info <- var.info %>%
  filter(tool %in% c("SCIPhIN", "SCIPhI", "Monovar") | (grepl("SIEVE|DelSIEVE", tool) & grepl(".+_stage_2$", tool_setup))) # grepl("SIEVE", tool))

fsa.var.info <- fsa.var.info %>%
  filter(tool %in% c("SCIPhIN", "SCIPhI", "Monovar") | (grepl("SIEVE|DelSIEVE", tool) & grepl(".+_stage_2$", tool_setup))) # grepl("SIEVE", tool))

ado.info <- ado.info %>%
  filter(grepl(".+_stage_2$", tool_setup))

# efficiency <- efficiency %>%
#   filter(tool != "Monovar")

# process tool names

# rename legends of tools if necessary and replace 'tool' column
tree.info$tool <- rename.tools(tree.info[c("tool", "snv_type", "tool_setup")])
# set colors
tree.info$tool <- as.factor(tree.info$tool)
tree.info.prettified <- prettify.colors(levels(tree.info$tool), FALSE)
tree.info$tool <- factor(tree.info$tool, levels = tree.info.prettified[["tool"]])

# rename legends of tools if necessary and replace 'tool' column
site.info$tool <- rename.tools(site.info[c("tool", "snv_type", "tool_setup")])
# set colors
site.info$tool <- as.factor(site.info$tool)
site.info.prettified <- prettify.colors(levels(site.info$tool), FALSE)
site.info$tool <- factor(site.info$tool, levels = site.info.prettified[["tool"]])

# rename legends of tools if necessary and replace 'tool' column
param.info$tool <- rename.tools(param.info[c("tool", "snv_type")])
# set colors
param.info$tool <- as.factor(param.info$tool)
param.info.prettified <- prettify.colors(levels(param.info$tool), FALSE)
param.info$tool <- factor(param.info$tool, levels = param.info.prettified[["tool"]])

# rename legends of tools if necessary and replace 'tool' column
var.info$tool <- rename.tools(var.info[c("tool", "tool_setup")])
fsa.var.info$tool <- rename.tools(fsa.var.info[c("tool", "tool_setup")])
# set colors
var.info$tool <- as.factor(var.info$tool)
var.info.prettified <- prettify.colors(levels(var.info$tool), FALSE)
var.info$tool <- factor(var.info$tool, levels = var.info.prettified[["tool"]])
fsa.var.info$tool <- as.factor(fsa.var.info$tool)
fsa.var.info.prettified <- prettify.colors(levels(fsa.var.info$tool), FALSE)
fsa.var.info$tool <- factor(fsa.var.info$tool, levels = fsa.var.info.prettified[["tool"]])

homo.var.prop <- var.info %>%
  select(
    cell_num, simulated_data, mutation_rate, coverage_mean, coverage_variance, dataset, tool,
    true_positive_homo_mu, false_positive_homo_mu, true_negative_homo_mu, false_negative_homo_mu, !!as.name(new.col.names[4])
  ) %>%
  mutate(
    true_homo = true_positive_homo_mu + false_negative_homo_mu,
    genotype_sum = true_positive_homo_mu + false_positive_homo_mu + true_negative_homo_mu + false_negative_homo_mu
  ) %>%
  mutate(
    true_homo_prop = (true_positive_homo_mu + false_negative_homo_mu) / genotype_sum
  )

# homo.var <- lapply(
#   split(homo.var.prop, homo.var.prop$simulated_data),
#   function(x) {
#     .x <- x %>% filter(tool != "Monovar")
#     ret <- list()
#     ret$info <- .x %>%
#       distinct(cell_num, simulated_data, mutation_rate, coverage_mean, coverage_variance, !!as.name(new.col.names[4]))
#
#     t1 <- .x %>%
#       select(dataset, tool, genotype_sum) %>%
#       arrange(tool, dataset)
#     ret$genotype_sum <- matrix(
#       t1$genotype_sum,
#       nrow = length(unique(t1$dataset)),
#       dimnames = list(
#         dataset = distinct(t1, dataset)$dataset,
#         tool = distinct(t1, tool)$tool
#       )
#     )
#
#     t2 <- .x %>%
#       select(dataset, tool, true_homo_prop) %>%
#       arrange(tool, dataset)
#     ret$true_homo_prop <- matrix(
#       t2$true_homo_prop,
#       nrow = length(unique(t2$dataset)),
#       dimnames = list(
#         dataset = distinct(t2, dataset)$dataset,
#         tool = distinct(t2, tool)$tool
#       )
#     )
#
#     num_tool <- dim(ret$genotype_sum)[2L]
#     if (num_tool > 1L) {
#       ret$info$identical_genotype_sum <- all(
#         vapply(
#           2L:num_tool,
#           function(y) identical(ret$genotype_sum[, 1L], ret$genotype_sum[, y]),
#           FUN.VALUE = logical(1L),
#           USE.NAMES = FALSE
#         )
#       )
#
#       ret$info$identical_true_homo_prop <- all(
#         vapply(
#           2L:num_tool,
#           function(y) identical(ret$true_homo_prop[, 1L], ret$true_homo_prop[, y]),
#           FUN.VALUE = logical(1L),
#           USE.NAMES = FALSE
#         )
#       )
#     }
#
#     ret$info$min_genotype_sum <- min(ret$genotype_sum[, 1L])
#     ret$info$max_genotype_sum <- max(ret$genotype_sum[, 1L])
#
#     ret$info$min_true_homo_prop <- min(ret$true_homo_prop[, 1L])
#     ret$info$max_true_homo_prop <- max(ret$true_homo_prop[, 1L])
#
#     return(ret)
#   }
# )
#
# homo.var.info <- bind_rows(
#   lapply(
#     homo.var,
#     function(x) x$info
#   )
# ) %>%
#   arrange(mutation_rate) %>%
#   group_by(mutation_rate) %>%
#   mutate(
#     global_min_genotype_sum = min(min_genotype_sum),
#     global_max_genotype_sum = max(max_genotype_sum),
#     global_min_true_homo_prop = min(min_true_homo_prop),
#     global_max_true_homo_prop = max(max_true_homo_prop)
#   )

# parallel.mut.prop <- parallel.mut %>%
#   group_by(mutation_rate) %>%
#   mutate(
#     min_prop_parallel_mu_sites = min(prop_parallel_mu_sites),
#     max_prop_parallel_mu_sites = max(prop_parallel_mu_sites),
#     min_prop_parallel_mu_entries = min(prop_parallel_mu_entries),
#     max_prop_parallel_mu_entries = max(prop_parallel_mu_entries)
#   ) %>%
#   distinct(
#     mutation_rate,
#     min_prop_parallel_mu_sites,
#     max_prop_parallel_mu_sites,
#     min_prop_parallel_mu_entries,
#     max_prop_parallel_mu_entries
#   ) %>%
#   arrange(mutation_rate)

# rename legends of tools if necessary and replace 'tool' column
ado.info$tool <- rename.tools(ado.info[c("tool", "tool_setup")])
# set colors
ado.info$tool <- as.factor(ado.info$tool)
ado.info.prettified <- prettify.colors(levels(ado.info$tool), FALSE)
ado.info$tool <- factor(ado.info$tool, levels = ado.info.prettified[["tool"]])

# # rename legends of tools if necessary and replace 'tool' column
# allelic.info$tool <- rename.tools(allelic.info[c("tool", "tool_setup")])
# # set colors
# allelic.info$tool <- as.factor(allelic.info$tool)
# allelic.info.prettified <- prettify.colors(levels(allelic.info$tool), FALSE)
# allelic.info$tool <- factor(allelic.info$tool, levels = allelic.info.prettified[["tool"]])

# efficiency$tool <- as.factor(efficiency$tool)
# efficiency.prettified <- prettify.colors(levels(efficiency$tool_setup))
# efficiency$tool <- factor(efficiency$tool, levels = efficiency.prettified[["tool"]])

common.legend <- list(
  tool = c(tree.info.prettified[[1]][1:2], var.info.prettified[[1]]),
  color = c(tree.info.prettified[[2]][1:2], var.info.prettified[[2]]),
  fill = c(tree.info.prettified[[3]][1:2], var.info.prettified[[3]])
)


# data_quality_legend <- list(
#   quality = c(
#     "High mean\nLow variance",
#     "High mean\nMedium variance",
#     "Low mean\nHigh variance"
#   ),
#   color = common.legend$color[1L:3L],
#   fill = common.legend$fill[1L:3L]
# )
# names(data_quality_legend$color) <- c("2", "10", "20")
# names(data_quality_legend$fill) <- c("2", "10", "20")

# set no legends for single image
no.legend <- theme(legend.position = "hidden")


# Tree objs --------------------

# Tree: normalized RF distance
normalized.rf.dist.plot <-
  ggplot(
    tree.info,
    aes(
      x = data_type,
      y = normalized_RF_distance
    )
  ) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    y = "Normalised Robinson-Foulds distance",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# Tree: branch score difference
branch.score.diff.plot <-
  ggplot(
    tree.info[!is.na(tree.info$rooted_branch_score_difference), ],
    aes(
      x = data_type,
      y = rooted_branch_score_difference
    )
  ) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 0.5) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    y = "Branch score distance",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# # Parameter objs --------------------
#
# # Parameters: ado rate
# ado.rate.plot <- ggplot(
#   param.info,
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = ado_rate
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = true.ado.rate[1],
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "ADO rate",
#     x = "CNA sites in input data"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, "mm"),
#     strip.text.x = strip.text.x,
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # Parameters: sequencing error rate
# eff.seq.err.rate.plot <- ggplot(
#   param.info[!is.na(param.info$eff_seq_err_rate) & !param.info$tool %in% c("SCIPhI", "SCIPhIN"), ],
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = eff_seq_err_rate
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = true.eff.seq.err.rate,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Effective sequencing error rate",
#     x = "CNA sites in input data"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, "mm"),
#     strip.text.x = strip.text.x,
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # Parameters: wildtype overdispersion
# wildtype.overdispersion.plot <- ggplot(
#   param.info[!is.na(param.info$wild_overdispersion), ],
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = wild_overdispersion
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = true.wildtype.overdispersion,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Wildtype overdispersion",
#     x = "CNA sites in input data"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, "mm"),
#     strip.text.x = strip.text.x,
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # Parameters: alternative overdispersion
# alt.overdispersion.plot <- ggplot(
#   param.info[!is.na(param.info$alternative_overdispersion), ],
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = alternative_overdispersion
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = true.alt.overdispersion,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Alternative overdispersion",
#     x = "CNA sites in input data"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, "mm"),
#     strip.text.x = strip.text.x,
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # Parameters: mean of allelic coverage
# mean.allelic.cov.plot <- ggplot(
#   param.info[!is.na(param.info$allelic_seq_cov), ],
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = allelic_seq_cov
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = 20,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Mean of allelic coverage",
#     x = "CNA sites in input data"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, "mm"),
#     strip.text.x = strip.text.x,
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # Parameters: variance of allelic coverage
# var.allelic.cov.plot <- ggplot(
#   param.info[!is.na(param.info$allelic_seq_cov_raw_var), ],
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = allelic_seq_cov_raw_var
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = 10,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Variance of allelic coverage",
#     x = "CNA sites in input data"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, "mm"),
#     strip.text.x = strip.text.x,
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # Parameters: relative deletion rate
# rela.deletion.rate.plot <- ggplot(
#   param.info[!is.na(param.info$deletion_rate), ],
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = deletion_rate
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = 0.25,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Relative deletion rate",
#     x = "CNA sites in input data"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, "mm"),
#     strip.text.x = strip.text.x,
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # Mutation calling objs --------------------
#
# # Variant calling: mutation plot
# # recall
# recall.plot <- ggplot(
#   var.info,# %>% filter(prop_true_hetero_mu > prop.threshold),
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = recall
#   )
# ) +
#
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     # title = "Mutations",
#     y = "Recall",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     plot.title = plot.title,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.x = element_blank(),
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(1, "mm"),
#     legend.key.size = unit(4, "mm"),
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.spacing = unit(-2, "mm"),
#     strip.text.x = strip.text.x,
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # precision
# precision.plot <- ggplot(
#   var.info,# %>% filter(prop_true_hetero_mu > prop.threshold),
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = precision
#   )
# ) +
#
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     # title = "Mutations",
#     y = "Precision",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     plot.title = plot.title,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.x = element_blank(),
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(1, "mm"),
#     legend.key.size = unit(4, "mm"),
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.spacing = unit(-2, "mm"),
#     strip.text.x = element_blank(),
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # f1 score
# f1.score.plot <- ggplot(
#   var.info,# %>% filter(prop_true_hetero_mu > prop.threshold),
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = f1_score
#   )
# ) +
#
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     title = "All types of mutations",
#     y = "F1 score",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     plot.title = plot.title,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     axis.title.y = element_blank(),
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(1, "mm"),
#     legend.key.size = unit(4, "mm"),
#     legend.position = "bottom",
#     # legend.box = "vertical",
#     # legend.spacing = unit(-2, "mm"),
#     strip.text.x = element_blank(),
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   ) +
#   guides(
#     fill = guide_legend(order = 1),
#     col = guide_legend(order = 2)
#   )
#
# # fallout
# fallout.plot <- ggplot(
#   var.info,# %>% filter(prop_true_hetero_mu > prop.threshold),
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = fall_out
#   )
# ) +
#
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   geom_hline(
#     yintercept = 0,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     # title = "Mutations",
#     y = "False positive rate",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     plot.title = plot.title,
#     axis.text = axis.text,
#     axis.title = axis.title,
#     # axis.text.x = element_blank(),
#     # axis.ticks.x = element_blank(),
#     # axis.title.x = element_blank(),
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(1, "mm"),
#     legend.key.size = unit(4, "mm"),
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.spacing = unit(-2, "mm"),
#     strip.text.x = element_blank(),
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )

# Heterozygous mutation calling objs --------------------

# Variant calling: heterozygous mutation plot
# recall
hetero.recall.plot <- ggplot(
  var.info[!grepl("^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>% filter(prop_true_hetero_mu > prop.threshold),
  aes(
    x = data_type,
    y = recall_hetero_mu
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0.5, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Single mutants",
    y = "Recall",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# precision
hetero.precision.plot <- ggplot(
  var.info[!grepl("^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>% filter(prop_true_hetero_mu > prop.threshold),
  aes(
    x = data_type,
    y = precision_hetero_mu
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0.5, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Single mutants",
    y = "Precision",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.title.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, "mm"),
    strip.text.x = element_blank(),
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# f1 score
hetero.f1.score.plot <- ggplot(
  var.info[!grepl("^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>% filter(prop_true_hetero_mu > prop.threshold),
  aes(
    x = data_type,
    y = f1_score_hetero_mu
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0.5, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Single mutants",
    y = "F1 score",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.position = "bottom",
    strip.text.x = element_blank(),
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# fallout
hetero.fallout.plot <- ggplot(
  var.info[!grepl("^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>% filter(prop_true_hetero_mu > prop.threshold),
  aes(
    x = data_type,
    y = fall_out_hetero_mu
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 0.04) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Single mutants",
    y = "False positive rate",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(-2, "mm"),
    legend.key.size = unit(4, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# Homozygous mutation calling objs --------------------

# Variant calling: homozygous mutation plot
# recall
homo.recall.plot <- ggplot(
  fsa.var.info[!grepl("^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>% filter(prop_true_homo_mu > prop.threshold),
  # fsa.var.info[!is.na(fsa.var.info$recall_homo_mu) & fsa.var.info$tool != "SCIPhIN",],
  aes(
    x = data_type,
    y = recall_homo_mu
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0.5, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Double mutants",
    y = "Recall",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# precision
homo.precision.plot <- ggplot(
  fsa.var.info[!grepl("^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>% filter(prop_true_homo_mu > prop.threshold),
  # fsa.var.info[!is.na(fsa.var.info$precision_homo_mu) & fsa.var.info$tool != "SCIPhI",],
  aes(
    x = data_type,
    y = precision_homo_mu
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0.5, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Double mutants",
    y = "Precision",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, "mm"),
    strip.text.x = element_blank(),
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# f1 score
homo.f1.score.plot <- ggplot(
  fsa.var.info[!grepl("^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>% filter(prop_true_homo_mu > prop.threshold),
  # fsa.var.info[!is.na(fsa.var.info$f1_score_homo_mu) & fsa.var.info$tool != "SCIPhI",],
  aes(
    x = data_type,
    y = f1_score_homo_mu
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0.5, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Double mutants",
    y = "F1 score",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "none",
    legend.key.size = unit(4, "mm"),
    strip.text.x = element_blank(),
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# fallout
homo.fallout.plot <- ggplot(
  fsa.var.info[!grepl("^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>% filter(prop_true_homo_mu > prop.threshold),
  # fsa.var.info[!is.na(fsa.var.info$fall_out_homo_mu) & fsa.var.info$tool != "SCIPhI",],
  aes(
    x = data_type,
    y = fall_out_homo_mu
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 0.002) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Double mutants",
    y = "False positive rate",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(-2, "mm"),
    legend.key.size = unit(4, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# Source analysis -----------------------
# Source analysis of false positives for heterozygous mutations: originally homozygous reference
# true.homo.ref.as.hetero.mu.in.called.pos.plot <- ggplot(
#   var.info,
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = prop_true_homo_ref_as_hetero_mu_in_called_pos
#   )
# ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       fill = tool,
#       color = tool
#     )
#   ) +
#   scale_fill_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]],
#     values = var.info.prettified[[3]]
#   ) +
#   scale_color_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]],
#     values = var.info.prettified[[2]]
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes(color = tool),
#     size = dot.size,
#     alpha = 0.5
#   ) +
#   facet_grid(
#     ~coverage_variance ~ ~cell_num,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var,
#       cell_num = cell.num.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   geom_hline(
#     yintercept = 0,
#     linetype="dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Proportion of true wildtype genotype in predicted single mutant genotype",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = element_text(size = 7),
#     axis.text = element_text(size = 7),
#     legend.text = element_text(size = 6),
#     legend.title = element_text(size = 6),
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, 'mm'),
# strip.text.x = strip.text.x,
#     strip.text.y = element_blank(),
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # Source analysis of false positives for heterozygous mutations: originally homozygous mutation
# true.homo.mu.as.hetero.mu.in.called.pos.plot <- ggplot(
#   var.info,
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = prop_true_homo_mu_as_hetero_mu_in_called_pos
#   )
# ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       fill = tool,
#       color = tool
#     )
#   ) +
#   scale_fill_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]],
#     values = var.info.prettified[[3]]
#   ) +
#   scale_color_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]],
#     values = var.info.prettified[[2]]
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes(color = tool),
#     size = dot.size,
#     alpha = 0.5
#   ) +
#   facet_grid(
#     ~coverage_variance ~ ~cell_num,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var,
#       cell_num = cell.num.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   geom_hline(
#     yintercept = 0,
#     linetype="dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Proportion of true double mutant genotypes in predicted single mutant genotype",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = element_text(size = 7),
#     axis.text = element_text(size = 7),
#     legend.text = element_text(size = 6),
#     legend.title = element_text(size = 6),
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, 'mm'),
# strip.text.x = strip.text.x,
# strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )

# Variant calling: deletion plot
# recall
# del.recall.plot <- ggplot(
#   var.info[which(!is.na(var.info$recall_del) & var.info$tool != "SIEVE"),],
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = recall_del
#   )
# ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       fill = tool,
#       color = tool
#     )
#   ) +
#   scale_fill_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]], #common.legend[[1]],
#     values = var.info.prettified[[3]], #common.legend[[3]]
#   ) +
#   scale_color_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]], #common.legend[[1]],
#     values = var.info.prettified[[2]], #common.legend[[2]]
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes(color = tool),
#     size = dot.size,
#     alpha = 0.5
#   ) +
#   facet_grid(
#     ~coverage_variance ~ ~rela_deletion_rate,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var,
#       rela_deletion_rate = rela.del.rate.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   # ylim(0.5, 1.0) +
#   geom_hline(
#     yintercept = 1,
#     linetype="dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Recall",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.x = element_blank(),
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, 'mm'),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank(),
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # precision
# del.precision.plot <- ggplot(
#   var.info[which(!is.na(var.info$precision_del) & var.info$tool != "SIEVE"),],
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = precision_del
#   )
# ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       fill = tool,
#       color = tool
#     )
#   ) +
#   scale_fill_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]], #common.legend[[1]],
#     values = var.info.prettified[[3]], #common.legend[[3]]
#   ) +
#   scale_color_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]], #common.legend[[1]],
#     values = var.info.prettified[[2]], #common.legend[[2]]
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes(color = tool),
#     size = dot.size,
#     alpha = 0.5
#   ) +
#   facet_grid(
#     ~coverage_variance ~ ~rela_deletion_rate,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var,
#       rela_deletion_rate = rela.del.rate.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   # ylim(0.55, 1.0) +
#   geom_hline(
#     yintercept = 1,
#     linetype="dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Precision",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, 'mm'),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank(),
#     strip.background = element_blank(),
#     panel.grid = element_blank(),
#     plot.margin = unit(c(0, 0, 0, 0),"mm")
#   )
#
# # f1 score
# del.f1.score.plot <- ggplot(
#   var.info[which(!is.na(var.info$f1_score_del) & var.info$tool != "SIEVE"),],
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = f1_score_del
#   )
# ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       fill = tool,
#       color = tool
#     )
#   ) +
#   scale_fill_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]], #common.legend[[1]],
#     values = var.info.prettified[[3]], #common.legend[[3]]
#   ) +
#   scale_color_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]], #common.legend[[1]],
#     values = var.info.prettified[[2]], #common.legend[[2]]
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes(color = tool),
#     size = dot.size,
#     alpha = 0.5
#   ) +
#   facet_grid(
#     ~coverage_variance ~ ~rela_deletion_rate,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var,
#       rela_deletion_rate = rela.del.rate.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   # ylim(0.65, 1.0) +
#   geom_hline(
#     yintercept = 1,
#     linetype="dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "F1 score",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.key.size = unit(4, 'mm'),
#     legend.position = "bottom",
#     strip.text.x = element_blank(),
# strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # fallout
# del.fallout.plot <- ggplot(
#   var.info[which(!is.na(var.info$fall_out_del) & var.info$tool != "SIEVE"),],
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = fall_out_del
#   )
# ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       fill = tool,
#       color = tool
#     )
#   ) +
#   scale_fill_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]],
#     values = var.info.prettified[[3]]
#   ) +
#   scale_color_manual(
#     name = "Method",
#     breaks = var.info.prettified[[1]],
#     values = var.info.prettified[[2]]
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes(color = tool),
#     size = dot.size,
#     alpha = 0.5
#   ) +
#   facet_grid(
#     ~coverage_variance ~ ~rela_deletion_rate,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var,
#       rela_deletion_rate = rela.del.rate.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   # ylim(0, 0.15) +
#   geom_hline(
#     yintercept = 0,
#     linetype="dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "False positive rate",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     axis.text = axis.text,
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.x = element_blank(),
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, 'mm'),
#     strip.text.x = element_blank(),
# strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )

# Alternative left deletion calling objs --------------------

# Variant calling: single deletion plot for 1/-
# recall
del.alt.left.recall.plot <- ggplot(
  fsa.var.info[!grepl("^sieve|^monovar$|^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>%
    filter(prop_true_del_alt_left > prop.threshold),
  aes(
    x = data_type,
    y = recall_del_alt_left
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Alternative-remaining single deletion",
    y = "Recall",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(1, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.position = "bottom",
    legend.box = "horizontal",
    # legend.spacing = unit(-2, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# precision
del.alt.left.precision.plot <- ggplot(
  fsa.var.info[!grepl("^sieve|^monovar$|^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>%
    filter(prop_true_del_alt_left > prop.threshold),
  aes(
    x = data_type,
    y = precision_del_alt_left
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Alternative-remaining single deletion",
    y = "Precision",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(1, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.position = "bottom",
    legend.box = "horizontal",
    # legend.spacing = unit(-2, "mm"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# f1 score
del.alt.left.f1.score.plot <- ggplot(
  fsa.var.info[!grepl("^sieve|^monovar$|^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>%
    filter(prop_true_del_alt_left > prop.threshold),
  aes(
    x = data_type,
    y = f1_score_del_alt_left
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Alternative-remaining single deletion",
    y = "F1 score",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(1, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(-2, "mm"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# fallout
del.alt.left.fallout.plot <- ggplot(
  fsa.var.info[!grepl("^sieve|^monovar$|^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>%
    filter(prop_true_del_alt_left > prop.threshold),
  aes(
    x = data_type,
    y = fall_out_del_alt_left
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  # ylim(0, 0.02) +
  scale_y_continuous(
    breaks = seq(0, 0.02, by = 0.01),
    limits = c(0, 0.02)#,
    # labels = scientific
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Alternative-remaining single deletion",
    y = "False positive rate",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(1, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(-2, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# Reference left deletion calling objs --------------------

# Variant calling: single deletion plot for 0/-
# recall
del.ref.left.recall.plot <- ggplot(
  fsa.var.info[!grepl("^sieve|^monovar$|^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>%
    filter(prop_true_del_ref_left > prop.threshold),
  aes(
    x = data_type,
    y = recall_del_ref_left
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Reference-remaining single deletion",
    y = "Recall",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.position = "bottom",
    legend.box = "horizontal",
    # legend.spacing = unit(-2, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# precision
del.ref.left.precision.plot <- ggplot(
  fsa.var.info[!grepl("^sieve|^monovar$|^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>%
    filter(prop_true_del_ref_left > prop.threshold),
  aes(
    x = data_type,
    y = precision_del_ref_left
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Reference-remaining single deletion",
    y = "Precision",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.position = "bottom",
    legend.box = "horizontal",
    # legend.spacing = unit(-2, "mm"),
    strip.text.x = element_blank(),
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# f1 score
del.ref.left.f1.score.plot <- ggplot(
  fsa.var.info[!grepl("^sieve|^monovar$|^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>%
    filter(prop_true_del_ref_left > prop.threshold),
  aes(
    x = data_type,
    y = f1_score_del_ref_left
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Reference-remaining single deletion",
    y = "F1 score",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(-2, "mm"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# fallout
del.ref.left.fallout.plot <- ggplot(
  fsa.var.info[!grepl("^sieve|^monovar$|^sciphi$|^sciphin$", fsa.var.info$tool, ignore.case = TRUE), ] %>%
    filter(prop_true_del_ref_left > prop.threshold),
  aes(
    x = data_type,
    y = fall_out_del_ref_left
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 0.3) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "Reference-remaining single deletion",
    y = "False positive rate",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(-2, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# Double deletion calling objs --------------------

# Variant calling: double deletion -
# recall
# del.all.recall.plot <- ggplot(
#   fsa.var.info[grepl("^delsieve", fsa.var.info$tool, ignore.case = TRUE), ] %>%
#     filter(as.numeric(as.character(!!as.name(new.col.names[2]))) > 1e-5) %>%
#     filter(prop_true_all_del > prop.threshold),
#   aes(
#     x = !!as.name(new.col.names[5]),
#     y = recall_all_del
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.7,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = doublet_rate,
#       fill = doublet_rate
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = doublet_rate#,
#       # fill = doublet_rate
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   scale_fill_manual(
#     name = "Doublet rate",
#     values = doublet.rate.colors
#   ) +
#   facet_wrap(
#     as.name(new.col.names[6]),
#     labeller = labeller(
#       .cols = cna.sites.prop.labels
#     ),
#     scales = "free_x"
#   ) +
#   scale_x_discrete(
#     breaks = levels(tree.info[[new.col.names[5]]]),
#     label = cna.status.labels
#   ) +
#   # ylim(0.5, 1.0) +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     title = "Double deletions",
#     y = "Recall",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     plot.title = plot.title,
#     axis.text = axis.text,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.key.size = unit(4, "mm"),
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     # legend.spacing = unit(-2, "mm"),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank(),
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # precision
# del.all.precision.plot <- ggplot(
#   fsa.var.info[grepl("^delsieve", fsa.var.info$tool, ignore.case = TRUE), ] %>%
#     filter(as.numeric(as.character(!!as.name(new.col.names[2]))) > 1e-5) %>%
#     filter(prop_true_all_del > prop.threshold),
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = precision_all_del
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = tool,
#       fill = !!as.name(new.col.names[4])
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = tool,
#       fill = !!as.name(new.col.names[4])
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Method",
#     breaks = tree.info.prettified[[1]], # common.legend[[1]],
#     values = tree.info.prettified[[2]], # common.legend[[2]]
#   ) +
#   scale_fill_manual(
#     name = "Simulated dropout type",
#     breaks = c("ADO", "LDO"),
#     values = c("white", bg.color)
#   ) +
#   # scale_fill_manual(
#   #   name = "Method",
#   #   breaks = tree.info.prettified[[1]], #common.legend[[1]],
#   #   values = tree.info.prettified[[3]], #common.legend[[3]]
#   # ) +
#   facet_grid(
#     ~coverage_variance ~ ~rela_deletion_rate,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var,
#       rela_deletion_rate = rela.del.rate.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
#     label = scientific,
#     drop = TRUE #FALSE
#   ) +
#   # ylim(0.55, 1.0) +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     title = "Double deletions",
#     y = "Precision",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     plot.title = plot.title,
#     axis.text = axis.text,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.key.size = unit(4, "mm"),
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     # legend.spacing = unit(-2, "mm"),
#     strip.text.x = element_blank(),
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   ) +
#   guides(
#     fill = guide_legend(order = 1),
#     col = guide_legend(order = 2)
#   )
#
# # f1 score
# del.all.f1.score.plot <- ggplot(
#   fsa.var.info[grepl("^delsieve", fsa.var.info$tool, ignore.case = TRUE), ] %>%
#     filter(as.numeric(as.character(!!as.name(new.col.names[2]))) > 1e-5) %>%
#     filter(prop_true_all_del > prop.threshold),
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = f1_score_all_del
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = tool,
#       fill = !!as.name(new.col.names[4])
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = tool,
#       fill = !!as.name(new.col.names[4])
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Method",
#     breaks = tree.info.prettified[[1]], # common.legend[[1]],
#     values = tree.info.prettified[[2]], # common.legend[[2]]
#   ) +
#   scale_fill_manual(
#     name = "Simulated dropout type",
#     breaks = c("ADO", "LDO"),
#     values = c("white", bg.color)
#   ) +
#   # scale_fill_manual(
#   #   name = "Method",
#   #   breaks = tree.info.prettified[[1]], #common.legend[[1]],
#   #   values = tree.info.prettified[[3]], #common.legend[[3]]
#   # ) +
#   facet_grid(
#     ~coverage_variance ~ ~rela_deletion_rate,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var,
#       rela_deletion_rate = rela.del.rate.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
#     label = scientific,
#     drop = TRUE #FALSE
#   ) +
#   ylim(0, 1.0) +
#   geom_hline(
#     yintercept = 1,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     title = "Double deletions",
#     y = "F1 score",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     plot.title = plot.title,
#     axis.text = axis.text,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.key.size = unit(4, "mm"),
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.spacing = unit(-2, "mm"),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank(),
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   ) +
#   guides(
#     fill = guide_legend(order = 1),
#     col = guide_legend(order = 2)
#   )
#
# # fallout
# del.all.fallout.plot <- ggplot(
#   fsa.var.info[grepl("^delsieve", fsa.var.info$tool, ignore.case = TRUE), ] %>%
#     filter(as.numeric(as.character(!!as.name(new.col.names[2]))) > 1e-5) %>%
#     filter(prop_true_all_del > prop.threshold),
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = fall_out_all_del
#   )
# ) +
#   box.func(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 1.0,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       color = tool,
#       fill = !!as.name(new.col.names[4])
#     ),
#     position = position_dodge(width = dodge.width)
#   ) +
#   jitter.func(
#     position = position_jitterdodge(
#       jitter.width = jitter.width,
#       dodge.width = dodge.width
#     ),
#     aes(
#       color = tool,
#       fill = !!as.name(new.col.names[4])
#     ),
#     size = dot.size,
#     alpha = jitter.alpha
#   ) +
#   scale_color_manual(
#     name = "Method",
#     breaks = tree.info.prettified[[1]], # common.legend[[1]],
#     values = tree.info.prettified[[2]], # common.legend[[2]]
#   ) +
#   scale_fill_manual(
#     name = "Simulated dropout type",
#     breaks = c("ADO", "LDO"),
#     values = c("white", bg.color)
#   ) +
#   # scale_fill_manual(
#   #   name = "Method",
#   #   breaks = tree.info.prettified[[1]], #common.legend[[1]],
#   #   values = tree.info.prettified[[3]], #common.legend[[3]]
#   # ) +
#   facet_grid(
#     ~coverage_variance ~ ~rela_deletion_rate,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var,
#       rela_deletion_rate = rela.del.rate.labels
#     )
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(var.info[[new.col.names[2]]])),
#     label = scientific,
#     drop = TRUE #FALSE
#   ) +
#   ylim(0, 0.015) +
#   geom_hline(
#     yintercept = 0,
#     linetype = "dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     title = "Double deletions",
#     y = "False positive rate",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = text,
#     plot.title = plot.title,
#     axis.text = axis.text,
#     legend.text = legend.text,
#     legend.title = legend.title,
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.key.size = unit(4, "mm"),
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.spacing = unit(-2, "mm"),
#     strip.text.x = element_blank(),
#     strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   ) +
#   guides(
#     fill = guide_legend(order = 1),
#     col = guide_legend(order = 2)
#   )

# Single ADO calling objs --------------------

# Single ADO calling: mutation plot
# recall
single.ado.recall.plot <- ggplot(
  ado.info[!is.na(ado.info$recall_single_ado), ],
  aes(
    x = data_type,
    y = recall_single_ado
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "ADO",
    y = "Recall",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(-2, "mm"),
    legend.key.size = unit(4, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# precision
single.ado.precision.plot <- ggplot(
  ado.info[!is.na(ado.info$precision_single_ado), ],
  aes(
    x = data_type,
    y = precision_single_ado
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "ADO",
    y = "Precision",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.title.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(-2, "mm"),
    legend.key.size = unit(4, "mm"),
    strip.text.x = element_blank(),
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# f1 score
single.ado.f1.score.plot <- ggplot(
  ado.info[!is.na(ado.info$f1_score_single_ado), ],
  aes(
    x = data_type,
    y = f1_score_single_ado
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 1) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "ADO",
    y = "F1 score",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.x = axis.text.x,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(-2, "mm"),
    legend.key.size = unit(4, "mm"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# fallout
single.ado.fallout.plot <- ggplot(
  ado.info[!is.na(ado.info$fall_out_single_ado), ],
  aes(
    x = data_type,
    y = fall_out_single_ado
  )
) +
  box.func(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = cnv_rate,
      fill = doublet_rate
    ),
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      fill = doublet_rate,
      shape = cnv_rate
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_color_manual(
    name = "CNA rate",
    values = cna.rate.colors
  ) +
  scale_fill_manual(
    name = "Doublet rate",
    values = doublet.rate.colors
  ) +
  scale_shape_manual(
    guide = "none",
    values = c(16, 16, 16)
  ) +
  facet_wrap(
    as.name(new.col.names[7]),
    labeller = labeller(
      .cols = cna.doublet.status
    ),
    scales = "free_x"
  ) +
  scale_x_discrete(
    breaks = levels(tree.info$data_type),
    label = cna.status.labels
  ) +
  ylim(0, 0.1) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray43",
    linewidth = 1 / .pt
  ) +
  labs(
    # title = "ADO",
    y = "False positive rate",
    x = "CNA sites in input data"
  ) +
  theme_bw() +
  theme(
    text = text,
    plot.title = plot.title,
    axis.text.y = axis.text.y,
    axis.title = axis.title,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = legend.text,
    legend.title = legend.title,
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing = unit(-2, "mm"),
    legend.key.size = unit(4, "mm"),
    strip.text.x = strip.text.x,
    strip.text.y = strip.text.y,
    strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    col = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# Locus ADO calling objs --------------------
# Locus ADO calling: mutation plot
# if ("locus_dropout" %in% ado.modes) {
#   # recall
#   locus.ado.recall.plot <- ggplot(
#     ado.info[grepl("LDO", ado.info$simulated_dropout_type, ignore.case = TRUE) & !is.na(ado.info$recall_locus_ado), ],
#     aes(
#       x = !!as.name(new.col.names[2]),
#       y = recall_locus_ado
#     )
#   ) +
#     box.func(
#       lwd = 1 / .pt,
#       fatten = 1.2 / .pt,
#       alpha = 1.0,
#       outlier.alpha = 1.0,
#       outlier.size = dot.size,
#       outlier.shape = NA,
#       aes(
#         color = tool,
#         fill = !!as.name(new.col.names[4])
#       ),
#       position = position_dodge(width = dodge.width)
#     ) +
#     scale_fill_manual(
#       name = "Simulated dropout type",
#       breaks = c("ADO", "LDO"),
#       values = c("white", bg.color)
#     ) +
#     scale_color_manual(
#       name = "Method",
#       breaks = ado.info.prettified[[1]],
#       values = ado.info.prettified[[2]]
#     ) +
#     jitter.func(
#       position = position_jitterdodge(
#         jitter.width = jitter.width,
#         dodge.width = dodge.width
#       ),
#       aes(
#         color = tool,
#         fill = !!as.name(new.col.names[4])
#       ),
#       size = dot.size,
#       alpha = jitter.alpha
#     ) +
#     facet_grid(
#       ~coverage_variance ~ ~rela_deletion_rate,
#       labeller = labeller(
#         coverage_variance = allelic.raw.var,
#         rela_deletion_rate = rela.del.rate.labels
#       )
#     ) +
#     scale_x_discrete(
#       breaks = as.numeric(levels(ado.info[[new.col.names[2]]])),
#       label = scientific
#     ) +
#     ylim(0, 1) +
#     geom_hline(
#       yintercept = 1,
#       linetype = "dashed",
#       color = "gray43",
#       linewidth = 1 / .pt
#     ) +
#     labs(
#       title = "LDO",
#       y = "Recall",
#       x = "Mutation rate"
#     ) +
#     theme_bw() +
#     theme(
#       text = text,
#       plot.title = plot.title,
#       axis.text = axis.text,
#       legend.text = legend.text,
#       legend.title = legend.title,
#       legend.title.align = 0.5,
#       legend.box.spacing = unit(0.3, "mm"),
#       legend.position = "bottom",
#       legend.box = "vertical",
#       legend.spacing = unit(-2, "mm"),
#       legend.key.size = unit(4, "mm"),
#       strip.text.x = element_blank(),
#       strip.text.y = element_blank(),
#       strip.background = element_blank(),
#       panel.grid = element_blank()
#     ) +
#     guides(
#       fill = guide_legend(order = 1),
#       col = guide_legend(order = 2)
#     )
#
#   # precision
#   locus.ado.precision.plot <- ggplot(
#     ado.info[grepl("LDO", ado.info$simulated_dropout_type, ignore.case = TRUE) & !is.na(ado.info$precision_locus_ado), ],
#     aes(
#       x = !!as.name(new.col.names[2]),
#       y = precision_locus_ado
#     )
#   ) +
#     box.func(
#       lwd = 1 / .pt,
#       fatten = 1.2 / .pt,
#       alpha = 1.0,
#       outlier.alpha = 1.0,
#       outlier.size = dot.size,
#       outlier.shape = NA,
#       aes(
#         color = tool,
#         fill = !!as.name(new.col.names[4])
#       ),
#       position = position_dodge(width = dodge.width)
#     ) +
#     scale_fill_manual(
#       name = "Simulated dropout type",
#       breaks = c("ADO", "LDO"),
#       values = c("white", bg.color)
#     ) +
#     scale_color_manual(
#       name = "Method",
#       breaks = ado.info.prettified[[1]],
#       values = ado.info.prettified[[2]]
#     ) +
#     jitter.func(
#       position = position_jitterdodge(
#         jitter.width = jitter.width,
#         dodge.width = dodge.width
#       ),
#       aes(
#         color = tool,
#         fill = !!as.name(new.col.names[4])
#       ),
#       size = dot.size,
#       alpha = jitter.alpha
#     ) +
#     facet_grid(
#       ~coverage_variance ~ ~rela_deletion_rate,
#       labeller = labeller(
#         coverage_variance = allelic.raw.var,
#         rela_deletion_rate = rela.del.rate.labels
#       ),
#       drop = FALSE
#     ) +
#     scale_x_discrete(
#       breaks = as.numeric(levels(ado.info[[new.col.names[2]]])),
#       label = scientific
#     ) +
#     ylim(0, 1) +
#     geom_hline(
#       yintercept = 1,
#       linetype = "dashed",
#       color = "gray43",
#       linewidth = 1 / .pt
#     ) +
#     labs(
#       title = "LDO",
#       y = "Precision",
#       x = "Mutation rate"
#     ) +
#     theme_bw() +
#     theme(
#       text = text,
#       plot.title = plot.title,
#       axis.text = axis.text,
#       legend.text = legend.text,
#       legend.title = legend.title,
#       legend.title.align = 0.5,
#       legend.box.spacing = unit(0.3, "mm"),
#       legend.position = "bottom",
#       legend.box = "vertical",
#       legend.spacing = unit(-2, "mm"),
#       legend.key.size = unit(4, "mm"),
#       strip.text.x = element_blank(),
#       strip.text.y = strip.text.y,
#       strip.background = element_blank(),
#       panel.grid = element_blank()
#     ) +
#     guides(
#       fill = guide_legend(order = 1),
#       col = guide_legend(order = 2)
#     )
#
#   # f1 score
#   locus.ado.f1.score.plot <- ggplot(
#     ado.info[grepl("LDO", ado.info$simulated_dropout_type, ignore.case = TRUE) & !is.na(ado.info$f1_score_locus_ado), ],
#     aes(
#       x = !!as.name(new.col.names[2]),
#       y = f1_score_locus_ado
#     )
#   ) +
#     box.func(
#       lwd = 1 / .pt,
#       fatten = 1.2 / .pt,
#       alpha = 1.0,
#       outlier.alpha = 1.0,
#       outlier.size = dot.size,
#       outlier.shape = NA,
#       aes(
#         color = tool,
#         fill = !!as.name(new.col.names[4])
#       ),
#       position = position_dodge(width = dodge.width)
#     ) +
#     scale_fill_manual(
#       name = "Simulated dropout type",
#       breaks = c("ADO", "LDO"),
#       values = c("white", bg.color)
#     ) +
#     scale_color_manual(
#       name = "Method",
#       breaks = ado.info.prettified[[1]],
#       values = ado.info.prettified[[2]]
#     ) +
#     jitter.func(
#       position = position_jitterdodge(
#         jitter.width = jitter.width,
#         dodge.width = dodge.width
#       ),
#       aes(
#         color = tool,
#         fill = !!as.name(new.col.names[4])
#       ),
#       size = dot.size,
#       alpha = jitter.alpha
#     ) +
#     facet_grid(
#       ~coverage_variance ~ ~rela_deletion_rate,
#       labeller = labeller(
#         coverage_variance = allelic.raw.var,
#         rela_deletion_rate = rela.del.rate.labels
#       ),
#       drop = FALSE
#     ) +
#     scale_x_discrete(
#       breaks = as.numeric(levels(ado.info[[new.col.names[2]]])),
#       label = scientific
#     ) +
#     ylim(0, 1) +
#     geom_hline(
#       yintercept = 1,
#       linetype = "dashed",
#       color = "gray43",
#       linewidth = 1 / .pt
#     ) +
#     labs(
#       title = "LDO",
#       y = "F1 score",
#       x = "Mutation rate"
#     ) +
#     theme_bw() +
#     theme(
#       text = text,
#       plot.title = plot.title,
#       axis.text = axis.text,
#       legend.text = legend.text,
#       legend.title = legend.title,
#       axis.text.y = element_blank(),
#       axis.title.y = element_blank(),
#       axis.ticks.y = element_blank(),
#       legend.title.align = 0.5,
#       legend.box.spacing = unit(0.3, "mm"),
#       legend.position = "bottom",
#       legend.box = "vertical",
#       legend.spacing = unit(-2, "mm"),
#       legend.key.size = unit(4, "mm"),
#       strip.text.x = strip.text.x,
#       strip.text.y = strip.text.y,
#       strip.background = element_blank(),
#       panel.grid = element_blank()
#     ) +
#     guides(
#       fill = guide_legend(order = 1),
#       col = guide_legend(order = 2)
#     )
#
#   # fallout
#   locus.ado.fallout.plot <- ggplot(
#     ado.info[grepl("LDO", ado.info$simulated_dropout_type, ignore.case = TRUE) & !is.na(ado.info$fall_out_locus_ado), ],
#     aes(
#       x = !!as.name(new.col.names[2]),
#       y = fall_out_locus_ado
#     )
#   ) +
#     box.func(
#       lwd = 1 / .pt,
#       fatten = 1.2 / .pt,
#       alpha = 1.0,
#       outlier.alpha = 1.0,
#       outlier.size = dot.size,
#       outlier.shape = NA,
#       aes(
#         color = tool,
#         fill = !!as.name(new.col.names[4])
#       ),
#       position = position_dodge(width = dodge.width)
#     ) +
#     scale_fill_manual(
#       name = "Simulated dropout type",
#       breaks = c("ADO", "LDO"),
#       values = c("white", bg.color)
#     ) +
#     scale_color_manual(
#       name = "Method",
#       breaks = ado.info.prettified[[1]],
#       values = ado.info.prettified[[2]]
#     ) +
#     jitter.func(
#       position = position_jitterdodge(
#         jitter.width = jitter.width,
#         dodge.width = dodge.width
#       ),
#       aes(
#         color = tool,
#         fill = !!as.name(new.col.names[4])
#       ),
#       size = dot.size,
#       alpha = jitter.alpha
#     ) +
#     facet_grid(
#       ~coverage_variance ~ ~rela_deletion_rate,
#       labeller = labeller(
#         coverage_variance = allelic.raw.var,
#         rela_deletion_rate = rela.del.rate.labels
#       )
#     ) +
#     scale_x_discrete(
#       breaks = as.numeric(levels(ado.info[[new.col.names[2]]])),
#       label = scientific
#     ) +
#     ylim(0, 0.2) +
#     geom_hline(
#       yintercept = 0,
#       linetype = "dashed",
#       color = "gray43",
#       linewidth = 1 / .pt
#     ) +
#     labs(
#       title = "LDO",
#       y = "False positive rate",
#       x = "Mutation rate"
#     ) +
#     theme_bw() +
#     theme(
#       text = text,
#       plot.title = plot.title,
#       axis.text = axis.text,
#       legend.text = legend.text,
#       legend.title = legend.title,
#       # axis.text.x = element_blank(),
#       # axis.ticks.x = element_blank(),
#       # axis.title.x = element_blank(),
#       legend.title.align = 0.5,
#       legend.box.spacing = unit(0.3, "mm"),
#       legend.position = "bottom",
#       legend.box = "vertical",
#       legend.spacing = unit(-2, "mm"),
#       legend.key.size = unit(4, "mm"),
#       strip.text.x = element_blank(),
#       strip.text.y = strip.text.y,
#       strip.background = element_blank(),
#       panel.grid = element_blank()
#     ) +
#     guides(
#       fill = guide_legend(order = 1),
#       col = guide_legend(order = 2)
#     )
# }

# Candidate variant sites calling ----------------

# recall
# candidate.recall.plot <- ggplot(
#   site.info,
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = recall
#   )
# ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.5,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       fill = cell_num,
#       color = cell_num
#     )
#   ) +
#   scale_fill_brewer(
#     name = "Number of cells",
#     palette = "Set2"
#   ) +
#   scale_color_brewer(
#     name = "Number of cells",
#     palette = "Set2"
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes(color = cell_num),
#     size = dot.size,
#     alpha = 1
#   ) +
#   facet_wrap(
#     ~coverage_variance,
#     ncol = 1,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var
#     ),
#     strip.position = "right"
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(site.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   ylim(0, 1.0) +
#   geom_hline(
#     yintercept = 1,
#     linetype="dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Recall",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = element_text(size = 7),
#     axis.text = element_text(size = 7),
#     legend.text = element_text(size = 6),
#     legend.title = element_text(size = 6),
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, 'mm'),
# strip.text.x = strip.text.x,
#     strip.text.y = element_blank(),
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # precision
# candidate.precision.plot <- ggplot(
#   site.info,
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = precision
#   )
# ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.5,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       fill = cell_num,
#       color = cell_num
#     )
#   ) +
#   scale_fill_brewer(
#     name = "Number of cells",
#     palette = "Set2"
#   ) +
#   scale_color_brewer(
#     name = "Number of cells",
#     palette = "Set2"
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes(color = cell_num),
#     size = dot.size,
#     alpha = 1
#   ) +
#   facet_wrap(
#     ~coverage_variance,
#     ncol = 1,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var
#     ),
#     strip.position = "right"
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(site.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   ylim(0, 1.0) +
#   geom_hline(
#     yintercept = 1,
#     linetype="dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "Precision",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = element_text(size = 7),
#     axis.text = element_text(size = 7),
#     legend.text = element_text(size = 6),
#     legend.title = element_text(size = 6),
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, 'mm'),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank(),
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # f1 score
# candidate.f1.score.plot <- ggplot(
#   site.info,
#   aes(
#     x = !!as.name(new.col.names[2]),
#     y = f1_score
#   )
# ) +
#   geom_boxplot(
#     lwd = 1 / .pt,
#     fatten = 1.2 / .pt,
#     alpha = 0.5,
#     outlier.alpha = 1.0,
#     outlier.size = dot.size,
#     outlier.shape = NA,
#     aes(
#       fill = cell_num,
#       color = cell_num
#     )
#   ) +
#   scale_fill_brewer(
#     name = "Number of cells",
#     palette = "Set2"
#   ) +
#   scale_color_brewer(
#     name = "Number of cells",
#     palette = "Set2"
#   ) +
#   geom_point(
#     position = position_jitterdodge(),
#     aes(color = cell_num),
#     size = dot.size,
#     alpha = 1
#   ) +
#   facet_wrap(
#     ~coverage_variance,
#     ncol = 1,
#     labeller = labeller(
#       coverage_variance = allelic.raw.var
#     ),
#     strip.position = "right"
#   ) +
#   scale_x_discrete(
#     breaks = as.numeric(levels(site.info[[new.col.names[2]]])),
#     label = scientific
#   ) +
#   ylim(0, 1.0) +
#   geom_hline(
#     yintercept = 1,
#     linetype="dashed",
#     color = "gray43",
#     linewidth = 1 / .pt
#   ) +
#   labs(
#     y = "F1 score",
#     x = "Mutation rate"
#   ) +
#   theme_bw() +
#   theme(
#     text = element_text(size = 7),
#     axis.text = element_text(size = 7),
#     legend.text = element_text(size = 6),
#     legend.title = element_text(size = 6),
#     legend.title.align = 0.5,
#     legend.box.spacing = unit(0.3, "mm"),
#     legend.position = "bottom",
#     legend.key.size = unit(4, 'mm'),
#     strip.text.x = element_blank(),
# strip.text.y = strip.text.y,
#     strip.background = element_blank(),
#     panel.grid = element_blank()
#   )

# Efficiency benchmarking ---------------
efficiency.time.plot <-
  ggplot(
    efficiency,
    aes(
      x = scenario,
      y = mean_m
    )
  ) +
  geom_boxplot(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    aes(
      color = tool_setup,
      fill = tool_setup
    ),
    position = position_dodge(width = dodge.width)
  ) +
  scale_fill_manual(
    name = "Method",
    values = efficiency.colors
  ) +
  scale_color_manual(
    name = "Method",
    values = efficiency.colors
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      color = tool_setup
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_x_discrete(
    breaks = c(1, 2),
    label = efficiency.labels
  ) +
  ylim(0, 120) +
  labs(
    y = "Run time (minutes)",
    x = "Scenario"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text.x = strip.text.x,
    strip.text.y = strip.text.y,
    # strip.background = element_blank(),
    panel.grid = element_blank()
  )

efficiency.num.sites.plot <-
  ggplot(
    efficiency %>%
      distinct(
        scenario,
        num_candidate_var_sites
      ),
    aes(
      x = scenario,
      y = num_candidate_var_sites
    )
  ) +
  geom_boxplot(
    width = box.width,
    lwd = box.linewidth,
    fatten = 1.2 / .pt,
    alpha = 0.7,
    outlier.alpha = 1.0,
    outlier.size = dot.size,
    outlier.shape = NA,
    position = position_dodge(width = dodge.width)
  ) +
  jitter.func(
    position = position_jitterdodge(
      jitter.width = jitter.width,
      dodge.width = dodge.width
    ),
    aes(
      color = scenario
    ),
    size = dot.size,
    alpha = jitter.alpha
  ) +
  scale_x_discrete(
    breaks = c(1, 2),
    label = efficiency.labels
  ) +
  scale_color_manual(
    values = c("black", "black")
  ) +
  ylim(0, 2000) +
  labs(
    y = "Number of candidate variant sites",
    x = "Scenario"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.title.align = 0.5,
    legend.box.spacing = unit(0.3, "mm"),
    legend.position = "bottom",
    legend.key.size = unit(4, 'mm'),
    strip.text.x = strip.text.x,
    strip.text.y = strip.text.y,
    # strip.background = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(
    color = "none"
  )

# Exploratory plots -------------

# blank.plot <- ggplot() +
#   geom_blank(aes(1,1)) +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank(),
#     axis.line = element_blank()
#   )

# tree
# tree.plot <- ggarrange(
#   branch.score.diff.plot,
#   normalized.rf.dist.plot,
#   ncol = 2L,
#   align = "h",
#   labels = c("a", "b"),
#   font.label = list(face = "bold", size = 7),
#   common.legend = TRUE,
#   legend = "bottom",
#   legend.grob = get_legend(normalized.rf.dist.plot)
# )
#
# ggsave(
#   here(output.prefix, "tree.pdf"),
#   plot = tree.plot,
#   device = "pdf",
#   width = 183,
#   height = 150,
#   units = "mm",
#   dpi = 300
# )

# parameters
# param.plot <- ggarrange(
#   ggarrange(
#     eff.seq.err.rate.plot,
#     ado.rate.plot,
#     nrow = 1L,
#     align = "h",
#     labels = c("a", "b"),
#     font.label = list(face = "bold", size = 8),
#     legend = "none"
#   ),
#   ggarrange(
#     wildtype.overdispersion.plot,
#     alt.overdispersion.plot,
#     nrow = 1L,
#     align = "h",
#     labels = c("c", "d"),
#     font.label = list(face = "bold", size = 8),
#     legend = "none"
#   ),
#   ggarrange(
#     rela.deletion.rate.plot,
#     # blank.plot,
#     nrow = 1L,
#     labels = "e",
#     font.label = list(face = "bold", size = 8),
#     legend = "none"
#   ),
#   ncol = 2L,
#   # align = "hv",
#   # heights = c(3, 2),
#   common.legend = TRUE,
#   legend = "bottom",
#   legend.grob = get_legend(eff.seq.err.rate.plot)
# )

# param.plot <- ggarrange(
#   eff.seq.err.rate.plot,
#   NULL,
#   rela.deletion.rate.plot,
#   NULL,
#   NULL,
#   NULL,
#   single.ado.rate.plot,
#   NULL,
#   locus.dropout.rate.plot,
#   NULL,
#   NULL,
#   NULL,
#   wildtype.overdispersion.plot,
#   NULL,
#   alt.overdispersion.plot,
#   NULL,
#   NULL,
#   NULL,
#   mean.allelic.cov.plot,
#   NULL,
#   var.allelic.cov.plot,
#   labels = c("a", "", "b", "", "", "", "c", "", "d", "", "", "", "e", "", "f", "", "", "", "g", "", "h"),
#   font.label = list(face = "bold", size = 8),
#   nrow = 7L,
#   ncol = 3L,
#   align = "hv",
#   widths = c(1, -0.09, 1),
#   heights = c(1, -0.11, 1, -0.11, 1, -0.11, 1),
#   common.legend = TRUE,
#   legend = "bottom",
#   legend.grob = get_legend(single.ado.rate.plot)
# )
#
# ggsave(
#   here(output.prefix, "params.pdf"),
#   plot = param.plot,
#   device = "pdf",
#   width = 183,
#   height = 247,
#   units = "mm",
#   dpi = 300
# )

# single mutant genotype calling
# single.mutant.plot <- ggarrange(
#   ggarrange(
#     hetero.recall.plot,
#     hetero.precision.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("a", "b"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ggarrange(
#     hetero.fallout.plot,
#     hetero.f1.score.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("c", "d"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ncol = 2,
#   align = "h",
#   common.legend = TRUE,
#   legend = "bottom",
#   legend.grob = get_legend(hetero.recall.plot)
# )
#
# ggsave(
#   here(output.prefix, "single_mutant.pdf"),
#   plot = single.mutant.plot,
#   device = "pdf",
#   width = 183,
#   height = 150,
#   units = "mm",
#   dpi = 300
# )

# double mutant genotype calling
# double.mutant.plot <- ggarrange(
#   ggarrange(
#     homo.recall.plot,
#     homo.precision.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("a", "b"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ggarrange(
#     homo.fallout.plot,
#     homo.f1.score.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("c", "d"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ncol = 2,
#   align = "h",
#   common.legend = TRUE,
#   legend = "bottom",
#   legend.grob = get_legend(homo.recall.plot)
# )
#
# ggsave(
#   here(output.prefix, "double_mutant.pdf"),
#   plot = double.mutant.plot,
#   device = "pdf",
#   width = 183,
#   height = 150,
#   units = "mm",
#   dpi = 300
# )

# # deletion calling
# del.plot <- ggarrange(
#   ggarrange(
#     del.recall.plot,
#     del.precision.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("a", "b"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ggarrange(
#     del.fallout.plot,
#     del.f1.score.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("c", "d"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ncol = 2,
#   align = "h",
#   common.legend = TRUE,
#   legend = "bottom",
#   legend.grob = get_legend(del.recall.plot)
# )
#
# ggsave(
#   paste0(output.prefix, "del.pdf"),
#   plot = del.plot,
#   device = "pdf",
#   width = 183,
#   height = 150,
#   units = "mm",
#   dpi = 300
# )

# single deletion calling: 1/-
# del.alt.left.plot <- ggarrange(
#   ggarrange(
#     del.alt.left.recall.plot,
#     del.alt.left.precision.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("a", "b"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ggarrange(
#     del.alt.left.fallout.plot,
#     del.alt.left.f1.score.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("c", "d"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ncol = 2,
#   align = "h",
#   common.legend = TRUE,
#   legend = "bottom",
#   legend.grob = get_legend(del.alt.left.recall.plot)
# )
#
# ggsave(
#   here(output.prefix, "del_alt_left.pdf"),
#   plot = del.alt.left.plot,
#   device = "pdf",
#   width = 183,
#   height = 150,
#   units = "mm",
#   dpi = 300
# )

# single deletion calling: 0/-
# del.ref.left.plot <- ggarrange(
#   ggarrange(
#     del.ref.left.recall.plot,
#     del.ref.left.precision.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("a", "b"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ggarrange(
#     del.ref.left.fallout.plot,
#     del.ref.left.f1.score.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("c", "d"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ncol = 2,
#   align = "h",
#   common.legend = TRUE,
#   legend = "bottom",
#   legend.grob = get_legend(del.ref.left.recall.plot)
# )
#
# ggsave(
#   here(output.prefix, "del_ref_left.pdf"),
#   plot = del.ref.left.plot,
#   device = "pdf",
#   width = 183,
#   height = 150,
#   units = "mm",
#   dpi = 300
# )

# double deletion calling: -
# del.all.plot <- ggarrange(
#   ggarrange(
#     del.all.recall.plot,
#     del.all.precision.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("a", "b"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ggarrange(
#     del.all.fallout.plot,
#     del.all.f1.score.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("c", "d"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ncol = 2,
#   align = "h",
#   common.legend = TRUE,
#   legend = "bottom",
#   legend.grob = get_legend(del.all.recall.plot)
# )
#
# ggsave(
#   here(output.prefix, "del_all.pdf"),
#   plot = del.all.plot,
#   device = "pdf",
#   width = 183,
#   height = 150,
#   units = "mm",
#   dpi = 300
# )

# Single ADO calling
# single.ado.plot <- ggarrange(
#   ggarrange(
#     single.ado.recall.plot,
#     single.ado.precision.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("a", "b"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ggarrange(
#     single.ado.fallout.plot,
#     single.ado.f1.score.plot,
#     nrow = 2,
#     align = "v",
#     labels = c("c", "d"),
#     font.label = list(face = "bold", size = 7),
#     legend = "none"
#   ),
#   ncol = 2,
#   align = "h",
#   common.legend = TRUE,
#   legend = "bottom",
#   legend.grob = get_legend(single.ado.recall.plot)
# )
#
# ggsave(
#   here(output.prefix, "single_ado.pdf"),
#   plot = single.ado.plot,
#   device = "pdf",
#   width = 183,
#   height = 150,
#   units = "mm",
#   dpi = 300
# )

# Locus ADO calling
# if ("locus_dropout" %in% ado.modes) {
#   locus.ado.plot <- ggarrange(
#     ggarrange(
#       locus.ado.recall.plot,
#       locus.ado.precision.plot,
#       nrow = 2,
#       align = "v",
#       labels = c("a", "b"),
#       font.label = list(face = "bold", size = 7),
#       legend = "none"
#     ),
#     ggarrange(
#       locus.ado.fallout.plot,
#       locus.ado.f1.score.plot,
#       nrow = 2,
#       align = "v",
#       labels = c("c", "d"),
#       font.label = list(face = "bold", size = 7),
#       legend = "none"
#     ),
#     ncol = 2,
#     align = "h",
#     common.legend = TRUE,
#     legend = "bottom",
#     legend.grob = get_legend(locus.ado.recall.plot)
#   )
#
#   ggsave(
#     here(output.prefix, "locus_ado.pdf"),
#     plot = locus.ado.plot,
#     device = "pdf",
#     width = 183,
#     height = 150,
#     units = "mm",
#     dpi = 300
#   )
# }

# Figure tree -------------------------------------
fig.2 <- ggarrange(
  adjust_facet_length(branch.score.diff.plot + no.legend),
  adjust_facet_length(normalized.rf.dist.plot + no.legend),
  labels = c("a", "b"),
  font.label = list(face = "bold", size = 7),
  nrow = 1,
  widths = c(1, 1),
  align = "h",
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(normalized.rf.dist.plot)
)

ggsave(
  here(output.prefix, "cna", "tree.pdf"),
  plot = fig.2,
  device = "pdf",
  width = 183,
  height = 100,
  units = "mm",
  dpi = 300
)

# Figure alt_left -------------------------------------
fig.3 <- ggarrange(
  ggarrange(
    adjust_facet_length(del.alt.left.recall.plot + no.legend),
    adjust_facet_length(del.alt.left.precision.plot + no.legend),
    labels = c("a", "b"),
    font.label = list(face = "bold", size = 7),
    ncol = 1,
    align = "v",
    common.legend = TRUE,
    legend = "none"
  ),
  ggarrange(
    adjust_facet_length(del.alt.left.fallout.plot + no.legend),
    adjust_facet_length(del.alt.left.f1.score.plot + no.legend),
    labels = c("c", "d"),
    font.label = list(face = "bold", size = 7),
    ncol = 1,
    align = "v",
    common.legend = TRUE,
    legend = "none"
  ),
  nrow = 1,
  align = "h",
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(del.alt.left.recall.plot)
)

ggsave(
  here(output.prefix, "cna", "alt_left.pdf"),
  plot = fig.3,
  device = "pdf",
  width = 183,
  height = 140,
  units = "mm",
  dpi = 300
)

# Figure ref_left-------------------------------------
fig.s2 <- ggarrange(
  ggarrange(
    adjust_facet_length(del.ref.left.recall.plot + no.legend),
    adjust_facet_length(del.ref.left.precision.plot + no.legend),
    labels = c("a", "b"),
    font.label = list(face = "bold", size = 7),
    ncol = 1,
    align = "v",
    common.legend = TRUE,
    legend = "none"
  ),
  ggarrange(
    adjust_facet_length(del.ref.left.fallout.plot + no.legend),
    adjust_facet_length(del.ref.left.f1.score.plot + no.legend),
    labels = c("c", "d"),
    font.label = list(face = "bold", size = 7),
    ncol = 1,
    align = "v",
    common.legend = TRUE,
    legend = "none"
  ),
  nrow = 1,
  align = "h",
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(del.ref.left.recall.plot)
)

ggsave(
  here(output.prefix, "cna", "ref_left.pdf"),
  plot = fig.s2,
  device = "pdf",
  width = 183,
  height = 140,
  units = "mm",
  dpi = 300
)

# Figure hetero_mut-------------------------------------
fig.s3 <- ggarrange(
  ggarrange(
    adjust_facet_length(hetero.recall.plot + no.legend),
    adjust_facet_length(hetero.precision.plot + no.legend),
    labels = c("a", "b"),
    font.label = list(face = "bold", size = 7),
    ncol = 1,
    align = "v",
    common.legend = TRUE,
    legend = "none"
  ),
  ggarrange(
    adjust_facet_length(hetero.fallout.plot + no.legend),
    adjust_facet_length(hetero.f1.score.plot + no.legend),
    labels = c("c", "d"),
    font.label = list(face = "bold", size = 7),
    ncol = 1,
    align = "v",
    common.legend = TRUE,
    legend = "none"
  ),
  nrow = 1,
  align = "h",
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(hetero.recall.plot)
)

ggsave(
  here(output.prefix, "cna", "hetero_mut.pdf"),
  plot = fig.s3,
  device = "pdf",
  width = 183,
  height = 140,
  units = "mm",
  dpi = 300
)

# Figure homo_mut-------------------------------------
fig.s4 <- ggarrange(
  ggarrange(
    adjust_facet_length(homo.recall.plot + no.legend),
    adjust_facet_length(homo.precision.plot + no.legend),
    labels = c("a", "b"),
    font.label = list(face = "bold", size = 7),
    ncol = 1,
    align = "v",
    common.legend = TRUE,
    legend = "none"
  ),
  ggarrange(
    adjust_facet_length(homo.fallout.plot + no.legend),
    adjust_facet_length(homo.f1.score.plot + no.legend),
    labels = c("c", "d"),
    font.label = list(face = "bold", size = 7),
    ncol = 1,
    align = "v",
    common.legend = TRUE,
    legend = "none"
  ),
  nrow = 1,
  align = "h",
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(homo.recall.plot)
)

ggsave(
  here(output.prefix, "cna", "homo_mut.pdf"),
  plot = fig.s4,
  device = "pdf",
  width = 183,
  height = 140,
  units = "mm",
  dpi = 300
)

# Figure single_ado-------------------------------------
fig.s5 <- ggarrange(
  ggarrange(
    adjust_facet_length(single.ado.recall.plot + no.legend),
    adjust_facet_length(single.ado.precision.plot + no.legend),
    labels = c("a", "b"),
    font.label = list(face = "bold", size = 7),
    ncol = 1,
    align = "v",
    common.legend = TRUE,
    legend = "none"
  ),
  ggarrange(
    adjust_facet_length(single.ado.fallout.plot + no.legend),
    adjust_facet_length(single.ado.f1.score.plot + no.legend),
    labels = c("c", "d"),
    font.label = list(face = "bold", size = 7),
    ncol = 1,
    align = "v",
    common.legend = TRUE,
    legend = "none"
  ),
  nrow = 1,
  align = "h",
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(homo.recall.plot)
)

ggsave(
  here(output.prefix, "cna", "single_ado.pdf"),
  plot = fig.s5,
  device = "pdf",
  width = 183,
  height = 140,
  units = "mm",
  dpi = 300
)

# Figure efficiency -------------------------------------
fig.efficiency <- ggarrange(
  efficiency.num.sites.plot,
  efficiency.time.plot + no.legend,
  labels = c("a", "b"),
  font.label = list(face = "bold", size = 7),
  nrow = 1,
  widths = c(1, 2),
  align = "h",
  common.legend = TRUE,
  legend = "bottom",
  legend.grob = get_legend(efficiency.time.plot)
)

ggsave(
  here(output.prefix, "cna", "efficiency.pdf"),
  plot = fig.efficiency,
  device = "pdf",
  width = 183,
  height = 100,
  units = "mm",
  dpi = 300
)
