repository <- "https://stat.ethz.ch/CRAN/"

if (!"stringr" %in% installed.packages()) {
  install.packages("stringr", repos = repository)
}
library(stringr)


prettify.colors <- function(tool.names, reorder.sieve.delsieve = TRUE) {
  sorted.tool.names <- character()
  sorted.tool.fills <- character()
  sorted.tool.colors <- character()

  tool.sifit.names <- character()
  tool.sifit.fills <- character()
  tool.sifit.colors <- character()

  tool.monovar.names <- character()
  tool.monovar.fills <- character()
  tool.monovar.colors <- character()

  tool.sciphi.names <- character()
  tool.sciphi.fills <- character()
  tool.sciphi.colors <- character()

  tool.sieve.names <- character()
  tool.sieve.fills <- character()
  tool.sieve.colors <- character()

  tool.delsieve.names <- character()
  tool.delsieve.fills <- character()
  tool.delsieve.colors <- character()

  # c(
  #   "#0173b2", "#de8f05", "#cc78bc", "#949494", "#ca9161"
  # )

  predefined.colors.sifit <- c(
    "#048a64",
    "#029e73"
    # "#008169",
    # "#10b698"
  )
  predefined.colors.monovar <- c(
    "#d775b9",
    "#fbafe4"
    # "#00dcb5",
    # "#78f9e2"
  )
  predefined.colors.sciphi <- c(
    "#149af2",
    "#72c2f7"
    # "#1f77ab",
    # "#56b4e9"
    # "#cfad22",
    # "#ffdf61",
    # "#9bcf22",
    # "#ceec89"
  )
  predefined.colors.sieve <- c(
    "#3c5488",
    "#7c93c5",
    "#149af2",
    "#72c2f7"
    # "#3c5488",
    # "#7c93c5",
    # "#a9ab13",
    # "#eef214",
    # "#1348ab",
    # "#4982ed"
  )
  predefined.colors.delsieve <- c(
    "#b95000",
    "#d55e00",
    "#e64b35",
    "#f09285"
    # "#990000",
    # "#ed6e6e",
    # "#9400e6",
    # "#c873f6",
    # "#66267d",
    # "#8c5f9c"
  )

  # sifit, monovar, sciphi, sieve, delsieve
  color.index <- c(1, 1, 1, 1, 1)
  for (tool in tool.names) {
    if (grepl("^sifit", tool, ignore.case = TRUE)) {
      tool.sifit.names <- c(tool.sifit.names, tool)
      tool.sifit.colors <- c(tool.sifit.colors, predefined.colors.sifit[color.index[1]])
      tool.sifit.fills <- c(tool.sifit.fills, predefined.colors.sifit[color.index[1] + 1])
      color.index[1] <- color.index[1] + 2
    } else if (grepl("^monovar", tool, ignore.case = TRUE)) {
      tool.monovar.names <- c(tool.monovar.names, tool)
      tool.monovar.colors <- c(tool.monovar.colors, predefined.colors.monovar[color.index[2]])
      tool.monovar.fills <- c(tool.monovar.fills, predefined.colors.monovar[color.index[2] + 1])
      color.index[2] <- color.index[2] + 2
    } else if (grepl("^sciphi", tool, ignore.case = TRUE)) {
      tool.sciphi.names <- c(tool.sciphi.names, tool)
      tool.sciphi.colors <- c(tool.sciphi.colors, predefined.colors.sciphi[color.index[3]])
      tool.sciphi.fills <- c(tool.sciphi.fills, predefined.colors.sciphi[color.index[3] + 1])
      color.index[3] <- color.index[3] + 2
    } else if (grepl("^delsieve", tool, ignore.case = TRUE)) {
      tool.delsieve.names <- c(tool.delsieve.names, tool)
      tool.delsieve.colors <- c(tool.delsieve.colors, predefined.colors.delsieve[color.index[5]])
      tool.delsieve.fills <- c(tool.delsieve.fills, predefined.colors.delsieve[color.index[5] + 1])
      color.index[5] <- color.index[5] + 2
    } else if (grepl("^sieve", tool, ignore.case = TRUE)) {
      tool.sieve.names <- c(tool.sieve.names, tool)
      tool.sieve.colors <- c(tool.sieve.colors, predefined.colors.sieve[color.index[4]])
      tool.sieve.fills <- c(tool.sieve.fills, predefined.colors.sieve[color.index[4] + 1])
      color.index[4] <- color.index[4] + 2
    }
  }

  if (reorder.sieve.delsieve) {
    adjusted.tool.sieve.names <- adjust_order(
      tool.sieve.names,
      c("^sieve_uncorrected$", "^sieve_neg_.*?", "^sieve$", "^sieve_pos_.*?"),
      ignore.case = TRUE,
      perl = TRUE
    )
    tool.sieve.names <- tool.sieve.names[adjusted.tool.sieve.names]
    tool.sieve.colors <- tool.sieve.colors[adjusted.tool.sieve.names]
    tool.sieve.fills <- tool.sieve.fills[adjusted.tool.sieve.names]

    adjusted.tool.delsieve.names <- adjust_order(
      tool.delsieve.names,
      c("^delsieve_uncorrected$", "^delsieve_neg_.*?", "^delsieve$", "^delsieve_pos_.*?"),
      ignore.case = TRUE,
      perl = TRUE
    )
    tool.delsieve.names <- tool.delsieve.names[adjusted.tool.delsieve.names]
    tool.delsieve.colors <- tool.delsieve.colors[adjusted.tool.delsieve.names]
    tool.delsieve.fills <- tool.delsieve.fills[adjusted.tool.delsieve.names]
  }

  sorted.tool.names <- c(tool.sifit.names, tool.monovar.names, tool.sciphi.names, tool.sieve.names, tool.delsieve.names)
  sorted.tool.colors <- c(tool.sifit.colors, tool.monovar.colors, tool.sciphi.colors, tool.sieve.colors, tool.delsieve.colors)
  sorted.tool.fills <- c(tool.sifit.fills, tool.monovar.fills, tool.sciphi.fills, tool.sieve.fills, tool.delsieve.fills)

  names(sorted.tool.colors) <- sorted.tool.names
  names(sorted.tool.fills) <- sorted.tool.names

  return(list(
    tool = sorted.tool.names,
    color = sorted.tool.colors,
    fill = sorted.tool.fills
  ))
}


adjust_order <- function(names, orderred_patterns, ...) {
  vapply(
    orderred_patterns,
    function(x, ...) grep(x, ...),
    integer(1),
    names,
    ...,
    USE.NAMES = FALSE
  )
}


rename.sieve <- function(tool.name, tool.setup, more.info = FALSE) {
  if (!is.na(tool.setup)) {
    if (grepl(".+?-.+?_nm$", tool.setup)) {
      return(paste(tool.name, "neg_manipulated", sep = "_"))
    } else if (grepl(".+?-.+?_pm$", tool.setup)) {
      return(paste(tool.name, "pos_manipulated", sep = "_"))
    } else if (grepl("-", tool.setup, fixed = TRUE)) {
      # return(paste(tool.name, "stage2", sep = "_"))
      if (more.info) {
        info <- str_split_i(str_split_i(tool.setup, "-", 1), "_", 3)

        if (info == "ld")
          return(paste(tool.name, "LDO", sep = "_"))
        else if (info == "sa")
          return(paste(tool.name, "ADO", sep = "_"))
        else
          return(tool.name)
      } else {
        return(tool.name)
      }
    } else if (grepl("_muindel$|_mudel$|_mu$", tool.setup, ignore.case = TRUE, fixed = FALSE)) {
      return(paste(tool.name, "combined", sep = "_"))
    } else if (grepl("_muindel_|_mudel_|_mu_", tool.setup, ignore.case = TRUE, fixed = FALSE)) {
      return(paste(tool.name, "uncorrected", sep = "_"))
    } else {
      stop(paste("Unhandled condition encountered:", tool.name, "and", tool.setup, sep = " "))
    }
  }
}


rename.tool <- function(tool.name, snv.type, tool.setup) {
  renamed.tool.setup <- ""
  if (!is.null(tool.setup) && !is.na(tool.setup)) {
    if (grepl("true_parameters$", tool.setup, ignore.case = TRUE)) {
      renamed.tool.setup <- "tp"
    } else if (grepl("true_parameters_zero$", tool.setup, ignore.case = TRUE)) {
      renamed.tool.setup <- "tp0"
    } else if (grepl("^ep$", tool.setup, ignore.case = TRUE)) {
      renamed.tool.setup <- "EP"
    } else if (grepl("^pl$", tool.setup, ignore.case = TRUE)) {
      renamed.tool.setup <- "PL"
    } else {
      stop(paste("Unhandled condition encountered:", tool.name, ",", snv.type, ",", tool.setup, sep = " "))
    }
  }

  renamed.snv.type <- ""
  if (!is.null(snv.type) && !is.na(snv.type)) {
    if (grepl("true_monovar_snvs$", snv.type, ignore.case = TRUE)) {
      renamed.snv.type <- "mnv"
    } else if (grepl("sieve_snvs$", snv.type, ignore.case = TRUE)) {
      renamed.snv.type <- "sv"
    } else {
      stop(paste("Unhandled condition encountered:", tool.name, ",", snv.type, ",", tool.setup, sep = " "))
    }
  }

  if (renamed.tool.setup != "" && renamed.snv.type != "") {
    return(paste(tool.name, renamed.snv.type, renamed.tool.setup, sep = "_"))
  } else if (renamed.tool.setup != "") {
    return(paste(tool.name, renamed.tool.setup, sep = "_"))
  } else if (renamed.snv.type != "") {
    return(paste(tool.name, renamed.snv.type, sep = "_"))
  } else {
    return(tool.name)
  }
}


rename.tools.internal <- function(row.data, rename.sieve, rename.delsieve, rename.sifit.monovar, more.info = FALSE) {
  if ((rename.sieve > 0 || rename.delsieve > 0) && grepl("sieve", row.data["tool"], ignore.case = TRUE)) {
    return(
      rename.sieve(
        row.data["tool"],
        row.data["tool_setup"],
        more.info = more.info
      )
    )
    # } else if (rename.cellphy > 0 && grepl("cellphy", row.data["tool"], ignore.case = TRUE)) {
    #   return(
    #     rename.tool.helper(
    #       rename.cellphy,
    #       row.data
    #       )
    #     )
  } else if (rename.sifit.monovar > 0 && grepl("sifit", row.data["tool"], ignore.case = TRUE)) {
    return(
      rename.tool.helper(
        rename.sifit.monovar,
        row.data
      )
    )
  } else if (rename.sifit.monovar > 0 && grepl("monovar", row.data["tool"], ignore.case = TRUE)) {
    return(
      rename.tool.helper(
        rename.sifit.monovar,
        row.data
      )
    )
  } else {
    return(row.data["tool"])
  }
}


rename.tool.helper <- function(rename.level, row.data) {
  if (rename.level == 0) {
    return(row.data["tool"])
  } else if (rename.level == 1) {
    return(
      rename.tool(
        row.data["tool"],
        row.data["snv_type"],
        NULL
      )
    )
  } else if (rename.level == 2) {
    return(
      rename.tool(
        row.data["tool"],
        NULL,
        row.data["tool_setup"]
      )
    )
  } else if (rename.level == 3) {
    return(
      rename.tool(
        row.data["tool"],
        row.data["snv_type"],
        row.data["tool_setup"]
      )
    )
  }
}


rename.tools <- function(data, more.info = FALSE) {
  rename.sieve <- determine.rename.level(data[data$tool == "SIEVE", ])
  rename.delsieve <- determine.rename.level(data[data$tool == "DelSIEVE", ])
  rename.sifit.monovar <- determine.rename.level(data[data$tool == "SiFit", ]) + determine.rename.level(data[data$tool == "Monovar", ])

  return(
    apply(
      data, 1, rename.tools.internal,
      rename.sieve = rename.sieve,
      rename.delsieve = rename.delsieve,
      rename.sifit.monovar = rename.sifit.monovar,
      more.info = more.info
    )
  )
}


# 0: do not rename
# 1: rename with "snv_type"
# 2: rename with "tool_setup"
# 3: rename with both "snv_type" and "tool_setup"
determine.rename.level <- function(data) {
  if (nrow(data) == 0) {
    return(0)
  }

  length.snv.type <- length(levels(factor(data$snv_type)))
  length.tool.setup <- length(levels(factor(data$tool_setup)))

  if (length.snv.type <= 1 && length.tool.setup <= 1) {
    return(0)
  } else if (length.snv.type > 1 && length.tool.setup > 1) {
    return(3)
  } else if (length.snv.type <= 1 && length.tool.setup > 1) {
    return(2)
  } else if (length.snv.type > 1 && length.tool.setup <= 1) {
    return(1)
  }
}
