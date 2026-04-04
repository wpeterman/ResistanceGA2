# Internal helpers for ResistanceGA2 result objects -----------------------

resga_add_class <- function(x, class_name) {
  current <- class(x)

  if (is.null(current)) {
    class(x) <- class_name
  } else {
    class(x) <- unique(c(class_name, current))
  }

  x
}

resga_pick_column <- function(x, candidates) {
  hit <- intersect(candidates, names(x))
  if (length(hit) == 0) {
    return(NULL)
  }

  hit[[1]]
}

resga_first_numeric_column <- function(x) {
  numeric_cols <- names(x)[vapply(x, is.numeric, logical(1))]
  if (length(numeric_cols) == 0) {
    return(NULL)
  }

  numeric_cols[[1]]
}

resga_label_column <- function(x) {
  resga_pick_column(x, c("Surface", "surface", "model", "name"))
}

resga_direction <- function(metric) {
  if (is.null(metric)) {
    return("ascending")
  }

  if (metric %in% c("AIC", "AICc", "avg.AIC", "avg.AICc", "avg.rank",
                    "delta", "delta.AICc", "RMSE", "avg.RMSE")) {
    return("ascending")
  }

  "descending"
}

resga_rank_table <- function(x, preferred) {
  x <- resga_plain_table(x)
  metric <- resga_pick_column(x, preferred)

  if (is.null(metric)) {
    metric <- resga_first_numeric_column(x)
  }

  direction <- resga_direction(metric)

  if (!is.null(metric)) {
    ord <- order(x[[metric]], decreasing = identical(direction, "descending"),
                 na.last = TRUE)
    x <- x[ord, , drop = FALSE]
  }

  list(
    table = x,
    metric = metric,
    direction = direction,
    label = resga_label_column(x)
  )
}

resga_head_table <- function(x, n = 6L) {
  utils::head(resga_plain_table(x), n = max(1L, as.integer(n)))
}

resga_plain_table <- function(x) {
  if (!is.data.frame(x)) {
    return(x)
  }

  class(x) <- setdiff(
    class(x),
    grep("^(resga_|summary\\.resga_)", class(x), value = TRUE)
  )

  if (length(class(x)) == 0) {
    class(x) <- "data.frame"
  }

  x
}

resga_minutes <- function(seconds) {
  as.numeric(seconds) / 60
}

resga_barplot <- function(values,
                          labels,
                          main,
                          xlab,
                          col = "steelblue",
                          ...) {
  graphics::barplot(
    rev(values),
    names.arg = rev(labels),
    horiz = TRUE,
    las = 1,
    col = col,
    border = NA,
    main = main,
    xlab = xlab,
    ...
  )
}

resga_is_replicated_all_comb <- function(x) {
  length(x) > 0 &&
    !is.null(names(x)) &&
    all(startsWith(names(x), "rep_")) &&
    all(vapply(x, function(run) {
      is.list(run) && "summary.table" %in% names(run)
    }, logical(1)))
}

resga_all_comb_runs <- function(x) {
  if (resga_is_replicated_all_comb(x)) {
    return(x)
  }

  list(rep_1 = x)
}

#' S3 methods for ResistanceGA2 result objects
#'
#' Lightweight S3 methods that provide concise printing, summaries, and quick
#' ranking plots for optimization and bootstrap outputs returned by
#' \code{\link{Resist.boot}}, \code{\link{SS_optim}}, \code{\link{MS_optim}},
#' and \code{\link{all_comb}}.
#'
#' @param x A ResistanceGA2 result object.
#' @param object A ResistanceGA2 result object.
#' @param n Number of rows to display in the printed summary.
#' @param metric Ranking metric to visualize.
#' @param top_n Maximum number of rows to display in a plot.
#' @param type Plot type for multisurface optimization results.
#' @param main Optional plot title.
#' @param ... Additional arguments passed to plotting methods.
#' @rawNamespace S3method(plot,resga_all_comb)
#' @rawNamespace S3method(plot,resga_bootstrap)
#' @rawNamespace S3method(plot,resga_ms_optim)
#' @rawNamespace S3method(plot,resga_ss_optim)
#' @rawNamespace S3method(print,resga_all_comb)
#' @rawNamespace S3method(print,resga_bootstrap)
#' @rawNamespace S3method(print,resga_ms_optim)
#' @rawNamespace S3method(print,resga_ss_optim)
#' @rawNamespace S3method(print,summary.resga_all_comb)
#' @rawNamespace S3method(print,summary.resga_bootstrap)
#' @rawNamespace S3method(print,summary.resga_ms_optim)
#' @rawNamespace S3method(print,summary.resga_ss_optim)
#' @rawNamespace S3method(summary,resga_all_comb)
#' @rawNamespace S3method(summary,resga_bootstrap)
#' @rawNamespace S3method(summary,resga_ms_optim)
#' @rawNamespace S3method(summary,resga_ss_optim)
#'
#' @name ResistanceGA2-result-methods
NULL

#' @rdname ResistanceGA2-result-methods
#' @method print resga_bootstrap
#' @export
print.resga_bootstrap <- function(x, n = 6L, ...) {
  print(summary(x, n = n, ...))
  invisible(x)
}

#' @rdname ResistanceGA2-result-methods
#' @method summary resga_bootstrap
#' @export
summary.resga_bootstrap <- function(object, n = 6L, ...) {
  ranked <- resga_rank_table(
    object,
    preferred = c("avg.rank", "avg.AICc", "avg.weight", "Percent.top", "support")
  )

  best_surface <- NA_character_
  if (!is.null(ranked$label) && nrow(ranked$table) > 0) {
    best_surface <- as.character(ranked$table[[ranked$label]][1])
  }

  structure(
    list(
      n_models = nrow(object),
      best_surface = best_surface,
      metric = ranked$metric,
      table = ranked$table,
      n = as.integer(n)
    ),
    class = "summary.resga_bootstrap"
  )
}

#' @rdname ResistanceGA2-result-methods
#' @method print summary.resga_bootstrap
#' @export
print.summary.resga_bootstrap <- function(x, ...) {
  cat("ResistanceGA2 bootstrap results\n")
  cat("  models evaluated:", x$n_models, "\n")

  if (!is.na(x$best_surface)) {
    cat("  top-supported model:", x$best_surface, "\n")
  }

  if (!is.null(x$metric)) {
    cat("  ranked by:", x$metric, "\n")
  }

  cat("\n")
  print(resga_head_table(x$table, x$n), row.names = FALSE)
  invisible(x)
}

#' @rdname ResistanceGA2-result-methods
#' @method plot resga_bootstrap
#' @export
plot.resga_bootstrap <- function(x,
                                 metric = c("avg.weight", "Percent.top", "avg.rank"),
                                 top_n = 10L,
                                 main = NULL,
                                 ...) {
  ranked <- resga_rank_table(x, preferred = metric)

  if (is.null(ranked$metric) || is.null(ranked$label)) {
    stop("No plottable ranking columns were found in this bootstrap result.")
  }

  top <- resga_head_table(ranked$table, top_n)
  if (is.null(main)) {
    main <- "Bootstrap model support"
  }

  xlab <- ranked$metric
  if (identical(ranked$direction, "ascending")) {
    xlab <- paste0(xlab, " (lower is better)")
  }

  resga_barplot(
    values = top[[ranked$metric]],
    labels = as.character(top[[ranked$label]]),
    main = main,
    xlab = xlab,
    ...
  )

  invisible(top)
}

#' @rdname ResistanceGA2-result-methods
#' @method print resga_ss_optim
#' @export
print.resga_ss_optim <- function(x, n = 6L, ...) {
  print(summary(x, n = n, ...))
  invisible(x)
}

#' @rdname ResistanceGA2-result-methods
#' @method summary resga_ss_optim
#' @export
summary.resga_ss_optim <- function(object, n = 6L, ...) {
  ranked <- resga_rank_table(
    object$AICc,
    preferred = c("AICc", "AIC", "R2m", "LL")
  )

  best_surface <- NA_character_
  best_aicc <- NA_real_

  if (!is.null(ranked$label) && nrow(ranked$table) > 0) {
    best_surface <- as.character(ranked$table[[ranked$label]][1])
  }
  if ("AICc" %in% names(ranked$table) && nrow(ranked$table) > 0) {
    best_aicc <- ranked$table$AICc[1]
  }

  structure(
    list(
      n_models = if (is.null(object$AICc)) 0L else nrow(object$AICc),
      best_surface = best_surface,
      best_aicc = best_aicc,
      run_time = object$Run.Time,
      categorical_surfaces = if (is.null(object$CategoricalResults)) 0L else nrow(object$CategoricalResults),
      continuous_surfaces = if (is.null(object$ContinuousResults)) 0L else nrow(object$ContinuousResults),
      table = ranked$table,
      n = as.integer(n)
    ),
    class = "summary.resga_ss_optim"
  )
}

#' @rdname ResistanceGA2-result-methods
#' @method print summary.resga_ss_optim
#' @export
print.summary.resga_ss_optim <- function(x, ...) {
  cat("ResistanceGA2 single-surface optimization\n")
  cat("  models evaluated:", x$n_models, "\n")
  cat("  categorical surfaces:", x$categorical_surfaces, "\n")
  cat("  continuous surfaces:", x$continuous_surfaces, "\n")

  if (!is.na(x$best_surface)) {
    cat("  best surface:", x$best_surface, "\n")
  }

  if (!is.na(x$best_aicc)) {
    cat("  best AICc:", signif(x$best_aicc, 6), "\n")
  }

  if (!is.null(x$run_time) && length(x$run_time) == 1L && is.finite(x$run_time)) {
    cat("  run time (min):", signif(resga_minutes(x$run_time), 6), "\n")
  }

  cat("\n")
  print(resga_head_table(x$table, x$n), row.names = FALSE)
  invisible(x)
}

#' @rdname ResistanceGA2-result-methods
#' @method plot resga_ss_optim
#' @export
plot.resga_ss_optim <- function(x,
                                metric = c("AICc", "delta.AICc", "R2m", "LL"),
                                top_n = 10L,
                                main = NULL,
                                ...) {
  tab <- x$AICc

  if (!"delta.AICc" %in% names(tab) && "AICc" %in% names(tab)) {
    tab$delta.AICc <- tab$AICc - min(tab$AICc, na.rm = TRUE)
  }

  ranked <- resga_rank_table(tab, preferred = metric)

  if (is.null(ranked$metric) || is.null(ranked$label)) {
    stop("No plottable ranking columns were found in this SS_optim result.")
  }

  top <- resga_head_table(ranked$table, top_n)
  if (is.null(main)) {
    main <- "Single-surface model ranking"
  }

  xlab <- ranked$metric
  if (identical(ranked$direction, "ascending")) {
    xlab <- paste0(xlab, " (lower is better)")
  }

  resga_barplot(
    values = top[[ranked$metric]],
    labels = as.character(top[[ranked$label]]),
    main = main,
    xlab = xlab,
    ...
  )

  invisible(top)
}

#' @rdname ResistanceGA2-result-methods
#' @method print resga_ms_optim
#' @export
print.resga_ms_optim <- function(x, n = 6L, ...) {
  print(summary(x, n = n, ...))
  invisible(x)
}

#' @rdname ResistanceGA2-result-methods
#' @method summary resga_ms_optim
#' @export
summary.resga_ms_optim <- function(object, n = 6L, ...) {
  ranked <- resga_rank_table(
    object$AICc.tab,
    preferred = c("AICc", "AIC", "R2m", "LL")
  )
  contribution <- resga_rank_table(
    object$percent.contribution,
    preferred = c("mean", "weight")
  )$table

  best_surface <- NA_character_
  if (!is.null(ranked$label) && nrow(ranked$table) > 0) {
    best_surface <- as.character(ranked$table[[ranked$label]][1])
  }

  structure(
    list(
      best_surface = best_surface,
      fit_table = ranked$table,
      contribution_table = contribution,
      n = as.integer(n)
    ),
    class = "summary.resga_ms_optim"
  )
}

#' @rdname ResistanceGA2-result-methods
#' @method print summary.resga_ms_optim
#' @export
print.summary.resga_ms_optim <- function(x, ...) {
  cat("ResistanceGA2 multisurface optimization\n")

  if (!is.na(x$best_surface)) {
    cat("  best surface:", x$best_surface, "\n")
  }

  cat("\nBest-fit table\n")
  print(resga_head_table(x$fit_table, 1L), row.names = FALSE)

  if (!is.null(x$contribution_table) && nrow(x$contribution_table) > 0) {
    cat("\nPercent contribution\n")
    print(resga_head_table(x$contribution_table, x$n), row.names = FALSE)
  }

  invisible(x)
}

#' @rdname ResistanceGA2-result-methods
#' @method plot resga_ms_optim
#' @export
plot.resga_ms_optim <- function(x,
                                type = c("contribution", "fit"),
                                metric = c("mean", "AICc", "R2m", "LL"),
                                top_n = 10L,
                                main = NULL,
                                ...) {
  type <- match.arg(type)

  if (identical(type, "contribution")) {
    tab <- resga_rank_table(
      x$percent.contribution,
      preferred = c("mean", "weight")
    )

    if (is.null(tab$metric) || is.null(tab$label)) {
      stop("No contribution columns were found in this MS_optim result.")
    }

    top <- resga_head_table(tab$table, top_n)
    if (is.null(main)) {
      main <- "Multisurface percent contribution"
    }

    resga_barplot(
      values = top[[tab$metric]],
      labels = as.character(top[[tab$label]]),
      main = main,
      xlab = tab$metric,
      col = "darkolivegreen3",
      ...
    )

    return(invisible(top))
  }

  ranked <- resga_rank_table(x$AICc.tab, preferred = metric)

  if (is.null(ranked$metric) || is.null(ranked$label)) {
    stop("No plottable fit columns were found in this MS_optim result.")
  }

  top <- resga_head_table(ranked$table, top_n)
  if (is.null(main)) {
    main <- "Multisurface model fit"
  }

  xlab <- ranked$metric
  if (identical(ranked$direction, "ascending")) {
    xlab <- paste0(xlab, " (lower is better)")
  }

  resga_barplot(
    values = top[[ranked$metric]],
    labels = as.character(top[[ranked$label]]),
    main = main,
    xlab = xlab,
    col = "darkolivegreen3",
    ...
  )

  invisible(top)
}

#' @rdname ResistanceGA2-result-methods
#' @method print resga_all_comb
#' @export
print.resga_all_comb <- function(x, n = 6L, ...) {
  print(summary(x, n = n, ...))
  invisible(x)
}

#' @rdname ResistanceGA2-result-methods
#' @method summary resga_all_comb
#' @export
summary.resga_all_comb <- function(object, n = 6L, ...) {
  runs <- resga_all_comb_runs(object)

  best_by_run <- lapply(seq_along(runs), function(i) {
    run_name <- names(runs)[i]
    run <- runs[[i]]
    summary_ranked <- resga_rank_table(
      run$summary.table,
      preferred = c("AICc", "delta.AICc", "weight")
    )
    boot_ranked <- NULL
    if (is.data.frame(run$boot.results)) {
      boot_ranked <- resga_rank_table(
        run$boot.results,
        preferred = c("avg.rank", "avg.AICc", "avg.weight", "Percent.top", "support")
      )
    }

    data.frame(
      replicate = run_name,
      best_model = if (is.null(summary_ranked$label) || nrow(summary_ranked$table) == 0) {
        NA_character_
      } else {
        as.character(summary_ranked$table[[summary_ranked$label]][1])
      },
      best_bootstrap_model = if (is.null(boot_ranked) || is.null(boot_ranked$label) ||
                                 nrow(boot_ranked$table) == 0) {
        NA_character_
      } else {
        as.character(boot_ranked$table[[boot_ranked$label]][1])
      },
      stringsAsFactors = FALSE
    )
  })

  best_by_run <- do.call(rbind, best_by_run)
  consensus_tab <- sort(table(best_by_run$best_model), decreasing = TRUE)
  consensus <- data.frame(
    best_model = names(consensus_tab),
    n_runs = as.integer(consensus_tab),
    stringsAsFactors = FALSE
  )

  single_run <- !resga_is_replicated_all_comb(object)
  summary_table <- NULL
  bootstrap_table <- NULL

  if (single_run) {
    summary_table <- resga_rank_table(
      object$summary.table,
      preferred = c("AICc", "delta.AICc", "weight")
    )$table
    if (is.data.frame(object$boot.results)) {
      bootstrap_table <- resga_rank_table(
        object$boot.results,
        preferred = c("avg.rank", "avg.AICc", "avg.weight", "Percent.top", "support")
      )$table
    }
  }

  structure(
    list(
      replicated = !single_run,
      n_runs = length(runs),
      best_by_run = best_by_run,
      consensus = consensus,
      summary_table = summary_table,
      bootstrap_table = bootstrap_table,
      n = as.integer(n)
    ),
    class = "summary.resga_all_comb"
  )
}

#' @rdname ResistanceGA2-result-methods
#' @method print summary.resga_all_comb
#' @export
print.summary.resga_all_comb <- function(x, ...) {
  cat("ResistanceGA2 all-combinations analysis\n")
  cat("  replicate runs:", x$n_runs, "\n")

  if (!x$replicated && !is.null(x$summary_table)) {
    cat("\nCandidate models\n")
    print(resga_head_table(x$summary_table, x$n), row.names = FALSE)

    if (!is.null(x$bootstrap_table)) {
      cat("\nBootstrap support\n")
      print(resga_head_table(x$bootstrap_table, x$n), row.names = FALSE)
    }

    return(invisible(x))
  }

  cat("\nBest model by replicate\n")
  print(x$best_by_run, row.names = FALSE)

  if (nrow(x$consensus) > 0) {
    cat("\nConsensus across replicates\n")
    print(x$consensus, row.names = FALSE)
  }

  invisible(x)
}

#' @rdname ResistanceGA2-result-methods
#' @method plot resga_all_comb
#' @export
plot.resga_all_comb <- function(x,
                                metric = c("avg.weight", "Percent.top", "AICc"),
                                top_n = 10L,
                                main = NULL,
                                ...) {
  if (!resga_is_replicated_all_comb(x)) {
    if (is.data.frame(x$boot.results)) {
      boot_ranked <- resga_rank_table(
        x$boot.results,
        preferred = metric
      )

      if (!is.null(boot_ranked$metric) && !is.null(boot_ranked$label)) {
        top <- resga_head_table(boot_ranked$table, top_n)
        if (is.null(main)) {
          main <- "All-combinations bootstrap support"
        }

        xlab <- boot_ranked$metric
        if (identical(boot_ranked$direction, "ascending")) {
          xlab <- paste0(xlab, " (lower is better)")
        }

        resga_barplot(
          values = top[[boot_ranked$metric]],
          labels = as.character(top[[boot_ranked$label]]),
          main = main,
          xlab = xlab,
          col = "tan3",
          ...
        )

        return(invisible(top))
      }
    }

    ranked <- resga_rank_table(
      x$summary.table,
      preferred = c("AICc", "delta.AICc", "weight")
    )

    if (is.null(ranked$metric) || is.null(ranked$label)) {
      stop("No plottable ranking columns were found in this all_comb result.")
    }

    top <- resga_head_table(ranked$table, top_n)
    if (is.null(main)) {
      main <- "All-combinations model ranking"
    }

    xlab <- ranked$metric
    if (identical(ranked$direction, "ascending")) {
      xlab <- paste0(xlab, " (lower is better)")
    }

    resga_barplot(
      values = top[[ranked$metric]],
      labels = as.character(top[[ranked$label]]),
      main = main,
      xlab = xlab,
      col = "tan3",
      ...
    )

    return(invisible(top))
  }

  runs <- summary(x)$consensus
  if (nrow(runs) == 0) {
    stop("No replicate consensus information is available to plot.")
  }

  if (is.null(main)) {
    main <- "Best model frequency across replicates"
  }

  resga_barplot(
    values = runs$n_runs,
    labels = as.character(runs$best_model),
    main = main,
    xlab = "Number of replicate wins",
    col = "tan3",
    ...
  )

  invisible(runs)
}
