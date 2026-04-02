# Internal helpers for optimization result tables ---------------------------

.rga_null_model_formula <- function(formula,
                                    data,
                                    fallback = NULL) {
  if (!is.null(formula)) {
    return(stats::update(formula, . ~ . - cd))
  }

  .mlpe_formula_from_data(
    data,
    response = "gd",
    predictor = "1",
    fallback = fallback
  )
}

.rga_fixed_effect_label <- function(fit.mod) {
  fixed_formula <- lme4::nobars(stats::formula(fit.mod))
  terms_object <- stats::terms(fixed_formula)
  fixed_terms <- attr(terms_object, "term.labels", exact = TRUE)

  if (length(fixed_terms) == 0L) {
    return("1")
  }

  paste(fixed_terms, collapse = " + ")
}

.rga_surface_k <- function(surface,
                           fit.mod,
                           GA.inputs) {
  fixed_k <- length(lme4::fixef(fit.mod))

  surface_idx <- match(surface, GA.inputs$layer.names)
  if (is.na(surface_idx)) {
    return(fixed_k)
  }

  surface_parms <- GA.inputs$parm.type$n.parm[[surface_idx]]

  if (GA.inputs$k.value == 1) {
    return(fixed_k)
  }

  if (GA.inputs$k.value == 2) {
    return(surface_parms + fixed_k - 1L)
  }

  if (GA.inputs$k.value == 3) {
    return(surface_parms + length(GA.inputs$layer.names) + fixed_k - 1L)
  }

  1L + fixed_k - 1L
}

.rga_multisurface_k <- function(fit.mod, GA.inputs) {
  fixed_k <- length(lme4::fixef(fit.mod))

  if (GA.inputs$k.value == 1) {
    return(fixed_k)
  }

  if (GA.inputs$k.value == 2) {
    return(sum(GA.inputs$parm.type$n.parm) + fixed_k - 1L)
  }

  if (GA.inputs$k.value == 3) {
    return(
      sum(GA.inputs$parm.type$n.parm) +
        length(GA.inputs$layer.names) +
        fixed_k - 1L
    )
  }

  length(GA.inputs$layer.names) + fixed_k - 1L
}

.rga_update_result_table_models <- function(result_table,
                                            fit_list,
                                            GA.inputs,
                                            n_pops) {
  if (is.null(result_table) || nrow(result_table) == 0L ||
      !"Surface" %in% names(result_table)) {
    return(result_table)
  }

  fixed_effects <- rep(NA_character_, nrow(result_table))

  for (i in seq_len(nrow(result_table))) {
    surface <- as.character(result_table$Surface[[i]])
    fit.mod <- fit_list[[surface]]

    if (is.null(fit.mod)) {
      next
    }

    fixed_effects[[i]] <- .rga_fixed_effect_label(fit.mod)
    result_table$k[[i]] <- .rga_surface_k(surface, fit.mod, GA.inputs)

    if ("LL" %in% names(result_table) && is.finite(result_table$LL[[i]])) {
      result_table$AIC[[i]] <- (-2 * result_table$LL[[i]]) +
        (2 * result_table$k[[i]])
      result_table$AICc[[i]] <- result_table$AIC[[i]] +
        (((2 * result_table$k[[i]]) * (result_table$k[[i]] + 1)) /
           max((n_pops - result_table$k[[i]] - 1), 1))
    }
  }

  result_table$Fixed.Effects <- fixed_effects

  front <- c(
    "Surface",
    "Fixed.Effects",
    paste0("obj.func_", GA.inputs$method),
    "k",
    "AIC",
    "AICc",
    "R2m",
    "R2c",
    "LL"
  )
  keep_cols <- intersect(front, names(result_table))
  other_cols <- setdiff(names(result_table), keep_cols)
  result_table[, c(keep_cols, other_cols), drop = FALSE]
}
