
#' Summarizing paar objects
#' @inheritParams base::summary
#' @returns A data.frame with the following columns:
#' \itemize{
#'  \item \code{condition} a character vector with the final condition.
#'  \item \code{n} a numeric vector with the number of rows for each condition.
#'  \item \code{percentage} a numeric vector with the percentage of rows for each condition.
#'  }
#'
#' @export
summary.paar <- function(object, ...) {
  condition <- object$condition
  condition[is.na(condition)] <- 'normal point'
  fr_condition_table <- table(condition)
  fr_condition_table <- sort(fr_condition_table, decreasing = TRUE)
  pr_condition_table <- prop.table(fr_condition_table) * 100
  pr_f_condition_table <-
    paste0(signif(pr_condition_table, 2), '%')
  # f_condition_table <- paste0(fr_condition_table, ' (',
  #                             pr_f_condition_table, ')')
  # names(f_condition_table) <- names(fr_condition_table)

  df_condition <-
    data.frame(
      'condition' = names(fr_condition_table),
      cbind('n' = fr_condition_table,
            'percentage' = pr_condition_table)
    )
  rownames(df_condition) <- NULL

  # message(f_condition_table)
  # f_condition_table
  df_condition <- structure(df_condition,
                            class = c('summary.paar', 'data.frame'))
  df_condition
}


#' Print summarized paar object
#' @inheritParams base::print
#' @export

print.summary.paar <- function(x, digits, ...) {
  fr_condition_table <- x[['n']]
  pr_f_condition_table <- x[['percentage']]
  pr_f_condition_table <-
    paste0(signif(pr_f_condition_table, 2), '%')
  f_condition_table <- paste0(fr_condition_table, ' (',
                              pr_f_condition_table, ')')
  names(f_condition_table) <- x[['condition']]
  print(f_condition_table,
        quote = FALSE,
        right = TRUE)
}

#' Print paar objects
#' @inheritParams base::print
#' @param n an integer vector specifying maximum number of rows or
#'   elements to print.
#' @export
print.paar <- function(x, n = 3, ...) {
  p_removed <- sum(!is.na(x$condition))/length(x$condition)*100
  p_removed_f <- paste0(signif(p_removed, 2), '%')
  n_normal <- sum(is.na(x$condition))
  cat('Depurated data has', n_normal,'rows.\n')
  cat('The process removed',
      p_removed_f,
      'of original data.\n\n')
  cat('$depurated_data\n')
  print(x$depurated_data, n = n)
  cat('\n\n')
  cat('$condition\n')
  cat('vector of length', paste0(length(x$condition),'.'),
      'First', n, 'elements:\n')
  print(utils::head(x$condition, n = n))

  invisible(x)

}



#' Bind outlier condition to an object.
#' @name bind
#' @return \code{cbind} called with m.
#' @param ... objects to bind.
#' @inheritParams base::cbind
#' @export
#'
cbind.paar = function(..., deparse.level = 1) {
  dots <- list(...)
  dots <- dots[!sapply(dots, is.null)]
  classes <- sapply(dots, function(x) inherits(x, 'paar'))
  if (any(unlist(classes))) {
    paar_object <- dots[classes]
    no_paar_object <- dots[!classes]
    if (length(paar_object) > 1) {
      stop('Only one paar object is allowed.', .call = FALSE)
    }
    if (length(no_paar_object) == 0) {
      stop('No object to bind.', .call = FALSE)
    }
    if (length(no_paar_object) > 1) {
      stop('Only one object to bind is allowed.', .call = FALSE)
    }
    dots[classes][[1]] <- dots[classes][[1]][['condition']]
    if (is.null(names(dots))) {names(dots) <- rep('', length(dots))}
    if (names(dots)[classes] == "") {
      names(dots)[classes] <- 'condition'
    }
    do.call(base::cbind, dots)
  } else {
    stop('No paar object found.')
  }

}
