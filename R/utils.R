
#' Summarizing paar objects
#' @inheritParams base::summary
#' @export
summary.paar <- function(object, ...) {
  condition <- object$condition
  condition[is.na(condition)] <- 'normal point'
  fr_condition_table <- table(condition)
  fr_condition_table <- sort(fr_condition_table, decreasing = TRUE)
  pr_condition_table <- prop.table(fr_condition_table) * 100
  pr_f_condition_table <-
    paste0(signif(pr_condition_table, 2), '%')
  f_condition_table <- paste0(fr_condition_table, ' (',
                              pr_f_condition_table, ')')
  names(f_condition_table) <- names(fr_condition_table)

  df_condition <-
    data.frame(
      'condition' = names(f_condition_table),
      cbind('n' = fr_condition_table,
            'percentage' = pr_condition_table)
    )
  rownames(df_condition) <- NULL

  print.default(f_condition_table,
                print.gap = 3L, quote = FALSE)

  invisible(df_condition)
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
