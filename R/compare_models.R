#' Compare models
#'
#' @param ... models to be compared. Must be models created by serosv. If models' names are not provided, indices will be used instead for the `model` column in the returned data.frame.
#'
#'
#' @return
#' a data.frame of 4 columns
#'   \item{model}{name or index of the model}
#'   \item{type}{model type of the given model (a serosv model name)}
#'   \item{AIC}{AIC value for the model (lower value indicates better fit)}
#'   \item{BIC}{BIC value for the model (lower value indicates better fit)}
#'
#' @importFrom magrittr %>%
#' @importFrom purrr imap_dfr
#' @importFrom stringr str_detect
#'
#' @export
compare_models <- function(...){
  list(...) %>%
    imap_dfr(~ {
      # return error if input contains non-serosv models
      if(!all(str_detect(class(.x), "_model"))) {
        stop("Inputs must be serosv models")
      }

      data.frame(
        model = .y,
        type = class(.x),
        AIC = AIC(.x$info),
        BIC = BIC(.x$info)
      )
    })
}


