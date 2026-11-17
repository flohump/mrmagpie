#' @title correctBendingTheCurve2
#' @description Sets missing and negative values of the BTC2 protocol layers to zero.
#' @param x magpie object provided by the read function
#' @param subtype see \code{\link{readBendingTheCurve2}}
#' @return magpie object
#' @author Florian Humpenoeder
#' @examples
#' \dontrun{
#' readSource("BendingTheCurve2", subtype = "pools", convert = "onlycorrect")
#' }
#'
correctBendingTheCurve2 <- function(x, subtype) {
  if (subtype == "targets" && anyNA(x)) {
    stop("Missing restoration targets")
  }
  x <- toolConditionalReplace(x, conditions = c("is.na()", "<0"), replaceby = 0)
  return(x)
}
