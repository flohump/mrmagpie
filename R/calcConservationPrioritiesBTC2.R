#' @title calcConservationPrioritiesBTC2
#' @description Conservation priority areas of \code{mrland::calcConservationPriorities} plus the
#' Bending the Curve 2.0 option "BTC2" (\code{\link{calcConservationBTC2}}).
#' @param cells "lpjcell" (67420 cells)
#' @return List with a magpie object (Mha) at 0.5 degree land cells
#' @author Florian Humpenoeder
#' @examples
#' \dontrun{
#' calcOutput("ConservationPrioritiesBTC2", aggregate = FALSE)
#' }
#'
calcConservationPrioritiesBTC2 <- function(cells = "lpjcell") {
  x <- calcOutput("ConservationPriorities", nclasses = "seven", cells = cells, aggregate = FALSE)
  btc2 <- calcOutput("ConservationBTC2", cells = cells, aggregate = FALSE)
  x <- mbind(x, btc2[getItems(x, dim = 1), , getItems(x, dim = 3.2)])
  return(list(x = x, weight = NULL, unit = "Mha",
              description = "Conservation priority areas, including the BTC2 protocol option",
              isocountries = FALSE))
}
