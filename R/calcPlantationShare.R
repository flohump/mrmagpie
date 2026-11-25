#' @title calcPlantationShare
#' @description Share of plantation forest in the planted forest of the 1995 land initialisation: the FRA 2025
#' country share (mrland \code{calcPlantedForest}) on every cell, weighted by the cell's planted forest area.
#' Clusters without planted forest get the land-area-weighted share.
#'
#' @param cells "lpjcell" for 67420 cells
#' @return magpie object in cellular resolution with the planted forest area as weight
#' @author Florian Humpenoeder
#' @seealso \code{\link{calcPlantationClimateClass}}
#'
#' @examples
#' \dontrun{
#' calcOutput("PlantationShare", aggregate = "cluster")
#' }
#' @importFrom magclass getItems setYears dimSums setNames
#' @importFrom madrat calcOutput
calcPlantationShare <- function(cells = "lpjcell") {
  landIni <- calcOutput("LanduseInitialisation", nclasses = "seven", aggregate = FALSE, cellular = TRUE,
                        cells = cells, input_magpie = TRUE)[, 1995, ]
  plantShare <- calcOutput("PlantedForest", aggregate = FALSE)
  share <- as.vector(plantShare)[match(getItems(landIni, dim = "iso", full = TRUE), getItems(plantShare, dim = 1))]
  share[is.na(share)] <- 0
  planted <- setYears(landIni[, , "forestry"], NULL)
  x <- setNames(planted, NULL)
  x[, , ] <- share
  landArea <- setYears(dimSums(landIni, dim = 3), NULL)
  # the tiny land-area term decides only in clusters without planted forest
  weight <- planted + 10^-15 * landArea

  return(list(x = x,
              weight = weight,
              unit = "share",
              min = 0,
              max = 1,
              description = "Share of plantation forest in planted forest (FRA 2025 country shares)",
              isocountries = FALSE))
}
