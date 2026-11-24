#' @title calcPlantationClimateClass
#' @description Koeppen-Geiger climate class shares of the timber plantation estate. The cell's class shares are
#' weighted by its plantation area in 1995 (forestry area of the land initialisation times the country's share of
#' plantation forest in planted forest, FRA 2025); clusters without plantations get land-area-weighted shares.
#'
#' @param cells "lpjcell" for 67420 cells
#' @return magpie object in cellular resolution with the plantation area as weight
#' @author Florian Humpenoeder
#'
#' @examples
#' \dontrun{
#' calcOutput("PlantationClimateClass", aggregate = "cluster")
#' }
#' @importFrom magclass getItems setYears dimSums
#' @importFrom madrat calcOutput
calcPlantationClimateClass <- function(cells = "lpjcell") {
  koeppen <- calcOutput("ClimateClass", datasource = "koeppen", years = "y2001", aggregate = FALSE)
  landIni <- calcOutput("LanduseInitialisation", nclasses = "seven", aggregate = FALSE, cellular = TRUE,
                        cells = cells, input_magpie = TRUE)[, 1995, ]
  plantShare <- calcOutput("PlantedForest", aggregate = FALSE)
  iso <- getItems(landIni, dim = "iso", full = TRUE)
  share <- as.vector(plantShare)[match(iso, getItems(plantShare, dim = 1))]
  share[is.na(share)] <- 0
  plantArea <- setYears(landIni[, , "forestry"], NULL) * share
  landArea <- setYears(dimSums(landIni, dim = 3), NULL)
  # the tiny land-area term decides only in clusters without plantations
  weight <- plantArea + 10^-15 * landArea

  return(list(x = setYears(koeppen[getItems(landIni, dim = 1), , ], NULL),
              weight = weight,
              unit = "share",
              description = "Koeppen-Geiger climate class shares of the timber plantation estate",
              isocountries = FALSE))
}
