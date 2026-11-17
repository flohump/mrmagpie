#' @title calcConservationBTC2
#' @description Land proposed for protection by 2030 in the Bending the Curve 2.0 (BTC2) protocol,
#' in addition to the WDPA protected area of 2020, by MAgPIE land type.
#'
#' The protocol gives the proposed area by IUCN habitat class. Each class is assigned to the land types
#' of the BTC2 legend, in proportion to the unprotected land of these types in the cell in 2020.
#' Area that finds no such land is substituted by other unprotected land of the cell: natural land
#' first, then pasture, forestry and cropland. Built-up classes are substituted in the same way.
#' The result is capped by the unprotected land of the cell.
#'
#' @param cells "lpjcell" (67420 cells)
#' @return List with a magpie object (Mha) at 0.5 degree land cells
#' @author Florian Humpenoeder
#' @seealso \code{\link{readBendingTheCurve2}}
#' @examples
#' \dontrun{
#' calcOutput("ConservationBTC2", aggregate = FALSE)
#' }
#'
calcConservationBTC2 <- function(cells = "lpjcell") {
  if (cells != "lpjcell") {
    stop("Only cells = 'lpjcell' is supported")
  }
  consv <- readSource("BendingTheCurve2", subtype = "consv2030", convert = "onlycorrect")

  land <- calcOutput("LanduseInitialisation", aggregate = FALSE, cellular = TRUE, nclasses = "seven",
                     input_magpie = TRUE, cells = cells, years = "y2020")
  wdpa <- calcOutput("ProtectedAreaBaseline", aggregate = FALSE, nclasses = "seven",
                     cells = cells, magpie_input = TRUE)[, "y2020", "WDPA"]
  types <- c("crop", "past", "forestry", "primforest", "secdforest", "other")
  avl <- as.array(collapseDim(land[, , types] - collapseDim(wdpa[, , types], dim = 3.1), dim = 2))[, 1, ]
  avl[avl < 0] <- 0
  consv <- as.array(collapseDim(consv, dim = 2))[rownames(avl), 1, ]

  # IUCN classes to land types of the BTC2 legend (protocol hybrid legend; Plantations as forestry)
  classTypes <- list(Forest = c("primforest", "secdforest", "forestry"),
                     Savanna = c("other", "past", "primforest", "secdforest"),
                     Shrubland = "other", Grassland = c("past", "other"), Wetlands_inland = "other",
                     Rocky_Areas = "other", Desert = "other", Arable_Land = "crop", Pastureland = "past",
                     Plantations = "forestry")
  if (!all(names(classTypes) %in% colnames(consv))) {
    stop("IUCN classes missing in the BTC2 conservation layer")
  }
  substitution <- list(c("primforest", "secdforest", "other"), "past", "forestry", "crop")

  # allocate area to land types in proportion to their remaining unprotected land, at most that land
  allocate <- function(area, typeCols) {
    avlSum <- rowSums(avl[, typeCols, drop = FALSE])
    alloc <- pmin(area, avlSum)
    share <- avl[, typeCols, drop = FALSE] / ifelse(avlSum > 0, avlSum, 1)
    out <- matrix(0, nrow(avl), ncol(avl), dimnames = dimnames(avl))
    out[, typeCols] <- alloc * share
    return(out)
  }

  out <- avl * 0
  shortfall <- rowSums(consv[, setdiff(colnames(consv), names(classTypes)), drop = FALSE])
  # single-type classes first, so that mixed classes take the remaining land
  ordered <- names(classTypes)[order(lengths(classTypes))]
  for (cl in ordered) {
    a <- allocate(consv[, cl], classTypes[[cl]])
    out <- out + a
    avl <- avl - a
    shortfall <- shortfall + consv[, cl] - rowSums(a)
  }
  for (typeCols in substitution) {
    a <- allocate(shortfall, typeCols)
    out <- out + a
    avl <- avl - a
    shortfall <- shortfall - rowSums(a)
  }

  x <- new.magpie(rownames(out), NULL, paste0("BTC2.", c(types, "urban")), fill = 0,
                  sets = c("x.y.iso", "t", "consv.land"))
  x[, , paste0("BTC2.", types)] <- out
  x <- x[, , paste0("BTC2.", c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other"))]

  message("BTC2 conservation: ", round(sum(consv)), " Mha proposed, ", round(sum(out)), " Mha allocated, ",
          round(sum(consv) - sum(out)), " Mha without unprotected land")
  return(list(x = x, weight = NULL, unit = "Mha",
              description = "Land proposed for protection by 2030 in addition to WDPA 2020 (BTC2 protocol)",
              isocountries = FALSE))
}
