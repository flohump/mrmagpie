#' @title calcReconversionTargetsBTC2
#' @description Cumulative reconversion targets of the Bending the Curve 2.0 (BTC2) protocol for
#' agricultural land (pool C, "ag") and managed forest (pool D, "forest"), by milestone:
#' m1 = 2030 target, m2 = cumulative 2050 target.
#'
#' The regional targets are distributed to the cells of each AgMIP region in descending order of the
#' restoration priority (RRp) until the target is met. A cell provides at most its pool area, and at
#' most its 2020 cropland and pasture (pool C) or its 2020 secondary forest outside WDPA and the BTC2
#' conservation layer (pool D), since protected forest is not harvested in MAgPIE. A warning lists
#' regions where this land is smaller than the target.
#'
#' @param cells "lpjcell" (67420 cells)
#' @return List with a magpie object (Mha) at 0.5 degree land cells
#' @author Florian Humpenoeder
#' @seealso \code{\link{readBendingTheCurve2}}
#' @examples
#' \dontrun{
#' calcOutput("ReconversionTargetsBTC2", aggregate = FALSE)
#' }
#'
calcReconversionTargetsBTC2 <- function(cells = "lpjcell") {
  if (cells != "lpjcell") {
    stop("Only cells = 'lpjcell' is supported")
  }
  pools <- collapseDim(readSource("BendingTheCurve2", subtype = "pools", convert = "onlycorrect"), dim = 2)
  rrp <- collapseDim(readSource("BendingTheCurve2", subtype = "rrp", convert = "onlycorrect"), dim = 2)
  targets <- readSource("BendingTheCurve2", subtype = "targets", convert = "onlycorrect")
  land <- collapseDim(calcOutput("LanduseInitialisation", aggregate = FALSE, cellular = TRUE, nclasses = "seven",
                                 input_magpie = TRUE, cells = cells, years = "y2020"), dim = 2)
  land <- land[getItems(pools, dim = 1), , ]
  wdpa <- calcOutput("ProtectedAreaBaseline", aggregate = FALSE, nclasses = "seven",
                     cells = cells, magpie_input = TRUE)[getItems(pools, dim = 1), "y2020", "WDPA.secdforest"]
  btc2 <- calcOutput("ConservationBTC2", aggregate = FALSE, cells = cells)[getItems(pools, dim = 1), , "BTC2.secdforest"]
  secdFree <- pmax(0, as.vector(land[, , "secdforest"]) - as.vector(wdpa) - as.vector(btc2))

  capacity <- list(ag = pmin(as.vector(pools[, , "poolC"]), as.vector(dimSums(land[, , c("crop", "past")], dim = 3))),
                   forest = pmin(as.vector(pools[, , "poolD"]), secdFree))
  poolOf <- c(ag = "poolC", forest = "poolD")

  mapping <- toolGetMapping("regionmappingAgMIP.csv", type = "regional", where = "mappingfolder")
  region <- mapping$RegionCode[match(getItems(pools, dim = "iso", full = TRUE), mapping$CountryCode)]
  if (anyNA(region) || !setequal(unique(region), getItems(targets, dim = 1))) {
    stop("AgMIP regions of the cells do not match the regions of the targets")
  }
  priority <- as.vector(rrp)

  # fill cells by descending priority; the same ranking for both milestones keeps m1 within m2
  fill <- function(area, prio, target) {
    ord <- order(prio, decreasing = TRUE)
    before <- c(0, utils::head(cumsum(area[ord]), -1))
    out <- numeric(length(area))
    out[ord] <- pmin(area[ord], pmax(0, target - before))
    return(out)
  }

  x <- new.magpie(getItems(pools, dim = 1), NULL, c("ag.m1", "ag.m2", "forest.m1", "forest.m2"), fill = 0,
                  sets = c("x.y.iso", "t", "pool.milestone"))
  short <- NULL
  for (r in getItems(targets, dim = 1)) {
    i <- which(region == r)
    for (p in names(poolOf)) {
      for (m in c("m1", "m2")) {
        target <- as.vector(targets[r, c(m1 = "y2030", m2 = "y2050")[m], poolOf[p]])
        x[i, , paste(p, m, sep = ".")] <- fill(capacity[[p]][i], priority[i], target)
      }
      gap <- as.vector(targets[r, "y2050", poolOf[p]]) - sum(x[i, , paste0(p, ".m2")])
      if (gap > 1e-6) {
        short <- c(short, paste0(r, " ", p, " ", round(gap, 2), " Mha"))
      }
    }
  }
  if (!is.null(short)) {
    warning("Land smaller than the BTC2 reconversion target in: ", paste(short, collapse = ", "))
  }
  return(list(x = x, weight = NULL, unit = "Mha",
              description = "Cumulative BTC2 reconversion targets of pools C (ag) and D (forest) since 2025",
              isocountries = FALSE))
}
