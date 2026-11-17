#' @title readBendingTheCurve2
#' @description Reads the conservation and restoration layers of the Bending the Curve 2.0 (BTC2)
#' modelling protocol (run 3, November 2026).
#' @param subtype
#' \itemize{
#'   \item "consv2030": area by IUCN habitat class, unprotected in 2020 and proposed as protected
#'   by 2030 (Jung et al. expansion priorities, WDPA Nov 2024), in Mha
#'   \item "consv2020": area by IUCN habitat class protected in 2020 (WDPA Nov 2024), in Mha
#'   \item "pools": degraded natural land (A), degraded agricultural land (B), agricultural land (C)
#'   and managed forest (D), in Mha
#'   \item "rrp": region normalized range size rarity (restoration priority), 0 to 1, mean over the
#'   0.1 degree cells of each 0.5 degree cell
#'   \item "fnf": share of the cell in forest biomes of the RESOLVE ecoregions
#'   \item "targets": restoration targets by AgMIP region and pool for 2030 and cumulative 2050, in Mha
#' }
#' @return magpie object at 0.5 degree land cells (lpjcell), or by AgMIP region for "targets"
#' @author Florian Humpenoeder
#' @examples
#' \dontrun{
#' readSource("BendingTheCurve2", subtype = "consv2030", convert = "onlycorrect")
#' }
#'
readBendingTheCurve2 <- function(subtype) {
  map <- toolGetMappingCoord2Country(pretty = TRUE)

  # raster layers (shares of pixel area, or values) at the land cells
  toCells <- function(r) {
    out <- as.magpie(as.matrix(terra::extract(r, map[c("lon", "lat")], ID = FALSE)), spatial = 1)
    dimnames(out) <- list("x.y.iso" = paste(map$coords, map$iso, sep = "."), "t" = NULL, "data" = names(r))
    return(out)
  }
  # pixel shares to Mha
  toArea <- function(r) {
    return(r * terra::cellSize(r[[1]], unit = "ha") / 1e6)
  }

  if (subtype %in% c("consv2030", "consv2020")) {
    folder <- c(consv2030 = "LC_protected_2030", consv2020 = "LC_protected_2020")[subtype]
    files <- list.files(file.path("conservation", folder), pattern = "\\.tif$", full.names = TRUE)
    r <- terra::rast(files)
    names(r) <- gsub(" ", "_", sub(".*__[0-9]+_(.*)__ver.*", "\\1", basename(files)))
    out <- toCells(toArea(r))
  } else if (subtype == "pools") {
    files <- c(A = "PoolA_LUH3_EII48.tif", B = "PoolB_LUH3.tif", C = "PoolC_LUH3.tif", D = "PoolD_LUH3.tif")
    r <- terra::rast(file.path("restoration", "PoolA-D", files))
    names(r) <- paste0("pool", names(files))
    out <- toCells(toArea(r))
  } else if (subtype == "rrp") {
    r <- terra::rast(file.path("restoration", "RRp", "RRp_2025_regionnormalized.tif"))
    r <- terra::aggregate(r, fact = 5, fun = "mean", na.rm = TRUE)
    names(r) <- "rrp"
    out <- toCells(r)
  } else if (subtype == "fnf") {
    r <- terra::rast(file.path("restoration", "forest-nonforest", "FNF.tif"))
    r <- terra::aggregate(r, fact = round(0.5 / terra::res(r)[1]), fun = "mean", na.rm = TRUE)
    names(r) <- "forest"
    out <- toCells(r)
  } else if (subtype == "targets") {
    x <- as.data.frame(readxl::read_excel(file.path("restoration", "AgMIP-restoration-targets-261027.xlsx"),
                                          col_names = FALSE, .name_repair = "minimal"))
    header <- as.character(x[3, ])
    rows <- x[[1]] %in% c("ANZ", "BRA", "CAN", "CHN", "EUR", "FSU", "IND", "MEN", "OAS", "OSA", "SEA", "SSA", "USA")
    out <- NULL
    for (pool in c("A", "B", "C", "D")) {
      col <- grep(paste0("^Target Pool ", pool, ":"), header)
      if (length(col) != 1) {
        stop("Target column of pool ", pool, " not found")
      }
      # columns: 2030, increment 2030-2050, cumulative 2050; Gha to Mha
      vals <- cbind(as.numeric(x[rows, col]), as.numeric(x[rows, col + 2])) * 1000
      out <- mbind(out, new.magpie(x[rows, 1], c("y2030", "y2050"), paste0("pool", pool), fill = vals))
    }
    getSets(out) <- c("region", "year", "data")
  } else {
    stop("Unknown subtype")
  }
  return(out)
}
