#' @title calcAgeClassDistribution
#' @description Forest area in 15 age classes. With \code{dataset = "GFAD"} (default) the
#' distribution comes directly from the Global Forest Age Dataset (GFAD v1.1, Poulter et al.
#' 2019): per-cell class fractions times cell area, summed over the four forest types.
#' With \code{dataset = "GAMI"} the *age distribution* comes from GAMI v2.1 (Besnard et al.
#' 2024) while the per-cell forest *area* is taken from the GFAD product. GAMI supplies only
#' the within-forest age shares (12 classes, summing to 1); keeping the GFAD forest-area map
#' as the spatial weight conserves each cell's total forest area and isolates the pure
#' age-structure change. GAMI's 12 twenty-year classes are mapped onto the 15 GFAD ten-year
#' classes (each split 50/50 into two adjacent classes; all classes >=140 yr, i.e. GAMI
#' `140-160 ... >299`, collapse into `class15`, the mature bin routed to `acx` in module 28).
#' Where GAMI reports no forest but GFAD does, the GFAD age shape is retained (area-conserving).
#'
#' For both datasets, forest of the 1995 land initialisation that the age map does not cover
#' (dry woodland and sparse boreal forest outside the GFAD forest mask, forest cleared
#' between 1995 and the 2010 age map, and with \code{backcast} the land of stands younger
#' than 15 yr in 2010) is added to \code{class15}, so that the age structure
#' describes the forest MAgPIE starts with, with the uncovered part treated as mature.
#' @param cells lpjcell for 67420 cells
#' @param dataset "GFAD" (default) or "GAMI"
#' @param gamiYear GAMI reference year when \code{dataset = "GAMI"}: "y2010" (default, matches
#'   the GFAD vintage) or "y2020"
#' @param youngCorrection with \code{dataset = "GAMI"}: move young forest that does not regrow to the
#'   older classes, see \code{.youngCorrection} (default TRUE)
#' @param backcast with \code{dataset = "GAMI"}: move the ages of the age map back from \code{gamiYear} to
#'   MAgPIE's start year 1995, see \code{.backcast} (default TRUE)
#'
#' @return magpie object in cluster resolution
#' @author Abhijeet Mishra, Felicitas Beier, Florian Humpenoeder
#'
#' @examples
#' \dontrun{
#' calcOutput("AgeClassDistribution", aggregate = FALSE)
#' calcOutput("AgeClassDistribution", dataset = "GAMI", aggregate = FALSE)
#' }

calcAgeClassDistribution <- function(cells = "lpjcell", dataset = "GFAD", gamiYear = "y2010",
                                     youngCorrection = TRUE, backcast = TRUE) {

  if (cells == "magpiecell") {
    stop("magpiecell is no longer supported")
  }

  if (dataset == "GFAD") {

    out <- .gfadAreaByClass()
    descr <- "Forest area in 15 age classes based on the Global Forest Age Dataset (GFAD) from Poulter et al. 2019"

  } else if (dataset == "GAMI") {

    # within-forest age shares from GAMI v2.1 (native 0.5-degree grid, 12 classes, sum to 1 where forest exists)
    gami <- readSource("GAMI", convert = "onlycorrect")

    # project GAMI onto the 67420 lpjcell grid and attach ISO codes (done here, not in the read function)
    coordMap <- mstools::toolGetMappingCoord2Country()
    gami <- gami[coordMap$coords, , ]
    getItems(gami, dim = 1, raw = TRUE) <- paste(coordMap$coords, coordMap$iso, sep = ".")
    getSets(gami) <- c("x", "y", "iso", "year", "ageClass")

    if (!gamiYear %in% getItems(gami, 2)) {
      stop("calcAgeClassDistribution: gamiYear must be one of ", paste(getItems(gami, 2), collapse = ", "))
    }
    shares <- collapseNames(gami[, gamiYear, ])

    # per-cell forest AREA from the GFAD product (same forest-area map -> isolates the age change)
    gfad  <- .gfadAreaByClass()
    cellIt <- getItems(gfad, 1)
    gfadM  <- matrix(as.numeric(gfad[cellIt, , ]), nrow = length(cellIt), ncol = 15,
                     dimnames = list(cellIt, paste0("class", 1:15)))
    forestArea <- rowSums(gfadM)
    sharesM <- matrix(as.numeric(shares[cellIt, , ]), nrow = length(cellIt), ncol = 12,
                      dimnames = list(cellIt, NULL))

    # remap 12 GAMI 20-yr classes -> 15 GFAD 10-yr classes
    remap <- matrix(0, nrow = 12, ncol = 15)
    for (g in 1:7) {
      remap[g, 2 * g - 1] <- 0.5
      remap[g, 2 * g]     <- 0.5
    }
    remap[8:12, 15] <- 1                                   # all GAMI classes >=140 yr -> class15 (acx)
    shape15 <- sharesM %*% remap                           # rowSums = 1 (forest) or 0

    # fallback: keep the GFAD age shape where GAMI has no forest but GFAD does (area-conserving)
    gfadShape <- gfadM / pmax(forestArea, 1e-12)
    useGami   <- rowSums(sharesM) > 1e-6
    shape15[!useGami, ] <- gfadShape[!useGami, ]
    moved <- 0
    if (youngCorrection) {
      yc <- .youngCorrection(shape15, gami, cellIt, useGami, forestArea)
      shape15 <- yc$shape
      moved <- yc$moved
    }
    if (backcast) {
      shape15 <- .backcast(shape15, as.integer(sub("y", "", gamiYear)) - 1995)
    }
    shape15 <- .toOld(shape15, moved)

    areaByClass <- shape15 * forestArea                    # per-cell total forest area == GFAD (conserved)

    out <- gfad                                            # reuse the GFAD magpie skeleton (cells x 1 x 15)
    out[cellIt, , ] <- areaByClass

    descr <- paste0("Forest area in 15 age classes: age distribution from GAMI v2.1 (Besnard et al. 2024, ",
                    gamiYear, ") on the GFAD-derived per-cell forest-area weight")

  } else {
    stop("calcAgeClassDistribution: 'dataset' must be 'GFAD' or 'GAMI'")
  }

  # forest of the land initialisation (MAgPIE's start year 1995) beyond the age map goes to the oldest class
  cellIt <- getItems(out, 1)
  land   <- calcOutput("LanduseInitialisation", nclasses = "seven", cellular = TRUE, cells = cells,
                       input_magpie = TRUE, years = "y1995", aggregate = FALSE)
  forest <- as.numeric(dimSums(land[cellIt, , c("primforest", "secdforest", "forestry")], dim = 3))
  excess <- pmax(forest - as.numeric(dimSums(out, dim = 3)), 0)
  out[cellIt, , "class15"] <- out[cellIt, , "class15"] + excess
  descr <- paste0(descr, "; forest of the 1995 land initialisation beyond the age map in class15")

  return(list(x = out,
              weight = NULL,
              unit = "Mha",
              description = descr,
              isocountries = FALSE))
}

# GFAD forest area per cell and age class (Mha): per-cell class fractions times cell area, summed over forest types
.gfadAreaByClass <- function() {
  poulterDataset <- readSource("GFAD", convert = "onlycorrect")

  mapping <- mstools::toolGetMappingCoord2Country()
  cb <- toolGetMapping("LPJ_CellBelongingsToCountries.csv",
                       type = "cell", where = "mrcommons")
  # 111e3 m ~ length of one degree of latitude; (111e3 * 0.5)^2 = area of a 0.5-degree cell
  # at the equator (m2); cos(lat) corrects for the convergence of meridians towards the poles
  cellArea <- (111e3 * 0.5) * (111e3 * 0.5) * cos(cb$lat / 180 * pi)
  cellArea <- as.magpie(cellArea, spatial = 1) /  1e10 # convert to Mha
  getItems(cellArea, dim = 1, raw = TRUE) <- paste(mapping$coords, mapping$iso, sep = ".")
  getSets(cellArea) <- c("x", "y", "iso", "year", "data")

  out <- poulterDataset * cellArea
  out <- dimSums(out, dim = 3.1)
  getNames(out) <- gsub(pattern = "X", replacement = "class", x = getNames(out))
  return(out)
}

# GAMI predicts the 2010 and 2020 ages separately from biomass, so low-biomass forest is young in both years.
# Forest aged 0-20 in 2020 beyond the ageing 2010 stands and ten years of stand-replacing loss (Sims et al. 2025)
# does not regrow; young forest of natural-disturbance origin (30 yr of wildfire and other natural disturbances)
# regrows outside land-use change, as in bookkeeping; fallow of shifting cultivation (ten years of its clearing) is cut
# again before it matures. All three leave the classes under 30 yr; returns the reduced shape and the moved share.
.youngCorrection <- function(shape15, gami, cellIt, useGami, forestArea) {
  iso  <- sub("^.*\\.", "", cellIt)
  loss <- dimSums(readSource("GFWLossByDriver")[, 2011:2020, ], dim = 2) / 10           # Mha per yr
  countryForest <- tapply(forestArea, iso, sum)
  rate <- function(drivers) {                                                          # share of forest per yr
    x <- dimSums(loss[, , drivers], dim = 3)
    r <- stats::setNames(as.numeric(x), getItems(x, 1))[iso] / pmax(countryForest[iso], 1e-6)
    as.numeric(pmin(ifelse(is.na(r), 0, r), 0.05)) * useGami                          # tapply index leaves a 1-d array
  }
  natural   <- rate(c("wildfire", "other_natural_disturbances"))
  replacing <- rate(c("wildfire", "other_natural_disturbances", "shifting_cultivation", "logging"))
  fallow    <- 10 * rate("shifting_cultivation")
  y0 <- function(year) as.numeric(gami[cellIt, year, "gami1"])                          # GAMI class 0-20
  persistent <- pmax(y0("y2020") - 0.5 * y0("y2010") - 10 * replacing, 0) * useGami

  young <- shape15[, 1:3, drop = FALSE]
  under20 <- rowSums(young[, 1:2, drop = FALSE])
  p <- pmin(persistent, under20)
  s <- pmin(fallow, under20 - p)
  young[, 1:2] <- young[, 1:2] * ifelse(under20 > 0, 1 - (p + s) / pmax(under20, 1e-12), 1)
  under30 <- rowSums(young)
  f <- pmin(30 * natural, under30)
  young <- young * ifelse(under30 > 0, 1 - f / pmax(under30, 1e-12), 1)
  return(list(shape = cbind(young, shape15[, 4:15, drop = FALSE]), moved = p + s + f))
}

# The age map observes the GAMI year; MAgPIE starts in 1995. Stands younger than the gap did not exist yet and are
# dropped; the 1995 forest on their land is added to class15 by the land-initialisation step. class15 stays.
.backcast <- function(shape15, years) {
  if (years < 0 || years %% 5 != 0) stop("calcAgeClassDistribution: backcast needs a multiple of 5 years")
  if (years == 0) return(shape15)
  halves  <- shape15[, rep(1:14, each = 2), drop = FALSE] / 2                         # 5-yr halves under 140 yr
  shifted <- cbind(halves[, -seq_len(years / 5), drop = FALSE], matrix(0, nrow(halves), years / 5))
  out <- cbind(shifted[, seq(1, 27, 2), drop = FALSE] + shifted[, seq(2, 28, 2), drop = FALSE], shape15[, 15])
  colnames(out) <- colnames(shape15)
  return(out)
}

# share moved out of the young classes goes to the classes of 30 yr and older in proportion (class15 if none)
.toOld <- function(shape15, moved) {
  old <- shape15[, 4:15, drop = FALSE]
  oldShare <- old / pmax(rowSums(old), 1e-12)
  oldShare[rowSums(old) <= 0, 12] <- 1
  shape15[, 4:15] <- old + moved * oldShare
  return(shape15)
}
