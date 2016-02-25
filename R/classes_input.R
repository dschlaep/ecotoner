## Class definitions

#' An S4-class describing raster grid information which all of the grids must correspond to
#' 
#' @slot crs A \link[sp::CRS]{CRS} object.
#' @slot res_m A numeric value. The cell size of the grids in meters. The cells must be square.
#' @slot longlat A logical. A flag indicating whether the grids' coordinates represent geographic longitude/latitude (TRUE) or not (FALSE).
#' @slot origin A numeric vector of length two. The origin of the grids.
#' @slot rotation A logical. A flag indicating whether the grids are rotated (TRUE) or not (FALSE).
#' 
GridInfo <- setClass("GridInfo",
						slots = list(crs = structure("CRS", package = "sp"),
									res_m = "numeric",
									longlat = "logical",
									origin = "numeric",
									rotation = "logical"
								),
						prototype = list(crs = sp::CRS(NA_character_), res_m = NA_real_, longlat = NA, origin = c(NA_real_, NA_real_), rotation = NA)
					)

# TODO(drs): validify that units of slot 'res_m' are correct

#' An S4-class to represent the information of the input raster grids
#'
#' @slot grid_env A \code{\linkS4class{RasterLayer}} object. The environmental grid representing the gradient to investigate, e.g., elevation
#' @slot grid_veg A \code{\linkS4class{RasterLayer}} object. The vegetation grid representing the two vegetation (types) that make up the ecotone/ecological boundary of consideration
#' @slot grid_flow A \code{\linkS4class{RasterLayer}} object. A grid representing the flow path of the environmental grid. Used only if \code{transect_type == 3}
#' @slot grid_abut A \code{\linkS4class{RasterLayer}} object. A grid marking the cells where neighboring cells are of the opposite vegetation (type).
#' @slot grid_aspect_mean A \code{\linkS4class{RasterLayer}} object. A grid representing the mean smoothed aspect of the environmental grid. Used only if \code{transect_type == 4}
#' @slot grid_aspect_sd A \code{\linkS4class{RasterLayer}} object. A grid representing the standard deviation of the smoothed aspect of the environmental grid. Used only if \code{transect_type == 4}
#' @slot specs_grid A \code{\linkS4class{GridInfo}} object.
#' @slot df_veg A \code{data.frame}. A (RAT) table of 'grid_veg' in which the \code{type@ids} contained.
#' @export
## Can extend and create with new("EcotonerGrids", ...). You can use EcotonerGrids() inside package, but not exported. See http://r-pkgs.had.co.nz/namespace.html
EcotonerGrids <- setClass("EcotonerGrids", 
							slots = list(grid_env = structure("RasterLayer", package = "raster"),
										grid_veg = structure("RasterLayer", package = "raster"),
										grid_flow = structure("RasterLayer", package = "raster"),
										grid_abut = structure("RasterLayer", package = "raster"),
										grid_aspect_mean = structure("RasterLayer", package = "raster"),
										grid_aspect_sd = structure("RasterLayer", package = "raster"),
										specs_grid = "GridInfo",
										df_veg = "data.frame"
									),
							prototype = list(df_veg = data.frame())
							)

setValidity("EcotonerGrids", function(object) {
	check_grid <- function(grid, specs) {
		if (!identical(grid, new(structure("RasterLayer", package = "raster")))) {
			test_res <- isTRUE(all.equal(raster::res(grid),  rep(specs@res_m, 2)))
			test_crs <- try(raster::compareCRS(raster::crs(grid), specs@crs), silent = TRUE)
			test_crs <- if (inherits(test_crs, "try-error")) FALSE else test_crs
			test_longlat <- identical(raster::isLonLat(grid), specs@longlat)
			test_origin <- isTRUE(all.equal(raster::origin(grid),  specs@origin))
			test_rotation <- identical(raster::rotated(grid), specs@rotation)
			
			if (all(test_res, test_longlat, test_crs, test_origin, test_rotation)) {
				TRUE		
			} else {
				paste("grid does not correspond in:", paste(c(
						if (!test_res) "cell resolution" else NULL,
						if (!test_longlat) "long/lat status" else NULL,
						if (!test_crs) "CRS" else NULL,
						if (!test_origin) "origin" else NULL,
						if (!test_rotation) "rotation" else NULL),
						sep = ", ")
					)
			}
		} else {
			TRUE
		}
	}
	
	grid_list <- list(grid_env = object@grid_env,
						grid_veg = object@grid_veg,
						grid_flow = object@grid_flow,
						grid_abut = object@grid_abut,
						grid_aspect_mean = object@grid_aspect_mean,
						grid_aspect_sd = object@grid_aspect_sd)
	
	vtests <- lapply(grid_list, check_grid, specs = object@specs_grid)
	id_good_grids <- sapply(vtests, function(x) is.logical(x) && x)
	
	if (all(id_good_grids)) {
		TRUE
	} else {
		invalid <- paste(names(vtests)[!id_good_grids], vtests[!id_good_grids]) # add grid name to invalid text
		paste(invalid, collapse = "\n")
	}
})


## Accessor generic function definitions
setGeneric("crs", signature = "x", function(x) standardGeneric("crs"))
setGeneric("res_m", signature = "x", function(x) standardGeneric("res_m"))
setGeneric("longlat", signature = "x", function(x) standardGeneric("longlat"))
setGeneric("origin", signature = "x", function(x) standardGeneric("origin"))
setGeneric("rotation", signature = "x", function(x) standardGeneric("rotation"))

#' @export
setGeneric("grid_env", signature = "x", function(x) standardGeneric("grid_env"))
#' @export
setGeneric("grid_veg", signature = "x", function(x) standardGeneric("grid_veg"))
#' @export
setGeneric("grid_flow", signature = "x", function(x) standardGeneric("grid_flow"))
#' @export
setGeneric("grid_abut", signature = "x", function(x) standardGeneric("grid_abut"))
#' @export
setGeneric("grid_aspect_mean", signature = "x", function(x) standardGeneric("grid_aspect_mean"))
#' @export
setGeneric("grid_aspect_sd", signature = "x", function(x) standardGeneric("grid_aspect_sd"))
#' @export
setGeneric("specs_grid", signature = "x", function(x) standardGeneric("specs_grid"))
#' @export
setGeneric("df_veg", signature = "x", function(x) standardGeneric("df_veg"))


## Accessor method definitions
setMethod("crs", "GridInfo", function(x) slot(x, "crs"))
setMethod("res_m", "GridInfo", function(x) slot(x, "res_m"))
setMethod("longlat", "GridInfo", function(x) slot(x, "longlat"))
setMethod("origin", "GridInfo", function(x) slot(x, "origin"))
setMethod("rotation", "GridInfo", function(x) slot(x, "rotation"))

setMethod("grid_env", "EcotonerGrids", function(x) slot(x, "grid_env"))
setMethod("grid_veg", "EcotonerGrids", function(x) slot(x, "grid_veg"))
setMethod("grid_flow", "EcotonerGrids", function(x) slot(x, "grid_flow"))
setMethod("grid_abut", "EcotonerGrids", function(x) slot(x, "grid_abut"))
setMethod("grid_aspect_mean", "EcotonerGrids", function(x) slot(x, "grid_aspect_mean"))
setMethod("grid_aspect_sd", "EcotonerGrids", function(x) slot(x, "grid_aspect_sd"))
setMethod("specs_grid", "EcotonerGrids", function(x) slot(x, "specs_grid"))
setMethod("df_veg", "EcotonerGrids", function(x) slot(x, "df_veg"))


## Replacer generic function definitions
setGeneric("crs<-", signature = "x", function(x, value) standardGeneric("crs<-"))
setGeneric("res_m<-", signature = "x", function(x, value) standardGeneric("res_m<-"))
setGeneric("longlat<-", signature = "x", function(x, value) standardGeneric("longlat<-"))
setGeneric("origin<-", signature = "x", function(x, value) standardGeneric("origin<-"))
setGeneric("rotation<-", signature = "x", function(x, value) standardGeneric("rotation<-"))

#' @export
setGeneric("grid_env<-", signature = "x", function(x, value) standardGeneric("grid_env<-"))
#' @export
setGeneric("grid_veg<-", signature = "x", function(x, value) standardGeneric("grid_veg<-"))
#' @export
setGeneric("grid_flow<-", signature = "x", function(x, value) standardGeneric("grid_flow<-"))
#' @export
setGeneric("grid_abut<-", signature = "x", function(x, value) standardGeneric("grid_abut<-"))
#' @export
setGeneric("grid_aspect_mean<-", signature = "x", function(x, value) standardGeneric("grid_aspect_mean<-"))
#' @export
setGeneric("grid_aspect_sd<-", signature = "x", function(x, value) standardGeneric("grid_aspect_sd<-"))
#' @export
setGeneric("specs_grid<-", signature = "x", function(x, value) standardGeneric("specs_grid<-"))
#' @export
setGeneric("df_veg<-", signature = "x", function(x, value) standardGeneric("df_veg<-"))


## Replacer method definitions
setReplaceMethod("crs", "GridInfo", function(x, value) initialize(x, crs = value))
setReplaceMethod("res_m", "GridInfo", function(x, value) initialize(x, res_m = value))
setReplaceMethod("longlat", "GridInfo", function(x, value) initialize(x, longlat = value))
setReplaceMethod("origin", "GridInfo", function(x, value) initialize(x, origin = value))
setReplaceMethod("rotation", "GridInfo", function(x, value) initialize(x, rotation = value))

set_grid <- function(x, name, value) {
	if (inherits(value, "RasterLayer")) {
		crs(specs_grid(x)) <- raster::crs(value)
		res_m(specs_grid(x)) <- mean(raster::res(value))
		longlat(specs_grid(x)) <- raster::isLonLat(value)
		origin(specs_grid(x)) <- raster::origin(value)
		rotation(specs_grid(x)) <- raster::rotated(value)
	}	
	slot(x, name) <- value
	validObject(x)
	x
}

setReplaceMethod("grid_env", "EcotonerGrids", function(x, value) set_grid(x, "grid_env", value))
setReplaceMethod("grid_veg", "EcotonerGrids", function(x, value) set_grid(x, "grid_veg", value))
setReplaceMethod("grid_flow", "EcotonerGrids", function(x, value) set_grid(x, "grid_flow", value))
setReplaceMethod("grid_abut", "EcotonerGrids", function(x, value) set_grid(x, "grid_abut", value))
setReplaceMethod("grid_aspect_mean", "EcotonerGrids", function(x, value) set_grid(x, "grid_aspect_mean", value))
setReplaceMethod("grid_aspect_sd", "EcotonerGrids", function(x, value) set_grid(x, "grid_aspect_sd", value))
setReplaceMethod("specs_grid", "EcotonerGrids", function(x, value) initialize(x, specs_grid = value))
setReplaceMethod("df_veg", "EcotonerGrids", function(x, value) initialize(x, df_veg = value))


## Generic function definitions
#' Checks if a raster grid has more than one cell.
#' 
#' @return A logical value.
#' @export
setGeneric("valid_grid", signature = "x", function(x) standardGeneric("valid_grid"))

## Method definitions
setMethod("valid_grid", structure("RasterLayer", package = "raster"), function(x) raster::ncell(x) > 1)
