## Class definitions

#' An S4-class describing the vegetation (type) of each side of the ecotone/ecological boundary
#'
#' @slot ids An integer (vector). The cell values of 'grid_veg' that belong to this vegetation (type).
#' @slot field A character string. The column name of 'df_veg' that was used to define this vegetation (type) [currently, not implemented].
#' @slot end_to_left A logical. If \code{TRUE}, then this vegetation (type) is on the 'left' side of the ecotone/ecological boundary.
#'
TypeInfo <- setClass("TypeInfo",
						slots = list(ids = "integer",
									field = "character",
									end_to_left = "logical"
								),
						prototype = list(ids = NA_integer_, field = NA_character_, end_to_left = NA)
					)


setClass("EcotonerPath", slots = list(path = "character"), prototype = list(path = NA_character_))

setValidity("EcotonerPath", function(object) {
	if (is.na(object@path) || file.exists(object@path)) {
		TRUE
	} else {
		TRUE
		##TODO(drs): how to determine if a path could be correct, but currently doesn't exist, e.g., volume is not mounted; file is not yet created
		# paste(object@path, "does not exist and/or is not a suitable path")
	}
})

setClass("EcotonerFile", contains = "EcotonerPath", prototype = list(path = NA_character_))

setValidity("EcotonerFile", function(object) {
	
	if (is.na(object@path)) {
		TRUE
	} else {
		types <- tolower(c("rds", "RData", "csv"))
		fextension <- tolower(tail(strsplit(basename(object@path), split = ".", fixed = TRUE)[[1]], n = 1))

		if (fextension %in% types) {
			TRUE
		} else {
			paste("The file", object@path, "is not of one of the accepted types:", paste(types, collapse = ", "))
		}
	}
})




#' An S4-class to represent the settings of an 'ecotoner' project
#'
#' @slot transect_type An integer. This flag selects the method to locate a transect, i.e., index to \code{get("transect_types", envir = etr_vars)}. See \code{\link{package?ecotoner}}. Default is \code{4} corresponding to the method by 'sufficient area of homogeneous aspect'.
#' @slot transect_type_desc A character string. A string describing the method selected by \code{transect_type}. One of \code{get("transect_types", envir = etr_vars)}.
#' @slot candidate_THAs_N An integer. The number of candidate transects that should be created if \code{transect_type == 4}. Default is \code{20}.
#' @slot transect_azimuth A numeric. The orientation of the environmental gradient within a transect. Defaults to \code{3 / 2 * pi}, i.e., "left is low & right is high".
#' @slot cores_N An integer. This project is attempted to be run on the local machine in parallel with \code{cores_N} cores, but not with more than the actual number of cores of this machine.
#' @slot rng_seed0 An integer. This is the seed to set the random number generator if used.
#' @slot rng_seed1 An integer or \code{NULL}. The 'serial' seed. This is the seed used to set the random number generator before parallel streams are initiated. If \code{reproducible} is \code{TRUE}, rng_seed1 corresponds to rng_seed0 else to \code{NULL}.
#' @slot rng_seed2 An integer vector or \code{NULL}. The 'parallel' seed. This slot captures the state of the random number generator after parallel streams of the kind of 'L'Ecuyer-CMRG' have been set up. This state will be used to re-set the generator if \code{reproducible} is \code{TRUE}.
#' @slot reproducible A logical. If \code{TRUE} then this flag sets the seed of the random number generator at the start.
#' @slot reseed A logical. If \code{TRUE} then this flag re-sets the random number generator each time a set of random numbers is generated.
#' @slot transect_N An integer. The number of search points to generated to initiate a transect search in the neighborhood of size(s) \code{neighborhoods}.
#' @slot inhibit_searchpoints A logical. If \code{TRUE} then search points to initiate a transect search are sampled by a random simple sequential inhibition process; if \code{FALSE}, then a random (stratified) process is used to sample search points.
#' @slot neighborhoods An integer vector. The number(s) of cells per side of a square representing local grid window(s) around a search point within which the transect methods searches for a suitable transect. Neighborhood(s) must be odd. Note: Below a neighborhood of 333 cells, the density function may become unreliable.
#' @slot neighbors_N An integer. The number of \code{neighborhoods}.
#' @slot stepsHullEdge An integer vector. The step length(s), in grid cells, of the hull edge method by Gastner et al.
#' @slot clumpAreaClasses_m2 A numeric vector of length \code{2}. The first number is the limit between small and medium patches, and the second number is the limit between medium and large patches.
#' @slot res_rotate_grid An integer. This number is used in the call to extract_tband_grids(); it represents variable 'fact' in raster::disaggregate(). The rotation of the grid from geographical to transect orientation requires a high resolution in cells because not all resulting cells would otherwise be represented well after the rotation. Instead, they would need to be interpolated which is tricky with categorical data such as vegetation cover. Hence, the results, e.g., patches, flowpaths, are relatively sensitive to this factor. \code{2} produced only a coarse representation of the rotated grids; \code{4} improved the resolution after rotation; \code{6} and \code{8} produced mostly similar results, but only \code{8} was satisfactorily. Default is thus \code{8}.
#' @slot Aspect_SDof201Mean_Maximum A numeric. If \code{transect_type == 4}, then this option limits the search area for line-starting points to areas of low standard variation of the mean aspect; default is \code{75 * pi / 180 = 1.31}
#' @slot Aspect_DeviationAlongCandidateLines_Maximum A numeric. If \code{transect_type == 4}, then this option sets the limit by which the mean aspect along candidate lines may differ from the mean aspect at the starting point differs; default is \code{30 * pi / 180 = 0.52}
#' @slot vegDefiningDensityTransect_low A numeric. Determines the start of an ecotone transect by setting a critical vegetation density which should not be reached during the first run of sufficiently low values; default is \code{0.05}
#' @slot vegDefiningDensityTransect_high A numeric. Determines the start of an ecotone transect by setting a critical vegetation density which should be reached during the first run of sufficiently high values; default is \code{0.75}
#' @slot vegDefiningDensityTransectExtended_min A numeric. Determines the location of an ecotone transect by setting a critical vegetation density which should be maintained; default is \code{0.5}
#' @slot bandTransect_width_cellN An integer. The width, in grid cell numbers, of the extracted band transect, half of which is on either side of the line transect. If \code{transect_type == 4}, then this number must match the smoothing window width of \code{\linkS4class{EcotonerGrids@grid_aspect_mean}} and \code{\linkS4class{EcotonerGrids@grid_aspect_sd}}.
#' @slot min_slope_with_aspect A numeric. The minimum gradient, above which an aspect can be meaningfully measured/calculated. If \code{transect_type == 4}, then this number must match the used during production of \code{\linkS4class{EcotonerGrids@grid_aspect_mean}} and \code{\linkS4class{EcotonerGrids@grid_aspect_sd}}. Default is \code{2 * pi / 180 = 0.035} based for elevational slopes.

#' @slot type_excl A \code{\linkS4class{TypeInfo}} object. This 'type' is to be excluded from 'grid_veg'.
#' @slot type_veg1 A \code{\linkS4class{TypeInfo}} object. This 'type' is one of the two vegetation (types) from 'grid_veg'.
#' @slot type_veg2 A \code{\linkS4class{TypeInfo}} object. This 'type' is one of the two vegetation (types) from 'grid_veg'.

#' @slot dir_prj An 'EcotonerPath' object, i.e., a character string representing a path. The path to the project directory. 
#' @slot dir_big An 'EcotonerPath' object or \code{NA}. If not \code{NA}, then the path where the transect data will be stored. 
#' @slot dir_init An 'EcotonerPath' object. The path to the folder where the search points and other data are stored before transects can be searched for and extracted. This path is automatically generated.
#' @slot dir_out An 'EcotonerPath' object. The path that contains the outputs.  This path is automatically generated from dir_prj or from dir_big if not NA.
#' @slot dir_out_fig An 'EcotonerPath' object. The path that contains the folders where the transect figures are stored.  This path is automatically generated.
#' @slot dir_out_dat An 'EcotonerPath' object. The path that contains the folders where the transect data are stored.  This path is automatically generated.
#' @slot dir_env An 'EcotonerPath' object. The path to where the raster grid representing the environmental gradient 'grid_env' is stored. 
#' @slot dir_veg An 'EcotonerPath' object. The path to where the raster grid representing the vegetation 'grid_veg' is stored. 
#' @slot dir_flow An 'EcotonerPath' object. The path to where the raster grid representing the flow paths 'grid_flow' is stored (required with transect_type == 3). 
#' @slot dir_abut An 'EcotonerPath' object. The path to where the raster grid representing the abutting cells of vegetation types 1 and 2 'grid_abut' is stored (required with transect_type == 4). 
#' @slot dir_aspect_mean An 'EcotonerPath' object. The path to where the raster grid representing the mean smoothed aspect 'grid_aspect_mean' is stored (required with transect_type == 4). 
#' @slot dir_aspect_sd An 'EcotonerPath' object. The path to where the raster grid representing the standard deviation of the smoothed aspect 'grid_aspect_sd' is stored (required with transect_type == 4). 
#' @slot file_searchpoints An 'EcotonerFile' object. The entire file path to where the search points are/will be stored.
#' @slot file_initwindow An 'EcotonerFile' object. The entire file path to where the spatstat::owin object for the search point generation is/will be stored.
#' @slot file_etsummary An 'EcotonerFile' object. The entire file path to where the table of the transect summary are stored.
#' @slot file_etsummary_temp An 'EcotonerFile' object. The entire file path to where temporary table of the transect summary are stored.
#' @slot file_etmeasure_base An 'EcotonerFile' object. The entire file path (minus a flag indicating which measure type) to where temporary tables of the transect measures are stored.
#'
#' @export
## Can extend and create with new("EcotonerSettings", ...). You can use EcotonerSettings() inside package, but not exported. See http://r-pkgs.had.co.nz/namespace.html
EcotonerSettings <- setClass("EcotonerSettings", 
							slots = list(transect_type = "integer",
										 transect_type_desc = "character",
										 
										 candidate_THAs_N = "integer",
										 transect_azimuth = "numeric",
										 cores_N = "integer",
										 rng_seed0 = "integer",
										 rng_seed1 = "ANY",
										 rng_seed2 = "ANY",
										 reproducible = "logical",
										 reseed = "logical",
										 transect_N = "integer",
										 inhibit_searchpoints = "logical",
										 neighborhoods = "integer",
										 neighbors_N = "integer",
										 stepsHullEdge = "integer",
										 clumpAreaClasses_m2 = "numeric",
										 
										 res_rotate_grid = "integer",
										 Aspect_SDof201Mean_Maximum = "numeric",
										 Aspect_DeviationAlongCandidateLines_Maximum = "numeric",

										 vegDefiningDensityTransect_low = "numeric",
										 vegDefiningDensityTransect_high = "numeric",
										 vegDefiningDensityTransectExtended_min = "numeric",
										 
										 bandTransect_width_cellN = "integer",
										 min_slope_with_aspect = "numeric",
										
										 type_excl = "TypeInfo",
										 type_veg1 = "TypeInfo",
										 type_veg2 = "TypeInfo",
										 
										 dir_prj = "EcotonerPath",
										 dir_big = "EcotonerPath",
										 dir_init = "EcotonerPath",
										 dir_out = "EcotonerPath",
										 dir_out_fig = "EcotonerPath",
										 dir_out_dat = "EcotonerPath",
										 
										 dir_env = "EcotonerPath",
										 dir_veg = "EcotonerPath",
										 dir_flow = "EcotonerPath",
										 dir_abut = "EcotonerPath",
										 dir_aspect_mean = "EcotonerPath",
										 dir_aspect_sd = "EcotonerPath",
										 										 
										 file_searchpoints = "EcotonerFile",
										 file_initwindow = "EcotonerFile",
										 file_etsummary = "EcotonerFile",
										 file_etsummary_temp = "EcotonerFile",
										 file_etmeasure_base = "EcotonerFile"
										),
							prototype = list(transect_type = 4L,
											 transect_type_desc = "HomogeneousAspect",
											 vegDefiningDensityTransect_low = 0.05,
											 vegDefiningDensityTransect_high = 0.75,
											 vegDefiningDensityTransectExtended_min = 0.5,
											 Aspect_SDof201Mean_Maximum = 75 * pi/180,
											 Aspect_DeviationAlongCandidateLines_Maximum = 30 * pi/180,
											 bandTransect_width_cellN = 200L,
											 res_rotate_grid = 8L,
											 min_slope_with_aspect = 2 * pi / 180,
											 transect_azimuth = 3 / 2 * pi,
											 candidate_THAs_N = 20L,
											 cores_N = 1L,

											 rng_seed0 = 123L,
											 reproducible = TRUE,
											 reseed = FALSE,
											 transect_N = 30L,
											 inhibit_searchpoints = FALSE,
											 neighborhoods = c(1667L, 999L),
											 neighbors_N = 2L,
											 stepsHullEdge = c(1L, 3L),
											 clumpAreaClasses_m2 = c(1e4, 1e6)
										)
							)


setValidity("EcotonerSettings", function(object) {
	
	tests <- list()
	tests[["transect_type"]] <- all(length(object@transect_type) == 1,
									object@transect_type %in% seq_along(get("transect_types", envir = etr_vars)))
	tests[["transect_type_desc"]] <- all(length(object@transect_type_desc) == 1,
										object@transect_type_desc %in% get("transect_types", envir = etr_vars))
	tests[["candidate_THAs_N"]] <- all(length(object@candidate_THAs_N) == 1,
										object@candidate_THAs_N > 0)
	tests[["transect_azimuth"]] <- all(length(object@transect_azimuth) == 1,
										object@transect_azimuth >= 0, object@transect_azimuth < (2 * pi))
	tests[["cores_N"]] <- all(length(object@cores_N) == 1,
							object@cores_N %in% (seq_len(parallel::detectCores()) - 1))
	tests[["rng_seed0"]] <- all(length(object@rng_seed0) == 1,
							object@rng_seed0 > 0)
	tests[["transect_N"]] <- all(length(object@transect_N) == 1,
									object@transect_N > 0)
	tests[["neighborhoods"]] <- all(object@neighborhoods > 0,
									is_odd(object@neighborhoods))
	tests[["neighbors_N"]] <- all(length(object@neighbors_N) == 1,
									object@neighbors_N == length(object@neighborhoods))
	tests[["stepsHullEdge"]] <- all(object@stepsHullEdge > 0, sort(object@stepsHullEdge) == sort(unique(object@stepsHullEdge)))
	tests[["clumpAreaClasses_m2"]] <- all(length(object@clumpAreaClasses_m2) == 2,
											object@clumpAreaClasses_m2 > 0,
											object@clumpAreaClasses_m2[1] < object@clumpAreaClasses_m2[2])
	tests[["res_rotate_grid"]] <- all(length(object@res_rotate_grid) == 1,
										object@res_rotate_grid > 0)
	tests[["Aspect_SDof201Mean_Maximum"]] <- all(length(object@Aspect_SDof201Mean_Maximum) == 1,
												object@Aspect_SDof201Mean_Maximum > 0)
	tests[["Aspect_DeviationAlongCandidateLines_Maximum"]] <- all(length(object@Aspect_DeviationAlongCandidateLines_Maximum) == 1,
																	object@Aspect_DeviationAlongCandidateLines_Maximum > 0)
	tests[["vegDefiningDensityTransect_low"]] <- all(length(object@vegDefiningDensityTransect_low) == 1,
													object@vegDefiningDensityTransect_low > 0,
													object@vegDefiningDensityTransect_low < 1,
													object@vegDefiningDensityTransect_low < object@vegDefiningDensityTransect_high)
	tests[["vegDefiningDensityTransect_high"]] <- all(length(object@vegDefiningDensityTransect_high) == 1,
														object@vegDefiningDensityTransect_high > 0,
														object@vegDefiningDensityTransect_high < 1,
														object@vegDefiningDensityTransect_high > object@vegDefiningDensityTransect_low)
	tests[["vegDefiningDensityTransectExtended_min"]] <- all(length(object@vegDefiningDensityTransectExtended_min) == 1,
																object@vegDefiningDensityTransectExtended_min > 0,
																object@vegDefiningDensityTransectExtended_min < 1,
																object@vegDefiningDensityTransectExtended_min < object@vegDefiningDensityTransect_high,
																object@vegDefiningDensityTransectExtended_min > object@vegDefiningDensityTransect_low)
	tests[["bandTransect_width_cellN"]] <- all(length(object@bandTransect_width_cellN) == 1, object@bandTransect_width_cellN >= 200L)
	tests[["min_slope_with_aspect"]] <- all(length(object@min_slope_with_aspect) == 1, object@min_slope_with_aspect > 0)
	
	ibad <- !unlist(tests)
	if (any(ibad)) {
		paste("The following slot(s) is/are set incorrectly:", names(ibad)[ibad], collapse = ", ")
	} else {
		TRUE
	}

})


## Generic definitions
#' @export
setGeneric("transect_type", signature = "x", function(x) standardGeneric("transect_type"))
#' @export
setGeneric("transect_type<-", signature = "x", function(x, value) standardGeneric("transect_type<-"))
#' @export
setGeneric("transect_type_desc", signature = "x", function(x) standardGeneric("transect_type_desc"))

#' @export
setGeneric("candidate_THAs_N", signature = "x", function(x) standardGeneric("candidate_THAs_N"))
#' @export
setGeneric("candidate_THAs_N<-", signature = "x", function(x, value) standardGeneric("candidate_THAs_N<-"))

#' @export
setGeneric("transect_azimuth", signature = "x", function(x) standardGeneric("transect_azimuth"))
#' @export
setGeneric("transect_azimuth<-", signature = "x", function(x, value) standardGeneric("transect_azimuth<-"))

#' @export
setGeneric("cores_N", signature = "x", function(x) standardGeneric("cores_N"))
#' @export
setGeneric("cores_N<-", signature = "x", function(x, value) standardGeneric("cores_N<-"))

setGeneric("init_sseed", signature = "x", function(x) standardGeneric("init_sseed"))
#' @export
setGeneric("init_pseed", signature = "x", function(x) standardGeneric("init_pseed"))
#' @export
setGeneric("get_seed", signature = "x", function(x) standardGeneric("get_seed"))
#' @export
setGeneric("set_seed<-", signature = "x", function(x, value) standardGeneric("set_seed<-"))
#' @export
setGeneric("get_sseed", signature = "x", function(x) standardGeneric("get_sseed"))
#' @export
setGeneric("get_pseed", signature = "x", function(x) standardGeneric("get_pseed"))

#' @export
setGeneric("reproducible", signature = "x", function(x) standardGeneric("reproducible"))
#' @export
setGeneric("reproducible<-", signature = "x", function(x, value) standardGeneric("reproducible<-"))

#' @export
setGeneric("reseed", signature = "x", function(x) standardGeneric("reseed"))
#' @export
setGeneric("reseed<-", signature = "x", function(x, value) standardGeneric("reseed<-"))

#' @export
setGeneric("transect_N", signature = "x", function(x) standardGeneric("transect_N"))
#' @export
setGeneric("transect_N<-", signature = "x", function(x, value) standardGeneric("transect_N<-"))

#' @export
setGeneric("inhibit_searchpoints", signature = "x", function(x) standardGeneric("inhibit_searchpoints"))
#' @export
setGeneric("inhibit_searchpoints<-", signature = "x", function(x, value) standardGeneric("inhibit_searchpoints<-"))

#' @export
setGeneric("neighborhoods", signature = "x", function(x) standardGeneric("neighborhoods"))
#' @export
setGeneric("neighborhoods<-", signature = "x", function(x, value) standardGeneric("neighborhoods<-"))
#' @export
setGeneric("neighbors_N", signature = "x", function(x) standardGeneric("neighbors_N"))

#' @export
setGeneric("stepsHullEdge", signature = "x", function(x) standardGeneric("stepsHullEdge"))
#' @export
setGeneric("stepsHullEdge<-", signature = "x", function(x, value) standardGeneric("stepsHullEdge<-"))

#' @export
setGeneric("clumpAreaClasses_m2", signature = "x", function(x) standardGeneric("clumpAreaClasses_m2"))
#' @export
setGeneric("clumpAreaClasses_m2<-", signature = "x", function(x, value) standardGeneric("clumpAreaClasses_m2<-"))

#' @export
setGeneric("res_rotate_grid", signature = "x", function(x) standardGeneric("res_rotate_grid"))
#' @export
setGeneric("res_rotate_grid<-", signature = "x", function(x, value) standardGeneric("res_rotate_grid<-"))

#' @export
setGeneric("Aspect_SDof201Mean_Maximum", signature = "x", function(x) standardGeneric("Aspect_SDof201Mean_Maximum"))
#' @export
setGeneric("Aspect_SDof201Mean_Maximum<-", signature = "x", function(x, value) standardGeneric("Aspect_SDof201Mean_Maximum<-"))

#' @export
setGeneric("Aspect_DeviationAlongCandidateLines_Maximum", signature = "x", function(x) standardGeneric("Aspect_DeviationAlongCandidateLines_Maximum"))
#' @export
setGeneric("Aspect_DeviationAlongCandidateLines_Maximum<-", signature = "x", function(x, value) standardGeneric("Aspect_DeviationAlongCandidateLines_Maximum<-"))

#' @export
setGeneric("vegDefiningDensityTransect_low", signature = "x", function(x) standardGeneric("vegDefiningDensityTransect_low"))
#' @export
setGeneric("vegDefiningDensityTransect_low<-", signature = "x", function(x, value) standardGeneric("vegDefiningDensityTransect_low<-"))

#' @export
setGeneric("vegDefiningDensityTransect_high", signature = "x", function(x) standardGeneric("vegDefiningDensityTransect_high"))
#' @export
setGeneric("vegDefiningDensityTransect_high<-", signature = "x", function(x, value) standardGeneric("vegDefiningDensityTransect_high<-"))

#' @export
setGeneric("vegDefiningDensityTransectExtended_min", signature = "x", function(x) standardGeneric("vegDefiningDensityTransectExtended_min"))
#' @export
setGeneric("vegDefiningDensityTransectExtended_min<-", signature = "x", function(x, value) standardGeneric("vegDefiningDensityTransectExtended_min<-"))

#' @export
setGeneric("bandTransect_width_cellN", signature = "x", function(x) standardGeneric("bandTransect_width_cellN"))
#' @export
setGeneric("bandTransect_width_cellN<-", signature = "x", function(x, value) standardGeneric("bandTransect_width_cellN<-"))

#' @export
setGeneric("min_slope_with_aspect", signature = "x", function(x) standardGeneric("min_slope_with_aspect"))
#' @export
setGeneric("min_slope_with_aspect<-", signature = "x", function(x, value) standardGeneric("min_slope_with_aspect<-"))

#' @export
setGeneric("type_ids", signature = "x", function(x) standardGeneric("type_ids"))
#' @export
setGeneric("type_ids<-", signature = "x", function(x, value) standardGeneric("type_ids<-"))
#' @export
setGeneric("type_field", signature = "x", function(x) standardGeneric("type_field"))
#' @export
setGeneric("type_field<-", signature = "x", function(x, value) standardGeneric("type_field<-"))
#' @export
setGeneric("end_to_left", signature = "x", function(x) standardGeneric("end_to_left"))
#' @export
setGeneric("end_to_left<-", signature = "x", function(x, value) standardGeneric("end_to_left<-"))
#' @export
setGeneric("type_excl", signature = "x", function(x) standardGeneric("type_excl"))
#' @export
setGeneric("type_excl<-", signature = "x", function(x, value) standardGeneric("type_excl<-"))
#' @export
setGeneric("type_veg1", signature = "x", function(x) standardGeneric("type_veg1"))
#' @export
setGeneric("type_veg1<-", signature = "x", function(x, value) standardGeneric("type_veg1<-"))
#' @export
setGeneric("type_veg2", signature = "x", function(x) standardGeneric("type_veg2"))
#' @export
setGeneric("type_veg2<-", signature = "x", function(x, value) standardGeneric("type_veg2<-"))



#' @export
setGeneric("dir_prj", signature = "x", function(x) standardGeneric("dir_prj"))
#' @export
setGeneric("dir_prj<-", signature = "x", function(x, value) standardGeneric("dir_prj<-"))
#' @export
setGeneric("dir_init", signature = "x", function(x) standardGeneric("dir_init"))
#' @export
setGeneric("dir_out", signature = "x", function(x) standardGeneric("dir_out"))
#' @export
setGeneric("dir_out_dat", signature = "x", function(x) standardGeneric("dir_out_dat"))
#' @export
setGeneric("dir_out_fig", signature = "x", function(x) standardGeneric("dir_out_fig"))
#' @export
setGeneric("dir_big", signature = "x", function(x) standardGeneric("dir_big"))
#' @export
setGeneric("dir_big<-", signature = "x", function(x, value) standardGeneric("dir_big<-"))
#' @export
setGeneric("dir_env", signature = "x", function(x) standardGeneric("dir_env"))
#' @export
setGeneric("dir_env<-", signature = "x", function(x, value) standardGeneric("dir_env<-"))
#' @export
setGeneric("dir_veg", signature = "x", function(x) standardGeneric("dir_veg"))
#' @export
setGeneric("dir_veg<-", signature = "x", function(x, value) standardGeneric("dir_veg<-"))
#' @export
setGeneric("dir_flow", signature = "x", function(x) standardGeneric("dir_flow"))
#' @export
setGeneric("dir_flow<-", signature = "x", function(x, value) standardGeneric("dir_flow<-"))
#' @export
setGeneric("dir_abut", signature = "x", function(x) standardGeneric("dir_abut"))
#' @export
setGeneric("dir_abut<-", signature = "x", function(x, value) standardGeneric("dir_abut<-"))
#' @export
setGeneric("dir_aspect_mean", signature = "x", function(x) standardGeneric("dir_aspect_mean"))
#' @export
setGeneric("dir_aspect_mean<-", signature = "x", function(x, value) standardGeneric("dir_aspect_mean<-"))
#' @export
setGeneric("dir_aspect_sd", signature = "x", function(x) standardGeneric("dir_aspect_sd"))
#' @export
setGeneric("dir_aspect_sd<-", signature = "x", function(x, value) standardGeneric("dir_aspect_sd<-"))
#' @export
setGeneric("file_searchpoints", signature = "x", function(x) standardGeneric("file_searchpoints"))
#' @export
setGeneric("file_searchpoints<-", signature = "x", function(x, value) standardGeneric("file_searchpoints<-"))
#' @export
setGeneric("file_initwindow", signature = "x", function(x) standardGeneric("file_initwindow"))
#' @export
setGeneric("file_initwindow<-", signature = "x", function(x, value) standardGeneric("file_initwindow<-"))
#' @export
setGeneric("file_etsummary", signature = "x", function(x) standardGeneric("file_etsummary"))
#' @export
setGeneric("file_etsummary<-", signature = "x", function(x, value) standardGeneric("file_etsummary<-"))
#' @export
setGeneric("file_etsummary_temp", signature = "x", function(x) standardGeneric("file_etsummary_temp"))
#' @export
setGeneric("file_etmeasure_base<-", signature = "x", function(x, value) standardGeneric("file_etmeasure_base<-"))
#' @export
setGeneric("file_etmeasure_base", signature = "x", function(x, value) standardGeneric("file_etmeasure_base"))



## Method definitions
setMethod("transect_type", "EcotonerSettings", function(x) slot(x, "transect_type"))
setReplaceMethod("transect_type", "EcotonerSettings", function(x, value) {x@transect_type <- as.integer(value); x@transect_type_desc <- get("transect_types", envir = etr_vars)[x@transect_type]; x})
setMethod("transect_type_desc", "EcotonerSettings", function(x) slot(x, "transect_type_desc"))

setMethod("candidate_THAs_N", "EcotonerSettings", function(x) slot(x, "candidate_THAs_N"))
setReplaceMethod("candidate_THAs_N", "EcotonerSettings", function(x, value) initialize(x, candidate_THAs_N = as.integer(value)))

setMethod("transect_azimuth", "EcotonerSettings", function(x) slot(x, "transect_azimuth"))
setReplaceMethod("transect_azimuth", "EcotonerSettings", function(x, value) initialize(x, transect_azimuth = value %% (2 * pi)))

setMethod("cores_N", "EcotonerSettings", function(x) slot(x, "cores_N"))
setReplaceMethod("cores_N", "EcotonerSettings", function(x, value) initialize(x, cores_N = as.integer(value)))

setMethod("init_sseed", "EcotonerSettings", function(x) initialize(x, rng_seed1 = if (x@reproducible) x@rng_seed0 else NULL))
setMethod("init_pseed", "EcotonerSettings", function(x) initialize(init_sseed(x), rng_seed2 = if (is.null(x@rng_seed1)) NULL else .Random.seed))
setMethod("get_seed", "EcotonerSettings", function(x) slot(x, "rng_seed0"))
setReplaceMethod("set_seed", "EcotonerSettings", function(x, value) init_sseed(initialize(x, rng_seed0 = as.integer(value))))
setMethod("get_sseed", "EcotonerSettings", function(x) slot(x, "rng_seed1"))
setMethod("get_pseed", "EcotonerSettings", function(x) slot(x, "rng_seed2"))

setMethod("reproducible", "EcotonerSettings", function(x) slot(x, "reproducible"))
setReplaceMethod("reproducible", "EcotonerSettings", function(x, value) init_pseed(init_sseed(initialize(x, reproducible = value))))

setMethod("reseed", "EcotonerSettings", function(x) slot(x, "reseed"))
setReplaceMethod("reseed", "EcotonerSettings", function(x, value) initialize(x, reseed = value))

setMethod("transect_N", "EcotonerSettings", function(x) slot(x, "transect_N"))
setReplaceMethod("transect_N", "EcotonerSettings", function(x, value) initialize(x, transect_N = as.integer(value)))

setMethod("inhibit_searchpoints", "EcotonerSettings", function(x) slot(x, "inhibit_searchpoints"))
setReplaceMethod("inhibit_searchpoints", "EcotonerSettings", function(x, value) initialize(x, inhibit_searchpoints = value))

setMethod("neighborhoods", "EcotonerSettings", function(x) slot(x, "neighborhoods"))
setReplaceMethod("neighborhoods", "EcotonerSettings", function(x, value) {x@neighborhoods <- as.integer(value); x@neighbors_N <- length(x@neighborhoods); x})
setMethod("neighbors_N", "EcotonerSettings", function(x) slot(x, "neighbors_N"))

setMethod("stepsHullEdge", "EcotonerSettings", function(x) slot(x, "stepsHullEdge"))
setReplaceMethod("stepsHullEdge", "EcotonerSettings", function(x, value) initialize(x, stepsHullEdge = sort(unique(as.integer(value)))))

setMethod("clumpAreaClasses_m2", "EcotonerSettings", function(x) slot(x, "clumpAreaClasses_m2"))
setReplaceMethod("clumpAreaClasses_m2", "EcotonerSettings", function(x, value) initialize(x, clumpAreaClasses_m2 = value))

setMethod("res_rotate_grid", "EcotonerSettings", function(x) slot(x, "res_rotate_grid"))
setReplaceMethod("res_rotate_grid", "EcotonerSettings", function(x, value) initialize(x, res_rotate_grid = as.integer(value)))

setMethod("Aspect_SDof201Mean_Maximum", "EcotonerSettings", function(x) slot(x, "Aspect_SDof201Mean_Maximum"))
setReplaceMethod("Aspect_SDof201Mean_Maximum", "EcotonerSettings", function(x, value) initialize(x, Aspect_SDof201Mean_Maximum = value))

setMethod("Aspect_DeviationAlongCandidateLines_Maximum", "EcotonerSettings", function(x) slot(x, "Aspect_DeviationAlongCandidateLines_Maximum"))
setReplaceMethod("Aspect_DeviationAlongCandidateLines_Maximum", "EcotonerSettings", function(x, value) initialize(x, Aspect_DeviationAlongCandidateLines_Maximum = value))

setMethod("vegDefiningDensityTransect_low", "EcotonerSettings", function(x) slot(x, "vegDefiningDensityTransect_low"))
setReplaceMethod("vegDefiningDensityTransect_low", "EcotonerSettings", function(x, value) initialize(x, vegDefiningDensityTransect_low = value))

setMethod("vegDefiningDensityTransect_high", "EcotonerSettings", function(x) slot(x, "vegDefiningDensityTransect_high"))
setReplaceMethod("vegDefiningDensityTransect_high", "EcotonerSettings", function(x, value) initialize(x, vegDefiningDensityTransect_high = value))

setMethod("vegDefiningDensityTransectExtended_min", "EcotonerSettings", function(x) slot(x, "vegDefiningDensityTransectExtended_min"))
setReplaceMethod("vegDefiningDensityTransectExtended_min", "EcotonerSettings", function(x, value) initialize(x, vegDefiningDensityTransectExtended_min = value))

setMethod("bandTransect_width_cellN", "EcotonerSettings", function(x) slot(x, "bandTransect_width_cellN"))
setReplaceMethod("bandTransect_width_cellN", "EcotonerSettings", function(x, value) initialize(x, bandTransect_width_cellN = as.integer(value)))

setMethod("min_slope_with_aspect", "EcotonerSettings", function(x) slot(x, "min_slope_with_aspect"))
setReplaceMethod("min_slope_with_aspect", "EcotonerSettings", function(x, value) initialize(x, min_slope_with_aspect = value))


setMethod("type_ids", "TypeInfo", function(x) slot(x, "ids"))
setReplaceMethod("type_ids", "TypeInfo", function(x, value) initialize(x, ids = as.integer(value)))

setMethod("type_field", "TypeInfo", function(x) slot(x, "field"))
setReplaceMethod("type_field", "TypeInfo", function(x, value) initialize(x, field = value))

setMethod("end_to_left", "TypeInfo", function(x) slot(x, "end_to_left"))
setReplaceMethod("end_to_left", "TypeInfo", function(x, value) initialize(x, end_to_left = as.logical(value)))

setMethod("type_excl", "EcotonerSettings", function(x) slot(x, "type_excl"))
setReplaceMethod("type_excl", "EcotonerSettings", function(x, value) initialize(x, type_excl = value))

setMethod("type_veg1", "EcotonerSettings", function(x) slot(x, "type_veg1"))
setReplaceMethod("type_veg1", "EcotonerSettings", function(x, value) initialize(x, type_veg1 = value))

setMethod("type_veg2", "EcotonerSettings", function(x) slot(x, "type_veg2"))
setReplaceMethod("type_veg2", "EcotonerSettings", function(x, value) initialize(x, type_veg2 = value))



setMethod("dir_prj", "EcotonerSettings", function(x) slot(slot(x, "dir_prj"), "path"))
setReplaceMethod("dir_prj", "EcotonerSettings", function(x, value) verify_project_paths(initialize(x, dir_prj = new("EcotonerPath", path = value))))
setMethod("dir_init", "EcotonerSettings", function(x) slot(slot(x, "dir_init"), "path"))
setMethod("dir_out", "EcotonerSettings", function(x) slot(slot(x, "dir_out"), "path"))
setMethod("dir_out_dat", "EcotonerSettings", function(x) slot(slot(x, "dir_out_dat"), "path"))
setMethod("dir_out_fig", "EcotonerSettings", function(x) slot(slot(x, "dir_out_fig"), "path"))
setMethod("dir_big", "EcotonerSettings", function(x) slot(slot(x, "dir_big"), "path"))
setReplaceMethod("dir_big", "EcotonerSettings", function(x, value) verify_project_paths(initialize(x, dir_big = new("EcotonerPath", path = value))))
setMethod("dir_env", "EcotonerSettings", function(x) slot(slot(x, "dir_env"), "path"))
setReplaceMethod("dir_env", "EcotonerSettings", function(x, value) initialize(x, dir_env = new("EcotonerPath", path = value)))
setMethod("dir_veg", "EcotonerSettings", function(x) slot(slot(x, "dir_veg"), "path"))
setReplaceMethod("dir_veg", "EcotonerSettings", function(x, value) initialize(x, dir_veg = new("EcotonerPath", path = value)))
setMethod("dir_flow", "EcotonerSettings", function(x) slot(slot(x, "dir_flow"), "path"))
setReplaceMethod("dir_flow", "EcotonerSettings", function(x, value) initialize(x, dir_flow = new("EcotonerPath", path = value)))
setMethod("dir_abut", "EcotonerSettings", function(x) slot(slot(x, "dir_abut"), "path"))
setReplaceMethod("dir_abut", "EcotonerSettings", function(x, value) initialize(x, dir_abut = new("EcotonerPath", path = value)))
setMethod("dir_aspect_mean", "EcotonerSettings", function(x) slot(slot(x, "dir_aspect_mean"), "path"))
setReplaceMethod("dir_aspect_mean", "EcotonerSettings", function(x, value) initialize(x, dir_aspect_mean = new("EcotonerPath", path = value)))
setMethod("dir_aspect_sd", "EcotonerSettings", function(x) slot(slot(x, "dir_aspect_sd"), "path"))
setReplaceMethod("dir_aspect_sd", "EcotonerSettings", function(x, value) initialize(x, dir_aspect_sd = new("EcotonerPath", path = value)))

setMethod("file_searchpoints", "EcotonerSettings", function(x) slot(slot(x, "file_searchpoints"), "path"))
setReplaceMethod("file_searchpoints", "EcotonerSettings", function(x, value) initialize(x, file_searchpoints = new("EcotonerFile", path = file.path(dir_init(x), basename(value)))))
setMethod("file_initwindow", "EcotonerSettings", function(x) slot(slot(x, "file_initwindow"), "path"))
setReplaceMethod("file_initwindow", "EcotonerSettings", function(x, value) initialize(x, file_initwindow = new("EcotonerFile", path = value)))
setMethod("file_etsummary", "EcotonerSettings", function(x) slot(slot(x, "file_etsummary"), "path"))
setReplaceMethod("file_etsummary", "EcotonerSettings", function(x, value) {
			x <- initialize(x, file_etsummary = new("EcotonerFile", path = file.path(dir_out(x), basename(value))))

			temp <- strsplit(basename(file_etsummary(x)), split = ".", fixed = TRUE)[[1]]
			temp <- paste0(paste(head(temp, n = -1), collapse = ""), "_temp.", tail(temp, n = 1))
			initialize(x, file_etsummary_temp = new("EcotonerFile", path = file.path(dir_out(x), temp)))
		})
setMethod("file_etsummary_temp", "EcotonerSettings", function(x) slot(slot(x, "file_etsummary_temp"), "path"))
setReplaceMethod("file_etmeasure_base", "EcotonerSettings", function(x, value) initialize(x, file_etmeasure_base = new("EcotonerFile", path = file.path(dir_out(x), basename(value)))))
setMethod("file_etmeasure_base", "EcotonerSettings", function(x, value) {
			x <- slot(slot(x, "file_etmeasure_base"), "path")
			temp <- strsplit(basename(x), split = ".", fixed = TRUE)[[1]]
			temp <- paste0(paste(head(temp, n = -1), collapse = ""), "_", value, ".", tail(temp, n = 1))
			file.path(dirname(x), temp)
		})


## Paths
#' @export
setGeneric("verify_project_paths", signature = "x", function(x) standardGeneric("verify_project_paths"))
#' @export
setGeneric("verify_grid_paths", signature = "x", function(x) standardGeneric("verify_grid_paths"))

setMethod("verify_project_paths", "EcotonerSettings", function(x) {
	# Project paths
	attempt_out <- FALSE
	
	if (!is.na(x@dir_prj@path)) {
		x@dir_prj@path <- normalizePath(x@dir_prj@path, mustWork = TRUE)
		dir.create(x@dir_prj@path, showWarnings = FALSE)

		x@dir_init@path <- file.path(x@dir_prj@path, "1_Inits")
		dir.create(x@dir_init@path, showWarnings = FALSE)
		if (!is.na(file_searchpoints(x))) file_searchpoints(x) <- basename(file_searchpoints(x))
		
		x@dir_out@path <- file.path(x@dir_prj@path, "4_Output")
		
		attempt_out <- TRUE
	}

	if (!is.na(x@dir_big@path)) {
		temp <- try(normalizePath(x@dir_big@path, mustWork = TRUE), silent = TRUE)
		if (!inherits(temp, "try-error")) {
			x@dir_big@path <- temp
			dir.create(x@dir_big@path, showWarnings=FALSE)
			attempt_out <- TRUE
		} else {
			attempt_out <- FALSE ## assuming that x@dir_big@path is real, but currently unavailable
		}

		x@dir_out@path <- file.path(x@dir_big@path, "4_Output")
	}

	if (!is.na(dir_out(x))) {
		x@dir_out_fig@path <- file.path(dir_out(x), "Figures")
		x@dir_out_dat@path <- file.path(dir_out(x), "Transects")
		if (!is.na(file_etsummary(x))) file_etsummary(x) <- basename(file_etsummary(x))
	}
	
	if (attempt_out) {
		dir.create(x@dir_out@path, recursive = TRUE, showWarnings = FALSE)
		dir.create(x@dir_out_fig@path, showWarnings = FALSE)
		dir.create(x@dir_out_dat@path, showWarnings = FALSE)
	}
	
	x
})

setMethod("verify_grid_paths", "EcotonerSettings", function(x) {
	# GIS data
	if (!is.na(x@dir_env@path)) {
		x@dir_env@path <- normalizePath(x@dir_env@path, mustWork = NA)
	} else {
		warning("'ecotoner' requires a 'grid_env' to locate ecotone transects: currently, no path to such a grid is set")
	}

	if (!is.na(x@dir_veg@path)) {
		x@dir_veg@path <- normalizePath(x@dir_veg@path, mustWork = NA)
	} else {
		warning("'ecotoner' requires a 'grid_veg' to locate ecotone transects: currently, no path to such a grid is set")
	}

	if (!is.na(x@dir_flow@path)) {
		x@dir_flow@path <- normalizePath(x@dir_flow@path, mustWork = NA)
	} else if (x@transect_type == 3L) {
		warning("Option 'transect_type == 3' (flow paths) requires a 'grid_flow': currently, no path to such a grid is set")
	}
	
	if (!is.na(x@dir_abut@path)) {
		x@dir_abut@path <- normalizePath(x@dir_abut@path, mustWork = NA)
	} else if (x@transect_type == 4L) {
		warning("Option 'transect_type == 4' (homogenous aspect) requires a 'grid_abut': currently, no path to such a grid is set")
	}
	
	if (!is.na(x@dir_aspect_mean@path)) {
		x@dir_aspect_mean@path <- normalizePath(x@dir_aspect_mean@path, mustWork = NA)
	} else if (x@transect_type == 4L) {
		warning("Option 'transect_type == 4' (homogenous aspect) requires a 'grid_aspect_mean': currently, no path to such a grid is set")
	}
	
	if (!is.na(x@dir_aspect_sd@path)) {
		x@dir_aspect_sd@path <- normalizePath(x@dir_aspect_sd@path, mustWork = NA)
	} else if (x@transect_type == 4L) {
		warning("Option 'transect_type == 4' (homogenous aspect) requires a 'grid_aspect_sd': currently, no path to such a grid is set")
	}
	
	x
})



