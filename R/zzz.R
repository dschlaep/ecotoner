# Function to remove temporary raster files at end of session or unloading of package
#.my_finalizer <- function(var) unlink(x = file.path(tempdir(), "raster"), recursive = TRUE, force = TRUE)

#TODO(drs): I worked on some code to use with ‘reg.finalizer’ to delete all 
#	remaining temporary raster files when R session ends or the ‘ecotoner’
#	package is unloaded (see function ‘.my_finalizer’ and its registration
#	in ‘.onAttach’): as it currently is, every temporary raster file would
#	be deleted even if not originating from this package - bad!
#	Better idead: have every raster function create a specific file to dir_temp(esets), which
#	will be deleted either at completion of the transect location call or with a reg.finalizer

.onAttach <- function(libname, pkgname) {
	if (interactive()) {
    	meta <- packageDescription("ecotoner")
    	packageStartupMessage("Package 'ecotoner', ", meta$Version, " (", meta$Date, ") attached/loaded.")
	}
	
	invisible()
}


.onLoad <- function(libname, pgkname) {
	#--- Define options and set default values
	# based on chapter "When you do need side-effects" by Wickham, H. 2015. R packages. O'Reilly and Associates.
	op_old <- options()
	op_ecotoner <- list()
	toset <- !(names(op_ecotoner) %in% names(op_old))
	if (any(toset)) options(op_ecotoner[toset])

	#--- Define package level variables that should be hidden from package user and should not be changed
	assign("transect_types", c("HighestElevation", "SteepestSlope", "Flowpath", "HomogeneousAspect"), envir = etr_vars)
	assign("migtypes", c("AllMigration", "OnlyGoodMigration"), envir = etr_vars)
	
	#--- Set finalizer
	#reg.finalizer(e = parent.env(environment()), # register the finalizer on the enclosing environment of the .onLoad function == environment of the package 'ecotoner'
	#			  f = function (env) .my_finalizer(env$etr_vars), onexit = TRUE)
	
	invisible()
}

		

.onUnload <- function(libpath) {
	#--- Remove package options
	op_old <- options()
	op_ecotoner <- list()
	toset <- names(op_ecotoner) %in% names(op_old)
	if (any(toset)) options(op_ecotoner[toset])

	invisible()

}
