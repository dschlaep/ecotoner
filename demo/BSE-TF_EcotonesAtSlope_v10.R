#------------------------------------------------------------#
if (FALSE) { # TODO(drs): conversion package
	libraries  <- c("ecotoner", "devtools")
	l <- lapply(libraries, FUN=function(lib) stopifnot(require(lib, character.only=TRUE, quietly=FALSE)))
	
	dir_dev <- "~/Dropbox/Work_Stuff/2_Research/Software/GitHub_Projects"
	setwd(file.path(dir_dev, "ecotoner"))
	devtools::document()
	devtools::load_all()	
}						

#------------------------------------------------------------#
##------PROJECT NEW/CONTD

prj_status <- "contd"		# one of "new" and "contd"; if "contd" then bname_settings (and bname_grids if locate_transects) must exist on disk
bname_settings <- "20160229_0757_ecotoner_settings.rds"
bname_grids <- "20160229_0757_ecotoner_grids.rds"
time_stamp <- Sys.time()


##------INPUTS
do.debug <- FALSE
do.demo <- TRUE		# If TRUE, then code uses the example data from the ecotoner package, i.e., can be run without additional input data


actions <- c(	locate_transects = TRUE,	# calls the transect functions detect_ecotone_transects_from_searchpoint()
#				make_summary = FALSE,		# creates a table with the results from calling detect_ecotone_transects_from_searchpoint()
				measure_transects = TRUE,	# uses methods to extract information about the located transects
				make_map = FALSE			# draws a map with all transects
			)

interactions <- c(	verbose = TRUE,		# prints progress statements
					figures = TRUE,		# saves figures as pdf to disk
					save_interim = TRUE	# saves workspace at two time points in each function call
				)
			


#------------------------------------------------------------#
##------SETUP R PACKAGES
#libraries  <- c("ecotoner", "doMC", "rgeos", "rgdal", "raster", "igraph", "foreign", "geosphere", "maptools", "akima", "modeest", "boot", "simpleboot", "zoo", "gstat", "spdep", "circular", "RANN")
libraries  <- c("ecotoner", "parallel", "foreign")
l <- lapply(libraries, function(lib) stopifnot(require(lib, character.only = TRUE, quietly = FALSE)))

#options(warn=2, error=quote({dump.frames(to.file=TRUE); q()}))	#turns all warnings into errors, dumps all to a file, and quits
#NOTE: raster(<=v2.2.5)::overlay() [testing fun] and pdf() can cause an error that normally is caught with try() or similar, but code fails if error=q()



##------PROJECT SETTINGS
# Set path dir.prj according to your own project and system
dir.prj <- "~/Dropbox/Work_Stuff/2_Research/200907_UofWyoming_PostDoc/Product17_EcotoneGradients/3_EcotoneGradient/4_Intermountain_v6"
setwd(dir.prj)

if (prj_status == "new") {
	esets <- new("EcotonerSettings")

	transect_type(esets) <- 4
	transect_N(esets) <- 15000
	cores_N(esets) <- min(12, parallel::detectCores() - 2)
	reproducible(esets) <- TRUE
	reseed(esets) <- FALSE
	neighborhoods(esets) <- 1667
	stepsHullEdge(esets) <- c(1, 3)
	clumpAreaClasses_m2(esets) <- c(1e4, 1e6)


	# Paths
	dir_prj(esets) <- dir.prj
	dir_big(esets) <- "/Volumes/BookDuo_12TB/BigData/Product17_EcotoneGradients/3_EcotoneGradient/4_Intermountain_v6"

	fname_settings <- file.path(dir_init(esets), paste0(format(time_stamp, format = "%Y%m%d_%H%M"), "_ecotoner_settings.rds"))
	bname_searchpoints <- paste0("SearchPoints_N", transect_N(esets), "_Veg1and2Abut.rds")

	if (do.debug || do.demo) {
		interactions[seq_along(interactions)] <- TRUE
		transect_N(esets) <- 10
		reproducible(esets) <- TRUE
		reseed(esets) <- do.debug
		if (do.demo) {
			dir_big(esets) <- NA_character_
			fname_searchpoints <- paste0("SearchPoints_N", transect_N(esets), "_Veg1and2Abut.rds")
		}
	}

	# Files
	file_etsummary(esets) <- paste0("Table_transect_summary", if (do.debug) "_debug", ".csv")
	file_searchpoints(esets) <- file.path(dir_init(esets), bname_searchpoints)

} else {
	dir.init <- file.path(dir.prj, "1_Inits")
	fname_settings <- file.path(dir.init, bname_settings)
	
	if (file.exists(fname_settings)) {
		esets <- readRDS(file = fname_settings)
	} else {
		stop("Requested file ", fname_settings, " was not found: 'ecotoner' cannot continue")
	}
}



##------RASTER GRIDS
if (actions["locate_transects"]) {
	if (prj_status == "new") {
		fname_grids <- file.path(dir_init(esets), paste0(format(time_stamp, format = "%Y%m%d_%H%M"), "_ecotoner_grids.rds"))

		if (!do.demo) {
			# Paths to grids
			dir.gis <- "/Volumes/BookDuo_12TB/BigData/GIS/Data"
			dir_env(esets) <- file.path(dir.gis, "Topography", "NED_USA", "NED_1arcsec")
			dir.veg <- file.path(dir.gis, "SpeciesVegetation", "GAP_2011_v2")
			dir_veg(esets) <- file.path(dir.veg, "Nat_GAP_LandCover")
			dir_abut(esets) <- file.path(dir.veg, "Types")
			dir_aspect_mean(esets) <- dir_aspect_sd(esets) <- dir_flow(esets) <- file.path(dir_env(esets), "terrain")

			esets <- verify_grid_paths(esets)
	
			# Load grids
			egrids <- new("EcotonerGrids")

			#--- Vegetation
			fgap <- file.path(dir_veg(esets), "natgaplandcov_v2_1.img")
			gap <- raster::raster(fgap)
			raster::crs(gap) <- sp::CRS(paste(raster::crs(gap), "+datum=NAD83"))
			grid_veg(egrids) <- gap
			df_veg(egrids) <- foreign::read.dbf(file=paste0(fgap, ".vat.dbf")) # previously: gap.rat
			rm(gap)

			#--- Topography
			grid_env(egrids) <- raster::raster(file.path(dir_env(esets), "ned_1s_westernUS_AEANAD83.tif"))

			if (transect_type(esets) == 4) {
				#--- Abutting cells
				grid_abut(egrids) <- raster::raster(file.path(dir_abut(esets), "gapv2_bseABUTtf.tif"))

				#--- Smoothed aspect
				# These rasters were calculated (/Users/drschlaep/Documents/drschlaepfer/3_BigData/200907_UofWyoming_PostDoc/GISonE/Data/Topography/NED_USA/NED_1arcsec/3_Calc_Terrain.R) assuming:
				#	- min_slope_with_aspect == 2 * pi/180
				# 	- bandTransect_width_cellN == 200
				stopifnot(min_slope_with_aspect(esets) == 2 * pi/180, bandTransect_width_cellN(esets) == 200)
				grid_aspect_mean(egrids) <- raster::raster(file.path(dir_aspect_mean(esets), "aspect_201Mean_ned_1s_westernUS_AEANAD83.tif"))
				grid_aspect_sd(egrids) <- raster::raster(file.path(dir_aspect_sd(esets), "aspect_201SD_ned_1s_westernUS_AEANAD83.tif"))
			}


			if (transect_type(esets) == 3) {
				#--- Flowpath
				#	If cell x has one of the following values, then flow from x is onto that position with that value
				#	32	64	128
				#	16	x	1
				#	8	4	2
				ftemp1 <- file.path(dir_flow, "arcgis_flowpath", "flow.tif")
				ftemp2 <- filename = file.path(dir_flow, "flow_Rraster.tif")
				if (file.exists(ftemp1)) {
					grid_flow(egrids) <- raster::raster(ftemp1)	#Flowpath generated in ArcGIS 10
				} else if (file.exists(ftemp2)) {
					grid_flow(egrids) <- raster::raster(ftemp2)	#Flowpath generated with raster::terrain
				} else {
					message(paste0(Sys.time(), ": 'raster::terrain' will calculate a flow path grid; this may take a moment."))
					grid_flow(egrids) <- raster::terrain(grid_env(egrids), opt = 'flowdir', filename = ftemp2, overwrite = FALSE)
				}
			}
		} else {
			path_demo <- system.file("extdata", "raster", package = "ecotoner")

			dir_env(esets) <- path_demo
			dir_veg(esets) <- path_demo
			dir_abut(esets) <- path_demo
			dir_aspect_mean(esets) <- dir_aspect_sd(esets) <- dir_flow(esets) <- path_demo

			esets <- verify_grid_paths(esets)
	
			# Load grids
			egrids <- new("EcotonerGrids")

			grid_veg(egrids) <- raster::raster(file.path(dir_veg(esets), "veg_eg.grd"))
			
			df_veg(egrids) <- readRDS(file.path(path_demo, "gap_rat.rds"))

			grid_env(egrids) <- raster::raster(file.path(dir_env(esets), "elev_eg.grd"))

			grid_abut(egrids) <- raster::raster(file.path(dir_abut(esets), "abutt_eg.grd"))
			grid_aspect_mean(egrids) <- raster::raster(file.path(dir_aspect_mean(esets), "asp201Mean_eg.grd"))
			grid_aspect_sd(egrids) <- raster::raster(file.path(dir_aspect_sd(esets), "asp201SD_eg.grd"))
		}
		
	} else {
		fname_grids <- file.path(dir.init, bname_grids)

		if (file.exists(fname_grids)) {
			egrids <- readRDS(file = fname_grids)
		} else {
			stop("Requested file ", fname_grids, " was not found: 'ecotoner' cannot continue")
		}
	}

	stopifnot(valid_grid(grid_veg(egrids)), valid_grid(grid_env(egrids)))
	if (prj_status == "new") saveRDS(egrids, file = fname_grids)
}


##------VEGETATION SETTINGS
if (prj_status == "new") {

	rat_veg <- if (exists("egrids") && inherits(egrids, "EcotonerGrids")) {
					df_veg(egrids)
				} else {
					if (!do.demo) {
						fgap <- file.path(dir_veg(esets), "natgaplandcov_v2_1.img")
						foreign::read.dbf(file=paste0(fgap, ".vat.dbf"))
					} else {
						path_demo <- system.file("extdata", "raster", package = "ecotoner")
						load(file.path(path_demo, "meta_eg.RData"))
						gap.rat
					}
				}

	# Human impact on land cover:
	#CL - NVC_CLASS
	#07	- "Agricultural Vegetation": all
	#08 - "Developed & Other Human Use" (Mining and developed areas): all
	#10 - "Recently Disturbed or Modified": all (i.e., generic disturbance, logged or harvested, disturbed/successional)
	#	except 'recently burned' (LEVEL3: 8301, 8302, 8303, 8304)
	type_ids(type_excl(esets)) <- na.omit(rat_veg$Value[rat_veg$CL %in% c("07", "08", "10") & !(rat_veg$LEVEL3 %in% c(8301, 8302, 8303, 8304))])
	type_field(type_excl(esets)) <- ""

	# Big sagebrush ecosystems:
	#Value - gap.rat$LEVEL3	- gap.rat$ECOLSYS_LU:
	#489	- 5706	- Inter-Mountain Basins Big Sagebrush Shrubland
	#490	- 5307	- Inter-Mountain Basins Big Sagebrush Steppe
	#491	- 5308	- Inter-Mountain Basins Montane Sagebrush Steppe
	type_ids(type_veg1(esets)) <- c(489, 490, 491)
	type_field(type_veg1(esets)) <- ""
	end_to_left(type_veg1(esets)) <- FALSE

	# Temperate forests:
	type_ids(type_veg2(esets)) <- na.omit(rat_veg$Value[rat_veg$NVC_FORM %in% c("Warm Temperate Forest", "Cool Temperate Forest")]) # i.e., excluding boreal, tropical and flooded forests
	type_field(type_veg2(esets)) <- ""
	end_to_left(type_veg2(esets)) <- TRUE
}



#------------------------------------------------------------#
##------SETUP COMPUTER
if (any(actions)) {
	#---Print system information
	print(sessionInfo())
	if (interactions["verbose"]) {
		blas <- system2(command = file.path(Sys.getenv()[["R_HOME"]], "R"), args = "R CMD config BLAS_LIBS", stdout = TRUE)
		blas <- sub("-L/", "/", (strsplit(blas, split = " ")[[1]][1]))
		lapack <- system2(command = file.path(Sys.getenv()[["R_HOME"]], "R"), args = "R CMD config LAPACK_LIBS", stdout = TRUE)
		lapack <- sub("-L/", "/", (strsplit(lapack, split = " ")[[1]][1]))
		get_ls <- if(identical(blas, lapack)) list(blas) else list(blas, lapack)
		temp <- lapply(get_ls, FUN = function(x) print(system2(command = "ls", args = paste("-l", x), stdout = TRUE)))
		rm(temp, blas, lapack, get_ls)
	}
	
	#---Verify folder setup
	esets <- verify_project_paths(esets)

	#---Setup: parallel, RNG, project information
	cat(format(Sys.time(), format = ""), ": set up parallel cluster'\n", sep = "")
	if (interactive()) cores_N(esets) <- 1
	cl  <- if (.Platform$OS.type == "unix") {
				parallel::makeCluster(cores_N(esets), type = "FORK", outfile = "log_ecotone.txt")
			} else if (.Platform$OS.type == "windows") {
				message("Running this code in parallel on a Windows OS computer has not been tested.")
				parallel::makeCluster(cores_N(esets), type = "PSOCK", outfile = "log_ecotone.txt")
			} else {
				stop("Running this code on this type of platform is currently not implemented.")
			}
	parallel::clusterEvalQ(cl, library("ecotoner"))
	
	# In case of abort and unable to call on.exit(); after crash, load cl and close cluster properly
	ftemp_cl <- file.path(dir_init(esets), "ClusterData.RData")
	save(cl, file = ftemp_cl)
	
	# Set up random number generator	
	if (do.debug) message("TODO(drs): Does load balancing allow replicability when used with ('L'Ecuyer-CMRG', master: clusterSetRNGStream(); tasks: set.seed())?")
	old_RNGkind <- RNGkind()
	RNGkind("L'Ecuyer-CMRG")
	parallel::clusterSetRNGStream(cl, get_sseed(esets))
	if (prj_status == "new") esets <- init_pseed(esets)

	# Save project information to disk
	if (prj_status == "new") saveRDS(esets, file = fname_settings)
}



#------------------------------------------------------------#
#------ECOTONER: LOCATE TRANSECTS
if (actions["locate_transects"]) {
	cat(format(Sys.time(), format = ""), ": using 'ecotoner' to locate ecotone transects along elevational gradients between big sagebrush and temperate forest vegetation\n", sep = "")
	
	stopifnot(exists("esets"), exists("egrids"))
	
	#---Get search points to initiate transects
	cat(format(Sys.time(), format = ""), ": sample ", transect_N(esets), " 'initpoints'\n", sep = "")
	initpoints <- get_transect_search_points(N = transect_N(esets),
											grid_mask1 = if (valid_grid(grid_abut(egrids))) grid_abut(egrids) else grid_veg(egrids),
											grid_maskNA = grid_env(egrids),
											initfile = file_searchpoints(esets),
											seed = get_sseed(esets),
											verbose = interactions["verbose"])
	

	#---Loop through random points
	cat(format(Sys.time(), format = ""), ": send ", transect_N(esets), " calls to the function 'detect_ecotone_transects_from_searchpoint'\n", sep = "")

	if (do.debug) message("TODO(drs): organize output; separate action:make_summary")

.Last <- function() {#TODO(drs): remove
	sessionInfo()
	print(warnings())
}
	resultTransects <- parallel::parSapplyLB(cl, seq_len(transect_N(esets)), detect_ecotone_transects_from_searchpoint,
							initpoints = initpoints,
							ecotoner_settings = esets,
							ecotoner_grids = egrids, 
							do_interim = interactions["save_interim"],
							verbose = interactions["verbose"],
							do_figures = interactions["figures"])
	
	print(warnings())
	#---Results
	write.csv(resultTransects, file = file_etsummary(esets))
}



#------------------------------------------------------------#
#------ECOTONER: MEASURE TRANSECTS
if (actions["measure_transects"]) {
	cat(format(Sys.time(), format = ""), ": using 'ecotoner' to measure/observe ecotone transects along elevational gradients between big sagebrush and temperate forest vegetation\n", sep = "")

	stopifnot(exists("esets"))

	#---Loop through random points
	cat(format(Sys.time(), format = ""), ": send ", transect_N(esets), " calls to the function 'measure_ecotone_per_transect'\n", sep = "")

et_methods_choices <- c("Danz2012JVegSci_1D", "Danz2012JVegSci_2D", "Eppinga2013Ecography", "Gastner2010AmNat")

	XXX <- parallel::parSapplyLB(cl, seq_len(transect_N(esets)), measure_ecotone_per_transect,
							et_methods = c("Danz2012JVegSci_2D"),
							ecotoner_settings = esets,
							verbose = interactions["verbose"],
							do_figures = interactions["figures"])
	
	print(warnings())
stop("TODO(drs): not implemented")	
et_methods2 <- c("InterZoneLocation", "InterZonePatchDistr")
	XXX <- measure_ecotones_all_transects()

	#---Results
	write.csv(XXX, file = fXXX)
}



#------------------------------------------------------------#
#Draw overall map
if (actions["make_map"]){
	transect.files <- list.files(path=dir.big.trans, pattern="_etransect.RData", full.names=TRUE)
	
	if(length(transect.files) > 0){
		if(do.verbose) print(paste(Sys.time(), "start overall map with", length(transect.files), "transects"))
		
		gadm_usa <- readRDS(file.path("~/Dropbox/Work_Stuff/2_Research/200907_UofWyoming_PostDoc/Product17_EcotoneGradients/2_Data", "20160126_GADMv28_USA_adm1.rds")) # WGS84
		gadm_usa_aeanad83 <- sp::spTransform(gadm_usa, raster::crs(elev))

		pdf(width=14, height=14, file=file.path(dir.out, paste0("Map_Transect_All_v2.pdf")))
			plot(elev, col=gray(0:255/255))
		
			raster::image(calc(gap, fun=function(x) ifelse(x %in% bseRatValue, 1, NA)), col=adjustcolor("red", alpha.f = 0.3), add=TRUE)
			raster::image(calc(gap, fun=function(x) ifelse(x %in% tfRatValue, 1, NA)), col=adjustcolor("darkgreen", alpha.f = 0.3), add=TRUE)
		
			for(i in seq_along(transect.files)){
				if(do.verbose) print(paste("debug", i, "out of", length(transect.files)))
				res_neighbors <- list()
				load(transect.files[i])
				if(length(res_neighbors) > 0) for(b in seq_along(res_neighbors)){
					if(!is.null(res_neighbors[[b]]$etline)){
						points(res_neighbors[[b]]$etline$pts, type="l", lwd=2, col="yellow")
					}
				}
			}
			
			plot(gadm_usa_aeanad83, border = "orange", add = TRUE)
			
		dev.off()
	}
}



#------------------------------------------------------------#
	# Parallel clean-up
	parallel::stopCluster(cl)
	unlink(ftemp_cl)
	RNGkind(kind = old_RNGkind[1], normal.kind = old_RNGkind[2])

