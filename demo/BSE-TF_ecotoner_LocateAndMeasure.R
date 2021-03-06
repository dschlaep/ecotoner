#------------------------------------------------------------#
# An ecotoner project
#
# by Daniel R Schlaepfer, 2016
#------------------------------------------------------------#

#------------------------------------------------------------#
##------PROJECT NEW/CONTD
prj_status <- "contd"		# one of "new" and "contd"; if "contd" then bname_settings (and bname_grids if locate_transects) must exist on disk
bname_settings <- "myproject_ecotoner_settings.rds"
bname_grids <- "myproject_ecotoner_grids.rds"
time_stamp <- Sys.time()


##------INPUTS
do.debug <- FALSE
do.demo <- TRUE		# If TRUE, then code uses the example data from the ecotoner package, i.e., can be run without additional input data


actions <- c(	preprocess_data = FALSE,		# mosaic NED tiles; project NED to crs(GAP); determine BSE-TF abutting; calculate smoothed aspect
				sample_searchpoints = TRUE,		# samples locations to initiate the location of transects
				locate_transects = TRUE,		# call the transect functions detect_ecotone_transects_from_searchpoint()
				make_map = TRUE,				# draw a map of all transects
				measure_transects = TRUE		# use methods to extract information about the located transects
			)

interactions <- c(	verbose = TRUE,		# prints progress statements
					figures = TRUE,		# saves figures as pdf to disk
					save_interim = TRUE	# saves workspace at two time points in each function call
				)
			

#------------------------------------------------------------#
##------SETUP R PACKAGES
pkg_reqd  <- c("ecotoner", "parallel", "foreign")
l <- lapply(pkg_reqd, function(lib) stopifnot(require(lib, character.only = TRUE, quietly = FALSE)))



##------PROJECT SETTINGS
# Set path dir.prj according to your own project and system
dir.prj <- "~"
setwd(dir.prj)

if (prj_status == "new") {
	esets <- new("EcotonerSettings")

	transect_type(esets) <- 4
	searchpoints_N(esets) <- 20000
	inhibit_searchpoints(esets) <- TRUE
	candidate_THAs_N(esets) <- 50
	cores_N(esets) <- min(20, parallel::detectCores() - 1) # 0 or 1 will turn off parallelization
	reproducible(esets) <- TRUE		# If TRUE, then transects set their own unique random seed (this setup is reproducible even after re-starting a parallel function call)
	neighborhoods(esets) <- 1667
	stepsHullEdge(esets) <- c(1, 3, 10)
	clumpAreaClasses_m2(esets) <- c(1e4, 1e6)


	# Paths
	dir_prj(esets) <- dir.prj
	dir_big(esets) <- dir.prj

	fname_settings <- file.path(dir_init(esets), paste0(format(time_stamp, format = "%Y%m%d_%H%M"), "_ecotoner_settings.rds"))
	bname_searchpoints <- paste0("SearchPoints_",
		if (inhibit_searchpoints(esets)) paste0(inhibit_dist_m(esets, 1), "cells_inhibited") else "Poisson", "_",
		searchpoints_N(esets), "N_",
		"Veg1and2Abut.rds")

	if (do.debug || do.demo) {
		interactions[seq_along(interactions)] <- TRUE
		searchpoints_N(esets) <- 6

		if (do.demo) {
			inhibit_searchpoints(esets) <- FALSE
			neighborhoods(esets) <- 667
			candidate_THAs_N(esets) <- 20
			dir_big(esets) <- NA_character_
			bname_searchpoints <- paste0("SearchPoints_",
				if (inhibit_searchpoints(esets)) paste0(inhibit_dist_m(esets, 1), "cells_inhibited") else "Poisson", "_",
				searchpoints_N(esets), "N_",
				"Veg1and2Abut.rds")
		}
		
	}

	# Files
	file_etsummary(esets) <- paste0("Table_transect_summary", if (do.debug) "_debug", ".csv")
	file_searchpoints(esets) <- file.path(dir_init(esets), bname_searchpoints)
	file_etmeasure_base(esets) <- paste0("Table_transect_measure", if (do.debug) "_debug", ".csv")
	file_initwindow(esets) <- file.path(dir_init(esets), paste0("owin_Veg1and2Abut", if (do.demo) "_demo", ".rds"))

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
if (actions["locate_transects"] || actions["make_map"]) {
	if (prj_status == "new") {
		fname_grids <- file.path(dir_init(esets), paste0(format(time_stamp, format = "%Y%m%d_%H%M"), "_ecotoner_grids.rds"))

		# Set the paths to the grids files on disk
		if (!do.demo) {
			dir.gis <- "/Volumes/BookDuo_12TB/BigData/GIS/Data"
			dir_env(esets) <- file.path(dir.gis, "Topography", "NED_USA", "NED_1arcsec")
			dir.veg <- file.path(dir.gis, "SpeciesVegetation", "GAP_2011_v2")
			dir_veg(esets) <- file.path(dir.veg, "Nat_GAP_LandCover")
			dir_veg1(esets) <- file.path(dir.veg, "Types", "BigSagebrushEcosystems")
			dir_veg2(esets) <- file.path(dir.veg, "Types", "ForestWoodlands")
			dir_abut(esets) <- file.path(dir.veg, "Types")
			dir_aspect_mean(esets) <- dir_aspect_sd(esets) <- dir_flow(esets) <- file.path(dir_env(esets), "terrain")
		} else {
			path_demo <- system.file("extdata", "spatial", package = "ecotoner")

			dir_env(esets) <- path_demo
			dir_veg(esets) <- dir_veg1(esets) <- dir_veg2(esets) <- path_demo
			dir_abut(esets) <- path_demo
			dir_aspect_mean(esets) <- dir_aspect_sd(esets) <- dir_flow(esets) <- path_demo
		}
		
		# Check the paths and get a new EcotonerGrid object
		esets <- verify_grid_paths(esets)
		egrids <- new("EcotonerGrids")

		#--- Vegetation
		if (!do.demo) {
			fgap <- file.path(dir_veg(esets), "natgaplandcov_v2_W97.tif")
			gap <- raster::raster(fgap)
			raster::crs(gap) <- sp::CRS(paste(raster::crs(gap), "+datum=NAD83"))
			grid_veg(egrids) <- gap
			df_veg(egrids) <- foreign::read.dbf(file = file.path(dir_veg(esets), "natgaplandcov_v2_1.img.vat.dbf"))
			rm(gap)
		} else {
			grid_veg(egrids) <- raster::raster(file.path(dir_veg(esets), "veg_eg.grd"))
			df_veg(egrids) <- readRDS(file.path(dir_veg(esets), "gap_rat.rds"))
		}

		#--- Vegetation metadata
		if (all(dim(df_veg(egrids)) > 0)) {
			# Human impact on land cover:
			#CL - NVC_CLASS
			#07	- "Agricultural Vegetation": all
			#08 - "Developed & Other Human Use" (Mining and developed areas): all
			#10 - "Recently Disturbed or Modified": all (i.e., generic disturbance, logged or harvested, disturbed/successional)
			#	except 'recently burned' (LEVEL3: 8301, 8302, 8303, 8304)
			type_ids(type_excl(esets)) <- na.omit(df_veg(egrids)$Value[df_veg(egrids)$CL %in% c("07", "08", "10") & !(df_veg(egrids)$LEVEL3 %in% c(8301, 8302, 8303, 8304))])
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
			type_ids(type_veg2(esets)) <- na.omit(df_veg(egrids)$Value[df_veg(egrids)$NVC_FORM %in% c("Warm Temperate Forest", "Cool Temperate Forest")]) # i.e., excluding boreal, tropical and flooded forests
			type_field(type_veg2(esets)) <- ""
			end_to_left(type_veg2(esets)) <- TRUE
		}
		
		#--- Vegetation types
		fveg1 <- if (!do.demo) {
						file.path(dir_veg1(esets), "gapv2_bse_W97.tif")
					} else {
						file.path(dir_veg1(esets), "veg1_bse_eg.grd")
					}
		if (actions["preprocess_data"] && !file.exists(fveg1)) {
			grid_veg1(egrids) <- extract_vegetation(
				grid_veg(egrids),
				ids = type_ids(type_veg1(esets)),
				filename = fveg1,
				parallel_N = cores_N(esets),
				dataType = "INT1S",
				options = c("COMPRESS=LWZ", "TFW=YES", "TILED=YES"))
		} else {
			grid_veg1(egrids) <- raster::raster(fveg1)
		}

		fveg2 <- if (!do.demo) {
						file.path(dir_veg2(esets), "gapv2_tf_W97.tif")
					} else {
						file.path(dir_veg2(esets), "veg2_tf_eg.grd")
					}
		if (actions["preprocess_data"] && !file.exists(fveg2)) {
			grid_veg2(egrids) <- extract_vegetation(
				grid_veg(egrids),
				ids = type_ids(type_veg2(esets)),
				filename = fveg2,
				parallel_N = cores_N(esets),
				dataType = "INT1S",
				options = c("COMPRESS=LWZ", "TFW=YES", "TILED=YES"))
		} else {
			grid_veg2(egrids) <- raster::raster(fveg2)
		}

		
		#--- Abutting cells
		if (transect_type(esets) == 4) {
			fabut <- if (!do.demo) {
						file.path(dir_abut(esets), "gapv2_bseABUTtf_W97.tif")
					} else {
						file.path(dir_veg(esets), "abutt_eg.grd")
					}
			if (actions["preprocess_data"] && !file.exists(fabut)) {
				grid_abut(egrids) <- determine_abutters(
					grid_veg(egrids),
					grid_veg1 = grid_veg1,
					grid_veg2 = grid_veg2,
					filename = fabut,
					dataType = "INT1S",
					options = c("COMPRESS=LWZ", "TFW=YES", "TILED=YES"))
			} else {
				grid_abut(egrids) <- raster::raster(fabut)
			}
		}

		
		#--- Topography
		if (!do.demo) {
			fenv <- file.path(dir_env(esets), "NED_1arcsec_AEA83_W97.tif")
			if (actions["preprocess_data"] && !file.exists(fenv)) {
				fenv_temp <- file.path(dir_env(esets), "ned_1s_westernUS_GeogrNAD83.tif")
				fenv_temp2 <- sub("_GeogrNAD83.tif", "_AEA83.tif", fenv_temp)
				
				grid_env_temp <- mosaic_tiles(dir_tiles = file.path(dir_env(esets), "grid_tiles"), chunksize = 10L, fname_grid_ned = fenv)
				raster::crs(grid_env_temp) <- sp::CRS(paste(raster::crs(grid_env_temp, asText = TRUE), "+datum=NAD83 +vunits=m"))
				
				grid_env_temp2 <- project_raster(
					grid_from = grid_env_temp,
					fname_grid_to = fenv_temp2,
					res_to = raster::res(gap),
					crs_to = raster::crs(gap),
					parallel_N = cores_N(esets),
					chunksize = 1e+05, maxmemory = 1e+06,
					options = c("COMPRESS=NONE", "TFW=YES", "TILED=YES"))

				grid_env_temp3 <- raster::extend(raster::crop(grid_env_temp2, raster::extent(grid_veg(egrids))), raster::extent(grid_veg(egrids)))
				grid_env(egrids) <- raster::overlay(
					grid_env_temp3, grid_veg(egrids),
					fun = function(x, y) ifelse(is.na(x) | is.na(y), NA, x),
					filename = fenv,
					dataype = "FLT4S",
					options = c("COMPRESS=LWZ", "TFW=YES", "TILED=YES"))

				rm(fenv_temp, grid_env_temp, grid_env_temp2, grid_env_temp3)
			} else {
				grid_env(egrids) <- raster::raster(fenv)
			}
		} else {
			grid_env(egrids) <- raster::raster(file.path(dir_env(esets), "elev_eg.grd"))
		}
		
		
		#--- Smoothed aspect
		if (transect_type(esets) == 4) {
			fslope <- if (!do.demo) {
						file.path(dir_aspect_mean(esets), "slope_ned_1s_westernUS_AEANAD83.tif")
					} else {
						file.path(dir_aspect_mean(esets), "slope_eg.grd")
					}
			if (actions["preprocess_data"] && !file.exists(fslope)) {
				grid_slope <- terrain_slope(
					grid_env(egrids),
					filename = fslope,
					options = c("COMPRESS=LWZ", "TFW=YES", "TILED=YES"))
			} else {
				grid_slope <- raster::raster(fslope)
			}
			
			faspect <- if (!do.demo) {
						file.path(dir_aspect_mean(esets), "aspect_ned_1s_westernUS_AEANAD83.tif")
					} else {
						file.path(dir_aspect_mean(esets), "aspect_eg.grd")
					}
			if (actions["preprocess_data"] && !file.exists(faspect)) {
				grid_aspect <- terrain_aspect(
					grid_env(egrids),
					grid_slope,
					min_slope = min_slope_with_aspect(esets),
					parallel_N = cores_N(esets),
					filename = faspect,
					options = c("COMPRESS=LWZ", "TFW=YES", "TILED=YES"))
			} else {
				# This rasters was calculated with: min_slope_with_aspect == 2 * pi/180
				stopifnot(min_slope_with_aspect(esets) == 2 * pi/180)
				grid_aspect <- raster::raster(faspect)
			}
			
			fasp_m <- if (!do.demo) {
						file.path(dir_aspect_mean(esets), "aspect_201Mean_ned_1s_westernUS_AEANAD83.tif")
					} else {
						file.path(dir_aspect_mean(esets), "asp201Mean_eg.grd")
					}
			if (actions["preprocess_data"] && !file.exists(fasp_m)) {
				grid_aspect_mean(egrids) <- homogenous_aspect(
					grid_aspect,
					fun = "mean",
					window_N = bandTransect_width_cellN(esets),
					filename = fasp_m,
					options = c("COMPRESS=LWZ", "TFW=YES", "TILED=YES"))
			} else {
				# This rasters was calculated with: bandTransect_width_cellN == 200
				stopifnot(bandTransect_width_cellN(esets) == 200)
				grid_aspect_mean(egrids) <- raster::raster(fasp_m)
			}
			
			fasp_sd <- if (!do.demo) {
						file.path(dir_aspect_sd(esets), "aspect_201SD_ned_1s_westernUS_AEANAD83.tif")
					} else {
						file.path(dir_aspect_sd(esets), "asp201SD_eg.grd")
					}
			if (actions["preprocess_data"] && !file.exists(fasp_sd)) {
				grid_aspect_sd(egrids) <- homogenous_aspect(
					grid_aspect,
					fun = "sd",
					window_N = bandTransect_width_cellN(esets),
					filename = fasp_sd,
					options = c("COMPRESS=LWZ", "TFW=YES", "TILED=YES"))
			} else {
				# This rasters was calculated with: bandTransect_width_cellN == 200
				stopifnot(bandTransect_width_cellN(esets) == 200)
				grid_aspect_sd(egrids) <- raster::raster(fasp_sd)
			}
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
				grid_flow(egrids) <- raster::terrain(
					grid_env(egrids),
					opt = 'flowdir',
					filename = ftemp2,
					overwrite = FALSE)
			}
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
	if (prj_status == "new") unlink(file_timing_locate(esets))

	#---Setup: random number generator	
	RNGkind_old <- RNGkind()

	#---Setup: parallel and project information
	if (interactive() || (actions["make_map"] && sum(actions) == 1)) cores_N(esets) <- 1
	
	if (cores_N(esets) > 1) {
		cat(format(Sys.time(), format = ""), ": setting up parallel cluster with ", cores_N(esets), " workers\n", sep = "")
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
	}
	
	pfun <- fun_to_apply_foreach_transect(esets)

	# Save project information to disk
	if (prj_status == "new") saveRDS(esets, file = fname_settings)
}



#------------------------------------------------------------#
##------ECOTONER: SAMPLE SEARCH POINTS
if (actions["sample_searchpoints"] || actions["locate_transects"] || actions["make_map"]) {
	cat(format(Sys.time(), format = ""), ": using 'ecotoner' to sample search locations from where to initiate the location of ecotone transects\n", sep = "")

	cat(format(Sys.time(), format = ""), ": sampling ", searchpoints_N(esets), " search points/locations\n", sep = "")
	initpoints <-
		if (prj_status == "new" || !file.exists(file_searchpoints(esets))) {
			get_transect_search_points(
				N = searchpoints_N(esets),
				grid_mask1 = if (valid_grid(grid_abut(egrids))) grid_abut(egrids) else grid_veg(egrids),
				inhibit_dist = inhibit_dist_m(esets, res_m(specs_grid(egrids))), 
				mywindow = if (inhibit_searchpoints(esets) && file.exists(file_initwindow(esets))) readRDS(file_initwindow(esets)) else NULL,
				grid_maskNA = grid_env(egrids),
				initfile = file_searchpoints(esets),
				initwindowfile = file_initwindow(esets),
				seed = if (reproducible(esets)) get_global_seed(esets) else NULL,
				verbose = interactions["verbose"])
		} else {
			readRDS(file_searchpoints(esets))
		}
	
	if (!identical(transect_N(esets), length(initpoints))) {
		transect_N(esets) <- length(initpoints)
		saveRDS(esets, file = fname_settings)
	}
}



#------------------------------------------------------------#
#------ECOTONER: LOCATE TRANSECTS
if (actions["locate_transects"]) {
	cat(format(Sys.time(), format = ""), ": using 'ecotoner' to locate ecotone transects along elevational gradients between big sagebrush and temperate forest vegetation\n", sep = "")
	
	#---Loop through random points
	cat(format(Sys.time(), format = ""), ": sending ", transect_N(esets), " calls to the function 'detect_ecotone_transects_from_searchpoint'\n", sep = "")

	seeds_locate <- if (reproducible(esets)) {
							prepare_RNG_streams(N = N_of_location_calls(esets),
												iseed = get_global_seed(esets))
					} else NULL
	
	.Last <- function() raster::removeTmpFiles(h = 0)
	resultTransects <- pfun(seq_len(transect_N(esets)),
							detect_ecotone_transects_from_searchpoint,
							initpoints = initpoints,
							ecotoner_settings = esets,
							ecotoner_grids = egrids,
							seed_streams = seeds_locate,
							do_interim = interactions["save_interim"],
							verbose = interactions["verbose"],
							do_figures = interactions["figures"])
	raster::removeTmpFiles(h = 0)

	#---Results
	write.csv(resultTransects, file = file_etsummary(esets))
}


#------------------------------------------------------------#
#------DRAW A MAP WITH TRANSECT LOCATIONS
if (actions["make_map"]){
	cat(format(Sys.time(), format = ""), ": calculating vegetation meta-data\n", sep = "")
	
	meta_veg <- characterize_veg_data(
		ecotoner_settings = esets,
		ecotoner_grids = egrids,
		initpoints = initpoints,
		inhibit_dist = inhibit_dist_m(esets, res_m(specs_grid(egrids))))


	cat(format(Sys.time(), format = ""), ": drawing a map of the study area\n", sep = "")

	id_transect_examples <- 1
	ext <- raster::union(raster::extent(grid_env(egrids)), raster::extent(grid_veg(egrids)))

	# map
	pdf(width = 7, height = 6, file = file.path(dir_out(esets), paste0("Map_InputData_TransectLocations.pdf")))
	op_old <- par(mar = c(2.5, 2.5, 0.5, 0.5))
	
		raster::plot(grid_env(egrids), asp = 1, ext = ext, col = gray(25:255/255), maxpixels = 50000)
		raster::image(grid_veg(egrids), col = "gray", add = TRUE)

		raster::image(grid_env(egrids), col = gray(25:255/255), add = TRUE)
		
		map <- try(readRDS(file.path(dir_out(esets), "..", "..", "..", "2_Data", "20160126_GADMv28_USA_adm1_AEANAD83.rds")), silent = TRUE)
		if (inherits(map, "SpatialPolygons"))
			sp::plot(map, border = "white", add = TRUE)

		if (raster::ncell(grid_veg1(egrids)) > 1)
			raster::image(grid_veg1(egrids), col = adjustcolor("red", alpha.f = 0.3), add = TRUE)
		if (raster::ncell(grid_veg2(egrids)) > 1)
			raster::image(grid_veg2(egrids), col = adjustcolor("darkgreen", alpha.f = 0.3), add = TRUE)

		if (!is.null(meta_veg[["initpoints_buffering"]]))
			sp::plot(meta_veg[["initpoints_buffering"]], col = adjustcolor("royalblue1", alpha.f = 0.4), border = NA, add = TRUE)
		if (raster::ncell(grid_abut(egrids)) > 1) {
			raster::image(grid_abut(egrids), col = "yellow", add = TRUE)
			raster::image(grid_abut(egrids), col = "yellow", add = TRUE)
		}

		if (!exists("resultTransects"))
			resultTransects <- try(read.csv(file = file_etsummary(esets)), silent = TRUE)
		if (!inherits(resultTransects, "try-error") && nrow(resultTransects) > 0) {
			sp_start <- sp::SpatialPoints(coords = resultTransects[, c("TransectLinear_StartPoint_WGS84_Long", "TransectLinear_StartPoint_WGS84_Lat")], proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
			sp_start <- sp::coordinates(sp::spTransform(sp_start, crs(specs_grid(egrids))))
			sp_end <- sp::SpatialPoints(coords = resultTransects[, c("TransectLinear_EndPoint_WGS84_Long", "TransectLinear_EndPoint_WGS84_Lat")], proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
			sp_end <- sp::coordinates(sp::spTransform(sp_end, crs(specs_grid(egrids))))
	
			segments(x0 = sp_start[, 1], y0 = sp_start[, 2],
					x1 = sp_end[, 1], y1 = sp_end[, 2],
					col = "orange", lwd = 1.5)
			
			if (length(id_transect_examples) > 0)
				itemp <- which(resultTransects[, "Transect_ID"] == id_transect_examples)
			if (length(itemp) > 0) {
				points(x = apply(cbind(sp_start[itemp, 1], sp_end[itemp, 1]), 1, mean),
						y = apply(cbind(sp_start[itemp, 2], sp_end[itemp, 2]), 1, mean),
						col = "black", pch = 1, lwd = 2)
				segments(x0 = sp_start[itemp, 1], y0 = sp_start[itemp, 2],
						x1 = sp_end[itemp, 1], y1 = sp_end[itemp, 2],
						col = "green", lwd = 2)			
			}		
		}

	par(op_old)
	dev.off()

}


#------------------------------------------------------------#
#------ECOTONER: MEASURE TRANSECTS
if (actions["measure_transects"]) {
	cat(format(Sys.time(), format = ""), ": using 'ecotoner' to measure/observe ecotone transects along elevational gradients between big sagebrush and temperate forest vegetation\n", sep = "")

	stopifnot(exists("esets"))

	# measure options
	et_methods_options <- c("Danz2012JVegSci_1D", "Danz2012JVegSci_2D", "Eppinga2013Ecography", "Gastner2010AmNat")

	#---Loop through random points
	cat(format(Sys.time(), format = ""), ": sending ", transect_N(esets), " calls to the function 'measure_ecotone_per_transect'\n", sep = "")

	seeds_measure1 <- if (reproducible(esets)) {
							prepare_RNG_streams(N = N_of_measure_calls(esets),
												iseed = get_global_seed(esets))
					  } else NULL					

	measureTransects1 <- pfun(seq_len(transect_N(esets)),
								measure_ecotone_per_transect,
								et_methods = et_methods_options[4],
								ecotoner_settings = esets,
								seed_streams = seeds_measure1,
								verbose = interactions["verbose"],
								do_figures = interactions["figures"])



	warning("TODO(drs): code not further implemented")	
	if (FALSE) {
	et_methods2 <- c("InterZoneLocation", "InterZonePatchDistr")

	if (!exists("resultTransects")) resultTransects <- read.csv(file = file_etsummary(esets), header = TRUE, row.names = 1)

	measureTransects2 <- measure_ecotones_all_transects(pfun,
								seq_len(transect_N(esets)),
								et_methods = c("Danz2012JVegSci_global_2D"),
								et_desc = resultTransects,
								ecotoner_settings = esets,
								verbose = interactions["verbose"],
								do_figures = interactions["figures"])
	}
}

print(sessionInfo())



#------------------------------------------------------------#
	# Parallel clean-up
	parallel::stopCluster(cl)
	unlink(ftemp_cl)
	RNGkind(kind = RNGkind_old[1], normal.kind = RNGkind_old[2])
