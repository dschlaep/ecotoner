#' @export
measure_ecotone_per_transect <- function(i, et_methods, ecotoner_settings, seed_streams = NULL, verbose = TRUE, do_figures = TRUE) {
	if (verbose) {
		idh <- 1
		cat("'ecotoner' measuring: tr = ", i, "; start at ", format(t1 <- Sys.time(), format = ""), "\n", sep = "")
	}
	
	# Output containers
	iflag <- flag_itransect(ecotoner_settings, i)
	migtypes <- get("migtypes", envir = etr_vars)
	
	# Determine what to do
	do_measure <- do_new <- TRUE
	do_only_write <- FALSE
	
	if (file.exists(fname_etmeasured(ecotoner_settings, iflag))) {
		load(fname_etmeasured(ecotoner_settings, iflag)) #load: i, b, etmeas
		do_new <- FALSE
		itemp <- et_methods %in% names(etmeas)
#TODO(drs): check not only presence of methods, but also whether all neighborhoods and migration types have been completed
		if (all(itemp)) {
			do_measure <- FALSE  # all measured
			do_only_write <- TRUE
		} else {
			et_methods <- et_methods[!(itemp)]
		}
	}
	
	if (do_measure && file.exists(fname_etlocated(ecotoner_settings, iflag))) {
		load(fname_etlocated(ecotoner_settings, iflag)) #i, b, and etransect loaded
	} else {
		do_measure <- FALSE # no suitable transect located for search point i
	}


	# Measure ecotones
	if (do_measure) {
		if (do_new) {
			template_etobs <- list(	gETmeas = vector(mode = "list", length = neighbors_N(ecotoner_settings)), # list for 'gETmeas' of each neighborhood
									seeds = vector(mode = "list", length = neighbors_N(ecotoner_settings) * length(migtypes)),
									etable = data.frame(matrix(NA, nrow = 1, ncol = 0))) #container combining output table results with rows for each neighborhood and migtype
			temp <- vector("list", length = length(migtypes))
			names(temp) <- migtypes
			template_etobs$gETmeas <- lapply(template_etobs$gETmeas, function(x) temp)
			etmeas <- vector(mode = "list", length = length(et_methods))
			names(etmeas) <- et_methods
		}
		
		for (etm in et_methods) etmeas[[etm]] <- template_etobs

		dir_fig <- file.path(dir_out_fig(ecotoner_settings), iflag)
		dir_create(dir_fig)

		for (b in seq_len(neighbors_N(ecotoner_settings))) {
			flag_bfig <- flag_basename(ecotoner_settings, iflag, b)

#TODO(drs):	eB_InterZone <- eB_PatchSizeDist <- template_etobs	#ecological boundaries

			#----3. Loop through both versions of vegetation for applying the methods : Veg$AllMigration and Veg$OnlyGoodMigration
			for (im in seq_along(migtypes)) {
				copy_FromMig1_TF <- if (migtypes[im] == "AllMigration") FALSE else !etransect$etbands[[b]]$Veg[["OnlyGoodMigration"]]$diffs_amongMigTypes_TF

				if (verbose) cat("'ecotoner' measuring: tr = ", i, "; neigh = ", b, ": prog: ", idh <- idh + 1, "; mig-type: ", migtypes[im], "; method: ", etm, "\n", sep = "")
				
				# loop through measurement methods 'et_methods'
				for (etm in et_methods) {
					iseed <- (b - 1) * length(migtypes) + im
					etmeas[[etm]][["seeds"]][[iseed]] <- if (is.null(seed_streams)) NULL else if (inherits(seed_streams, "list")) seed_streams[[((i - 1) * neighbors_N(ecotoner_settings) + (b - 1)) * length(migtypes) + im]] else NA
					set_RNG_stream(etmeas[[etm]][["seeds"]][[iseed]])

					etmeas[[etm]] <- do.call(what = etm, args = list(i = i, b = b, migtype = migtypes[im], 
																		ecotoner_settings = ecotoner_settings,
																		etband = etransect[["etbands"]][[b]],
																		etmeasure = etmeas[[etm]],
																		copy_FromMig1_TF = copy_FromMig1_TF,
																		do_figures = do_figures, dir_fig = dir_fig, flag_bfig = flag_bfig,
																		seed = NA))

					# Write data to disk file
					index <- which(etmeas[[etm]]$etable[, "Neighbor_Cells"] == neighborhoods(ecotoner_settings)[b] &
									etmeas[[etm]]$etable[, "Migration_Type"] == migtypes[im])
					write_ecotoner_row(data_row = etmeas[[etm]]$etable[index, ],
										filename = file_etmeasure_base(ecotoner_settings, etm),
										tag_fun = 'measure_ecotone_per_transect',
										tag_id = paste0(i, " 'etmeas[[", etm, "]]$etable[", paste(index, collapse = "-"), ", ]'"))

					# Save data to RData disk file
					save(i, b, im, etmeas, file = fname_etmeasured(ecotoner_settings, iflag))
				}
			} # end loop through migtypes
		} # end loop through neighborhoods
	}

	# Write already measured information
	if (do_only_write) {
		for (b in seq_len(neighbors_N(ecotoner_settings))) for (etm in et_methods) {
			index <- which(etmeas[[etm]]$etable[, "Neighbor_Cells"] == neighborhoods(ecotoner_settings)[b])
			write_ecotoner_row(data_row = etmeas[[etm]]$etable[index, ],
								filename = file_etmeasure_base(ecotoner_settings, etm),
								tag_fun = 'measure_ecotone_per_transect',
								tag_id = paste0(i, " 'etmeas[[", etm, "]]$etable[", paste(index, collapse = "-"), ", ]'"))
		}
	}
	
		
	if (verbose) cat("'ecotoner' measuring: tr = ", i, " completed ", format(t2 <- Sys.time(), format = ""), " after ", round(difftime(t2, t1, units="hours"), 2), " hours\n", sep = "")

	i
}	


#' @export
measure_ecotones_all_transects <- function(pfun, X, et_methods, et_desc, ecotoner_settings, verbose, do_figures) {
	
	migtypes <- get("migtypes", envir = etr_vars)
	veg_types <- c("Veg1", "Veg2")
	
	# Load data of transects
#TODO(drs): Is 'transect_data2d' going to be so large that it needs to be a big data object such as ff::ffdf?
	transect_data2d <- pfun(X, get_transect_grids_as_df, etransect, et_desc, ecotoner_settings, migtypes, veg_types)
	colnames(transect_data2d)[3:8] <- c("env", "reps", paste(rep(paste("y", veg_types, sep = "_"), times = length(migtypes)), rep(migtypes, each = length(veg_types)), sep = "_"))
			
	m1 <- glm(y_Veg1_AllMigration ~ env * Transect_ID + 
									Transect_Azimuth_deg + reps,
				family = binomial, data = transect_data2d)
	X
}

