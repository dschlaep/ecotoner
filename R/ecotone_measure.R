#' @export
measure_ecotone_per_transect <- function(i, et_methods, ecotoner_settings, seed_streams = NULL, verbose = TRUE, do_figures = TRUE) {
	# Output containers
	iflag <- flag_itransect(ecotoner_settings, i)
	migtypes <- get("migtypes", envir = etr_vars)
	
	# Determine what to do
	what_measure <- array(TRUE, dim = c(length(et_methods), neighbors_N(ecotoner_settings), length(migtypes)),
							dimnames = list(et_methods, neighborhoods(esets), migtypes))
	
	if (file.exists(fname_etmeasured(ecotoner_settings, iflag))) {
		what_measure[] <- FALSE
		load(fname_etmeasured(ecotoner_settings, iflag)) #load: i, b, im, etmeas
		
		# Collect up-to-date method version numbers
		cur_versions <- lapply(et_methods, function(etm) getFromNamespace(paste0("version_", etm), "ecotoner")())
		names(cur_versions) <- et_methods
		
		# Determine measure methods status
		#	- do those that are not in etmeas yet
		started_etmeth <- et_methods %in% names(etmeas)
		what_measure[!started_etmeth, , ] <- TRUE
		#	- do those that are only list()
		started_etmeas <- names(etmeas) %in% et_methods
		if (any(started_etmeas)) {
			redo_etmeas <- sapply(etmeas[started_etmeas], function(x) length(x) == 0L)
			if (any(redo_etmeas)) {
				what_measure[et_methods %in% names(etmeas)[started_etmeas][redo_etmeas], , ] <- TRUE
			}
		}
		#	- check if some measure methods have only partially completed or were calculated by an outdated version
		if (any(started_etmeas)) {
			todo_neighsXmigs <- sapply(etmeas[started_etmeas], function(x)
									sapply(x$gETmeas, function(b)
										sapply(b, function(m)
													isTRUE(length(m) <= 1) ||
													is.null(m$meta) ||
													isTRUE(m$meta$version < cur_versions[[m$meta$method]]))),
									simplify = "array")
			todo_neighsXmigs <- aperm(todo_neighsXmigs, perm = rev(seq_along(dim(todo_neighsXmigs))))
			if (length(dim(todo_neighsXmigs)) == 3L &&
				identical(dim(todo_neighsXmigs), dim(what_measure[names(etmeas)[started_etmeas], , , drop = FALSE]))) {
				what_measure[names(etmeas)[started_etmeas], , ] <- todo_neighsXmigs
			} else {
				what_measure[names(etmeas)[started_etmeas], , ] <- TRUE
			}
		}
		
		do_new_etmeas <- FALSE
	} else {
		do_new_etmeas <- TRUE
	}

	do_measure <- any(what_measure)
	do_use_template <- apply(what_measure, 1, all)
	do_only_write <- all(!what_measure)

	if (verbose) {
		idh <- 1
		cat("'ecotoner' measuring: tr = ", i, "; start at ", format(t1 <- Sys.time(), format = ""), "\n", sep = "")
	}

	# Measure ecotones
	if (do_measure && file.exists(fname_etlocated(ecotoner_settings, iflag))) {
		load(fname_etlocated(ecotoner_settings, iflag)) #i, b, and etransect loaded

		# Data containers
		template_etobs <- list(	gETmeas = vector(mode = "list", length = neighbors_N(ecotoner_settings)), # list for 'gETmeas' of each neighborhood
								seeds = vector(mode = "list", length = neighbors_N(ecotoner_settings) * length(migtypes)),
								etable = data.frame(matrix(NA, nrow = 1, ncol = 0))) #container combining output table results with rows for each neighborhood and migtype
		temp <- vector("list", length = length(migtypes))
		names(temp) <- migtypes
		template_etobs$gETmeas <- lapply(template_etobs$gETmeas, function(x) temp)

		if (do_new_etmeas) {
			etmeas <- vector(mode = "list", length = length(et_methods))
			names(etmeas) <- et_methods
		}
		
		dir_fig <- file.path(dir_out_fig(ecotoner_settings), iflag)
		dir_create(dir_fig)
		
		for (b in seq_len(neighbors_N(ecotoner_settings))) {
			flag_bfig <- flag_basename(ecotoner_settings, iflag, b)

			#----3. Loop through both versions of vegetation for applying the methods : Veg$AllMigration and Veg$OnlyGoodMigration
			for (im in seq_along(migtypes)) {
				copy_FromMig1_TF <- if (migtypes[im] == "AllMigration") FALSE else !etransect$etbands[[b]]$Veg[["OnlyGoodMigration"]]$diffs_amongMigTypes_TF

				# loop through measurement methods 'et_methods'
				for (etm in et_methods) {
					
					# check whether there is anything to do
					if (what_measure[etm, b, im]) {
						if (verbose) cat("'ecotoner' measuring: tr = ", i, "; neigh = ", b, ": prog: ", idh <- idh + 1, "; mig-type: ", migtypes[im], "; method: ", etm, "\n", sep = "")
					
						if (is.null(etmeas[[etm]]) || (b == 1L && im == 1L && do_use_template[etm])) {
							etmeas[[etm]] <- template_etobs
						}
					
						iseed <- (b - 1) * length(migtypes) + im
						etmeas[[etm]][["seeds"]][[iseed]] <- if (is.null(seed_streams)) NULL else if (inherits(seed_streams, "list")) seed_streams[[((i - 1) * neighbors_N(ecotoner_settings) + (b - 1)) * length(migtypes) + im]] else NA
						set_RNG_stream(etmeas[[etm]][["seeds"]][[iseed]])

#trace(what = ecotoner:::calc_Eppinga2013_stats,
#	tracer = quote({print("calc_Eppinga2013_stats")
#					str(data_T17)
#					print(sum(is.na(data_T17)))}),
#	at = list(c(7, 3, 11, 3, 2)))
#
#trace(what = stats:::na.fail.default,
#	tracer = quote({print("na.fail.default")
#					print(sys.calls()[13:16])
#					str(object)
#					print(which(is.na(object)))}),
#	at = 3)
#
						temp <- try(do.call(what = etm,
											args = list(i = i, b = b, migtype = migtypes[im], 
														ecotoner_settings = ecotoner_settings,
														etband = etransect[["etbands"]][[b]],
														etmeasure = etmeas[[etm]],
														copy_FromMig1_TF = copy_FromMig1_TF,
														do_figures = do_figures,
														dir_fig = dir_fig,
														flag_bfig = flag_bfig,
														seed = NA)),
									silent = TRUE)
				
						if (inherits(temp, "try-error")) {
							warning("'measure_ecotone_per_transect': tr = ", i, "; neigh = ", b, "; mig-type: ", migtypes[im], "; method: ", etm, , ": ", temp, immediate. = TRUE)
						} else {
							etmeas[[etm]] <- temp

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
					}
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
	transect_data2d <- pfun(X, get_transect_grids_as_df, et_desc, ecotoner_settings, migtypes, veg_types)
	colnames(transect_data2d)[3:9] <- c("env", "cols", "rows", paste(rep(paste("y", veg_types, sep = "_"), times = length(migtypes)), rep(migtypes, each = length(veg_types)), sep = "_"))
			
	m1 <- glm(y_Veg1_AllMigration ~ env * Transect_ID + 
									Transect_Azimuth_deg + cols + rows,
				family = binomial, data = transect_data2d)
	X
}

