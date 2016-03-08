#' export
measure_ecotones_all_transects <- function(x) {return(NULL)}

#' export
measure_ecotone_per_transect <- function(i, ecotoner_settings, et_methods, verbose = TRUE, do_figures = TRUE) {
	if (verbose) {
		idh <- 1
		cat("'ecotoner' measuring: tr = ", i, "; start at ", format(t1 <- Sys.time(), format = ""), "\n", sep = "")
	}
	
	# Output containers
	iflag <- flag_itransect(i, ecotoner_settings)
	
	temp <- vector("list", length = length(get("migtypes", envir = etr_vars)))
	names(temp) <- get("migtypes", envir = etr_vars)
	template_etobs <- list(gETmeas = vector(mode = "list", length = neighbors_N(ecotoner_settings)), # list for 'gETmeas' of each neighborhood
					etable = data.frame(matrix(NA, nrow=1, ncol=0))) #container combining output table results from each neighborhood and migtype
	template_etobs$gETmeas <- lapply(template_etobs$gETmeas, function(x) temp)
	etmeas <- vector(mode = "list", length = length(et_methods))
	names(etmeas) <- et_methods
	
	do_measure <- TRUE
	if (file.exists(fname_etmeasured(iflag, ecotoner_settings))) {
		load(ftemp3) #load: i, b, etmeas
		itemp <- et_methods %in% names(etmeas)
#TODO(drs): check not only presence of methods, but also whether all neighborhoods and migration types have been completed
		if (all(itemp)) {
			do_measure <- FALSE  # all done
		} else {
			et_methods <- et_methods[!(itemp)]
		}
	}
	
	if (file.exists(fname_etlocated(iflag, ecotoner_settings))) {
		load(fname_etlocated(iflag, ecotoner_settings)) #i, b, and etransect loaded
	} else {
		do_measure <- FALSE # no suitable transect located for search point i
	}
	
	if (do_measure) {
		for (im in et_methods) etmeas[[im]] <- template_etobs

		dir_fig <- file.path(dir_out_fig(ecotoner_settings), iflag)
		dir.create(dir_fig, showWarnings = FALSE)

		seed <- if (reseed(ecotoner_settings)) get_pseed(ecotoner_settings) else NA

		for (b in seq_len(neighbors_N(ecotoner_settings))) {
			flag_bfig <- flag_basename(iflag, b, ecotoner_settings)

		
#TODO(drs):	eB_InterZone <- eB_PatchSizeDist <- template_etobs	#ecological boundaries

			#----3. Loop through both versions of vegetation for applying the methods : Veg$AllMigration and Veg$OnlyGoodMigration
			for (migtype in get("migtypes", envir = etr_vars)) {
				copy_FromMig1_TF <- if (migtype == "AllMigration") FALSE else !etransect$etbands[[b]]$Veg[["OnlyGoodMigration"]]$diffs_amongMigTypes_TF

				if (verbose) cat("'ecotoner' measuring: tr = ", i, "; neigh = ", b, ": prog: ", idh <- idh + 1, "; mig-type: ", migtype, "; method: ", etm, "\n", sep = "")
				
				# loop through measurement methods 'et_methods'
				for (etm in et_methods) {
					etmeas[[etm]] <- do.call(what = etm, args = list(i = i, b = b, migtype = migtype, 
																		ecotoner_settings = ecotoner_settings,
																		etband = etransect[["etbands"]][[b]],
																		etmeasure = etmeas[[etm]],
																		flag_bfig = flag_bfig,
																		copy_FromMig1_TF = copy_FromMig1_TF,
																		do_figures = do_figures,
																		seed = seed))

					#Temporarily save data to disk file
					if (all(sapply(etmeas[[etm]]$etable[b, ], function(x) class(x)) %in% c("numeric", "logical", "integer", "character"))) {
#TODO(drs): determine XXX based on measurement method
						appT <- file.exists(XXX)
						write.table(etmeas[[etm]]$etable[b, ], file = XXX,
									append = appT, sep = ",", dec = ".", qmethod = "double",
									row.names = FALSE, col.names = !appT)
					} else {
						warning("'measure_ecotone_per_transect': ", i, " 'etmeas[[", etm, "]]$etable[", b, ", ]' contains unsuitable elements, e.g., a 'list'", immediate. = TRUE)
					}

					# Save data to RData disk file
					save(i, b, migtype, etmeas, file=fname_etmeasured(iflag, ecotoner_settings))
				}
			} # end loop through migtypes
		} # end loop through neighborhoods
	}
		
	if (verbose) cat("'ecotoner' measuring: tr = ", i, " completed ", format(t2 <- Sys.time(), format = ""), " after ", round(difftime(t2, t1, units="hours"), 2), " hours\n", sep = "")

	i
}	
