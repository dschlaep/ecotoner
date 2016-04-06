#------------------------------------------------------------#
# Functions drawing figures

#' @export
plot_transect_density_profiles <- function(filename, distX1, distX2, elev, dens_Veg1, dens_BSE, dens_Veg2, dens_TF, start1_m, end1_m, start2_m, end2_m) {
	
	pdf(width = 7, height = 4.5, file = filename)
	layout(mat = 1:2, widths = 1, heights = c(8, 1)/9)
	par_old <- par(mar = c(3, 2.5, 0.1, 2.1), mgp = c(1.5, 0.5, 0), cex = cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)

	panel_transect_density_profiles(distX1, distX2, elev, dens_Veg1, dens_BSE, dens_Veg2, dens_TF, start1_m, end1_m, start2_m, end2_m, cex)

	par(mar = c(0, 2.5, 0, 0))
	plot.new()
	legend("left", inset = 0.05, bty = "n", cex = 2/4*cex, ncol = 3,
		legend = c("Elevation", "Transect region", "All BSE", "Majority BSE", "All TF", "Majority TF"),
		col = c("black", "black", "darkgoldenrod1", "red", "aquamarine", "darkgreen"),
		lwd = c(2, -1, 1, 1, 1, 1),
		pch = c(-1, 22, -1, -1, -1, -1),
		pt.cex = c(-1, 2, -1, -1, -1, -1),
		pt.bg = "gray")

	invisible()
}
	
panel_transect_density_profiles <- function(distX1, distX2, elev, dens_Veg1, dens_BSE, dens_Veg2, dens_TF, start1_m, end1_m, start2_m, end2_m, cex) {
	emin <- min(elev)
	emax <- max(elev)
	
	plot(distX1, elev, type = "l", lwd = 2, xlab = "Transect length (m)", ylab = "Elevation (m)", axes = FALSE)
	at1 <- c(axTicks(1), max(distX1))
	if ((temp <- tail(at1, n = 2))[1] > temp[2]) at1 <- c(at1[1:(length(at1)-2)], at1[length(at1)])
	axis(1, pos = emin, at = at1)
	axis(2, pos = 0, at = round(c(emin, axTicks(2), emax)))
	rect(xleft = start1_m, ybottom = emin, xright = end1_m, ytop = emax, col = adjustcolor(gray(0.6), alpha.f = 0.3), border = "black", lty = 2)
	rect(xleft = start2_m, ybottom = emin, xright = end2_m, ytop = emax, col = adjustcolor(gray(0.6), alpha.f = 0.5), border = "black", lty = 1)
	lines(distX1, transform_01_to_R(dens_BSE, emin, emax), col = "darkgoldenrod1")
	lines(distX2, transform_01_to_R(dens_Veg1, emin, emax), col = "red")
	lines(distX1, transform_01_to_R(dens_TF, emin, emax), col = "aquamarine")
	lines(distX2, transform_01_to_R(dens_Veg2, emin, emax), col = "darkgreen")
	axis(4, pos = max(at1), at = at2 <- seq(emin, emax, length = 11), labels = round(transform_R_to_01(at2, emin, emax), 1))
	mtext(text = "Vegetation density", side = 4, cex = cex, line = 1)
	
	invisible()
}

#' @export
map_transect <- function(filename, stline_pts, etline_pts, gB_Env, elevCropped, bse, tf, Veg1, Veg2, human) {
	pdf(width = 7, height = 7, file = filename)
	par_old <- par(mfrow = c(1, 1), mar = c(4, 2, 0.5, 0), mgp = c(2, 0.5, 0), cex = 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)

	etemp <- raster::extent(etline_pts)
	xex <- yex <- 1
	if ((temp <- etemp@ymax - etemp@ymin) < 0.1 && temp > 0) yex <- 0.1 / temp
	if ((temp <- etemp@xmax - etemp@xmin) < 0.1 && temp > 0) xex <- 0.1 / temp
	etemp <- 1.1 * c(xex, yex) * etemp
	etemp <- raster::union(etemp, raster::extent(elevCropped))
	
	raster::plot(elevCropped, ext = etemp, col = gray(0:255/255), bigplot = c(0.1, 0.89, 0.2, 0.99), smallplot = c(0.1, 0.89, 0.10, 0.12), horizontal = TRUE)

	raster::image(bse, ext = etemp, col = adjustcolor("darkgoldenrod1", alpha.f = 0.5), add = TRUE)
	raster::image(Veg1, ext = etemp, col = adjustcolor("red", alpha.f = 0.3), add = TRUE)
	raster::image(tf, ext = etemp, col = adjustcolor("aquamarine", alpha.f = 0.5), add = TRUE)
	raster::image(Veg2, ext = etemp, col = adjustcolor("darkgreen", alpha.f = 0.3), add = TRUE)
	
	if (raster::cellStats(human, 'sum') > 0) raster::image(human, ext = etemp, col = adjustcolor("blue1", alpha.f = 0.3), add = TRUE)

	points(stline_pts, pch = ".", col = "blue")
	points(etline_pts, pch = ".", col = "yellow")

	raster::plot(sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(sp::rbind.SpatialPoints(gB_Env$band.pts, gB_Env$band.pts[1, ]))), ID=0)), proj4string=raster::crs(gB_Env$elev$grid)), pch = 16, border = "yellow", add = TRUE)
	
	invisible()
}

#' @export
map_flowpaths <- function(filename, elev, bse, bse_remain, paths_bse, tf, tf_remain, paths_tf) {
	removed <- function(full, part) raster::overlay(full, part, fun=function(full, part) ifelse(!is.na(full) & is.na(part), 1, NA))

	pdf(width=7, height=7, file=filename)
	par_old <- par(mfrow=c(1, 1), mar=c(0, 0.1, 1, 1), mgp=c(2, 0.5, 0), cex=cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)

	ext1 <- raster::extent(elev)
	xlim <- c(-1000, ext1@xmax)
	ylim <- c(-1000, ext1@ymax)

	raster::image(elev, col=gray(100:255/255), xlim=xlim, ylim=ylim, main="", xlab="", ylab="", asp=1, axes=FALSE)
	atx <- c((atx <- axTicks(1))[atx >= 0 & atx < ext1@xmax], ext1@xmax)
	axis(1, pos=0, at=atx)
	axis(2, pos=0, at=c(0, 2000, 4000, 6000))
	text(x=ext1@xmax/2, y=-strheight("0", units="user", cex=cex)*(0.5+2), labels="Transect length (m)")
	text(x=-strwidth("0", units="user", cex=cex)*(0.5+2.5), y=ext1@ymax/2, labels="Transect width (m)", srt=90)

	temp <- removed(bse, bse_remain)
	if (raster::cellStats(temp, 'sum') > 0) raster::image(temp, col=adjustcolor("yellow", alpha.f = 0.3), maxpixel=raster::ncell(temp), add=TRUE)
	raster::image(bse_remain, col=adjustcolor("red", alpha.f = 0.3), maxpixel=raster::ncell(bse_remain), add=TRUE)
	temp <- removed(tf, tf_remain)
	if (raster::cellStats(temp, 'sum') > 0) raster::image(temp, col=adjustcolor("blue", alpha.f = 0.3), maxpixel=raster::ncell(temp), add=TRUE)
	raster::image(tf_remain, col=adjustcolor("darkgreen", alpha.f = 0.3), maxpixel=raster::ncell(tf_remain), add=TRUE)
	
	temp <- lapply(paths_bse, FUN=function(p) points(p, pch=".", type="l", col="darkred"))
	temp <- lapply(paths_tf, FUN=function(p) points(p, pch=".", col="darkgreen"))
	
	invisible()
}

#' @export
plot_interaction_zones <- function(filename, datLine, eB_Env, eB_Veg, datFit) {
	iz_types <- names(datFit)
	iz_N <- length(iz_types)
	cols <- rainbow(n=iz_N+1)
	
	pdf(width=7, height=8, file=filename)
	par_old <- par(mfrow=c(nr <- 2, 1), mar=c(0, 2.5, 0.5, 1.5), mgp=c(1.5, 0.5, 0), cex=cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)


	#Panel a: vegetation and edge densities
	#Veg1
	plot(eB_Env$DistAlongXaxis_m, eB_Veg$Veg1$density, col="red", lwd=1, lty=1, type="l", ylim=c(0, 1), xlab="Transect length (m)", ylab="Density", xaxs="i", axes=FALSE)
	axis(1, at=c(0, axTicks(1), max(eB_Env$DistAlongXaxis_m)), labels=FALSE, xpd=TRUE, pos=0)
	axis(2, pos=0)
	#Veg2
	lines(eB_Env$DistAlongXaxis_m, eB_Veg$Veg2$density, col="darkgreen", lwd=1, lty=1)
	#Human footprint
	hf_rle <- rle(eB_Env$human$density > 0)
	if ((lt <- sum(hf_rle$values)) > 0) {
		hf_istarts <- 1 + (itemp <- c(0, cumsum(hf_rle$lengths)))[hf_rle$values]
		hf_iends <- itemp[1 + which(hf_rle$values)]
		for (i in seq_len(lt)) {
			i_hf <- c(max(1, hf_istarts[i]-1), hf_istarts[i]:hf_iends[i], min(hf_iends[i]+1, length(eB_Env$human$density)))
			lines(eB_Env$DistAlongXaxis_m[i_hf], eB_Env$human$density[i_hf], col="blue1", lwd=1, lty=1)
		}
	}
	
	legend("topleft", bty="n", legend=c("BSE density", "TF density", "Human footprint"), col=c("red", "darkgreen", "blue1"), lty=c(1, 1, 1), cex=2/3)
	mtext(text="(a)", line=-0.5, cex=cex, adj=0.01)

	#Panel b: locations of interaction zones
	par(mar=c(2.5, 2.5, 0, 1.5))
	plot(eB_Env$DistAlongXaxis_m, NULL, xlim=c(0, max(eB_Env$DistAlongXaxis_m)), ylim=c(0, 1), xlab="Transect length (m)", ylab="Methods", xaxs="i", axes=FALSE)
	axis(1, at=c(0, axTicks(1), max(eB_Env$DistAlongXaxis_m)), xpd=TRUE, pos=0)
	
	for (iz in seq_len(iz_N)) {
		if (!is.na(datFit[[iz]]$start_m) && !is.na(datFit[[iz]]$end_m && datFit[[iz]]$start_m != datFit[[iz]]$end_m)) {
			arrows(x0=datFit[[iz]]$start_m, y0=iz/(iz_N+1), x1=datFit[[iz]]$end_m, y1=iz/(iz_N+1), col=cols[iz], length=0.03*par()$pin[1]/nr, lwd=2, code=3)
			text(x=datFit[[iz]]$start_m+(datFit[[iz]]$end_m-datFit[[iz]]$start_m)/2, y=iz/(iz_N+1)-strheight("0", units="user", cex=2/3*cex), labels=iz_types[iz], cex=2/3*cex)
		}
	}

	mtext(text="(b)", line=-1.5, cex=cex, adj=0.01)

	invisible()
}

