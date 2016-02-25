#------------------------------------------------------------#
# Functions drawing figures

#' @export
plot_transect_density_profiles <- function(filename, distX1, distX2, elev, dens_Veg1, dens_BSE, dens_Veg2, dens_TF, start1_m, end1_m, start2_m, end2_m) {
	
	pdf(width = 7, height = 4.5, file = filename)
	layout(mat = 1:2, widths = 1, heights = c(8, 1)/9)
	par_old <- par(mar = c(3, 2.5, 0.1, 2.1), mgp = c(1.5, 0.5, 0), cex = cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)
	
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
plot_Danz2012_abruptness <- function(filename, eB_Env, eB_Veg, datFit) {

	pdf(width = 7, height = 8, file = filename)
	par_old <- par(mfrow = c(nr <- 2, 1), mar = c(2.5, 2.5, 0.5, 2.5), mgp = c(1.5, 0.5, 0), cex = cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)


	#Panel (a): Danz et al. 2012 Fig. 3 and 4
	emin <- min(eB_Env$elev$YMeans_ForEachX)
	emax <- max(eB_Env$elev$YMeans_ForEachX)
	
	plot(eB_Env$DistAlongXaxis_m, eB_Env$elev$YMeans_ForEachX, type = "l", lwd = 2, xlab = "Transect length (m)", ylab = "Elevation (m)", axes = FALSE)
	at1 <- c(axTicks(1), max(eB_Env$DistAlongXaxis_m))
	if((temp <- tail(at1, n = 2))[1] > temp[2]) at1 <- c(at1[1:(length(at1)-2)], at1[length(at1)])
	axis(1, pos = emin, at = at1)
	axis(2, pos = 0, at = round(c(emin, axTicks(2), emax)))

	lines(eB_Env$DistAlongXaxis_m, transform_01_to_R(eB_Veg$Veg1$density, emin, emax), col = "red")
	lines(eB_Env$DistAlongXaxis_m, transform_01_to_R(eB_Veg$Veg2$density, emin, emax), col = "darkgreen")
	axis(4, pos = max(at1), at = at <- seq(emin, emax, length = 11), labels = round(transform_R_to_01(at, emin, emax), 1))
	mtext(text = "Vegetation density", side = 4, cex = cex, line = 1)
	mtext(text = "(a)", line = -1, cex = cex, adj = 0.05)
	legend("top", bty = "n", cex = 1/3*cex, legend = c("Elevation", "Big sagebrush ecosystem", "Temperate Forest", "Sigmoidal fit, center, & 90% range"), col = c("black", "red", "darkgreen", "black"), lty = c(1, 1, 1, 2), lwd = c(2, 1, 1, 2))

	#Veg vs Dist: BSE
	if (!inherits(datFit$Veg1VsDist$sigmoidal$m, "try-error")) {
		ypred <- datFit$Veg1VsDist$ymin + predict(datFit$Veg1VsDist$sigmoidal$m) * (datFit$Veg1VsDist$ymax - datFit$Veg1VsDist$ymin)
		lines(eB_Env$DistAlongXaxis_m, transform_01_to_R(ypred, emin, emax), lty = 2, lwd = 2, col = "red")
		line_y <- par()$usr[3] + 0.02*(par()$usr[4]-par()$usr[3])
		arr_x0 <- datFit$Veg1VsDist$sigmoidal$center_x
		arr_y0 <- transform_01_to_R(datFit$Veg1VsDist$ymin + 0.5 * (datFit$Veg1VsDist$ymax - datFit$Veg1VsDist$ymin), emin, emax)
		arrows(x0 = arr_x0, y0 = arr_y0, x1 = arr_x0, y1 = line_y, lwd = 2, col = "red", length = 0.03*par()$pin[1]/nr)
		seg_x0 <- datFit$Veg1VsDist$sigmoidal$x_5percOfy
		seg_y0 <- transform_01_to_R(datFit$Veg1VsDist$ymin + 0.05 * (datFit$Veg1VsDist$ymax - datFit$Veg1VsDist$ymin), emin, emax)
		seg_x1 <- datFit$Veg1VsDist$sigmoidal$x_95percOfy
		seg_y1 <- transform_01_to_R(datFit$Veg1VsDist$ymin + 0.95 * (datFit$Veg1VsDist$ymax - datFit$Veg1VsDist$ymin), emin, emax)
		segments(x0 = seg_x0, y0 = line_y, x1 = seg_x1, y1 = line_y, lwd = 2, col = "red")
		segments(x0 = seg_x0, y0 = seg_y0, x1 = seg_x0, y1 = line_y, lwd = 2, col = "red")
		segments(x0 = seg_x1, y0 = seg_y1, x1 = seg_x1, y1 = line_y, lwd = 2, col = "red")
	}
	#Veg vs Dist: TF
	if (!inherits(datFit$Veg2VsDist$sigmoidal$m, "try-error")) {
		ypred <- datFit$Veg2VsDist$ymin + predict(datFit$Veg2VsDist$sigmoidal$m) * (datFit$Veg2VsDist$ymax - datFit$Veg2VsDist$ymin)
		lines(eB_Env$DistAlongXaxis_m, transform_01_to_R(ypred, emin, emax), lty = 2, lwd = 2, col = "darkgreen")
		line_y <- par()$usr[3] + 0.03*(par()$usr[4]-par()$usr[3])
		arr_x0 <- datFit$Veg2VsDist$sigmoidal$center_x
		arr_y0 <- transform_01_to_R(datFit$Veg2VsDist$ymin + 0.5 * (datFit$Veg2VsDist$ymax - datFit$Veg2VsDist$ymin), emin, emax)
		arrows(x0 = arr_x0, y0 = arr_y0, x1 = arr_x0, y1 = line_y, lwd = 2, col = "darkgreen", length = 0.03*par()$pin[1]/nr)
		seg_x0 <- datFit$Veg2VsDist$sigmoidal$x_5percOfy
		seg_y0 <- transform_01_to_R(datFit$Veg2VsDist$ymin + 0.05 * (datFit$Veg2VsDist$ymax - datFit$Veg2VsDist$ymin), emin, emax)
		seg_x1 <- datFit$Veg2VsDist$sigmoidal$x_95percOfy
		seg_y1 <- transform_01_to_R(datFit$Veg2VsDist$ymin + 0.95 * (datFit$Veg2VsDist$ymax - datFit$Veg2VsDist$ymin), emin, emax)
		segments(x0 = seg_x0, y0 = line_y, x1 = seg_x1, y1 = line_y, lwd = 2, col = "darkgreen")
		segments(x0 = seg_x0, y0 = seg_y0, x1 = seg_x0, y1 = line_y, lwd = 2, col = "darkgreen")
		segments(x0 = seg_x1, y0 = seg_y1, x1 = seg_x1, y1 = line_y, lwd = 2, col = "darkgreen")
	}

	#Panel (b): Danz et al. 2012 Fig. 3 inset and Fig. 6
	#Veg vs Elev: BSE
	plot(sort(eB_Env$elev$YMeans_ForEachX), eB_Veg$Veg1$density[order(eB_Env$elev$YMeans_ForEachX)], xlim = c(min(eB_Env$elev$YMeans_ForEachX), max(eB_Env$elev$YMeans_ForEachX)), ylim = c(0, 1), type = "l", col = "red", xlab = "Elevation (m)", ylab = "Vegetation density", axes = FALSE)
	at2 <- c(emin, axTicks(1), emax)
	if ((temp <- tail(at2, n = 2))[1] > temp[2]) at2 <- c(at2[1:(length(at2)-2)], at2[length(at2)])
	if ((temp <- head(at2, n = 2))[1] > temp[2]) at2 <- c(at2[1], at2[3:length(at2)])
	axis(1, pos = 0, at = round(at2))
	axis(2, pos = emin)

	if (!inherits(datFit$Veg1VsElev$sigmoidal$m, "try-error")) {
		ypred <- datFit$Veg1VsElev$ymin + predict(datFit$Veg1VsElev$sigmoidal$m) * (datFit$Veg1VsElev$ymax - datFit$Veg1VsElev$ymin)
		lines(eB_Env$elev$YMeans_ForEachX, ypred, col = "red", lty = 2, lwd = 2)
		line_y <- par()$usr[3] + 0.02*(par()$usr[4]-par()$usr[3])
		arr_x0 <- datFit$Veg1VsElev$sigmoidal$center_x
		arr_y0 <- datFit$Veg1VsElev$ymin + 0.5 * (datFit$Veg1VsElev$ymax - datFit$Veg1VsElev$ymin)
		arrows(x0 = arr_x0, y0 = arr_y0, x1 = arr_x0, y1 = line_y, lwd = 2, col = "red", length = 0.03*par()$pin[1]/nr)
		seg_x0 <- datFit$Veg1VsElev$sigmoidal$x_5percOfy
		seg_y0 <- datFit$Veg1VsElev$ymin + 0.05 * (datFit$Veg1VsElev$ymax - datFit$Veg1VsElev$ymin)
		seg_x1 <- datFit$Veg1VsElev$sigmoidal$x_95percOfy
		seg_y1 <- datFit$Veg1VsElev$ymin + 0.95 * (datFit$Veg1VsElev$ymax - datFit$Veg1VsElev$ymin)
		segments(x0 = seg_x0, y0 = line_y, x1 = seg_x1, y1 = line_y, lwd = 2, col = "red")
		segments(x0 = seg_x0, y0 = seg_y0, x1 = seg_x0, y1 = line_y, lwd = 2, col = "red")
		segments(x0 = seg_x1, y0 = seg_y1, x1 = seg_x1, y1 = line_y, lwd = 2, col = "red")
	}
	
	#Veg vs Elev: TF
	lines(sort(eB_Env$elev$YMeans_ForEachX), eB_Veg$Veg2$density[order(eB_Env$elev$YMeans_ForEachX)], col = "darkgreen")
	
	if (!inherits(datFit$Veg2VsElev$sigmoidal$m, "try-error")) {
		ypred <- datFit$Veg2VsElev$ymin + predict(datFit$Veg2VsElev$sigmoidal$m) * (datFit$Veg2VsElev$ymax - datFit$Veg2VsElev$ymin)
		lines(eB_Env$elev$YMeans_ForEachX, ypred, col = "darkgreen", lty = 2, lwd = 2)
		line_y <- par()$usr[3] + 0.03*(par()$usr[4]-par()$usr[3])
		arr_x0 <- datFit$Veg2VsElev$sigmoidal$center_x
		arr_y0 <- datFit$Veg2VsElev$ymin + 0.5 * (datFit$Veg2VsElev$ymax - datFit$Veg2VsElev$ymin)
		arrows(x0 = arr_x0, y0 = arr_y0, x1 = arr_x0, y1 = line_y, lwd = 2, col = "darkgreen", length = 0.03*par()$pin[1]/nr)
		seg_x0 <- datFit$Veg2VsElev$sigmoidal$x_5percOfy
		seg_y0 <- datFit$Veg2VsElev$ymin + 0.05 * (datFit$Veg2VsElev$ymax - datFit$Veg2VsElev$ymin)
		seg_x1 <- datFit$Veg2VsElev$sigmoidal$x_95percOfy
		seg_y1 <- datFit$Veg2VsElev$ymin + 0.95 * (datFit$Veg2VsElev$ymax - datFit$Veg2VsElev$ymin)
		segments(x0 = seg_x0, y0 = line_y, x1 = seg_x1, y1 = line_y, lwd = 2, col = "darkgreen")
		segments(x0 = seg_x0, y0 = seg_y0, x1 = seg_x0, y1 = line_y, lwd = 2, col = "darkgreen")
		segments(x0 = seg_x1, y0 = seg_y1, x1 = seg_x1, y1 = line_y, lwd = 2, col = "darkgreen")
	}

	mtext(text = "(b)", line = -0.5, cex = cex, adj = 0.05)
		
	invisible()
}

#' @export
map_front_runners_Eppinga2013 <- function(filename, eB_Env, eB_Veg, datFit) {
	pdf(width=7, height=7, file=filename)
	par_old <- par(mfrow=c(1, 1), mar=c(0, 0.1, 1, 1), mgp=c(2, 0.5, 0), cex=cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)

	ext1 <- raster::extent(eB_Env$elev$grid)
	xlim <- c(-1000, ext1@xmax)
	ylim <- c(-1000, ext1@ymax+1000)
	
	raster::image(eB_Env$elev$grid, col=gray(0:255/255), xlim=xlim, ylim=ylim, main="", xlab="", ylab="", asp=1, axes=FALSE)
	raster::image(eB_Veg$Veg1$grid, col=adjustcolor("red", alpha.f = 0.3), add=TRUE)
	raster::image(eB_Veg$Veg2$grid, col=adjustcolor("darkgreen", alpha.f = 0.3), add=TRUE)
	atx <- c((atx <- axTicks(1))[atx >= 0 & atx < ext1@xmax], ext1@xmax)
	axis(1, pos=0, at=atx)
	axis(2, pos=0, at=c(0, 2000, 4000, 6000))
	text(x=ext1@xmax/2, y=-strheight("0", units="user", cex=cex)*(0.5+2), labels="Transect length (m)")
	text(x=-strwidth("0", units="user", cex=cex)*(0.5+2.5), y=ext1@ymax/2, labels="Transect width (m)", srt=90)

	isnotna <- !is.na(datFit$adv_veg1$deltaFrontRunners_m)
	res_m <- raster::xres(eB_Env$elev$grid)
	points(x=((opt <- datFit$optim$pos_m) + datFit$adv_veg1$deltaFrontRunners_m)[isnotna], y=(ys <- (length(datFit$adv_veg1$deltaFrontRunners_m):1) * res_m)[isnotna], pch=".", col="magenta")
	isnotna <- !is.na(datFit$adv_veg2$deltaFrontRunners_m)
	points(x=(opt + datFit$adv_veg2$deltaFrontRunners_m)[isnotna], y=ys[isnotna], pch=".", col="green")
	segments(x0=opt, y0=ext1@ymin, x1=opt, y1=ext1@ymax, lwd=2, col="yellow")
	segments(x0=pos1 <- opt + datFit$adv_veg1$deltaFrontRunners_Mean_m, y0=ext1@ymin, x1=pos1, y1=ext1@ymax, lwd=2, col=adjustcolor("magenta", alpha.f = 0.7))
	segments(x0=pos2 <- opt + datFit$adv_veg2$deltaFrontRunners_Mean_m, y0=ext1@ymin, x1=pos2, y1=ext1@ymax, lwd=2, col=adjustcolor("green", alpha.f = 0.7))

	if (datFit$adv_stats$FrontsAdvBeyondOptBoundary) {
		p12 <- datFit$adv_stats$Front1Larger2_p < 0.05
		p21 <- datFit$adv_stats$Front2Larger1_p < 0.05
		if (p12) {
			arrows(x0=opt, y0=ext1@ymax+30, x1=pos1, y1=ext1@ymax+30, lwd=2, col="magenta")
			ptext <- paste0("adv(Veg1; BSE) > adv(Veg2; TF):\np < ", ifelse(datFit$adv_stats$Front1Larger2_p > 0.001, formatC(datFit$adv_stats$Front1Larger2_p, format="f", digits=3), "0.001"))
		}
		if (p21) {
			arrows(x0=opt, y0=ext1@ymax+30, x1=pos2, y1=ext1@ymax+30, lwd=2, col="green")
			ptext <- paste0("adv(Veg2; TF) > adv(Veg1; BSE):\np < ", ifelse(datFit$adv_stats$Front2Larger1_p > 0.001, formatC(datFit$adv_stats$Front2Larger1_p, format="f", digits=3), "0.001"))
		}
		if (!p12 & !p21) ptext <- paste0("adv(Veg1; BSE) = adv(Veg2; TF):\np > ", formatC(datFit$adv_stats$Front1Larger2_p, format="f", digits=3))
	} else {
		ptext <- paste0("adv() < optimal boundary")
	}
	text(x=ext1@xmax/2, y=ext1@ymax+strheight("0", units="user", cex=2/3*cex)*2, labels=ptext, cex=2/3)
	
	invisible()
}

#' @export
plot_Gastner2010_hulledge <- function(filename, eB_Env, eB_Veg, datFit) {
	pdf(width=7, height=10, file=filename)
	par_old <- par(mfrow=c(nr <- 2, 1), mar=c(0, 0.1, 1, 1), mgp=c(1.5, 0.5, 0), cex=cex <- 1.5)
	on.exit({par(par_old); dev.off()}, add = TRUE)


	#Panel a: map
	ext1 <- raster::extent(eB_Env$elev$grid)
	xlim <- c(-1000, ext1@xmax)
	ylim <- c(-1000, ext1@ymax)
	
	raster::image(eB_Env$elev$grid, col=gray(0:255/255), xlim=xlim, ylim=ylim, main="", xlab="", ylab="", asp=1, axes=FALSE)
	raster::image(eB_Veg$Veg1$grid, col=adjustcolor("red", alpha.f = 0.3), add=TRUE)
	raster::image(datFit$Veg1$grid_LargestPatch, col=adjustcolor("red", alpha.f = 0.3), add=TRUE)
	raster::image(eB_Veg$Veg2$grid, col=adjustcolor("darkgreen", alpha.f = 0.3), add=TRUE)
	raster::image(datFit$Veg2$grid_LargestPatch, col=adjustcolor("darkgreen", alpha.f = 0.3), add=TRUE)
	atx <- c((atx <- axTicks(1))[atx >= 0 & atx < ext1@xmax], ext1@xmax)
	axis(1, pos=0, at=atx)
	axis(2, pos=0, at=c(0, 2000, 4000, 6000))
	text(x=ext1@xmax/2, y=-strheight("0", units="user", cex=cex)*(0.5+2), labels="Transect length (m)", xpd=NA)
	text(x=-strwidth("0", units="user", cex=cex)*(0.5+2.5), y=ext1@ymax/2, labels="Transect width (m)", srt=90)

	lines(x=coordinates(datFit$Veg1$spLine_hullEdge)[[1]][[1]], col="orange")
	lines(x=coordinates(datFit$Veg2$spLine_hullEdge)[[1]][[1]], col="green")
	mtext(text="(a)", line=-1, cex=cex, adj=0.01)

	#Panel b: densities
	par(mar=c(2.5, 2.5, 1, 1.5))
	#Veg1
	plot(eB_Env$DistAlongXaxis_m, eB_Veg$Veg1$density, lty=1, col="red", type="l", xlim=c(0, ext1@xmax), ylim=c(0, 1), xlab="Transect length (m)", ylab="Density", xaxs="i", axes=FALSE)
	atx <- c((atx <- axTicks(1))[atx >= 0 & atx < ext1@xmax], ext1@xmax)
	axis(1, pos=0, at=atx)
	axis(2, pos=0)
	include <- datFit$Veg1$density > 0
	lines(eB_Env$DistAlongXaxis_m[include], datFit$Veg1$density[include], col="red", lwd=2)
	arrows(x0=pos <- datFit$Veg1$stats$position_m, y0=ytemp <- 1, x1=pos, y1=0, lwd=2, col="orange", length=0.03*par()$pin[1]/nr)
	arrows(x0=pos - datFit$Veg1$stats$width_m, y0=ytemp, x1=pos + datFit$Veg1$stats$width_m, y1=ytemp, lwd=2, col="orange", length=0.03*par()$pin[1]/nr, code=3)
	#Veg2
	lines(eB_Env$DistAlongXaxis_m, eB_Veg$Veg2$density, lty=1, col="darkgreen")
	include <- datFit$Veg2$density > 0
	lines(eB_Env$DistAlongXaxis_m[include], datFit$Veg2$density[include], col="darkgreen", lwd=2)
	arrows(x0=pos <- datFit$Veg2$stats$position_m, y0=ytemp-0.01, x1=pos, y1=0, lwd=2, col="green", length=0.03*par()$pin[1]/nr)
	arrows(x0=pos - datFit$Veg2$stats$width_m, y0=ytemp-0.01, x1=pos + datFit$Veg2$stats$width_m, y1=ytemp-0.01, lwd=2, col="green", length=0.03*par()$pin[1]/nr, code=3)

	mtext(text="(b)", line=-0.5, cex=cex, adj=0.01)

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

