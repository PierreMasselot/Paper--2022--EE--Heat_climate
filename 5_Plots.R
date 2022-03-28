#################################################################
#
#           AUDACE - Summer mortality prediction
#                        Figures
#
#################################################################

#------------------------
# Parameters
#------------------------

# Considered colors for heat percentiles
colvec <- c("black", rev(heat_hcl(length(heatper) + 1))[-1])
colclim <- 4

# Colors for areas
colmet <- c("forestgreen", "cornflowerblue")

#------------------------
# Map of regions
#------------------------

#----- Get map info
# Get Canadian boundaries
provs <- st_as_sf(states50)
canprovs <- subset(provs, admin == "Canada")
canbox <- st_bbox(canprovs)

# Get maps
cityboxes <- lapply(geolims, function(x){
  b <- st_bbox(x)
  names(b) <- c("left", "bottom", "right", "top")
  b
})
bgmap <- lapply(cityboxes, get_stamenmap,
  maptype = "terrain-background", force = T, zoom = 10)

# Compute areas of municipalities
areas <- lapply(geolims, st_area)
biggest <- lapply(areas, function(x)
  order(x, decreasing = T)[1:3])

#----- Plot with patchwork

# Canadian provinces plot
pcan <- ggplot(canprovs) + theme_void() + 
  geom_sf() + 
  geom_label(label = "CANADA", x = mean(canbox[c(1,3)]), 
    y = mean(canbox[c(2,4)]), size = 10) + 
  geom_sf_label(aes(label = ifelse(name == "Québec", name, NA)), size = 7) +
  geom_sf(data = geocenter, col = colmet, 
    size = 5, pch = 22, bg = NA, stroke = 2)

# Plot metropolitan areas
pmet <- list()
for (i in seq_along(region)){
  pmet[[i]] <- ggmap(bgmap[[i]]) + theme_void() + 
    theme(panel.border = element_rect(colour = colmet[i], 
      size = 5, fill = NA)) + 
    geom_sf(data = geolims[[i]], fill = NA, size = 1, inherit.aes = FALSE,
      col = grey(.2)) +
    geom_sf_label(data = geolims[[i]][biggest[[i]],], 
      mapping = aes(label = MUS_NM_MRC), inherit.aes = FALSE, size = 5) + 
    ggtitle(reg_lab[i]) + 
    theme(plot.title = element_text(colour = colmet[i], size = rel(3)))
}

# Put plot together
design <- "AAA
  B#C"
patch <- wrap_plots(A = pcan, B = pmet[[1]], C = pmet[[2]],
  design = design) + plot_layout(width = c(4, 1, 4), height = c(2, 1))

# Save
ggsave(patch, filename = "Figures/Figure1.pdf", 
  units = "in", height = 20, width = 15, device = "pdf")

#------------------------
# Index vs AF
#------------------------

# Choose which to plot
plotclim <- "AMO"
plotmonths <- 1:3
plotper <- "mmt"

# Years to highlight
highyears <- list(1983:1984, 1987:1988, 1998:1999, 2005, 2010:2013)

# AF matrix
afmat <- do.call(cbind, datatab[grep("AF", names(datatab))])

# Extract average of period for climate index
climplot <- scale(rowMeans(datatab[[plotclim]][,plotmonths])) 

# Plot AFs
x11()
layout(matrix(c(1:3, 0), 2, 2), width = c(4, 1))
par(mar = c(5, 4, 2, 1) + .1)
matplot(datatab$Year, afmat, col = colvec, type = "b",
  pch = rep(15:16, each = length(heatper) + 1), xlab = "Year", ylab = "AF")
for (i in seq_along(highyears)){
  rect(highyears[[i]][1] - .5, par("usr")[3], 
    tail(highyears[[i]], 1) + .5, par("usr")[4],
    border = NA, col = adjustcolor("grey", .5))
}

# plot climate index
plot(datatab$Year, climplot, type = "b", col = colclim, pch = 16, lwd = 2,
  xlab = "Year", ylab = plotclim)
for (i in seq_along(highyears)){
  rect(highyears[[i]][1] - .5, par("usr")[3], 
    tail(highyears[[i]], 1) + .5, par("usr")[4],
    border = NA, col = adjustcolor("grey", .5))
}

# Add legend
par(mar = c(5, 0, 4, 0) + .1)
plot.new()
lg <- legend("topleft", legend = c("MMT", heatper * 100), 
  fill = c(colvec), xpd = NA, bty = "n", title = "AF definition")
legend(lg$rect$left, lg$rect$top - lg$rect$h, reg_lab, 
  col = 1, pch = 15:16, xpd = NA, bty = "n")

# Save
dev.print(pdf, file = "Figures/Figure2.pdf")

#------------------------
# Second-stage functional regression result
#------------------------

#----- Parameters
# Check which indices should be plotted
sel_bl <- unique(unlist(sel_vars))

# Choose which to plot
plotclim <- "AMO"
plotind <- which(climate_labs == plotclim)

# Resolution of curve
curres <- 40
sgrid <- seq(clag[1] - 1, tail(clag, 1), length.out = 40)

#----- Plot

x11(width = 8)
# Loop on regions
layout(cbind(1:2, 3), widths = c(4, 1))
for (i in seq_along(region)){

  # Get the functional coefficient for each heat percentile
  betai <- sapply(reslist[[i]], function(x) {
    coef(x$main$final_mod, n1 = curres)$smterms[[plotind]]$value
  })
  
  # Get confidence intervals
  betaCIi <- sapply(reslist[[i]], function(x){
      bootcurves <- x$main$modCI$raw_results[[plotind + 1]]
      bootcurves <- split(bootcurves, 
        rep(1:B, each = length(attr(bootcurves, "x"))))
      bootcurves <- do.call(cbind, bootcurves)
      apply(bootcurves, 1, quantile, c(.025, .975))},
    simplify = "array")
  
  # Put them on comparable scale
  afinds <- grep(sprintf("AF_.*_%s", region[i]), names(datatab))
  sdafs <- sapply(datatab[afinds], sd)
  betai <- betai %*% diag(1 / sdafs)
  betaCIi[1,,] <- betaCIi[1,,] %*% diag(1 / sdafs)
  betaCIi[2,,] <- betaCIi[2,,] %*% diag(1 / sdafs)
  
  # Plot
  matplot(sgrid, betai, type = "l", lwd = 3, main = reg_lab[i],
    xlab = "Month", ylab = "AF change", xaxt = "n", ylim = c(-0.5, 0.5),
    col = colvec)
  for (j in seq_along(heatper)) {
    polygon(c(sgrid, rev(sgrid)), c(betaCIi[1,,j], rev(betaCIi[2,,j])),
      col = adjustcolor(colvec[j], .1), border = NA)
  }
  # matlines(sgrid, betai, col = colvec, lwd = 3)
  # matlines(sgrid, betaCIi[1,,], lty = 2,
  #   col = c("black", rev(heat_hcl(length(heatper)))))
  # matlines(sgrid, betaCIi[2,,], lty = 2,
  #   col = c("black", rev(heat_hcl(length(heatper)))))
  axis(1, at = c(datatab$s[1] - .5, datatab$s + .5), labels = F)
  axis(1, at = datatab$s, labels = mon_vec, cex.axis = .7, tick = F)
  abline(h = 0)
  abline(v = which(mon_vec == "Jan") - 1, lty = 2)
  text(tapply(datatab$s, year_vec, mean), par("usr")[3], 
    c("Same year", "Previous year"), pos = 3, cex = .8)
}
# Add legend
par(mar = c(5, 0, 4, 0))
plot.new()
legend("center", legend = c("MMT", heatper * 100), 
  lty = seq_len(length(heatper) + 1), 
  lwd = 2, col = colvec, title = "Percentile",
  xpd = NA, bty = "n", cex = 1.2)

# Save
dev.print(pdf, file = "Figures/Figure3.pdf")


#------------------------
# First stage DLNM result
#------------------------

yearsplot <- c(1981, 2018)
yearscol <- c(1,4)

#----- Overall relationship with selected heat cutpoints
# Common plot scales
ylim <- c(0.5, 3)
xlim <- range(sapply(dlist, "[[", "Tmean"))

x11()
# Plot in a grid layout
layout(cbind(1:2, 3), widths = c(4, 1))
for (i in seq_along(region)){
  plot(NA, xlim = xlim, ylim = ylim, xlab = "Temperature", 
    ylab = "RR", main = reg_lab[i])
  for (j in seq_along(yearsplot)){
    yj <- which(years == yearsplot[j]) 
    
    # select crosspred
    cp <- cplist[[i]][[yj]]
    
    # Plot overall relationship
    lines(cp, ptype = "overall", col = yearscol[j], lwd = 2, ci = "area",
      ci.arg = list(density = 15, col = adjustcolor(yearscol[j], .5), 
        angle = 45 * j))
    abline(v = mmt[yj,i], lty = 2, col = yearscol[j], lwd = 2)
  }
  abline(h = 1)
  abline(v = quantile(daily_dat$Tmean, heatper), lty = 2, lwd = 2,
    col = colvec[-1])
}
# Add legend
par(mar = c(5, 0, 4, 0))
plot.new()
h <- legend("topleft", legend = yearsplot, col = yearscol, lwd = 2, 
  title = "ERF", xpd = NA, bty = "n")
legend(h$rect$left, h$rect$top - h$rect$h, 
  legend = c(sprintf("MMT %i", yearsplot), heatper * 100), 
  lty = 2, lwd = 2, col = c(yearscol, colvec[-1]), title = "Percentile",
  xpd = NA, bty = "n")

# Save
dev.print(png, file = "Figures/SupFigure2.png", res = 100, unit = "in")

#------------------------
# Residuals analysis
#------------------------

# Extract residuals
resids <- do.call(cbind, lapply(reslist, sapply, function(x){
  x$final_mod$resid()
}))


#----- Residuals versus year
x11(width = 10)

layout(matrix(1:2, nrow = 2), heights = c(4,1))
par(mar = c(5, 4, 2, 4) + .1)
matplot(datatab$Year, resids, col = colvec, type = "b",
  pch = rep(15:16, each = length(heatper) + 1), 
  xlab = "Year", ylab = "Residual")
abline(h = 0)
# Add legend
par(mar = c(0, 4, 0, 4) + .1)
plot.new()
lg <- legend("left", legend = c("MMT", heatper * 100), 
  fill = c(colvec), xpd = NA, bty = "n", ncol = length(heatper) + 1)
legend(lg$rect$left + lg$rect$w, lg$rect$top, reg_lab, 
  col = rep(1, length(region)), pch = 15:17, xpd = NA, bty = "n",
  ncol = length(region) + 1)

# Save
dev.print(png, file = "Figures/SupFigure2.png", res = 100, unit = "in")

#----- Distribution of residuals

# Compute kernel density estimate
resid_kde <- apply(resids, 2, density)
kde_x <- sapply(resid_kde, "[[", "x")
kde_y <- sapply(resid_kde, "[[", "y")

# Plot KDE
x11(height = 10)

layout(cbind(1:2, 3), widths = c(4, 1))
for (i in seq_along(region)){
  ind <- (length(heatper) + 1) * (i - 1) + seq_len(length(heatper) + 1)
  matplot(kde_x[,ind], kde_y[,ind], type = "l", col = colvec, lwd = 3,
    xlab = "Residuals", ylab = "Density", main = reg_lab[i])
  abline(v = 0, lty = 2)
  abline(h = 0)
}
# Add legend
par(mar = c(5, 0, 4, 0))
plot.new()
legend("center", legend = c("MMT", heatper * 100), lty = seq_along(heatper), 
  lwd = 3, col = colvec, title = "Percentile",
  xpd = NA, bty = "n")

# Save
dev.print(png, file = "Figures/SupFigure3.png", res = 100, unit = "in")

#------------------------
# Plot climate indices
#------------------------

# Parameters
nma <- 3

# plot

x11(width = 15, height = 10)
par(mfrow = n2mfrow(ncl))

for(i in seq_len(ncl)){
  # Extract climate values
  cli <- subset(climate_dat, Y %in% years, climate_labs[i])[[1]]
  
  # Smooth
  smcli <- forecast::ma(cli, order = nma)
  
  # create dates
  clidates <- as.Date(sprintf("%i-%i-01", rep(years, each = 12), rep(1:12, ny)))
  
  # plot it
  plot(clidates, cli, type = "h", col = ifelse(cli > 0, 4, 2), 
    xlab = "Month", ylab = climate_labs[i], main = climate_labs[i])
  lines(clidates, smcli, lwd = 2)
  abline(h = 0)
}

dev.print(png, file = "Figures/SupFigure1.png", res = 100, unit = "in")