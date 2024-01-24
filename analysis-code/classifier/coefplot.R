plot_coeffs_S <- function(mlr_model, plot, title = "") {
    #coeffs <- sort(coefficients(mlr_model), decreasing = TRUE)  ### changed
    coeffs <- sort(coef(summary(mlr_model))[, "z value"], decreasing = TRUE)  ### changed to use z-values
    coeffs <- coeffs[names(coeffs) != "(Intercept)"]
    #print("coeffs", coeffs)
    #coeffs <- coeffs[2:length(coeffs)]
    pdf(plot)
        mp <- barplot(coeffs, horiz = F, col="#3F97D0", xaxt='n',main = title, ylab = "Weight", xlab = "Feature")
        lablist <- names(coeffs)
        #print(lablist)
        lablist <- gsub("_HUMAN.H11MO.0..", " motif", lablist)
        #print(lablist)
        text(mp, par("usr")[3], labels = lablist, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.6)
    dev.off()
}

# published on http://www.r-statistics.com/2010/07/visualization-of-regression-coefficients-in-r
# originally written by "<a href="http://statmath.wu.ac.at/~zeileis/">Achim Zeileis</a>"
# GPL-2

coefplot <- function(object, df = NULL, level = 0.95, parm = NULL,
					 labels = TRUE, xlab = "Coefficient confidence intervals", ylab = "",
					 xlim = NULL, ylim = NULL, title = "",
					 las = 1, lwd = 1, lty = c(1, 2), pch = 19, col = 1,
					 length = 0, angle = 30, code = 3, plotname = "test", ...)
{
	 cf <- coef(object)
	 se <- sqrt(diag(vcov(object)))
	 if(is.null(parm)) parm <- seq_along(cf)
	 if(is.numeric(parm) | is.logical(parm)) parm <- names(cf)[parm]
	 if(is.character(parm)) parm <- which(names(cf) %in% parm)
	 cf <- cf[parm]
	 se <- se[parm]
	 k <- length(cf)

	 if(is.null(df)) {
	   df <- if(identical(class(object), "lm")) df.residual(object) else 0
	 }

	 critval <- if(df > 0 & is.finite(df)) {
	   qt((1 - level)/2, df = df)
	 } else {
	   qnorm((1 - level)/2)
	 }
	 ci1 <- cf + critval * se
	 ci2 <- cf - critval * se

	 lwd <- rep(lwd, length.out = 2)
	 lty <- rep(lty, length.out = 2)
	 pch <- rep(pch, length.out = k)
	 col <- rep(col, length.out = k)

	 if(is.null(xlim)) xlim <- range(c(0, min(ci1), max(ci2)))
	 if(is.null(ylim)) ylim <- c(1 - 0.05 * k, 1.05 * k)

	 if(isTRUE(labels)) labels <- names(cf)
	 if(identical(labels, FALSE)) labels <- ""
	 labels <- rep(labels, length.out = k)

     pdf(plotname)
	 plot(0, 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
	   axes = FALSE, type = "n", las = las, main = title, ...)
	 arrows(ci1, 1:k, ci2, 1:k, lty = lty[1], lwd = lwd[1], col = col,
	   length = length, angle = angle, code = code)
	 points(cf, 1:k, pch = pch, col = col)
	 abline(v = 0, lty = lty[2], lwd = lwd[2])
	 axis(1)
	 axis(2, at = 1:k, labels = labels, las = las)
	 box()
     dev.off()
}
