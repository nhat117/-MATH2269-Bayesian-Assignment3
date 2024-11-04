HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

mode <- function(sampleVec) {
  dres <- density(sampleVec)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}


DbdaAcfPlot <- function(codaObjectPar) {
  
  nChain <- length(codaObjectPar)
  
  xMat <- NULL
  yMat <- NULL
  
  # Loop through chains
  for (cIdx in 1:nChain) {
    acfInfo <- acf(codaObjectPar[[cIdx]], plot = FALSE)
    xMat <- cbind(xMat, acfInfo$lag)
    yMat <- cbind(yMat, acfInfo$acf)
  }
  
  # Plot the ACF for each chain
  matplot(xMat, yMat, type = "o", pch = 20, ylim = c(0, 1),
          main = "Autocorrelation Plot", xlab = "Lag", ylab = "Autocorrelation")
  abline(h = 0, lty = "dashed")
  
  # Calculate and display the effective sample size (ESS)
  EffChnLngth <- effectiveSize(codaObjectPar)
  text(x = max(xMat), y = max(yMat), adj = c(1.0, 1.0), cex = 1.25,
       labels = paste("ESS =", round(EffChnLngth, 1)))
}


DbdaDensPlot <- function(codaObjectPar) {
  
  nChain = length(codaObjectPar)
  
  xMat = NULL
  yMat = NULL
  hdiLims = NULL
  
  # Loop over chains
  for (cIdx in 1:nChain) {
    densInfo = density(codaObjectPar[[cIdx]]) 
    xMat = cbind(xMat, densInfo$x)
    yMat = cbind(yMat, densInfo$y)
    hdiLims = cbind(hdiLims, HDIofMCMC(codaObjectPar[[cIdx]]))
  }
  
  # Create density plot with smooth solid lines
  matplot(xMat, yMat, type = "l", lty = 1, lwd = 2, col = 1:nChain,
          main = "Density Plots", xlab = "Param. Value", ylab = "Density")
  
  # Add baseline
  abline(h = 0)
  
  # Add 95% HDI markers
  points(hdiLims[1,], rep(0, nChain), col = 1:nChain, pch = "|")
  points(hdiLims[2,], rep(0, nChain), col = 1:nChain, pch = "|")
  text(mean(hdiLims), 0, "95% HDI", adj = c(0.5, -0.2))
  
  # Compute Monte Carlo Standard Error (MCSE)
  EffChnLngth = effectiveSize(codaObjectPar)
  MCSE = sd(as.matrix(codaObjectPar)) / sqrt(EffChnLngth)
  
  # Add text for MCSE
  text(max(xMat), max(yMat), adj = c(1.0, 1.0), cex = 1.25,
       paste("MCSE =\n", signif(MCSE, 3)))
}


mcmcDiagnostics = function(codaObjectPar, title="MCMC Diagnostics"){
  
  
  par(mfrow=c(2,2), oma = c(0, 0, 3, 0))
  
  # Traceplot
  coda::traceplot(codaObjectPar , main="Trace Plot" , ylab="Param. Value")
  
  # # Gelman.plot
  # coda::gelman.plot( codaObjectPar , main="Gelman Plot" , auto.layout=FALSE)
  
  # Autocorrelation
  DbdaAcfPlot(codaObjectPar)
  
  # Density
  DbdaDensPlot(codaObjectPar)
  
  mtext(title, outer = TRUE, line = -1, cex = 1.5)
  
  par(mfrow=c(1,1))
  
}


calculate_explained_variance <- function(data, target_column) {
  
  # Extract the target and features
  target <- data[[target_column]]
  features <- data[, !(names(data) %in% target_column)]
  bss <- numeric(ncol(features))
  wss <- numeric(ncol(features))
  
  # Calculate overall means and class means for each feature
  overall_means <- colMeans(features)
  class_means <- aggregate(features, by = list(class = target), FUN = mean)
  
  for (i in 1:ncol(features)) {
    feature_name <- colnames(features)[i]
    
    # Between-class variance (BSS)
    bss[i] <- sum(sapply(1:nrow(class_means), function(k) {
      n_k <- sum(target == class_means$class[k])
      (class_means[k, feature_name] - overall_means[feature_name])^2 * n_k
    }))
    
    # Within-class variance (WSS)
    wss[i] <- sum(sapply(unique(target), function(cls) {
      n_i <- sum(target == cls)
      var(features[target == cls, feature_name]) * (n_i - 1)
    }))
  }
  
  # Calculate discriminating power (BSS / WSS)
  discriminating_power <- bss / wss
  
  # Calculate explained variance (BSS / (BSS + WSS))
  explained_variance <- bss / (bss + wss)
  
  
  result <- data.frame(
    Feature = colnames(features),
    BSS = bss,
    WSS = wss,
    Discriminating_Power = discriminating_power,
    Explained_Variance = explained_variance
  )
  
  return(result)
}
