library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(stringr)
library(caret)
library(runjags)
library(coda)

# runJagsOut1: Noninformative - Nonimputed -------------------------------------------------------------
runJagsOut1 <- readRDS("~/Desktop/MATH2269 Applied Bayesian Statistics/Assignments/A3/Non-imputed/runJagsOut/runJagsOut1.rds")

runJagsOut1$timetaken

codaSamples1 <- as.mcmc.list( runJagsOut1 )

codadf1 <- codaSamples1 %>% 
  as.matrix(chains = TRUE) %>%
  as.data.frame() %>%
  select(-c(CHAIN, contains("train")))

summary1 <- codadf1 %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  group_by(parameter) %>%
  summarise(mean = mean(value),
            mode = mode(value),
            median = median(value),
            var = var(value),
            lower = HDIofMCMC(value, credMass = 0.95)[1],
            upper = HDIofMCMC(value, credMass = 0.95)[2],
            ESS = effectiveSize(value))






# runJagsOut2: Informative - Nonimputed -------------------------------------------------------------
runJagsOut2 <- readRDS("~/Desktop/MATH2269 Applied Bayesian Statistics/Assignments/A3/Non-imputed/runJagsOut/runJagsOut2.rds")

codaSamples2 <- as.mcmc.list( runJagsOut2 )

codadf2 <- codaSamples2 %>% 
  as.matrix(chains = TRUE) %>%
  as.data.frame() %>%
  select(-c(CHAIN, contains("train")))

summary2 <- codadf2 %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  group_by(parameter) %>%
  summarise(mean = mean(value),
            mode = mode(value),
            median = median(value),
            var = var(value),
            lower = HDIofMCMC(value, credMass = 0.95)[1],
            upper = HDIofMCMC(value, credMass = 0.95)[2],
            ESS = effectiveSize(value))

runJagsOut2$timetaken



# runJagsOut3: Noninformative - Imputed -------------------------------------------------------------
runJagsOut3 <- readRDS("~/Desktop/MATH2269 Applied Bayesian Statistics/Assignments/A3/Imputed/runJagsOut/runJagsOut1.rds")

runJagsOut3$timetaken
codaSamples3 <- as.mcmc.list( runJagsOut3 )

codadf3 <- codaSamples3 %>% 
  as.matrix(chains = TRUE) %>%
  as.data.frame() %>%
  select(-CHAIN)

summary3 <- codadf3 %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  group_by(parameter) %>%
  summarise(mean = mean(value),
            mode = mode(value),
            median = median(value),
            var = var(value),
            lower = HDIofMCMC(value, credMass = 0.95)[1],
            upper = HDIofMCMC(value, credMass = 0.95)[2],
            ESS = effectiveSize(value))



# runJagsOut4:Informative - Imputed -------------------------------------------------------------
runJagsOut4 <- readRDS("~/Desktop/MATH2269 Applied Bayesian Statistics/Assignments/A3/Imputed/runJagsOut/runJagsOut2.rds")

runJagsOut4$timetaken
codaSamples4 <- as.mcmc.list( runJagsOut4 )

codadf4 <- codaSamples4 %>% 
  as.matrix(chains = TRUE) %>%
  as.data.frame() %>%
  select(-CHAIN)

summary4 <- codadf4 %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  group_by(parameter) %>%
  summarise(mean = mean(value),
            mode = mode(value),
            median = median(value),
            var = var(value),
            lower = HDIofMCMC(value, credMass = 0.95)[1],
            upper = HDIofMCMC(value, credMass = 0.95)[2],
            ESS = effectiveSize(value))


# Helper ------------------------------------------------------------------

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



# Convergence -------------------------------------------------------------
parameters <- c("zbeta0","zbeta[1]", "zbeta[2]", "zbeta[3]", "zbeta[4]", "zbeta[5]", "guess")
mcmcDiagnostics(codaSamples1[, parameters[1]], title = "Intercept")
mcmcDiagnostics(codaSamples1[, parameters[2]], title = "DEROG")
mcmcDiagnostics(codaSamples1[, parameters[3]], title = "DELINQ")
mcmcDiagnostics(codaSamples1[, parameters[4]], title = "CLAGE")
mcmcDiagnostics(codaSamples1[, parameters[5]], title = "NINQ")
mcmcDiagnostics(codaSamples1[, parameters[6]], title = "DEBTINC")
mcmcDiagnostics(codaSamples1[, parameters[7]], title = "Guess")

codadf1 %>% 
  select(all_of(parameters), -c("zbeta0", "guess")) %>% 
  rename("DEROG" = "zbeta[1]",
         "DELINQ" = "zbeta[2]",
         "CLAGE" = "zbeta[3]",
         "NINQ" = "zbeta[4]",
         "DEBTINC" = "zbeta[5]") %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_interval(.width = c(.50, .80, .95, .99)) +
  stat_dotsinterval(quantiles = 90) +
  scale_color_brewer() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(title = "Model 1's Coefficient Significance",
       x = "",
       y = "")


parameters <- c("zbeta0","zbeta[1]", "zbeta[2]", "zbeta[3]", "zbeta[4]", "zbeta[5]", "guess")
mcmcDiagnostics(codaSamples2[, parameters[1]], title = "Intercept")
mcmcDiagnostics(codaSamples2[, parameters[2]], title = "DEROG")
mcmcDiagnostics(codaSamples2[, parameters[3]], title = "DELINQ")
mcmcDiagnostics(codaSamples2[, parameters[4]], title = "CLAGE")
mcmcDiagnostics(codaSamples2[, parameters[5]], title = "NINQ")
mcmcDiagnostics(codaSamples2[, parameters[6]], title = "DEBTINC")
mcmcDiagnostics(codaSamples2[, parameters[7]], title = "Guess")

codadf2 %>% 
  select(all_of(parameters), -c("zbeta0", "guess")) %>% 
  rename("DEROG" = "zbeta[1]",
         "DELINQ" = "zbeta[2]",
         "CLAGE" = "zbeta[3]",
         "NINQ" = "zbeta[4]",
         "DEBTINC" = "zbeta[5]") %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_interval(.width = c(.50, .80, .95, .99)) +
  stat_dotsinterval(quantiles = 90) +
  scale_color_brewer() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(title = "Model 2's Coefficient Significance",
       x = "",
       y = "")


parameters <- c("zbeta0","beta[1]", "beta[2]", "beta[3]", "zbeta[1]", "zbeta[2]", "guess")
mcmcDiagnostics(codaSamples3[, parameters[1]], title = "Intercept")
mcmcDiagnostics(codaSamples3[, parameters[2]], title = "DEROG")
mcmcDiagnostics(codaSamples3[, parameters[3]], title = "DELINQ")
mcmcDiagnostics(codaSamples3[, parameters[4]], title = "NINQ")
mcmcDiagnostics(codaSamples3[, parameters[5]], title = "CLAGE")
mcmcDiagnostics(codaSamples3[, parameters[6]], title = "DEBTINC")
mcmcDiagnostics(codaSamples3[, parameters[7]], title = "Guess")

codadf3 %>% 
  select(all_of(parameters), -c("zbeta0", "guess")) %>% 
  rename("DEROG" = "beta[1]",
         "DELINQ" = "beta[2]",
         "CLAGE" = "zbeta[1]",
         "NINQ" = "beta[3]",
         "DEBTINC" = "zbeta[2]") %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_interval(.width = c(.50, .80, .95, .99)) +
  stat_dotsinterval(quantiles = 90) +
  scale_color_brewer() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(title = "Model 3's Coefficient Significance",
       x = "",
       y = "")


parameters <- c("zbeta0","beta[1]", "beta[2]", "beta[3]", "zbeta[1]", "zbeta[2]", "guess")
mcmcDiagnostics(codaSamples4[, parameters[1]], title = "Intercept")
mcmcDiagnostics(codaSamples4[, parameters[2]], title = "DEROG")
mcmcDiagnostics(codaSamples4[, parameters[3]], title = "DELINQ")
mcmcDiagnostics(codaSamples4[, parameters[4]], title = "NINQ")
mcmcDiagnostics(codaSamples4[, parameters[5]], title = "CLAGE")
mcmcDiagnostics(codaSamples4[, parameters[6]], title = "DEBTINC")
mcmcDiagnostics(codaSamples4[, parameters[7]], title = "Guess")

codadf4 %>% 
  select(all_of(parameters), -c("zbeta0", "guess")) %>% 
  rename("DEROG" = "beta[1]",
         "DELINQ" = "beta[2]",
         "CLAGE" = "zbeta[1]",
         "NINQ" = "beta[3]",
         "DEBTINC" = "zbeta[2]") %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_interval(.width = c(.50, .80, .95, .99)) +
  stat_dotsinterval(quantiles = 90) +
  scale_color_brewer() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(title = "Model 4's Coefficient Significance",
       x = "",
       y = "")


