library(dplyr)
library(tidyr)
library(tidymodels)
library(runjags)

set.seed(123)
ratio <- 0.8


# data --------------------------------------------------------------------
predictors <- c("DEROG", "DELINQ", "CLAGE", "NINQ", "DEBTINC")
target <- "BAD"

dataset1 <- read.csv("Non-imputed/nonImputedHMEQ.csv") %>% 
  select(all_of(c(target, predictors))) %>% 
  drop_na()

dataset2 <- read.csv("Imputed/imputedHMEQ.csv") %>% 
  select(all_of(c(target, predictors))) %>%
  drop_na()


# Onto dataset 1 ----------------------------------------------------


# train-test split
data_split <- initial_split(dataset1, prop = ratio, strata = response)
train_data <- training(data_split)
test_data <- testing(data_split)


Xtrain <- train_data %>% select(-c(response)) %>% as.matrix()
XTrainSD <- apply(Xtrain, 2, sd)
XTrainMean <- apply(Xtrain, 2, mean)
scaledXTrain <- scale(Xtrain, center = XTrainMean, scale = XTrainSD)

scaledXTest <- test_data %>% select(-c(response)) %>% as.matrix() %>% scale(center = XTrainMean, scale = XTrainSD)

yTrain <- train_data %>% select(response) %>% as.matrix() %>% as.numeric()

# Frequentist Approach for Inits
frequentist_model <- glm(BAD ~ ., data = train_data, family = binomial(link = "logit"))

initsList <- list(
  zbeta0=coef(frequentist_model)["(Intercept)"],
  zbeta=coef(frequentist_model)[-1]
)



# datalist

dataList <- list(
  y = yTrain, 
  X = scaledXTrain, 
  xPred = scaledXTest, 
  Ntotal = nrow(scaledXTrain), 
  Nx = ncol(scaledXTrain), 
  Npred = nrow(scaledXTest), 
  xsd = XTrainSD, 
  xm = XTrainMean
)


# model

monitorList = c( "zbeta0" ,  "zbeta" , "guess", "beta0", "beta", "train", "pred" )

modelString = "
model {
  
  # Likelihood
  for ( i in 1:Ntotal ) {
    y[i] ~ dbern( theta[i] )
    theta[i] <- ( guess*(1/2) + (1.0-guess)*ilogit(zbeta0+sum(zbeta[1:Nx]*X[i,1:Nx])) )
  }
  
  # Priors
  zbeta0 ~ dnorm( 0 , 1.0E-6 )              #Intercept
  zbeta[1] ~ dnorm(  0 , 1.0E-6)            #DEROG
  zbeta[2] ~ dnorm(0.375, 1 / (0.10^2))     #CLAGE
  zbeta[3] ~ dnorm(  0 , 1.0E-6 )           #DELINQ
  zbeta[4] ~ dnorm(  0 , 1.0E-6 )           #NINQ
  zbeta[5] ~ dnorm(0.375, 1 / (0.10^2))     #DEBTINC
 
  guess ~ dbeta(1,19)
  
  # Training Accuracy
  for (i in 1:Ntotal){
  train[i] <- guess*(1/2) + (1.0-guess)*ilogit(zbeta0 + sum(zbeta[1:Nx] * X[i,1:Nx]))
  }
  
  # Predictions
  for (k in 1:Npred){
  pred[k] <- guess*(1/2) + (1.0-guess) *ilogit(zbeta0 + sum(zbeta[1:Nx] * xPred[k,1:Nx]))
  }
  
  # Scale back for Inference
  beta[1:Nx] <- zbeta[1:Nx] / xsd[1:Nx]
  beta0 <- zbeta0 - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )
}
" 

writeLines( modelString , con="TEMPmodel.txt" )


# Tuning

adaptSteps <- 100  

burnInSteps <- 5000

numSavedSteps <- 1000 

thinSteps <- 100

nChains <- 4

nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

runjagsMethod <- "parallel"


# Run

runJagsOut <- run.jags( method=runjagsMethod ,
                        model="TEMPmodel.txt" ,
                        monitor=monitorList ,
                        data=dataList ,
                        inits=initsList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps ,
                        summarise=FALSE ,
                        plots=FALSE )




# On dataset 2 ------------------------------------------------------------

# For dataset2
data_split <- initial_split(dataset2, prop = ratio, strata = response)
train_data <- training(data_split)
test_data <- testing(data_split)


Xtrain <- train_data %>% select(-c(response)) %>% as.matrix()

sdCLAGE <- Xtrain[,"CLAGE"] %>% sd()
sdDEBTINC <- Xtrain[,"DEBTINC"] %>% sd()
meanCLAGE <- Xtrain[,"CLAGE"] %>% mean()
meanDEBTINC <- Xtrain[,"DEBTINC"] %>% mean()
scaledXTrain <- Xtrain
scaledXTrain[,c("CLAGE", "DEBTINC")] <- Xtrain[,c("CLAGE", "DEBTINC")] %>% scale(center = c(meanCLAGE, meanDEBTINC), scale = c(sdCLAGE, sdDEBTINC))

scaledXTest <- test_data %>% select(-c(response)) %>% as.matrix()
scaledXTest[,c("CLAGE", "DEBTINC")] <- scaledXTest[,c("CLAGE", "DEBTINC")] %>% scale(center = c(meanCLAGE, meanDEBTINC), scale = c(sdCLAGE, sdDEBTINC))

yTrain <- train_data %>% select(response) %>% as.matrix() %>% as.numeric()
yTest <- test_data %>% select(response) %>% as.matrix() %>% as.numeric()


# Frequentist Approach for Inits
frequentist_model <- glm(BAD ~ ., data = train_data, family = binomial(link = "logit"))

coef(frequentist_model)

initsList <- list(
  zbeta0=coef(frequentist_model)["(Intercept)"],
  beta = coef(frequentist_model)[c("DEROG", "DELINQ", "NINQ")],
  zbeta=coef(frequentist_model)[c("CLAGE", "DEBTINC")],
  guess=0.3
)

# Datalist

dataList <- list(
  y = yTrain, 
  X = scaledXTrain, 
  Xpred = scaledXTest,
  Ntotal = nrow(scaledXTrain), 
  Npred = nrow(scaledXTest), 
  sd = c(sdCLAGE, sdDEBTINC),
  mean = c(meanCLAGE, meanDEBTINC)
)


# Model

monitorList = c( "zbeta0" ,  "zbeta", "beta",
                 "guess", "beta0", "pred", "betaDEBTINC", "betaCLAGE" )

modelString = "
model {
  
  # Likelihood
  for ( i in 1:Ntotal ) {
    y[i] ~ dbern( theta[i] )
    theta[i] <- ( guess*(1/2) + (1.0-guess)*ilogit(zbeta0+sum(beta[1:3]*X[i,1:3])+sum(zbeta[1:2]*X[i,4:5])) )
  }
  
  # Priors
  zbeta0 ~ dnorm( 0 , 1.0E-6 )              #Intercept
  beta[1] ~ dnorm(  0 , 1.0E-6 )            #DEROG
  beta[2] ~ dnorm( 0 , 1.0E-6 )             #DELINQ
  beta[3] ~ dnorm(  0 , 1.0E-6 )            #NINQ
  zbeta[1] ~ dnorm(  0.004201792 , 0.01 )           #CLAGE
  zbeta[2] ~ dnorm( 20.3454 ,  0.0002777778 )            #DEBTINC
 
  guess ~ dbeta(1,19)
  
  # Predictions
  for (k in 1:Npred){
  pred[k] <- ilogit(zbeta0+sum(beta[1:3]*Xpred[k,1:3])+sum(zbeta[1:2]*Xpred[k,4:5]))
  }
  
  # Scale back for Inference
  betaCLAGE <- zbeta[1] * sd[1]
  betaDEBTINC <- zbeta[2] * sd[2]
  beta0 <- zbeta0 - sum( zbeta[1:2] * mean / sd)
}
" 


# Tuning

adaptSteps <- 100  

burnInSteps <- 5000

numSavedSteps <- 500

thinSteps <- 150

nChains <- 4

nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

runjagsMethod <- "parallel"


# Run

runJagsOut <- run.jags( method=runjagsMethod ,
                        model="TEMPmodel.txt" ,
                        monitor=monitorList ,
                        data=dataList ,
                        inits=initsList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps ,
                        summarise=FALSE ,
                        plots=FALSE )





