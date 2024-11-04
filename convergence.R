library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(stringr)
library(caret)
library(runjags)
library(coda)

# runJagsOut1: Noninformative - Nonimputed -------------------------------------------------------------
runJagsOut1 <- readRDS("runJagsOut/runJagsOut1.rds")

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
runJagsOut2 <- readRDS("runJagsOut/runJagsOut2.rds")

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
runJagsOut3 <- readRDS("runJagsOut/runJagsOut1.rds")

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
runJagsOut4 <- readRDS("runJagsOut/runJagsOut2.rds")

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


