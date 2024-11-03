library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)
library(tidymodels)

setwd("~/Desktop/MATH2269_A3/submission")


# Data --------------------------------------------------------------------
data("HMEQ")
dataset1 <- read.csv("data/nonImputedHMEQ.csv")
dataset2 <- read.csv("data/ImputedHMEQ.csv")

# Explained Variance and Discrimination Power -----------------------------

HMEQ %>% 
  select(-c(JOB)) %>% 
  drop_na() %>%
  mutate(BAD = as.factor(BAD),
         REASON = as.numeric(as.factor(REASON))) %>% 
  calculate_explained_variance(target_column = target) %>% 
  arrange(desc(Explained_Variance))


# Wrangling ---------------------------------------------------------------

dataset1 %>%
  ggpairs(aes(color = as.factor(BAD)),
          lower = list(continuous = wrap("points", size = 1, alpha = 0.8)),
          upper = list(continuous = wrap("cor", size = 3)),
          diag = list(continuous = wrap("barDiag")))

dataset2 %>%
  ggpairs(aes(color = as.factor(BAD)),
          lower = list(continuous = wrap("points", size = 1, alpha = 0.8)),
          upper = list(continuous = wrap("cor", size = 3)),
          diag = list(continuous = wrap("barDiag")))


