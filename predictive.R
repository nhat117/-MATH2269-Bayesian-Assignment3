library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(stringr)
library(caret)
library(tidymodels)


# data --------------------------------------------------------------------
set.seed(123)

dataset1 <- read.csv("data/nonImputedHMEQ.csv")
dataset2 <- read.csv("data/ImputedHMEQ.csv")

data_test12 <- initial_split(dataset1, prop = ratio, strata = "BAD") %>% testing()
data_test34 <- initial_split(dataset2, prop = ratio, strata = "BAD") %>% testing()

testSummary1 <- summary1 %>% 
  filter(str_detect(parameter, "pred"))

testSummary2 <- summary2 %>%
  filter(str_detect(parameter, "pred"))

testSummary3 <- summary3 %>%
  filter(str_detect(parameter, "pred"))

testSummary4 <- summary4 %>%
  filter(str_detect(parameter, "pred"))

threshold_value <- 0.3
result1 <- data_test12 %>% 
  mutate(parameter = paste0("pred[", rownames(data_test12), "]")) %>% 
  select(parameter, BAD) %>% 
  left_join(testSummary1, by = "parameter") %>%
  mutate(BAD = as.factor(BAD),
         Predicted = as.factor(ifelse(mode > threshold_value, 1, 0)),
         Precision = case_when(
           BAD == 1 & Predicted == 1 ~ "True Positive",
           BAD == 0 & Predicted == 1 ~ "False Positive",
           BAD == 1 & Predicted == 0 ~ "False Negative",
           BAD == 0 & Predicted == 0 ~ "True Negative"
         )) %>%
  select(parameter, mode, BAD, Predicted, Precision, lower, upper)
confusionMatrix(result1$Predicted, result1$BAD)
true_values <- as.numeric(as.character(result1$BAD))
predicted_probabilities <- as.numeric(as.character(result1$Predicted))
ks.test(
  predicted_probabilities[true_values == 1],  # Probabilities for the positive class
  predicted_probabilities[true_values == 0]   # Probabilities for the negative class
)


threshold_value <- 0.3
result2 <- data_test12 %>%
  mutate(parameter = paste0("pred[", rownames(data_test12), "]")) %>% 
  select(parameter, BAD) %>% 
  left_join(testSummary2, by = "parameter") %>%
  mutate(BAD = as.factor(BAD),
         Predicted = as.factor(ifelse(mode > threshold_value, 1, 0)),
         Precision = case_when(
           BAD == 1 & Predicted == 1 ~ "True Positive",
           BAD == 0 & Predicted == 1 ~ "False Positive",
           BAD == 1 & Predicted == 0 ~ "False Negative",
           BAD == 0 & Predicted == 0 ~ "True Negative"
         )) %>%
  select(parameter, mode, BAD, Predicted, Precision, lower, upper)
confusionMatrix(result2$Predicted, result2$BAD)

true_values <- as.numeric(as.character(result2$BAD))
predicted_probabilities <- as.numeric(as.character(result2$Predicted))

ks.test(
  predicted_probabilities[true_values == 1],  # Probabilities for the positive class
  predicted_probabilities[true_values == 0]   # Probabilities for the negative class
)

threshold_value <- 0.5
result3 <- data_test34 %>%
  mutate(parameter = paste0("pred[", rownames(data_test34), "]")) %>% 
  select(parameter, BAD) %>% 
  left_join(testSummary3, by = "parameter") %>%
  mutate(BAD = as.factor(BAD),
         Predicted = as.factor(ifelse(mode > threshold_value, 1, 0)),
         Precision = case_when(
           BAD == 1 & Predicted == 1 ~ "True Positive",
           BAD == 0 & Predicted == 1 ~ "False Positive",
           BAD == 1 & Predicted == 0 ~ "False Negative",
           BAD == 0 & Predicted == 0 ~ "True Negative"
         )) %>%
  select(parameter, mode, BAD, Predicted, Precision, lower, upper)
confusionMatrix(result3$Predicted, result3$BAD)
true_values <- as.numeric(as.character(result3$BAD))

predicted_probabilities <- as.numeric(as.character(result3$Predicted))
ks.test(
  predicted_probabilities[true_values == 1],  # Probabilities for the positive class
  predicted_probabilities[true_values == 0]   # Probabilities for the negative class
)

threshold_value <- 0.3
result4 <- data_test34 %>%
  mutate(parameter = paste0("pred[", rownames(data_test34), "]")) %>% 
  select(parameter, BAD) %>% 
  left_join(testSummary4, by = "parameter") %>%
  mutate(BAD = as.factor(BAD),
         Predicted = as.factor(ifelse(mode > threshold_value, 1, 0)),
         Precision = case_when(
           BAD == 1 & Predicted == 1 ~ "True Positive",
           BAD == 0 & Predicted == 1 ~ "False Positive",
           BAD == 1 & Predicted == 0 ~ "False Negative",
           BAD == 0 & Predicted == 0 ~ "True Negative"
         )) %>%
  select(parameter, mode, BAD, Predicted, Precision, lower, upper)
confusionMatrix(result4$Predicted, result4$BAD)

true_values <- as.numeric(as.character(result4$BAD))
predicted_probabilities <- as.numeric(as.character(result4$Predicted))
ks.test(
  predicted_probabilities[true_values == 1],  # Probabilities for the positive class
  predicted_probabilities[true_values == 0]   # Probabilities for the negative class
)


# Model Performance and Threshold Sensitivity Analysis ------------------------------------------
library(readxl)
compare_data <- read_excel("model comparison.xlsx")
str(compare_data)

plot_data <- compare_data %>% 
  select(-c("Mcnemar's Test", "KS D-Statistic")) %>% 
  pivot_longer(cols = -c("Model", "Threshold"), names_to = "Metric", values_to = "Value") %>% 
  mutate(Value = Value * 100)

plot_data %>% 
  ggplot(aes(x = Value, y = Model, colour = factor(Threshold))) +
  facet_grid(~Metric, scales = "free") +
  geom_point() +
  labs(title = "Model Performance and Threshold Sensitivity Analysis", x = "Value", y = "", color = "Threshold")

# Model 4 at threshold = 0.3
# Model 3 at threshold  = 0.5
# Model 2 at threshold = 0.3
# Model 1 at threshold = 0.3


# Dicrimination Power and Stability -----------------------------------------------------

result1 <- result1 %>% mutate(model = "Model 1 - Threshold = 0.3 - KS D-Statistics = 0.31", )

result2 <- result2 %>% mutate(model = "Model 2 - Threshold = 0.3 - KS Dstatistics = 0.26")

result3 <- result3 %>% mutate(model = "Model 3 - Threshold = 0.5 - KS Dstatistics = 0.29")

result4 <- result4 %>% mutate(model = "Model 4 - Threshold = 0.3 - KS Dstatistics = 0.31")


result <- bind_rows(result1, result2, result3, result4)
str(result)

result %>% 
  ggplot(aes(x = parameter, y = mode, color = Precision)) +
  geom_point(alpha = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  facet_wrap(~ model, scales = "free") +
  labs(title = "Model Performance at Chosen Thresholds", y = "Predicted Probabilities", color = "Mode with 95% HDI") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")

result %>% 
  ggplot(aes(x = parameter, y = mode, color = Precision)) +
  geom_point(alpha = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  facet_wrap(~ model, scales = "free") +
  labs(title = "Model Performance at Chosen Thresholds", y = "Predicted Probabilities", color = "Mode with 95% HDI") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")


# Model Intepretation -------------------------------------------------
coefficientSamples <- codadf3 %>% 
  select(starts_with("beta"), "guess") %>% 
           rename("Baseline" = "guess" ,
                  "Intercept" = "beta0",
                  "DEROG" = "beta[1]",
                  "DELINQ" = "beta[2]",
                  "CLAGE" = "betaCLAGE",
                  "NINQ" = "beta[3]",
                  "DEBTINC" = "betaDEBTINC")

coefficients <- c("beta0", "beta[1]", "beta[2]", "beta[3]", "betaCLAGE", "betaDEBTINC", "guess")
coefficientsSummary <- summary3 %>% 
  filter(parameter %in% coefficients) %>% 
  mutate(parameter = case_when(
    parameter == "beta0" ~ "Intercept",
    parameter == "beta[1]" ~ "DEROG",
    parameter == "beta[2]" ~ "DELINQ",
    parameter == "beta[3]" ~ "NINQ",
    parameter == "betaCLAGE" ~ "CLAGE",
    parameter == "betaDEBTINC" ~ "DEBTINC",
    parameter == "guess" ~ "Baseline"
  ))


# Model Selection ---------------------------------------------------------
# Chosen model is model 3
coefficientSamples %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  left_join(coefficientsSummary, by = "parameter") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(fill = "skyblue", color = "black") +
  geom_text(aes(label = paste0("mode = ", round(mode, 3))),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3) +
  facet_wrap(~ parameter, scales = "free") +
  geom_segment(aes(x = lower, xend = upper, y = -0.05, yend = -0.05), colour = "black", size = 1.5) +
  geom_text(aes(x = (lower + upper) / 2, y = 40, label = "95% HDI"),
            size = 3, color = "black", hjust = 0.5) +
  geom_text(aes(x = lower, y = -30, label = round(lower, 3)),
            size = 3, color = "black", hjust = 1) +
  geom_text(aes(x = upper, y = -30, label = round(upper, 3)),
            size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  labs(title = "Parameter Interpretation",
       x = "Value",
       y = "")



