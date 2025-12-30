library(priogenapi)
library(tidyverse)
library(poolr)



key <- Sys.getenv("PRILIM_API_KEY")
data <- call_api(key, "results")

data1 <- data %>%
  dplyr::select(sample, rxn_type, MPR, MS, RAF, TtT) %>%
  mutate_at(c("MPR", "MS", "RAF", "TtT"), as.numeric)

threshold <- 5

data_sum <- data1 %>%
  group_by(sample, rxn_type) %>%
  summarize(
    mean_MPR = mean(MPR),
    median_MPR = median(MPR),
    positive = median_MPR > threshold,
  )

data2 <- left_join(data1, data_sum)


# Nano-QuIC probabilities

observation <- 5

nano_df <- data2 %>%
  filter(rxn_type == "Nano-QuIC")

# P(B|A)
num_pos <- nrow(nano_df[which(nano_df$positive),])
true_pos <- nrow(nano_df[which(nano_df$positive & nano_df$MPR > observation),])
# likelihood <- true_pos / num_pos
# 
# # P(A)
prior_prob <- length(nano_df$positive[which(nano_df$positive)]) / nrow(nano_df)
# 
# # P(B)
# evidence <- nrow(nano_df[which(nano_df$MPR > observation),]) / nrow(nano_df)
# 
# # P(A|B)
# posterior <- likelihood * prior_prob / evidence


likelihoods <- sapply(
  nano_df$MPR, 
  function(x) nrow(nano_df[which(nano_df$positive & nano_df$MPR > x),]) / num_pos
)

evidences <- sapply(
  nano_df$MPR,
  function(x) nrow(nano_df[which(nano_df$MPR > x),]) / nrow(nano_df)
)

nano_df_result <- nano_df %>%
  mutate(
    likelihood = likelihoods,
    prior_prob = prior_prob,
    evidence   = evidences,
    posterior  = (likelihood * prior_prob) / evidence
  )

nano_result_sum <- nano_df_result %>%
  group_by(sample) %>%
  summarize(
    mean_MPR   = median(MPR),
    likelihood = fisher(likelihood)$p,
    prior_prob = fisher(prior_prob)$p,
    evidence   = fisher(evidence)$p,
    posterior  = (likelihood * prior_prob) / evidence
  )




# mod <- glm(positive ~ MPR, data=nano_df, family="binomial")
# mod_sum <- summary(mod)
# 
# evidence <- predict.glm(mod, nano_df, type="response")









nano_df %>%
  mutate(positive = as.integer(positive)) %>%
  ggplot(aes(MPR, positive)) +
  geom_point() +
  geom_smooth(method="glm", method.args = list(family="binomial"))
