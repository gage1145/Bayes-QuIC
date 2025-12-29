library(priogenapi)
library(tidyverse)
library(MASS)



key <- Sys.getenv("PRILIM_API_KEY")
data <- call_api(key, "results")

data1 <- data %>%
  # filter(!(sample %in% c("N", "P"))) %>%
  select(sample, rxn_type, MPR, MS, RAF, TtT) %>%
  mutate(
    MPR_crossed = MPR > 3,
    MS_crossed = MS > 0.3
  ) %>%
  mutate_at(c("MPR", "MS", "RAF", "TtT"), as.numeric)

data_sum <- data1 %>%
  group_by(sample, rxn_type) %>%
  summarize(
    mean_MPR = mean(MPR),
    median_MPR = median(MPR),
    # mean_MS  = mean(MS),
    MPR_positive = median_MPR > 5,
    # MS_positive = mean(MS_crossed) > 0.25,
    # positive = as.numeric(MPR_positive & MS_positive)
  ) %>%
  mutate_at(c("MPR_positive"), as.numeric)

data1 <- left_join(data1, data_sum)


# Nano-QuIC probabilities
nano_sum <- filter(data_sum, rxn_type == "Nano-QuIC")
nano_num_pos <- length(nano_sum$MPR_positive[which(nano_sum$MPR_positive == 1)])
prior_prob <- nano_num_pos / nrow(nano_sum)

nano_df <- data1 %>%
  filter(rxn_type == "Nano-QuIC")

mod <- glm(MPR_positive ~ MPR, data=nano_df, family="binomial")
mod_sum <- summary(mod)
# plot(mod)

predictions <- predict.glm(mod, nano_df, type="response")

fit_ln <- fitdistrplus::fitdist(nano_df$MPR, "lnorm")

mu <- fit_ln$estimate["meanlog"]
sigma <- fit_ln$estimate["sdlog"]



nano_df <- nano_df %>%
  mutate(evidence = 1 - pnorm(log(MPR), mu, sigma))









nano_df %>%
  ggplot(aes(MPR, MPR_positive)) +
  geom_point() +
  geom_smooth(method="glm", method.args = list(family="binomial"))
