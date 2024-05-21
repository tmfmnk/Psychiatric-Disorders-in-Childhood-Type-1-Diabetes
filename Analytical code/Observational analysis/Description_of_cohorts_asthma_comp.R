#Libraries

library(tidyverse)

#Import data

load(file = "path/Data/Data_models/data_models_exact_matching_w_mortality_asthma_comp.RData")

#Descriptive statistics

data_models[[1]] %>%
 group_by(group) %>%
  summarise(total_n = formatC(round(n(), 2), format = "f", big.mark = ",", digits = 0),
            across(POHL,
                   ~ paste(formatC(sum(. == 1), format = "f", big.mark = ",", digits = 0),
                           paste0("(",
                                  formatC(round(sum(. == 1)/n() * 100, 2), format = "f", digits = 2),
                                  ")")),
                   .names = "males"),
           across(VEK,
                  ~ paste(formatC(round(mean(.), 2), format = "f", digits = 2),
                          paste0("(",
                                 formatC(round(sd(.), 2), format = "f", digits = 2),
                                 ")")),
                  .names = "age"),
           across(year_discharge,
                  ~ paste(median(.),
                          paste0("(",
                                 paste0(quantile(., 0.25, na.rm = TRUE),
                                        "-",
                                        quantile(., 0.75, na.rm = TRUE)),
                                 ")"))),
           across(month_discharge,
                  ~ paste(median(.),
                          paste0("(",
                                 paste0(quantile(., 0.25, na.rm = TRUE),
                                        "-",
                                        quantile(., 0.75, na.rm = TRUE)),
                                 ")")))) %>%
 ungroup() %>%
 pivot_longer(-group,
              names_to = "variable") %>%
 pivot_wider(names_from = "group",
             names_glue = "{group}_diabetes",
             values_from = "value") %>%
 write.csv(file = "path/Results/Description_of_cohorts_asthma_comp.csv",
           row.names = FALSE)

