#Libraries

library(tidyverse)

#Import tables

tab1 <- read.csv(file = "path/Results/Cox_ph_models_ex_matching_up_to_9y_w_mortality.csv")
tab2 <- read.csv(file = "path/Results/Cox_ph_models_ex_matching_w_mortality_asthma_comp.csv")
tab3 <- read.csv(file = "path/Results/Cox_ph_models_ex_matching_w_mortality_incident.csv")
tab4 <- read.csv(file = "path/Cox_ph_models_ex_matching_w_mortality_previous_hosp.csv")

#Outcomes

outcomes <- setNames(c("F1", "F10", "F1[1-9]", "F11", "F12", "F1[3-9]",
                       "F2", "F20", "F2[1-9]",
                       "F3", "F3[0-1]", "F3[2-3]", "F3[4-9]",
                       "F4", "F41", "F410", "F43", "F4[0, 2, 4-8]",
                       "F5", "F50", "F50[0-1]", "F50[2-3]", "F50[4-9]", "F5[1-9]",
                       "F6", "F60", "F6[1-9]"),
                     c("Substance use disorders", "Alcohol use disorders", "Drug use disorders",
                       "Opioid use disorders", "Cannabis use disorders", "Other non-alcohol substance use disorders",
                       "Psychotic disorders", "Schizophrenia", "Other psychotic disorders",
                       "Mood disorders", "Bipolar disorder", "Depression", "Other mood disorders",
                       "Anxiety disorders", "Other anxiety disorders", "Panic disorder", 
                       "Reaction to severe stress, and adjustment disorders", "All other anxiety disorders", 
                       "Behavioural syndromes", "Eating disorders", "Anorexia nervosa", 
                       "Bulimia nervosa", "Other eating disorders", "Other behavioural syndromes", 
                       "Personality disorders", "Specific personality disorders", "Other personality disorders"))

#Combine tables

bind_rows(tab1, tab2, tab3, tab4) %>%
 pivot_wider(names_from = "analysis",
             values_from = "estimate") %>%
 select(outcome, outcome_names, 3, 5, 6, 4) %>%
 arrange(match(outcome, outcomes)) %>%
 write.csv(file = "path/Results/Cox_ph_models_combined_senstivity_analyses.csv",
           row.names = FALSE)
