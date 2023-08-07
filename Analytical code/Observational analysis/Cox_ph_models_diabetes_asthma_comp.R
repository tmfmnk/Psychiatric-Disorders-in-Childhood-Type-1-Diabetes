#Libraries

library(tidyverse)
library(broom)
library(survival)
library(data.table)
library(patchwork)

#Import data

load(file = "path/Data/Data_models/data_models_exact_matching_w_mortality_asthma_comp.RData")

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

#Indices of datasets where both groups had at least 5 outcomes

ind <- map_lgl(.x = data_models,
               ~ .x %>%
                group_by(group) %>%
                summarise(across(ends_with("_binary"), sum)) %>%
                ungroup() %>%
                summarise(cond = all(across(ends_with("_binary")) >= 5)) %>%
                pull(cond))

#Stratified Cox proportional hazards models
#Fully adjusted

mfull_diabetes <- map2_dfr(.x = outcomes[ind],
                           .y = data_models[ind],
                           ~ tidy(coxph(as.formula(paste("Surv(", 
                                                         paste0("`", paste0(.x, "_time_to_event"), "`"),
                                                         ", ",
                                                         paste0("`", paste0(.x, "_binary"), "`"),
                                                         ")",
                                                         "~ group + VEK + as.factor(POHL) + month_discharge + year_discharge + strata(RODCIS2_exposed)")),
                                        data = .y),
                                  conf.int = TRUE,
                                  exponentiate = TRUE) %>%
                             mutate(outcome = .x) %>%
                             filter(term == "groupexposed"))
#Combine results

mfull_diabetes %>%
 transmute(outcome,
           estimate = paste(formatC(round(estimate, 2), format = "f", digits = 2),
                            paste0("(", 
                                   formatC(round(conf.low, 2), format = "f", digits = 2),
                                   "; ",
                                   formatC(round(conf.high, 2), format = "f", digits = 2),
                                   ")"))) %>%
 mutate(outcome_names = names(outcomes[match(outcome, outcomes)])) %>%
 select(outcome, outcome_names, estimate) %>%
 arrange(match(outcome, outcomes)) %>%
 mutate(analysis = "asthma comparison") %>%
 write.csv(file = "path/Results/Cox_ph_models_ex_matching_w_mortality_asthma_comp.csv",
           row.names = FALSE)

#Plotting the results from stratified Cox proportional hazards models
#Subplot with HRs

reg_plot <- mfull_diabetes %>%
 mutate(outcome_names = names(outcomes[match(outcome, outcomes)])) %>%
 arrange(match(outcome, outcomes)) %>%
 mutate(outcome_names = factor(outcome_names, unique(outcome_names)),
        outcome_names = fct_rev(outcome_names),
        color = rep(c("white", "gray95"), length.out = n())) %>%
 {ggplot(data = ., aes(x = estimate, y = outcome_names, xmin = conf.low, xmax = conf.high)) +
   geom_hline(aes(yintercept = outcome_names, color = color), size = 7) + 
   geom_pointrange(shape = 22, fill = "black") +
   geom_vline(xintercept = 1, linetype = 3) +
   theme_classic() +
   scale_colour_identity() +
   scale_x_log10(limits = c(0.1, 10), 
                 breaks = c(0.1, 0.3, 0.5, 0.7, 1, 1.25, 1.5, 2, 3, 5, 7, 10)) +
   xlab("aHR (95% CI)") +
   theme(axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.ticks.x = element_blank())}

#Subplot with counts
#Prepare data with counts

count_tab <- map_dfr(.x = data_models,
                     ~ .x %>%
                      group_by(group) %>%
                      summarise(across(ends_with("binary"), 
                                       ~ sub("_binary", "", cur_column()), 
                                       .names = "outcome"),
                                n = paste0(pull(across(ends_with("binary"), ~ sum(.))), "/", n())) %>%
                      ungroup() %>%
                      pivot_wider(names_from = group,
                                  values_from = n,
                                  names_prefix = "n_") %>%
                      mutate(outcome_names = names(outcomes[match(outcome, outcomes)]),
                             outcome_names = reorder(outcome_names, desc(match(outcome, outcomes))))) %>%
 left_join(mfull_diabetes %>%
            transmute(outcome,
                      estimate = paste(formatC(round(estimate, 2), format = "f", digits = 2),
                                       paste0("(", 
                                              formatC(round(conf.low, 2), format = "f", digits = 2),
                                              "; ",
                                              formatC(round(conf.high, 2), format = "f", digits = 2),
                                              ")"))),
           by = c("outcome")) %>%
 mutate(outcome_names = names(outcomes[match(outcome, outcomes)])) %>%
 arrange(match(outcome, outcomes)) %>%
 mutate(outcome_names = factor(outcome_names, unique(outcome_names)),
        outcome_names = fct_rev(outcome_names),
        outcome_names = fct_relabel(outcome_names, 
                                    ~ ifelse(.x %in% c("Substance use disorders", 
                                                       "Psychotic disorders", 
                                                       "Mood disorders", 
                                                       "Anxiety disorders", 
                                                       "Behavioural syndromes", 
                                                       "Personality disorders"),
                                             .x,
                                             paste0("       ", .x))),
        color = rep(c("white", "gray95"), length.out = n())) 

#Prepare subplot 

count_plot <- count_tab %>%
 pivot_longer(-c(outcome, outcome_names, color),
              values_transform = list(value = as.character),
              names_to = "variables",
              values_to = "values") %>%
 mutate(variables = factor(variables, 
                           levels = c("n_unexposed", "n_exposed", "estimate"))) %>%
 ggplot(aes(x = variables, y = outcome_names, label = values)) +
 geom_hline(aes(yintercept = outcome_names, color = color), size = 7) +
 geom_text(size = 3) +
 scale_x_discrete(position = "top", 
                  labels = c("Unexposed \nevents/total", "Exposed \nevents/total", "aHR \n(95% CI)")) +
 scale_colour_identity() +
 labs(y = NULL, x = NULL) +
 theme_classic() +
 theme(strip.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.border = element_blank(),
       axis.line = element_blank(),
       axis.text.y = element_text(hjust = 0),
       axis.ticks = element_blank(),
       axis.title = element_text(face = "bold"))

#Combine subplots

count_plot + reg_plot + plot_layout(widths = c(5, 10))

ggsave(file = "path/Results/HR_plot_asthma_comp.eps",
       device = "eps",
       width = 40,
       height = 20,
       units = "cm",
       dpi = 300)
