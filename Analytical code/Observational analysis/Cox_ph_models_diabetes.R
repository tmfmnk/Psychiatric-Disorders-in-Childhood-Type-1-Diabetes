#Libraries

library(tidyverse)
library(broom)
library(survival)
library(data.table)
library(EValue)
library(patchwork)
library(survminer)
library(ggpubr)
library(scales)

#Import data

load(file = "path/Data/Data_models/data_models_exact_matching_w_mortality.RData")

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

#Testing the proportionality assumption using Schoenfeld residuals

map2(.x = outcomes,
         .y = data_models,
         ~ cox.zph(coxph(as.formula(paste("Surv(", 
                                       paste0("`", paste0(.x, "_time_to_event"), "`"),
                                       ", ",
                                       paste0("`", paste0(.x, "_binary"), "`"),
                                       ")",
                                       "~ group + VEK + as.factor(POHL) + month_discharge + year_discharge + strata(RODCIS2_exposed)")),
                      data = .y)))

#Stratified Cox proportional hazards models
#Fully adjusted

mfull_diabetes <- map2_dfr(.x = outcomes,
                           .y = data_models,
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
 write.csv(file = "path/Results/Cox_ph_models_ex_matching_w_mortality.csv",
           row.names = FALSE)

#E-value

mfull_diabetes %>%
      transmute(outcome,
                e_value = pmap_dbl(across(c(estimate, conf.low, conf.high)), 
                                                     ~ as_tibble(evalues.HR(..1, ..2, ..3, rare = TRUE), 
                                                                 rownames = "name") %>%
                                                      filter(name == "E-values") %>%
                                                      pluck("point")),
                e_value = ifelse(data.table::between(1, conf.low, conf.high), NA_real_, e_value),
                e_value = formatC(round(e_value, 2), format = "f", digits = 2)) %>%
 mutate(outcome_names = names(outcomes[match(outcome, outcomes)])) %>%
 select(outcome, outcome_names, e_value) %>%
 arrange(match(outcome, outcomes)) %>%
 write.csv(file = "path/Results/E_values.csv",
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
   scale_x_log10(limits = c(0.25, 13), 
                 breaks = c(0.25, 0.4, 0.6, 0.8, 1, 1.25, 1.5, 2, 3, 5, 7, 10, 13)) +
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

ggsave(file = "path/Results/HR_plot.eps",
       device = "eps",
       width = 40,
       height = 20,
       units = "cm",
       dpi = 300)

#Survival probability plots
#Define custom plotting function

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title.x = element_text(face = "bold", size = 10),
          axis.title.y = element_text(face = "bold", size = 10, angle = 90),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10),
          axis.ticks.x = element_blank(),    
          axis.ticks.y = element_blank(),    
          legend.text = element_text(face = "bold", size = 10),
          legend.title = element_blank())
}

#Programmatic ploting
#Diabetes cohort
           
pmap(list(surv_fit(formula = set_names(map(.x = outcomes, 
                                           ~ as.formula(paste("Surv(", 
                                                              paste0("`", paste0(.x, "_time_to_event"), "`"),
                                                              ", ",
                                                              paste0("`", paste0(.x, "_binary"), "`"),
                                                              ")",
                                                              "~ group"))), 
                                       outcomes),
                   data = data_models,
                   match.fd = TRUE),
          outcomes,
          map_dbl(.x = data_models,
                  ~ .x %>%
                    summarise(across(ends_with("time_to_event"), max)) %>%
                   pull())),
     
     function(...) {
       
       p <- ggsurvplot(..1, 
                       fun = "event",
                       xlim = c(0, 9000),
                       ylim = c(0.0, 0.1),
                       break.x.by = ceiling((..3 + 1)/4),
                       conf.int = TRUE, 
                       conf.int.alpha = 0.3,
                       risk.table = TRUE,
                       cumcensor = TRUE,
                       cumevents = TRUE,
                       legend.labs = c("Individuals with childhood-onset type 1 diabetes", "Unexposed counterparts"),
                       title = paste0(..2, " outcome by childhood-onset type 1 diabetes status"),
                       xlab = "Time (days) since the index hospitalization",
                       ylab = "Cumulative event (95% CI)",
                       palette = c("blue", "red"),
                       legend.title = "",
                       ggtheme = custom_theme())
       
       p1 = p$plot
       p2 = p$table
       p3 = p$ncensor.plot
       p4 <- p$cumevents
       plots = cowplot::plot_grid(p1, p2, p3, p4, align = "v", ncol = 1, rel_heights = c(4, 1, 1, 1))
       
       ggsave(plot = plots,
              filename = paste0("Cumulative_event_plot_diabetes_", ..2, ".png"),
              path = "path/Results/Cumulative_event_plots/",
              device = "png",
              width = 12, 
              height = 7, 
              dpi = 300)
     }
)

#Models stratified by sex
#Males
#Indices of datasets where both groups had at least 5 outcomes

ind_males <- map_lgl(.x = data_models,
                     ~ .x %>%
                      filter(POHL == 1) %>%
                      group_by(group) %>%
                      summarise(across(ends_with("_binary"), sum)) %>%
                      ungroup() %>%
                      summarise(cond = all(across(ends_with("_binary")) >= 5)) %>%
                      pull(cond))

mfull_diabetes_males <- map2_dfr(.x = outcomes[ind_males],
                                 .y = data_models[ind_males],
                                 ~ tidy(coxph(as.formula(paste("Surv(", 
                                                               paste0("`", paste0(.x, "_time_to_event"), "`"),
                                                               ", ",
                                                               paste0("`", paste0(.x, "_binary"), "`"),
                                                               ")",
                                                               "~ group + VEK + month_discharge + year_discharge + strata(RODCIS2_exposed)")),
                                              subset = POHL == 1,
                                              data = .y),
                                        conf.int = TRUE,
                                        exponentiate = TRUE) %>%
                                  mutate(outcome = .x) %>%
                                  filter(term == "groupexposed"))


#Females
#Indices of datasets where both groups had at least 5 outcomes

ind_females <- map_lgl(.x = data_models,
                     ~ .x %>%
                      filter(POHL == 2) %>%
                      group_by(group) %>%
                      summarise(across(ends_with("_binary"), sum)) %>%
                      ungroup() %>%
                      summarise(cond = all(across(ends_with("_binary")) >= 5)) %>%
                      pull(cond))

mfull_diabetes_females <- map2_dfr(.x = outcomes[ind_females],
                                   .y = data_models[ind_females],
                                   ~ tidy(coxph(as.formula(paste("Surv(", 
                                                                 paste0("`", paste0(.x, "_time_to_event"), "`"),
                                                                 ", ",
                                                                 paste0("`", paste0(.x, "_binary"), "`"),
                                                                 ")",
                                                                 "~ group + VEK + month_discharge + year_discharge + strata(RODCIS2_exposed)")),
                                                subset = POHL == 2,
                                                data = .y),
                                          conf.int = TRUE,
                                          exponentiate = TRUE) %>%
                                    mutate(outcome = .x) %>%
                                    filter(term == "groupexposed"))

#Combine results

imap(mget(ls(pattern = "mfull_diabetes_males|mfull_diabetes_females")),
     ~ .x %>%
      transmute(outcome,
                !!.y := paste(formatC(round(estimate, 2), format = "f", digits = 2),
                              paste0("(", 
                                     formatC(round(conf.low, 2), format = "f", digits = 2),
                                     "; ",
                                     formatC(round(conf.high, 2), format = "f", digits = 2),
                                     ")")))) %>%
 reduce(left_join,
        by = "outcome") %>%
 mutate(outcome_names = names(outcomes[match(outcome, outcomes)])) %>%
 select(outcome, outcome_names, contains("diabetes")) %>%
 arrange(match(outcome, outcomes)) %>%
 write.csv(file = "path/Results/Cox_ph_models_ex_matching_w_mortality_sex_stratified.csv",
           row.names = FALSE)
