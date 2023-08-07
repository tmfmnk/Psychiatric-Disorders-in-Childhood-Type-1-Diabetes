#Libraries

library(data.table)
library(tidyverse)
library(purrr)
library(lubridate)

#Variable definition

#RODCIS2 = unique personal identifier
#DATPRI = date of admission
#DATUKO = date of discharge
#VEK = age
#POHL = sex
#NBYDL = region of residence
#ZDG = primary diagnosis
#DDG2 = secondary diagnosis
#DAUMR = date of death

#Import all hospitalizations from 1994 to 2017

hospitalizations_1994_2015 <- list.files(path = "path/", 
                                         full.names = TRUE) %>%
 .[str_detect(., "1[6-7].csv$", negate = TRUE)] %>%
 map_dfr(~ fread(.,
                 select = c("RODCIS2" = "character", 
                            "DATPRI" = "integer", 
                            "DATUKO" = "integer", 
                            "VEK"= "integer", 
                            "POHL" = "integer",
                            "NBYDL"= "character",
                            "ZDG" = "character",
                            "DDG2" = "character",
                            "DDG3" = "character",
                            "DDG4" = "character",
                            "DDG5" = "character"),
                 header = TRUE,
                 sep = ",",
                 dec = ".",
                 fill = TRUE,
                 encoding = "Latin-1",
                 nThread = 8))

hospitalizations_2016_2017 <- list.files(path = "path/", 
                                         full.names = TRUE) %>%
 .[str_detect(., "1[6-7].csv$", negate = FALSE)] %>%
 map_dfr(~ fread(.,
                 select = c("RODCIS" = "character", 
                            "DATPRI" = "integer", 
                            "DATUKO" = "integer", 
                            "VEK"= "integer", 
                            "POHL" = "integer",
                            "NBYDL"= "character",
                            "ZDG" = "character",
                            "DDG2" = "character"),
                 header = TRUE,
                 sep = ";",
                 dec = ".",
                 fill = TRUE,
                 encoding = "Latin-1",
                 nThread = 8)) %>%
 rename(RODCIS2 = RODCIS)

#Combine all data 

hospitalizations_1994_2017 <- hospitalizations_1994_2015 %>%
 bind_rows(hospitalizations_2016_2017)

#Remove partial data

rm(hospitalizations_1994_2015)
rm(hospitalizations_2016_2017)

#Import data on mortality

deaths_1994_2013 <- fread(file = "path/zem_1994_2013.csv")
deaths_2014 <- fread(file = "path/zem_2014.csv")
deaths_2015 <- fread(file = "path/zem_2015.csv")
deaths_2016 <- fread(file = "path/zem_2016.csv")
deaths_2017 <- fread(file = "path/zem_2017.csv")

#Unifying the format of mortality data

deaths_1994_2017 <- deaths_1994_2013 %>%
 transmute(RODCIS2 = RC,
           DAUMR = dmy(DAUMR),
           cause_of_death = trimws(DGP),
           external_cause_of_death = trimws(DGE)) %>%
 bind_rows(deaths_2014 %>%
            transmute(RODCIS2 = RC,
                      DAUMR = dmy(DAUMR),
                      cause_of_death = trimws(DGP),
                      external_cause_of_death = trimws(DGE)),
           deaths_2015 %>%
            transmute(RODCIS2 = RODCIS2,
                      DAUMR = ymd(DAUMR),
                      cause_of_death = trimws(DGP),
                      external_cause_of_death = trimws(DGE)),
           deaths_2016 %>%
            transmute(RODCIS2 = RCZEMAN2,
                      DAUMR = ymd(paste0(UMROK, str_pad(UMRMM, 2, pad = "0"), UMRDD)),
                      cause_of_death = trimws(DGUMR),
                      external_cause_of_death = trimws(DGUMR2)),
           deaths_2017 %>%
            transmute(RODCIS2 = RCZEMAN2,
                      DAUMR = ymd(paste0(UMROK, str_pad(UMRMM, 2, pad = "0"), UMRDD)),
                      cause_of_death = trimws(DGUMR),
                      external_cause_of_death = trimws(DGUMR2)))

#Remove partial data

rm(deaths_1994_2013)
rm(deaths_2014)
rm(deaths_2015)
rm(deaths_2016)
rm(deaths_2017)

#Excluding records with missing values on key variables
#Excluding records with invalid dates

hospitalizations_1994_2017 <- hospitalizations_1994_2017 %>%
 filter(rowSums(is.na(across(c(RODCIS2, DATPRI, DATUKO, VEK, POHL, NBYDL, ZDG)))) == 0) %>%
 filter(!is.na(ymd(DATPRI)) & !is.na(ymd(DATUKO))) 

#Excluding individuals with more than one date of death
#Excluding individuals with hospitalizations after death

hospitalizations_1994_2017 <- hospitalizations_1994_2017 %>%
 anti_join(hospitalizations_1994_2017 %>%
            inner_join(deaths_1994_2017, 
                       by = c("RODCIS2" = "RODCIS2")) %>%
            mutate(DATUKO = ymd(DATUKO)) %>%
            group_by(RODCIS2) %>%
            filter(DAUMR < max(DATUKO) | n_distinct(DAUMR) > 1) %>%
            ungroup(),
           by = c("RODCIS2" = "RODCIS2"))

#Excluding records with invalid date overlaps (discharge date after the admission date of another record)

hospitalizations_1994_2017 <- hospitalizations_1994_2017 %>%
 anti_join(map_dfr(.x = hospitalizations_1994_2017 %>%
                    group_split(split_ID = frank(RODCIS2, ties.method = "dense") %/% 10000),
                   ~ .x %>% 
                    select(RODCIS2,
                           DATPRI_index = DATPRI,
                           DATUKO_index = DATUKO) %>%
                    inner_join(.x %>%
                                select(RODCIS2,
                                       DATPRI_non_index = DATPRI,
                                       DATUKO_non_index = DATUKO), 
                               by = c("RODCIS2" = "RODCIS2")) %>%
                    filter(DATPRI_non_index < DATPRI_index & DATUKO_non_index > DATPRI_index) %>%
                    pivot_longer(-RODCIS2,
                                 names_to = c(".value", "type"), 
                                 names_pattern = "([^_]+)_(.*)")),
           by = c("RODCIS2" = "RODCIS2",
                  "DATPRI" = "DATPRI",
                  "DATUKO" = "DATUKO"))

#Establishing the exposed cohort
#Keeping hospitalizations between 1994 and 2007
#Keeping individuals aged 0-14
#Excluding individuals with residence outside of Czechia

hospitalizations_exp_1994_2007 <- hospitalizations_1994_2017 %>%
 filter(year(ymd(DATPRI)) >= 1994 & year(ymd(DATUKO)) <= 2007) %>%
 filter(VEK <= 14) %>%
 mutate(ZDG = trimws(ZDG),
        DDG2 = trimws(DDG2),
        DDG3 = trimws(DDG3),
        DDG4 = trimws(DDG4),
        DDG5 = trimws(DDG5)) %>%
 filter(grepl("^E10", ZDG) | 
         grepl("^E10", DDG2) | 
         grepl("^E10", DDG3) | 
         grepl("^E10", DDG4) | 
         grepl("^E10", DDG5))

#Selecting the first hospitalization in the examined time period

set.seed(123)
first_hospitalizations_exp_1994_2007 <- hospitalizations_exp_1994_2007 %>%
 group_by(RODCIS2, DATPRI, DATUKO) %>%
 filter(row_number() == sample(1:n(), 1)) %>%
 group_by(RODCIS2) %>%
 filter(DATPRI == min(DATPRI)) %>%
 filter(DATUKO == min(DATUKO)) %>%
 ungroup() %>%
 filter(!grepl("^99", NBYDL)) %>%
 filter(!(grepl("^F1|^F2|^F3|^F4|^F5|^F6", ZDG) | grepl("^F1|^F2|^F3|^F4|^F5|^F6", DDG2) | grepl("^F1|^F2|^F3|^F4|^F5|^F6", DDG3) | grepl("^F1|^F2|^F3|^F4|^F5|^F6", DDG4) | grepl("^F1|^F2|^F3|^F4|^F5|^F6", DDG5)))

rm(hospitalizations_exp_1994_2007)

#Establishing the unexposed cohort
#Keeping hospitalizations between 1994 and 2007
#Keeping individuals aged 0-14

hospitalizations_unexp_1994_2007 <- hospitalizations_1994_2017  %>%
 filter(year(ymd(DATPRI)) >= 1994 & year(ymd(DATUKO)) <= 2007) %>%
 filter(VEK <= 14) %>%
 mutate(ZDG = trimws(ZDG),
        DDG2 = trimws(DDG2),
        DDG3 = trimws(DDG3),
        DDG4 = trimws(DDG4),
        DDG5 = trimws(DDG5)) %>%
 mutate(cond = DATPRI < first_hospitalizations_exp_1994_2007$DATPRI[match(RODCIS2, first_hospitalizations_exp_1994_2007$RODCIS2)]) %>%
 filter(is.na(cond) | cond == TRUE) %>%
 select(-cond)

#Randomly selecting one record when multiple records with the same admission and discharge date are present

set.seed(123)
hospitalizations_unexp_1994_2007 <- hospitalizations_unexp_1994_2007 %>%
 group_by(RODCIS2, DATPRI, DATUKO) %>%
 filter(row_number() == sample(1:n(), 1)) %>%
 ungroup()

#Matching
#Randomly selecting 10 unexposed individuals per each exposed individual

set.seed(123)
matched_sampled_pairs <- first_hospitalizations_exp_1994_2007 %>%
 transmute(RODCIS2_exposed = RODCIS2,
           VEK,
           POHL,
           month_discharge = month(ymd(DATUKO)),
           year_discharge = year(ymd(DATUKO))) %>%
 inner_join(hospitalizations_unexp_1994_2007 %>%
             mutate(month_discharge = month(ymd(DATUKO)),
                    year_discharge = year(ymd(DATUKO))),
            by = c("VEK" = "VEK",
                   "POHL" = "POHL",
                   "month_discharge" = "month_discharge",
                   "year_discharge" = "year_discharge")) %>%
 filter(RODCIS2_exposed != RODCIS2) %>%
 filter(!grepl("^99", NBYDL)) %>%
 filter(!(grepl("^F1|^F2|^F3|^F4|^F5|^F6", ZDG) | grepl("^F1|^F2|^F3|^F4|^F5|^F6", DDG2) | grepl("^F1|^F2|^F3|^F4|^F5|^F6", DDG3) | grepl("^F1|^F2|^F3|^F4|^F5|^F6", DDG4) | grepl("^F1|^F2|^F3|^F4|^F5|^F6", DDG5))) %>%
 group_by(RODCIS2) %>%
 filter(row_number() %in% sample(1:n(), 1)) %>%
 group_by(RODCIS2_exposed) %>%
 filter(row_number() %in% sample(1:n(), 10)) %>%
 ungroup()

#Combine datasets

cohorts_baseline <- first_hospitalizations_exp_1994_2007 %>%
 mutate(RODCIS2_exposed = RODCIS2,
        group = "exposed") %>%
 mutate(month_discharge = month(ymd(DATUKO)),
        year_discharge = year(ymd(DATUKO))) %>%
 bind_rows(matched_sampled_pairs %>%
            mutate(group = "unexposed"))

#Establish the presence of mental disorders up to 2017

outcomes <- map(.x = c("^F1", "^F10", "^F1[1-9]", "^F11", "^F12", "^F1[3-9]",
                       "^F2", "^F20", "^F2[1-9]",
                       "^F3", "^F3[0-1]", "^F3[2-3]", "^F3[4-9]",
                       "^F4", "^F41", "^F410", "^F43", "^F4[0, 2, 4-8]",
                       "^F5", "^F50", "^F50[0-1]", "^F50[2-3]", "^F50[4-9]", "^F5[1-9]",
                       "^F6", "^F60", "^F6[1-9]"),
                ~ cohorts_baseline %>%
                 transmute(RODCIS2,
                           followup_start = ymd(DATUKO),
                           followup_end = ymd("2017-12-31")) %>%
                 inner_join(hospitalizations_1994_2017,
                            by = c("RODCIS2" = "RODCIS2")) %>%
                 mutate(hospitalization_start = ymd(DATPRI),
                        hospitalization_end = ymd(DATUKO)) %>%
                 filter(int_overlaps(interval(hospitalization_start, hospitalization_end),
                                     interval(followup_start, followup_end))) %>%
                 mutate(ZDG = trimws(ZDG),
                        DDG2 = trimws(DDG2),
                        DDG3 = trimws(DDG3),
                        DDG4 = trimws(DDG4),
                        DDG5 = trimws(DDG5)) %>%
                 mutate(!!sub("\\^", "", .x) := grepl(.x, ZDG) | grepl(.x, DDG2) | grepl(.x, DDG3) | grepl(.x, DDG4) | grepl(.x, DDG5)) %>%
                 mutate(across(any_of(sub("\\^", "", .x)), ~ replace(., . == FALSE, NA))) %>%
                 pivot_longer(any_of(sub("\\^", "", .x)),
                              names_to = "outcomes",
                              values_to = "binary",
                              values_drop_na = TRUE) %>%
                 group_by(RODCIS2, 
                          outcomes) %>%
                 summarise(first_hosp_start_date = min(hospitalization_start),
                           binary = first(binary)) %>%
                 ungroup() %>%
                 pivot_wider(names_from = outcomes,
                             names_glue = "{outcomes}_{.value}",
                             values_from = c("first_hosp_start_date",
                                             "binary"))) 

#Combine datasets

data_models <- map(.x = outcomes,
                   ~ cohorts_baseline %>%
                    left_join(.x,
                              by = c("RODCIS2" = "RODCIS2")) %>%
                    left_join(deaths_1994_2017,
                              by = c("RODCIS2" = "RODCIS2")) %>%
                    mutate(across(ends_with("_binary"), 
                                  ~ as.numeric(replace(., is.na(.), FALSE))),
                           across(ends_with("_first_hosp_start_date"), 
                                  ~ as.numeric(ifelse(is.na(.), pmin(DAUMR - ymd(DATUKO), ymd("2017-12-31") - ymd(DATUKO), na.rm = TRUE), . - ymd(DATUKO))),
                                  .names = '{stringr::str_remove(.col, "_first_hosp_start_date")}_time_to_event')) %>%
                    mutate(group = factor(group, levels = c("unexposed", "exposed"))))

#Export data

save(data_models, 
     file = "path/Data/Data_models/data_models_exact_matching_w_mortality.RData")
