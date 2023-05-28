library(tidyverse)

jobs <- as.data.frame(
    memisc::as.data.set(
        memisc::spss.portable.file(
            file.path("~/Google Drive/A.PAPERS/PIsensAnalysis",
                      "JOBSII_DataFolder",
                      "DS0011",
                      "02739-0011-Data.por")))
) %>%
    filter(!is.na(V9106), RISK01=="high") %>%

    mutate(treat      = ifelse(V9106=="EXP", 1, 0),
           complier   = 1 * (V9107=="SHOW"),

           age        = V9002,
           sex        = V9001,
           race       = V1424,
           edu        = V1531,
           marital    = V1407,
           hh.kids    = V1408,
           hh.income  = V1532,
           econ.hard  = V1503,
           occu       = V1401A,
           wks.unemp  = V1405,
           part.motiv = V9101,
           seek.motiv = V1535,
           seek.effi  = V1510,
           assertive  = V1522,
           depress1   = V1518,
           depress3   = V3518,
           work3      = V3461,
           earn3      = V3404,
           id         = V1
    ) %>%

    # recode work3 from 0 or NA to 1 if report positive earn3
    mutate(work3 = ifelse(!is.na(work3) & work3==1,
                          1,
                          ifelse(!is.na(earn3) & earn3>0,
                                 1,
                                 work3))) %>%
    mutate(earn3 = ifelse(!is.na(work3) & work3==0, 0, earn3)) %>%

    select(id, treat, complier,
           age, sex, race, edu, marital, hh.kids,
           hh.income, econ.hard, occu, wks.unemp,
           part.motiv, seek.motiv, seek.effi, assertive,
           depress1,
           work3, earn3, depress3)


px.fullX <- jobs %>%
    select(-c(complier, depress3, work3, earn3)) %>%
    drop_na() %>%
    pull(id)
px.complete <- jobs %>%
    select(-complier) %>%
    drop_na() %>%
    pull(id)

jobs <- jobs %>%
    mutate(fullX    = ifelse(id %in% px.fullX,    1, 0),
           complete = ifelse(id %in% px.complete, 1, 0))

rm(px.fullX, px.complete)


# analysis data
dat <- jobs %>%
    filter(complete==1) %>%
    select(-c(fullX, complete))

rm(jobs)


dat <- dat %>%
    # truncate 24 (12) extreme values of wks.unemp
    mutate(wks.unemp = ifelse(wks.unemp>52, 52, wks.unemp)) %>%

    # truncate 3 (1) extreme values of hh.kids
    mutate(hh.kids = ifelse(hh.kids>5, 5, hh.kids)) %>%

    # dichotomize race (collapsing small categories)
    mutate(race = ifelse(race=="WHITE", "white", "nonwhite"),
           race = factor(race)) %>%

    # collapse marital to three categories
    mutate(marital = ifelse(marital=="NEVMARR",
                            "nevmarr",
                            ifelse(marital=="MARRIED",
                                   "married",
                                   "divsepwid")),
           marital = factor(marital, levels = c("nevmarr",
                                                "married",
                                                "divsepwid"))) %>%

    mutate(sex = factor(sex, levels = c(0, 1), labels = c("male", "female")))

levels(dat$occu) <- c("professional",
                      "managerial",
                      "clerical",
                      "sales",
                      "crafts/foremen",
                      "operative",
                      "labor/service")


saveRDS(dat, file = "data-raw/dat.rds")
