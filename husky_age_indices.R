# ------------------------------------------------------------------------------
# Here comparison is based on each synis_id being allocated to a strata
#  I.e. the approach use up to and including year 2017
#  Note that the max year here is 2016, because husky::stratas_df has not been
#       updated
# This approach never gave the same result as husky - see comparison at the
#      bottom. Code kept here for historical reasons

library(tidyverse)
library(fjolst)    # for the time being
library(pax)
stratas <-
  husky::stratas_df %>%
  select(strata, area = rall.area)
SPECIES <- 1
lengthclass <- c(seq(4.5, 109.5, by = 5), 119.5, 139.5)
ind <- c(31931, 31932, 32131, 36731, 37031, 37131, 37132, 37231, 41431, 41531, 42231, 42232, 47431, 52331)
st1 <-
  bind_rows(husky::STODVAR) %>%
  filter(tognumer < 20 | index %in% ind) %>%
  mutate(region = ifelse(area %in% c(1, 9, 10), "south",
                         ifelse(area %in% c(2:8), "north", NA))) %>%
  filter(!is.na(region)) %>%
  select(synis.id, ar, towlength = toglengd, region, strata = newstrata)

st2 <-
  bind_rows(husky::STODVAR) %>%
  filter(area %in% 1:10) %>%
  bind_rows(lesa.stodvar(leidangur="A4-2001")) %>%
  mutate(region = ifelse(area %in% c(1, 9, 10), "south",
                         ifelse(area %in% c(2:8), "north", NA))) %>%
  filter(!is.na(region)) %>%
  select(synis.id, ar, towlength = toglengd, region, strata = newstrata)

x <- calc_age_indices(st_length = st1, st_ototliths = st2, species = 1)
x$aggr %>%
  mutate(n = round(n/1e6, 2)) %>%
  select(-n.cv) %>%
  filter(aldur %in% 1:11) %>% # just to fit things on the screen
  spread(aldur, n) %>%
  as.data.frame()

# Comparison with the vatican
attach("/net/hafkaldi/export/u2/reikn/R/SurveyWork/SMB/AgeIndex/codindexold.rdata")
x2 <-
  codindextable %>%
  filter(reg == "Tot") %>%
  select(year, age, n.hoski = fj, cv.hoski = cv) %>%
  as_tibble()
x$aggr %>%
  rename(year = ar, age = aldur, cv = n.cv) %>%
  left_join(x2) %>%
  mutate(p.n = 1- (n / 1e6) / n.hoski,
         p.cv = 1 - cv / cv.hoski) %>%
  select(year, age, p.n, p.cv) %>%
  gather(variable, value, p.n:p.cv) %>%
  ggplot(aes(year, value)) +
  geom_col() +
  facet_grid(age ~ variable, scale = "free_y") +
  labs(x = NULL, y = "1 - n.pax / n.hoski")
