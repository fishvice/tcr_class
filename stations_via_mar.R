library(mar)
STODVAR.all.B <-
  lesa_stodvar(con) %>%
  filter(synaflokkur == 30,
         # Two tows from 2009 are excluded, why?
         !synis_id %in% c(330573, 330931)) %>%
  mutate(index = reitur * 100 + tognumer,
         lon = ifelse(is.na(hift_v_lengd),
                      kastad_v_lengd,
                      (kastad_v_lengd + hift_v_lengd) / 2),
         lat = ifelse(is.na(hift_n_breidd),
                      kastad_n_breidd,
                      (kastad_n_breidd + hift_n_breidd) / 2)) %>%
  # the smb_index_strata only contain strata for tognumer <40
  left_join(tbl_mar(con, "ops$einarhj.smb_index_strata"),
            by = "index") %>%
  # Need to collect here because of steps below
  collect(n = Inf) %>%
  # Tognumer >40 are "variable stations" each year
  #  Bjarki: Could we use reitur-smareitur???
  mutate(oldstrata = husky::inside.strata(., stratalist = husky::ralllist,
                                          stratas = husky::STRATAS) %>%
           pull(newstrata),
         newstrata = husky::inside.strata(., stratalist = husky::smblist$nr,
                                          stratas = husky::NEWSTRATAS) %>%
           pull(newstrata)) %>%
  # Smá lagfæringar
  mutate(oldstrata = if_else(oldstrata %in% 23, 88, oldstrata),
         newstrata = if_else(index %in% 71741 & ar == 2009, 37, newstrata),
         newstrata = if_else(index %in% 61141 & ar == 2010, 42, newstrata)) %>%
  # Create only one column
  mutate(oldstrata = if_else(tognumer < 40, stdoldstrata, oldstrata),
         newstrata = if_else(tognumer < 40, stdnewstrata, newstrata)) %>%
  select(-stdoldstrata, -stdnewstrata)

# ------------------------------------------------------------------------------
# Comparison
load("/net/hafkaldi/export/u2/reikn/R/SurveyWork/SMB/Stations.rdata")

tidy <-
  STODVAR.all.B %>%
  select(synis_id, reitur, tognumer, index, oldstrata, newstrata) %>%
  arrange(synis_id) %>%
  mutate_if(is.double, as.integer)
hoski <-
  STODVAR.all %>%
  select(synis_id = synis.id, reitur, tognumer, index, oldstrata, newstrata) %>%
  arrange(synis_id) %>%
  as_tibble() %>%
  mutate_if(is.double, as.integer)
identical(tidy, hoski)
