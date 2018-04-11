# Survey indices  via mar
library(mar)

# ------------------------------------------------------------------------------
# Additional functions

#' Calculate overall cv from stratfied summary statistics
#'
#' @param m Mean value within strata
#' @param s Standard deviation within strata
#' @param area The area (e.g. survey strata area)
#' @param n number of samples within strata
#'
#' @export
#'
calc_cv <- function(m, s, area, n) {

  Mean = sum(m * area) / sum(area)
  Sum = sum(m * area)
  tmpsum = sum(m[!is.na(s)] * area[!is.na(s)])
  Calc.sdev = sqrt(sum(s[!is.na(s)]^2 * area[!is.na(s)]^2/  n[!is.na(s)])   / sum(area[!is.na(s)])^2)
  Sdev = Calc.sdev * Sum/tmpsum
  cv = Sdev/Mean

  return(cv)
}

# ------------------------------------------------------------------------------
# Constants
std.cv        <- 1             # The cv if only one station in a strata
std.towlength <- 4             # Standard tow length for SMB is 4 nautical miles
std.width     <- 17 / 1852     # Standard sweep width in nautical miles

min.towlength <- 2             # Minimum "acceptable" towlength
max.towlength <- 8             # Maximum "acceptable" towlength

# ------------------------------------------------------------------------------
# Variable stuff
Synaflokkur <- 30
Tognumer <- c(1:39, NA)
Species <- 1                      # Select a species
Length.min <- 5                   # Minimum length for indices calculation
Length.max <- 500                 # Maximum length for indices calculation

# ------------------------------------------------------------------------------
# Calculations

# A. Data gathering ------------------------------------------------------------

by.length <-
  # 1. get survey stations -----------------------------------------------------
  lesa_stodvar(con) %>%
  filter(synaflokkur == Synaflokkur) %>%
  mutate(index = reitur * 100 + tognumer) %>%
  select(synis_id, ar, index, reitur, tognumer, veidarfaeri, toglengd) %>%
  # 2. get length data ---------------------------------------------------------
  left_join(lesa_lengdir(con) %>%
            filter(tegund %in% Species,
                   lengd >= Length.min,
                   lengd < Length.max) %>%
            group_by(synis_id, tegund, lengd) %>%
            summarise(fjoldi = sum(fjoldi, na.rm = TRUE)),
          by = "synis_id") %>%
  ungroup() %>%
  # 0. A temporary fix, for zero stations --------------------------------------
  #     TODO: Find a more permanent solution so scripts works for more than
  #           one species (via group_by( ..., tegund))
  mutate(tegund = if_else(is.na(tegund), Species, tegund),
         lengd  = if_else(is.na(lengd), Length.min, lengd),
         fjoldi = if_else(is.na(fjoldi), 0, fjoldi)) %>%

  # 3. get count data ----------------------------------------------------------
  #    Note: If only counted at station this code approach fails
  left_join(lesa_numer(con) %>%
            mutate(r = ifelse(fj_talid==0 | is.na(fj_talid),
                                     1,
                                     1 + fj_talid / ifelse(fj_maelt == 0 | is.na(fj_maelt), 1, fj_maelt))) %>%
            select(synis_id, tegund, r),
            by = c("synis_id", "tegund")) %>%
  # 4. scale by counted --------------------------------------------------------
  mutate(N = r * fjoldi / 1e3) %>%   # units of thousand
  # 5. trim towlength ----------------------------------------------------------
  mutate(toglengd = if_else(toglengd > max.towlength, max.towlength, toglengd),
         toglengd = if_else(toglengd < min.towlength, min.towlength, toglengd)) %>%
  # 6.a standardize by towlength ------------------------------------------------
  mutate(N = N / toglengd * std.towlength) %>%      # standardize to per 4 miles
  # 6.b standardize to area swept
  #     this does not make much sense here because we already have this above
  #     need to pass this function further down in the code path
  mutate(N = N / if_else(veidarfaeri == 78, 1.25 * std.width, std.width)) %>%
  # 7. calculate_biomass from numbers, length and a and b ----------------------
  # 7.a get the length weight coefficients
  left_join(tbl_mar(con, "ops$einarhj.lwcoeff"),
            by = "tegund") %>%
  # 7.b use Newton's law if lwcoefficient for species not specified
  mutate(a = ifelse(is.na(a), 0.01, a),
         b = ifelse(is.na(b), 3.00, b),
         B  = ifelse(is.na(N), 0, N) * a * lengd^b / 1e3) %>%
  select(-a, -b) %>%
  # 8. Cumulative calculation
  #     Note: This step is a precursor for calculating things via
  #           abundance less than and biomass greater than
  arrange(tegund, synis_id, lengd) %>%
  group_by(tegund, synis_id) %>%
  mutate(cN = cumsum(N),
         cB = sum(B, na.rm = TRUE) - cumsum(B) + B) %>%
  ungroup()

# B. Summarise abundance and biomass by station --------------------------------
by.station <-
  by.length %>%
  # 8. summarise by station ----------------------------------------------------
  # NOTE: here is the first step where statistics by length is dropped
  #       some (minor) recoding above would be needed if one were to take things
  #       forward by each length class
  group_by(synis_id, index, reitur, tognumer, veidarfaeri, tegund, ar) %>%
  summarise(N = sum(N, na.rm = TRUE),
            B = sum(B, na.rm = TRUE)) %>%
  # Zero stations - THOUGHT THIS STEP SHOULD BE REDUNDANT
  mutate(N = ifelse(is.na(N), 0, N),
         B = ifelse(is.na(B), 0, B))


# C. Calculate mean and standard deviation of abundance and biomass of stations
#    within each strata and raise estimates by the area of the strata

by.strata <-
  by.station %>%
  # 9. filter stations ---------------------------------------------------------
  filter(tognumer %in% Tognumer) %>%
  # 10. summarise by strata ----------------------------------------------------
  # 10.a  Get the strata for each station
  left_join(tbl_mar(con, "ops$einarhj.smb_index_strata") %>%
            select(index, strata = stdoldstrata)) %>%
  # 10.b group by year and strata and calculate number of stations, mean and sd
  group_by(tegund, ar, strata) %>%
  summarise(sN  = n(),   # number of stations within strata
            n_m  = mean(N, na.rm = TRUE),
            n_d  = ifelse(n() == 1, mean(N, na.rm = TRUE) * std.cv, sd(N)),
            b_m  = mean(B, na.rm = TRUE),
            b_d  = ifelse(n() == 1, mean(B, na.rm = TRUE) * std.cv, sd(B))) %>%

  # 11. raise to strata area ---------------------------------------------------
  # 11.a get area of the strata
  left_join(tbl_mar(con, "ops$einarhj.oldstrataarea") %>%
            select(strata = oldstrata, area = rall.area) %>%
            #  area is above is in km2, here convert nm2
            mutate(area = area / 1.852^2)) %>%
  # 11.b do the strata raising
  mutate(n     = n_m  * area / std.towlength,
         b     = b_m  * area / std.towlength)

# D. Summarise data by year ----------------------------------------------------

by.year <-
  by.strata %>%
  # ----------------------------------------------------------------------------
  # Up to now we only have been operating within Oracle. I.e. sql-scripts via R.
  # Have to collect here because of calc_cv function in the year aggregate step
  # TODO: Fix that, do internally in Oracle
  collect(n = Inf) %>%
  # ----------------------------------------------------------------------------
  # some data fall outside strata, drop them - needs some double checking of code
  drop_na() %>%
  # 11. summarise by year ------------------------------------------------------
  group_by(tegund, ar) %>%
  summarise(n = sum(n, na.rm = TRUE),
            # A la Höski
            n.cv = calc_cv(n_m, n_d, area, sN),
            b = sum(b, na.rm = TRUE),
            # A la Höski
            b.cv = calc_cv(b_m, b_d, area, sN)) %>%
  ungroup()

# ------------------------------------------------------------------------------
# The results
glimpse(by.year)
by.year %>%
  ggplot(aes(ar + 3/12, b)) +
  geom_point() +
  geom_linerange(aes(ymin = b * (1 - b.cv),
                     ymax = b * (1 + b.cv))) +
  expand_limits(y = 0) +
  labs(x = NULL, y = NULL,
       title = "Spring survey biomass indices") +
  facet_wrap(~ tegund, scale = "free_y")

# ------------------------------------------------------------------------------
# Comparison with höski
attach("/net/hafkaldi/export/u2/reikn/R/SurveyWork/SMB/Allaggroldsmbindex.rdata")
vatican <-
  Allaggroldsmbindex %>%
  filter(species == Species,
         svaedi == "Heild",
         diurnal == 0,
         fixed == 0) %>%
  select(ar, lengd, b = bio.staerri, b.cv = cv.bio.staerri) %>%
  mutate(source = "Vatican") %>%
  as_tibble()
vatican <-
  vatican %>%
  filter(lengd == max(min(lengd), Length.min))
detach("file:/net/hafkaldi/export/u2/reikn/R/SurveyWork/SMB/Allaggroldsmbindex.rdata")

by.year %>%
  mutate(source = "tidy",
         ar = ar + 3/12) %>%
  bind_rows(vatican) %>%
  ggplot(aes(ar, b, colour = source)) +
  geom_point() +
  geom_pointrange(aes(ymin = b * (1 - b.cv),
                      ymax = b * (1 + b.cv))) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = c(0.25, 0.8)) +
  expand_limits(y = 0) +
  labs(x = NULL, y = NULL,
       title = "Biomass indices",
       subtitle = "Comparison of the orthodoxy (Vatican) and tidy")

