# ======================================================================================
#
# Mixed-effects logit model 
# Most probable mode switch procedure
#
# Mode use ~ f(travel impedance, individual's features, transport supply)
# School vs recreation 
#
# Young people - random effects 
# Travel impedance - travel distance OR travel time
# Individual's features - gender, age, settlement type
# Stacked set of modes per trip - random effects
#
# For reporting: Average Marginal Effects (AMEs) with 95% confidence intervals (CI)
# Additional: Odds Ratios (ORs) with 95% CI
#
# ======================================================================================

# ----------------------------
# Packages
# ----------------------------

pkgs <- c("aod", "dplyr", "patchwork", "MASS", "tidyr", "ggplot2", "car", "nnet", "lme4", "mclogit", "tidyverse", "marginaleffects", "blme", "posterior", "bayesplot", "marginaleffects")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(mclogit)
library(tidyverse)
library(marginaleffects)
library(posterior)
library(bayesplot)
library(nnet)
library(lme4)
library(blme)
library(car)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(scales)
library(patchwork)
library(aod)
library(purrr)


# ----------------------------
# Load data
# ----------------------------

url_csv <- "https://raw.githubusercontent.com/elaszkiewicz/research/main/database_stack.csv"

df <- read.csv(url_csv, stringsAsFactors = FALSE, sep = ";", dec = ".")

dim(dfa)


# ----------------------------------------------------------------------------------------------
# Recoding variables and NA's + assignation of mode-specific travel impedance measures
# ----------------------------------------------------------------------------------------------


modes <- c("TM_walk","TM_cycle","TM_rollerskates_skateboard","TM_scooter_other",
           "TM_publict","TM_car")


to_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", trimws(x))))

df3 <- df %>%
  mutate(
    across(c(D_walk, T_walk, D_cycle, T_cycle, D_publict, T_publict, D_car, T_car), to_num),

    TM_bin = ifelse(is.na(TM_bin), 0L, as.integer(TM_bin)),

    respondent = factor(respondent),
    trip_id    = factor(trip_id),
    trip_uid   = interaction(respondent, trip_id, drop = TRUE),

    # control extreme values 

    T_walk    = pmin(T_walk,    quantile(T_walk,    0.85, na.rm = TRUE)),
    D_walk    = pmin(D_walk,    quantile(D_walk,    0.85, na.rm = TRUE)),

    T_cycle   = pmin(T_cycle,   quantile(T_cycle,   0.85, na.rm = TRUE)),
    D_cycle   = pmin(D_cycle,   quantile(D_cycle,   0.85, na.rm = TRUE)),

    T_publict = pmin(T_publict, quantile(T_publict, 0.85, na.rm = TRUE)),
    D_publict = pmin(D_publict, quantile(D_publict, 0.85, na.rm = TRUE)),

    T_car     = pmin(T_car,     quantile(T_car,     0.85, na.rm = TRUE)),
    D_car     = pmin(D_car,     quantile(D_car,     0.85, na.rm = TRUE)),


    # replace NA for public transport by extreme values

    D_publict = ifelse(is.na(D_publict), 100000, as.integer(D_publict)),
    T_publict = ifelse(is.na(T_publict), 100000, as.integer(T_publict)),


    # recoding control variables

    gender = factor(gender),

    living_place = factor(living_place,
                          levels = c(1,2,3),
                          labels = c("Large city (Lodz)", "Medium-sized cities", "Small cities and rural areas")),

    age_group = case_when(
      age >= 7  & age <= 15 ~ "7–15",
      age >= 16 & age <= 18 ~ "16–18",
      TRUE ~ NA_character_
    ),
    age_group = factor(age_group, levels = c("7–15","16–18")),

    # spatial fixed effects for powiat

      spatial_fe = factor(ID_powiat),

    # trip purposes: school (utilitarian), recreation (outdoor)

	category1 = factor(
 	 ifelse(category == "S", "School", "Recreation"),
  			levels = c("School","Recreation")),

    # travel_mode labels (only if travel_mode is coded 1..6)

    travel_mode = factor(travel_mode, levels = 1:6, labels = modes),

    # assignation of time/distance for THIS row's observed mode

    T = case_when(
      travel_mode == "TM_walk"    ~ T_walk,
      travel_mode == "TM_cycle"   ~ T_cycle,
      travel_mode == "TM_rollerskates_skateboard"   ~ T_cycle,		
      travel_mode == "TM_publict" ~ T_publict,
      travel_mode == "TM_scooter_other" ~ T_car,				# assignation of car
      travel_mode == "TM_car"     ~ T_car,
      TRUE ~ NA_real_
    ),
    D = case_when(
      travel_mode == "TM_walk"    ~ D_walk,
      travel_mode == "TM_cycle"   ~ D_cycle,
      travel_mode == "TM_rollerskates_skateboard"   ~ D_cycle,		
      travel_mode == "TM_publict" ~ D_publict,
      travel_mode == "TM_scooter_other" ~ D_car,				# assignation of car
      travel_mode == "TM_car"     ~ D_car,
      TRUE ~ NA_real_
     ),

    # transport infrastructure measures - origin-level (scaled)

    TS_active_z0 = as.numeric(scale(SKRP_origin)),
    TS_PT_z0 = as.numeric(scale(KUKO_OIKM_origin)),
    TS_car_z0 = as.numeric(scale(SKJZ_origin)),

    # assignation of transport supply measures for THIS row's observed mode

    transport_supply0 = case_when(
      travel_mode == "TM_walk"    ~ TS_active_z0,
      travel_mode == "TM_cycle"   ~ TS_active_z0,
      travel_mode == "TM_rollerskates_skateboard"   ~ TS_active_z0,		
      travel_mode == "TM_publict" ~ TS_PT_z0,
      travel_mode == "TM_scooter_other" ~ TS_car_z0,			
      travel_mode == "TM_car"     ~ TS_car_z0,
      TRUE ~ NA_real_
     ),

    # transport infrastructure measures - NUTS4-level (scaled)

    TS_walk_z = as.numeric(scale(walk_lengper100km2_nuts4)),
    TS_active_z = as.numeric(scale(bike_lengper100km2_nuts4+walk_lengper100km2_nuts4)),
    TS_PT_z = as.numeric(scale(PT_countper100km2_nuts4)),
    TS_car_z = as.numeric(scale(road_lengper100km2_nuts4)),

    # assignation of transport supply measures for THIS row's observed mode

    transport_supply = case_when(
      travel_mode == "TM_walk"    ~ TS_walk_z,
      travel_mode == "TM_cycle"   ~ TS_active_z,
      travel_mode == "TM_rollerskates_skateboard"   ~ TS_active_z,		
      travel_mode == "TM_publict" ~ TS_PT_z,
      travel_mode == "TM_scooter_other" ~ TS_car_z,			
      travel_mode == "TM_car"     ~ TS_car_z,
      TRUE ~ NA_real_
     )
)

# remove missings

dim(df3)

df3 <- df3 %>%
  filter(
    !is.na(TM_bin),
    !is.na(travel_mode),
    !is.na(category),
    !is.na(respondent),
    !is.na(trip_uid),
    !is.na(gender),
    !is.na(age),
    !is.na(D_walk)
  )


dim(df3)



# preprocessing for stacked mode-use data : remove trips with no selected modes

df3 <- df3 %>%
  group_by(trip_uid) %>%
  mutate(
    n_modes = sum(TM_bin, na.rm = TRUE),
    multimodal = factor(as.integer(n_modes > 1))
  ) %>%
  ungroup() %>%
  dplyr::select(-n_modes)

valid_trips <- df3 %>%
  group_by(trip_id) %>%
  summarise(n_modes = sum(TM_bin)) %>%
  filter(n_modes > 0) %>%
  pull(trip_id)

df3 <- df3 %>% filter(trip_id %in% valid_trips)
df3 <- df3 %>% mutate(D_log = log1p(D))
df3 <- df3 %>% mutate(T_log = log1p(T))

dim(df3)


df3 <- df3 %>%
  dplyr::mutate(multimodal = factor(multimodal, levels = c(0,1),
                                    labels = c("Single mode","Multiple modes")))


# check data 

table(df3$travel_mode)
summary(df3$D_log)
summary(df3$D)
summary(df3$transport_supply)
summary(df3$transport_supply0)



# residualisation of transport infrastructure measure

res_TI <- lm(transport_supply ~ living_place, data = df3)
summary(res_TI)

# residuals

df3$transport_supply_resid <- residuals(res_TI)
df3$transport_supply_resid_z <- as.numeric(scale(df3$transport_supply_resid))

# correlation with original variable

cor(df3$transport_supply, df3$transport_supply_resid)





# ------------------------------------------------------------------------------------------------------------------
# RQ1: How do travel distance and time influence young people’s mode use for school and recreation trips?
# ------------------------------------------------------------------------------------------------------------------



# mixed-effect logit model



df3$travel_impedance <- df3$T_log		# change between travel distance and time


XA <- model.matrix(~ 0 + travel_mode +
         travel_mode:gender +
         travel_mode:age_group +
         travel_mode:travel_impedance +
         travel_mode:travel_impedance:living_place +
         travel_mode:category1 +
	   travel_mode:multimodal +
	   travel_mode:living_place +
         travel_mode:transport_supply_resid_z,
	   data = df3)

pA <- ncol(XA)

m_base <- bglmer(
  	   TM_bin ~ 0 + travel_mode +
         travel_mode:gender +
         travel_mode:age_group +
         travel_mode:travel_impedance+
         travel_mode:travel_impedance:living_place +
         travel_mode:category1 +
	   travel_mode:multimodal +
	   travel_mode:living_place +
         travel_mode:transport_supply_resid_z +
         (1|trip_id) +
         (1|respondent),
  data = df3, family = binomial(),
  fixef.prior = normal(cov = diag(2.5^2, pA)),
  control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(m_base)

################################################################################


### Wald tests for explanatory variables

b <- fixef(m_base)
V <- as.matrix(vcov(m_base))

aod::wald.test(b = b, Sigma = V, Terms = grep("gender", names(b)))
aod::wald.test(b = b, Sigma = V, Terms = grep("age_group", names(b)))
aod::wald.test(b = b, Sigma = V, Terms = grep("living_place", names(b)))
aod::wald.test(b = b, Sigma = V, Terms = grep("travel_impedance", names(b)))
aod::wald.test(b = b, Sigma = V, Terms = grep("travel_mode.*travel_impedance.*living_place", names(b)))
aod::wald.test(b = b, Sigma = V, Terms = grep("category1", names(b)))
aod::wald.test(b = b, Sigma = V, Terms = grep("multimodal", names(b)))
aod::wald.test(b = b, Sigma = V, Terms = grep("transport_supply_resid_z", names(b)))

### Odds ratios and 95% CI

library(dplyr)

b  <- fixef(m_base)
se <- sqrt(diag(vcov(m_base)))

results_table <- data.frame(
  Variable = names(b),
  OR = exp(b),
  OR_low = exp(b - 1.96*se),
  OR_high = exp(b + 1.96*se),
  p_value = 2*pnorm(abs(b/se), lower.tail = FALSE)
)

results_table


### Average Marginal Effects (AMEs) [Table 2]


ame_impedance <- avg_slopes(
  m_base,
  variables = "travel_impedance",
  by = c("travel_mode", "living_place"),

  type = "response"
)

ame_impedance2 <- avg_slopes(
  m_base,
  variables = "travel_impedance",
  by = c("travel_mode"),

  type = "response"
)


ame_ti <- avg_slopes(
  m_base,
  variables = "transport_supply_resid_z",
  by = "travel_mode",
  type = "response"
)


ame_gender <- avg_comparisons(
  m_base,
  variables = "gender",
  by = "travel_mode",
  type = "response"
)


ame_age <- avg_comparisons(
  m_base,
  variables = "age_group",
  by = "travel_mode",
  type = "response"
)


ame_purpose <- avg_comparisons(
  m_base,
  variables = "category1",
  by = "travel_mode",
  type = "response"
)


ame_multi <- avg_comparisons(
  m_base,
  variables = "multimodal",
  by = "travel_mode",
  type = "response"
)

ame_place <- avg_comparisons(
  m_base,
  variables = "living_place",
  by = "travel_mode",
  type = "response"
)


ame_all <- bind_rows(
  Gender = ame_gender,
  Age = ame_age,
  Settlement = ame_place,
  Purpose = ame_purpose,
  Multimodal = ame_multi,
  TravelImpedance = ame_impedance,
  TravelImpedanceLivingPlace = ame_impedance2,
  TransportInfrastructure = ame_ti,
  .id = "Variable"
)

write.csv(ame_all, "average_marginal_effects_D.csv", row.names = FALSE)		# travel distance or time


# ------------------------------------------------------------------------------------------------------------------------------------
# RQ2: At what distance and time thresholds does the most probable mode shift, and how do these thresholds vary by trip purpose?
# ------------------------------------------------------------------------------------------------------------------------------------


# -------------------------------
# Travel distance
# -------------------------------

# -------------------------------
# Travel distance
# -------------------------------

df3$travel_impedance <- df3$D_log

XA <- model.matrix(~ 0 + travel_mode +
  travel_mode:gender +
  travel_mode:age_group +
  travel_mode:travel_impedance +
  travel_mode:travel_impedance:living_place +
  travel_mode:category1 +
  travel_mode:multimodal +
  travel_mode:living_place +
  travel_mode:transport_supply_resid_z,
  data = df3
)

pA <- ncol(XA)

m_base <- bglmer(
  TM_bin ~ 0 + travel_mode +
    travel_mode:gender +
    travel_mode:age_group +
    travel_mode:travel_impedance +
    travel_mode:travel_impedance:living_place +
    travel_mode:category1 +
    travel_mode:multimodal +
    travel_mode:living_place +
    travel_mode:transport_supply_resid_z +
    (1 | trip_id) +
    (1 | respondent),
  data = df3,
  family = binomial(),
  fixef.prior = normal(cov = diag(2.5^2, pA)),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(m_base)

# 0) Categories + distance grid

cats_to_plot <- intersect(c("School", "Recreation"), levels(df3$category1))
stopifnot(length(cats_to_plot) == 2)

D_seq <- seq(0, 10000, length.out = 160)

# 1) Prediction grid

newdat <- expand.grid(
  travel_mode  = levels(df3$travel_mode),
  category1    = cats_to_plot,
  multimodal   = levels(df3$multimodal),
  gender       = levels(df3$gender),
  age_group    = levels(df3$age_group),
  living_place = levels(df3$living_place),
  D            = D_seq
) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    travel_mode  = factor(travel_mode,  levels = levels(df3$travel_mode)),
    category1    = factor(category1,    levels = levels(df3$category1)),
    multimodal   = factor(multimodal,   levels = levels(df3$multimodal)),
    gender       = factor(gender,       levels = levels(df3$gender)),
    age_group    = factor(age_group,    levels = levels(df3$age_group)),
    living_place = factor(living_place, levels = levels(df3$living_place)),

    travel_impedance = log1p(D),
    transport_supply_resid_z = 0,

    trip_id    = factor(df3$trip_id[1], levels = levels(df3$trip_id)),
    respondent = factor(df3$respondent[1], levels = levels(df3$respondent))
  )

# 2) Predict: fixed effects only

pred <- newdat %>%
  dplyr::mutate(
    eta  = as.numeric(predict(m_base, newdata = newdat, re.form = NA)),
    prob = plogis(eta)
  )

# 3) Wald CI for eta -> probability

beta <- fixef(m_base)
V    <- as.matrix(vcov(m_base))

form_fixed <- lme4::nobars(formula(m_base))
tt <- delete.response(terms(form_fixed))

X <- model.matrix(tt, pred)
X <- X[, names(beta), drop = FALSE]

se_eta <- sqrt(diag(X %*% V %*% t(X)))

pred <- pred %>%
  dplyr::mutate(
    lo = plogis(eta - 1.96 * se_eta),
    hi = plogis(eta + 1.96 * se_eta)
  )

# 4) Average over gender + age_group + multimodal
#    Keep living_place because it interacts with travel_impedance

pred_avg <- pred %>%
  dplyr::group_by(category1, living_place, travel_mode, D) %>%
  dplyr::summarise(
    eta  = mean(eta,  na.rm = TRUE),
    prob = mean(prob, na.rm = TRUE),
    lo   = mean(lo,   na.rm = TRUE),
    hi   = mean(hi,   na.rm = TRUE),
    .groups = "drop"
  )

# 5) Settlement-specific mode switches

best_path <- pred_avg %>%
  dplyr::group_by(category1, living_place, D) %>%
  dplyr::slice_max(eta, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(category1, living_place, D) %>%
  dplyr::group_by(category1, living_place) %>%
  dplyr::mutate(
    prev_mode = dplyr::lag(travel_mode),
    prev_D    = dplyr::lag(D)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(prev_mode), travel_mode != prev_mode) %>%
  dplyr::rename(to_mode = travel_mode, D2 = D) %>%
  dplyr::mutate(from_mode = prev_mode)

eta_lookup <- pred_avg %>%
  dplyr::select(category1, living_place, D, travel_mode, eta, prob)

switches_best <- best_path %>%
  dplyr::transmute(category1, living_place, prev_D, D2, from_mode, to_mode) %>%
  dplyr::left_join(
    eta_lookup %>%
      dplyr::rename(
        prev_D = D,
        from_mode = travel_mode,
        eta_from_D1 = eta,
        p_from_D1 = prob
      ),
    by = c("category1", "living_place", "prev_D", "from_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup %>%
      dplyr::rename(
        D2 = D,
        from_mode = travel_mode,
        eta_from_D2 = eta,
        p_from_D2 = prob
      ),
    by = c("category1", "living_place", "D2", "from_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup %>%
      dplyr::rename(
        prev_D = D,
        to_mode = travel_mode,
        eta_to_D1 = eta
      ),
    by = c("category1", "living_place", "prev_D", "to_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup %>%
      dplyr::rename(
        D2 = D,
        to_mode = travel_mode,
        eta_to_D2 = eta
      ),
    by = c("category1", "living_place", "D2", "to_mode")
  ) %>%
  dplyr::mutate(
    d1 = eta_from_D1 - eta_to_D1,
    d2 = eta_from_D2 - eta_to_D2,
    denom = d2 - d1,
    D_int = prev_D + (0 - d1) * (D2 - prev_D) / denom,
    p_int = p_from_D1 + (D_int - prev_D) * (p_from_D2 - p_from_D1) / (D2 - prev_D),
    label = paste0(round(D_int / 1000, 2), " km")
  ) %>%
  dplyr::filter(
    is.finite(D_int),
    abs(denom) > 1e-10,
    is.finite(p_int),
    D_int >= min(pred_avg$D),
    D_int <= max(pred_avg$D)
  ) %>%
  dplyr::select(category1, living_place, from_mode, to_mode, D_int, p_int, label)

print(switches_best)

# 6) Average over settlement type

pred_avg_settlement <- pred %>%
  dplyr::group_by(category1, travel_mode, D) %>%
  dplyr::summarise(
    eta  = mean(eta,  na.rm = TRUE),
    prob = mean(prob, na.rm = TRUE),
    lo   = mean(lo,   na.rm = TRUE),
    hi   = mean(hi,   na.rm = TRUE),
    living_place = "Averaged over settlement type",
    .groups = "drop"
  )

best_path_settlement <- pred_avg_settlement %>%
  dplyr::group_by(category1, D) %>%
  dplyr::slice_max(eta, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(category1, D) %>%
  dplyr::group_by(category1) %>%
  dplyr::mutate(
    prev_mode = dplyr::lag(travel_mode),
    prev_D    = dplyr::lag(D)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(prev_mode), travel_mode != prev_mode) %>%
  dplyr::rename(to_mode = travel_mode, D2 = D) %>%
  dplyr::mutate(from_mode = prev_mode)

eta_lookup_settlement <- pred_avg_settlement %>%
  dplyr::select(category1, D, travel_mode, eta, prob)

switches_best_settlement <- best_path_settlement %>%
  dplyr::transmute(category1, prev_D, D2, from_mode, to_mode) %>%
  dplyr::left_join(
    eta_lookup_settlement %>%
      dplyr::rename(
        prev_D = D,
        from_mode = travel_mode,
        eta_from_D1 = eta,
        p_from_D1 = prob
      ),
    by = c("category1", "prev_D", "from_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup_settlement %>%
      dplyr::rename(
        D2 = D,
        from_mode = travel_mode,
        eta_from_D2 = eta,
        p_from_D2 = prob
      ),
    by = c("category1", "D2", "from_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup_settlement %>%
      dplyr::rename(
        prev_D = D,
        to_mode = travel_mode,
        eta_to_D1 = eta
      ),
    by = c("category1", "prev_D", "to_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup_settlement %>%
      dplyr::rename(
        D2 = D,
        to_mode = travel_mode,
        eta_to_D2 = eta
      ),
    by = c("category1", "D2", "to_mode")
  ) %>%
  dplyr::mutate(
    d1 = eta_from_D1 - eta_to_D1,
    d2 = eta_from_D2 - eta_to_D2,
    denom = d2 - d1,
    D_int = prev_D + (0 - d1) * (D2 - prev_D) / denom,
    p_int = p_from_D1 + (D_int - prev_D) * (p_from_D2 - p_from_D1) / (D2 - prev_D),
    label = paste0(round(D_int / 1000, 2), " km")
  ) %>%
  dplyr::filter(
    is.finite(D_int),
    abs(denom) > 1e-10,
    is.finite(p_int),
    D_int >= min(pred_avg_settlement$D),
    D_int <= max(pred_avg_settlement$D)
  ) %>%
  dplyr::select(category1, from_mode, to_mode, D_int, p_int, label)

print(switches_best_settlement)

# 7) Label adjustment

mode_labs <- c(
  "TM_walk" = "Walking",
  "TM_cycle" = "Cycling",
  "TM_rollerskates_skateboard" = "Rollerskates and skateboard",
  "TM_scooter_other" = "Other modes",
  "TM_publict" = "Public transport",
  "TM_car" = "Car"
)

th_D  <- 500
th_p  <- 0.04
delta <- 0.06

switches_adj <- switches_best %>%
  dplyr::group_by(category1, living_place) %>%
  dplyr::arrange(D_int, .by_group = TRUE) %>%
  dplyr::mutate(
    new_cluster = dplyr::if_else(
      dplyr::row_number() == 1,
      TRUE,
      (abs(D_int - dplyr::lag(D_int)) > th_D) |
        (abs(p_int - dplyr::lag(p_int)) > th_p)
    ),
    cluster = cumsum(new_cluster),
    k = dplyr::row_number(),
    k_in_cluster = ave(k, cluster, FUN = seq_along),
    sign = ifelse(k_in_cluster %% 2 == 1, 1, -1),
    p_lab = pmin(0.98, pmax(0.02, p_int + sign * delta))
  ) %>%
  dplyr::ungroup()

switches_adj_settlement <- switches_best_settlement %>%
  dplyr::mutate(
    living_place = "Averaged over settlement type"
  ) %>%
  dplyr::group_by(category1, living_place) %>%
  dplyr::arrange(D_int, .by_group = TRUE) %>%
  dplyr::mutate(
    new_cluster = dplyr::if_else(
      dplyr::row_number() == 1,
      TRUE,
      (abs(D_int - dplyr::lag(D_int)) > th_D) |
        (abs(p_int - dplyr::lag(p_int)) > th_p)
    ),
    cluster = cumsum(new_cluster),
    k = dplyr::row_number(),
    k_in_cluster = ave(k, cluster, FUN = seq_along),
    sign = ifelse(k_in_cluster %% 2 == 1, 1, -1),
    p_lab = pmin(0.98, pmax(0.02, p_int + sign * delta))
  ) %>%
  dplyr::ungroup()

# 8) Combine settlement-specific and settlement-averaged results

pred_all <- dplyr::bind_rows(
  pred_avg,
  pred_avg_settlement
)

switches_all <- dplyr::bind_rows(
  switches_adj,
  switches_adj_settlement
)

pred_all$living_place <- factor(
  pred_all$living_place,
  levels = c(
    levels(df3$living_place),
    "Averaged over settlement type"
  )
)

switches_all$living_place <- factor(
  switches_all$living_place,
  levels = levels(pred_all$living_place)
)

# 9) Final combined plot

p_combined_distance <- ggplot(
  pred_all,
  aes(x = D, y = prob, color = travel_mode, fill = travel_mode)
) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  facet_grid(
    living_place ~ category1,
    drop = FALSE
  ) +
  geom_segment(
    data = switches_all,
    aes(
      x = D_int,
      xend = D_int,
      y = p_lab,
      yend = p_int
    ),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  geom_label(
    data = switches_all,
    aes(
      x = D_int,
      y = p_lab,
      label = label
    ),
    inherit.aes = FALSE,
    label.size = 0.2,
    alpha = 0.85,
    size = 3
  ) +
  scale_color_discrete(
    name = "Travel mode:",
    labels = mode_labs
  ) +
  scale_fill_discrete(
    name = "Travel mode:",
    labels = mode_labs
  ) +
  labs(
    x = "Travel distance [m]",
    y = "Predicted probability of mode use"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

p_combined_distance

ggsave(
  "most probable mode switch distance.png",
  p_combined_distance,
  width = 9.5,
  height = 8.5,
  dpi = 600
)


# --------------------------------
# Travel time
# --------------------------------


# Travel time plot [Figure 2]

cats_to_plot <- intersect(c("School", "Recreation"), levels(df3$category1))
stopifnot(length(cats_to_plot) == 2)

T_seq <- seq(0, 75, length.out = 160)

# 1) Prediction grid

newdat <- expand.grid(
  travel_mode  = levels(df3$travel_mode),
  category1    = cats_to_plot,
  multimodal   = levels(df3$multimodal),
  gender       = levels(df3$gender),
  age_group    = levels(df3$age_group),
  living_place = levels(df3$living_place),
  T            = T_seq
) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    travel_mode  = factor(travel_mode,  levels = levels(df3$travel_mode)),
    category1    = factor(category1,    levels = levels(df3$category1)),
    multimodal   = factor(multimodal,   levels = levels(df3$multimodal)),
    gender       = factor(gender,       levels = levels(df3$gender)),
    age_group    = factor(age_group,    levels = levels(df3$age_group)),
    living_place = factor(living_place, levels = levels(df3$living_place)),

    travel_impedance = log1p(T),
    transport_supply_resid_z = 0,

    trip_id    = factor(df3$trip_id[1], levels = levels(df3$trip_id)),
    respondent = factor(df3$respondent[1], levels = levels(df3$respondent))
  )

# 2) Predict: fixed effects only

pred <- newdat %>%
  dplyr::mutate(
    eta  = as.numeric(predict(m_base, newdata = newdat, re.form = NA)),
    prob = plogis(eta)
  )

# 3) Wald CI for eta -> probability

beta <- fixef(m_base)
V    <- as.matrix(vcov(m_base))

form_fixed <- lme4::nobars(formula(m_base))
tt <- delete.response(terms(form_fixed))

X <- model.matrix(tt, pred)
X <- X[, names(beta), drop = FALSE]

se_eta <- sqrt(diag(X %*% V %*% t(X)))

pred <- pred %>%
  dplyr::mutate(
    lo = plogis(eta - 1.96 * se_eta),
    hi = plogis(eta + 1.96 * se_eta)
  )

# 4) Average over gender + age_group + multimodal
#    Keep living_place because it interacts with travel_impedance

pred_avg <- pred %>%
  dplyr::group_by(category1, living_place, travel_mode, T) %>%
  dplyr::summarise(
    eta  = mean(eta,  na.rm = TRUE),
    prob = mean(prob, na.rm = TRUE),
    lo   = mean(lo,   na.rm = TRUE),
    hi   = mean(hi,   na.rm = TRUE),
    .groups = "drop"
  )

# 5) Settlement-specific mode switches

best_path <- pred_avg %>%
  dplyr::group_by(category1, living_place, T) %>%
  dplyr::slice_max(eta, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(category1, living_place, T) %>%
  dplyr::group_by(category1, living_place) %>%
  dplyr::mutate(
    prev_mode = dplyr::lag(travel_mode),
    prev_T    = dplyr::lag(T)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(prev_mode), travel_mode != prev_mode) %>%
  dplyr::rename(to_mode = travel_mode, T2 = T) %>%
  dplyr::mutate(from_mode = prev_mode)

eta_lookup <- pred_avg %>%
  dplyr::select(category1, living_place, T, travel_mode, eta, prob)

switches_best <- best_path %>%
  dplyr::transmute(category1, living_place, prev_T, T2, from_mode, to_mode) %>%
  dplyr::left_join(
    eta_lookup %>%
      dplyr::rename(
        prev_T = T,
        from_mode = travel_mode,
        eta_from_T1 = eta,
        p_from_T1 = prob
      ),
    by = c("category1", "living_place", "prev_T", "from_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup %>%
      dplyr::rename(
        T2 = T,
        from_mode = travel_mode,
        eta_from_T2 = eta,
        p_from_T2 = prob
      ),
    by = c("category1", "living_place", "T2", "from_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup %>%
      dplyr::rename(
        prev_T = T,
        to_mode = travel_mode,
        eta_to_T1 = eta
      ),
    by = c("category1", "living_place", "prev_T", "to_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup %>%
      dplyr::rename(
        T2 = T,
        to_mode = travel_mode,
        eta_to_T2 = eta
      ),
    by = c("category1", "living_place", "T2", "to_mode")
  ) %>%
  dplyr::mutate(
    d1 = eta_from_T1 - eta_to_T1,
    d2 = eta_from_T2 - eta_to_T2,
    denom = d2 - d1,
    T_int = prev_T + (0 - d1) * (T2 - prev_T) / denom,
    p_int = p_from_T1 + (T_int - prev_T) * (p_from_T2 - p_from_T1) / (T2 - prev_T),
    label = paste0(round(T_int, 1), " min")
  ) %>%
  dplyr::filter(
    is.finite(T_int),
    abs(denom) > 1e-10,
    is.finite(p_int),
    T_int >= min(pred_avg$T),
    T_int <= max(pred_avg$T)
  ) %>%
  dplyr::select(category1, living_place, from_mode, to_mode, T_int, p_int, label)

print(switches_best)

# 6) Average over settlement type

pred_avg_settlement <- pred %>%
  dplyr::group_by(category1, travel_mode, T) %>%
  dplyr::summarise(
    eta  = mean(eta,  na.rm = TRUE),
    prob = mean(prob, na.rm = TRUE),
    lo   = mean(lo,   na.rm = TRUE),
    hi   = mean(hi,   na.rm = TRUE),
    living_place = "Averaged over settlement type",
    .groups = "drop"
  )

best_path_settlement <- pred_avg_settlement %>%
  dplyr::group_by(category1, T) %>%
  dplyr::slice_max(eta, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(category1, T) %>%
  dplyr::group_by(category1) %>%
  dplyr::mutate(
    prev_mode = dplyr::lag(travel_mode),
    prev_T    = dplyr::lag(T)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(prev_mode), travel_mode != prev_mode) %>%
  dplyr::rename(to_mode = travel_mode, T2 = T) %>%
  dplyr::mutate(from_mode = prev_mode)

eta_lookup_settlement <- pred_avg_settlement %>%
  dplyr::select(category1, T, travel_mode, eta, prob)

switches_best_settlement <- best_path_settlement %>%
  dplyr::transmute(category1, prev_T, T2, from_mode, to_mode) %>%
  dplyr::left_join(
    eta_lookup_settlement %>%
      dplyr::rename(
        prev_T = T,
        from_mode = travel_mode,
        eta_from_T1 = eta,
        p_from_T1 = prob
      ),
    by = c("category1", "prev_T", "from_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup_settlement %>%
      dplyr::rename(
        T2 = T,
        from_mode = travel_mode,
        eta_from_T2 = eta,
        p_from_T2 = prob
      ),
    by = c("category1", "T2", "from_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup_settlement %>%
      dplyr::rename(
        prev_T = T,
        to_mode = travel_mode,
        eta_to_T1 = eta
      ),
    by = c("category1", "prev_T", "to_mode")
  ) %>%
  dplyr::left_join(
    eta_lookup_settlement %>%
      dplyr::rename(
        T2 = T,
        to_mode = travel_mode,
        eta_to_T2 = eta
      ),
    by = c("category1", "T2", "to_mode")
  ) %>%
  dplyr::mutate(
    d1 = eta_from_T1 - eta_to_T1,
    d2 = eta_from_T2 - eta_to_T2,
    denom = d2 - d1,
    T_int = prev_T + (0 - d1) * (T2 - prev_T) / denom,
    p_int = p_from_T1 + (T_int - prev_T) * (p_from_T2 - p_from_T1) / (T2 - prev_T),
    label = paste0(round(T_int, 1), " min")
  ) %>%
  dplyr::filter(
    is.finite(T_int),
    abs(denom) > 1e-10,
    is.finite(p_int),
    T_int >= min(pred_avg_settlement$T),
    T_int <= max(pred_avg_settlement$T)
  ) %>%
  dplyr::select(category1, from_mode, to_mode, T_int, p_int, label)

print(switches_best_settlement)

# 7) Label adjustment

mode_labs <- c(
  "TM_walk" = "Walking",
  "TM_cycle" = "Cycling",
  "TM_rollerskates_skateboard" = "Rollerskates and skateboard",
  "TM_scooter_other" = "Other modes",
  "TM_publict" = "Public transport",
  "TM_car" = "Car"
)

th_T  <- 5
th_p  <- 0.04
delta <- 0.06

switches_adj <- switches_best %>%
  dplyr::group_by(category1, living_place) %>%
  dplyr::arrange(T_int, .by_group = TRUE) %>%
  dplyr::mutate(
    new_cluster = dplyr::if_else(
      dplyr::row_number() == 1,
      TRUE,
      (abs(T_int - dplyr::lag(T_int)) > th_T) |
        (abs(p_int - dplyr::lag(p_int)) > th_p)
    ),
    cluster = cumsum(new_cluster),
    k = dplyr::row_number(),
    k_in_cluster = ave(k, cluster, FUN = seq_along),
    sign = ifelse(k_in_cluster %% 2 == 1, 1, -1),
    p_lab = pmin(0.98, pmax(0.02, p_int + sign * delta))
  ) %>%
  dplyr::ungroup()

switches_adj_settlement <- switches_best_settlement %>%
  dplyr::mutate(
    living_place = "Averaged over settlement type"
  ) %>%
  dplyr::group_by(category1, living_place) %>%
  dplyr::arrange(T_int, .by_group = TRUE) %>%
  dplyr::mutate(
    new_cluster = dplyr::if_else(
      dplyr::row_number() == 1,
      TRUE,
      (abs(T_int - dplyr::lag(T_int)) > th_T) |
        (abs(p_int - dplyr::lag(p_int)) > th_p)
    ),
    cluster = cumsum(new_cluster),
    k = dplyr::row_number(),
    k_in_cluster = ave(k, cluster, FUN = seq_along),
    sign = ifelse(k_in_cluster %% 2 == 1, 1, -1),
    p_lab = pmin(0.98, pmax(0.02, p_int + sign * delta))
  ) %>%
  dplyr::ungroup()

# 8) Combine settlement-specific and settlement-averaged results

pred_all <- dplyr::bind_rows(
  pred_avg,
  pred_avg_settlement
)

switches_all <- dplyr::bind_rows(
  switches_adj,
  switches_adj_settlement
)

pred_all$living_place <- factor(
  pred_all$living_place,
  levels = c(
    levels(df3$living_place),
    "Averaged over settlement type"
  )
)

switches_all$living_place <- factor(
  switches_all$living_place,
  levels = levels(pred_all$living_place)
)

# 9) Final combined plot

p_combined <- ggplot(
  pred_all,
  aes(x = T, y = prob, color = travel_mode, fill = travel_mode)
) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  facet_grid(
    living_place ~ category1,
    drop = FALSE
  ) +
  geom_segment(
    data = switches_all,
    aes(
      x = T_int,
      xend = T_int,
      y = p_lab,
      yend = p_int
    ),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  geom_label(
    data = switches_all,
    aes(
      x = T_int,
      y = p_lab,
      label = label
    ),
    inherit.aes = FALSE,
    label.size = 0.2,
    alpha = 0.85,
    size = 3
  ) +
  scale_color_discrete(
    name = "Travel mode:",
    labels = mode_labs
  ) +
  scale_fill_discrete(
    name = "Travel mode:",
    labels = mode_labs
  ) +
  labs(
    x = "Travel time [min]",
    y = "Predicted probability of mode use"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

p_combined


ggsave("most probable mode shift time.png", p_combined, width = 9.5, height = 5.5, dpi = 800)


# ---------------------------------------------------------
# Most probable mode switches for any grouping variables
# ---------------------------------------------------------

best_switch_table <- function(pred_avg, group_vars, x_var = "T") {

  join_vars <- c(group_vars, x_var)
  mode_var <- "travel_mode"

  best_path <- pred_avg %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(join_vars))) %>%
    dplyr::slice_max(eta, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(c(group_vars, x_var)))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::mutate(
      prev_mode = dplyr::lag(.data[[mode_var]]),
      prev_x    = dplyr::lag(.data[[x_var]])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(prev_mode), .data[[mode_var]] != prev_mode) %>%
    dplyr::rename(
      to_mode = travel_mode,
      x2 = dplyr::all_of(x_var)
    ) %>%
    dplyr::mutate(from_mode = prev_mode)

  eta_lookup <- pred_avg %>%
    dplyr::select(
      dplyr::all_of(group_vars),
      dplyr::all_of(x_var),
      travel_mode,
      eta,
      prob
    )

  switches <- best_path %>%
    dplyr::transmute(
      dplyr::across(dplyr::all_of(group_vars)),
      prev_x,
      x2,
      from_mode,
      to_mode
    ) %>%
    dplyr::left_join(
      eta_lookup %>%
        dplyr::rename(
          prev_x = dplyr::all_of(x_var),
          from_mode = travel_mode,
          eta_from_x1 = eta,
          p_from_x1 = prob
        ),
      by = c(group_vars, "prev_x", "from_mode")
    ) %>%
    dplyr::left_join(
      eta_lookup %>%
        dplyr::rename(
          x2 = dplyr::all_of(x_var),
          from_mode = travel_mode,
          eta_from_x2 = eta,
          p_from_x2 = prob
        ),
      by = c(group_vars, "x2", "from_mode")
    ) %>%
    dplyr::left_join(
      eta_lookup %>%
        dplyr::rename(
          prev_x = dplyr::all_of(x_var),
          to_mode = travel_mode,
          eta_to_x1 = eta
        ),
      by = c(group_vars, "prev_x", "to_mode")
    ) %>%
    dplyr::left_join(
      eta_lookup %>%
        dplyr::rename(
          x2 = dplyr::all_of(x_var),
          to_mode = travel_mode,
          eta_to_x2 = eta
        ),
      by = c(group_vars, "x2", "to_mode")
    ) %>%
    dplyr::mutate(
      d1 = eta_from_x1 - eta_to_x1,
      d2 = eta_from_x2 - eta_to_x2,
      denom = d2 - d1,
      x_int = prev_x + (0 - d1) * (x2 - prev_x) / denom,
      p_int = p_from_x1 + (x_int - prev_x) * (p_from_x2 - p_from_x1) / (x2 - prev_x),
      label = paste0(round(x_int / 1000, 2))
    ) %>%
    dplyr::filter(
      is.finite(x_int),
      abs(denom) > 1e-10,
      is.finite(p_int),
      x_int >= min(pred_avg[[x_var]], na.rm = TRUE),
      x_int <= max(pred_avg[[x_var]], na.rm = TRUE)
    ) %>%
    dplyr::select(
      dplyr::all_of(group_vars),
      from_mode,
      to_mode,
      x_int,
      p_int,
      label
    )

  switches
}

# ----------------------------------------------------------------------------------------
# Predictions averaged for each subgroup -- either distance or time, depending on m_base
# ----------------------------------------------------------------------------------------

# 1) Switches within category1 x gender
pred_gender <- pred %>%
  dplyr::group_by(category1, gender, travel_mode, T) %>%
  dplyr::summarise(
    eta  = mean(eta,  na.rm = TRUE),
    prob = mean(prob, na.rm = TRUE),
    .groups = "drop"
  )

switch_gender <- best_switch_table(
  pred_avg = pred_gender,
  group_vars = c("category1",  "gender"),
  x_var = "T"
)

print(switch_gender)


# 2) Switches within category1 x age_group
pred_age <- pred %>%
  dplyr::group_by(category1, age_group, travel_mode, T) %>%
  dplyr::summarise(
    eta  = mean(eta,  na.rm = TRUE),
    prob = mean(prob, na.rm = TRUE),
    .groups = "drop"
  )

switch_age <- best_switch_table(
  pred_avg = pred_age,
  group_vars = c("category1", "age_group"),
  x_var = "T"
)

print(switch_age)


# 3) Switches within category1 x settlement type
pred_place <- pred %>%
  dplyr::group_by(category1, living_place, travel_mode, T) %>%
  dplyr::summarise(
    eta  = mean(eta,  na.rm = TRUE),
    prob = mean(prob, na.rm = TRUE),
    .groups = "drop"
  )

switch_place <- best_switch_table(
  pred_avg = pred_place,
  group_vars = c("category1", "living_place"),
  x_var = "T"
)

print(switch_place)



