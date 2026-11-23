# ======================================================================================
#
# Mixed-effects logit models of young people's mode use (stacked mode-level data)
# Probability crossover points in predicted mode-use probabilities
#
#======================================================================================
#
# Mode use ~ f(travel impedance, young people's characteristics, settlement type,
#              transport infrastructure, trip purpose)
#
# Random intercepts: respondent, origin-destination pair (trip_id)
# Travel impedance: travel time OR travel distance 
# Trip purpose: school, local recreation (OR_N1, OR_N2), citywide recreation (OR_C1, OR_C2)
#
#
# Comments/assumptions:
#
# 1) Travel impedance specified as a natural cubic spline of log(1 + x) (3 df, interior knots
#    at the terciles); log-linear and linear specifications are compared in 5).
# 2) Origin-destination pairs without a public transport connection are handled with an indicator.
# 3) Probability crossover points between pairs of modes. 
#    Curves are averaged on the probability scale over the regional population composition.
# 4) Simulation-based 95% intervals for crossovers, pointwise bands for the curves,
#    data-support checks (observations near each crossover). 
# 5) Linear / log / spline impedance specifications compared with AIC, BIC, CV Brier score and log-loss
# 6) Trip purpose x impedance interaction: coefficients with intervals, joint Wald test,
#    second differences in AMEs, AIC, BIC, CV Brier score and log-loss.
#
#
# ======================================================================================


# --------------------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------------------

QUICK_TEST <- FALSE      # TRUE = fast run to check that the code works 

SET <- list(
  data_url   = "https://raw.githubusercontent.com/elaszkiewicz/research/main/database_stack.csv",
  out_dir    = "revision_categories3new",
  cap_q      = 0.85,                         											# winsorisation quantile 
  pt_missing = "indicator",  
  axis       = "od",         
  n_sim      = 1000,                         											# coefficient draws for crossover intervals
  run_cv     = TRUE,                        											# cross-validated diagnostics 
  k_folds    = 5,
  run_boot   = FALSE,                       										      # respondent-level cluster bootstrap
  n_boot     = 300,
  sesoi      = 0.05,                        										      # smallest effect size of interest
  seed       = 12345,
  fast_all   = FALSE,       			   											# TRUE only for testing: no Hessian anywhere (SEs approximate)
  weighting  = "population",
  three_purposes = TRUE,     # TRUE: school / local / citywide recreation; FALSE: school / recreation
  population_file = "https://raw.githubusercontent.com/elaszkiewicz/research/main/population_shares.csv"    # Regional population shares by gender x age group x settlement type
)

if (QUICK_TEST) {
  SET$n_sim <- 100
  SET$k_folds <- 2
  SET$n_boot <- 2
  SET$fast_all <- TRUE
}

MEASURES <- list(
  time = list(
    name = "time", var = "T", od_var = "T_walk",
    unit = 10, unit_label = "+10 min",
    grid = seq(1, 120, length.out = 160), band = 5,   # 0-75 min
    fmt = function(x) sprintf("%.0f min", x),
    axis_label = c(od   = "Walking time between origin and destination [min]",
                   same = "Travel time by each mode [min]")
  ),
  distance = list(
    name = "distance", var = "D", od_var = "D_walk",
    unit = 1000, unit_label = "+1 km",
    grid = seq(50, 11500, length.out = 160), band = 500,   # 0-10 km
    fmt = function(x) sprintf("%.1f km", x / 1000),
    axis_label = c(od   = "Walking distance between origin and destination [m]",
                   same = "Travel distance by each mode [m]")
  )
)

MODES <- c("TM_walk", "TM_cycle", "TM_rollerskates_skateboard", "TM_scooter_other",
           "TM_publict", "TM_car")

PURPOSES <- if (SET$three_purposes) c("School", "Local recreation", "Citywide recreation") else
  c("School", "Recreation")

MODE_LABS <- c(
  TM_walk = "Walking",
  TM_cycle = "Cycling",
  TM_rollerskates_skateboard = "Roller skates & skateboard",
  TM_scooter_other = "Other modes",
  TM_publict = "Public transport",
  TM_car = "Car"
)

# Routed Google Maps impedance assigned to each model category

IMP_SOURCE <- c(
  TM_walk = "walk",
  TM_cycle = "cycle",
  TM_rollerskates_skateboard = "cycle",
  TM_scooter_other = "car",
  TM_publict = "publict",
  TM_car = "car"
)

# Pairs of modes for which crossovers are reported

PAIRS <- list(
  c("TM_walk", "TM_car"),
  c("TM_walk", "TM_publict"),
  c("TM_walk", "TM_cycle"),
  c("TM_cycle", "TM_car"),
  c("TM_cycle", "TM_publict"),
  c("TM_publict", "TM_car")
)


# --------------------------------------------------------------------------------------
# Packages
# --------------------------------------------------------------------------------------

pkgs <- c("lme4", "blme", "marginaleffects", "dplyr", "tidyr", "ggplot2", "MASS")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(lme4)
  library(blme)
  library(marginaleffects)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(splines)
})

dir.create(SET$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(name) file.path(SET$out_dir, name)
save_csv <- function(x, name) write.csv(x, out(name), row.names = FALSE)


# --------------------------------------------------------------------------------------
# Data preparation 
# --------------------------------------------------------------------------------------

prepare_data <- function(raw, cap_q = SET$cap_q, pt_missing = SET$pt_missing) {

  imp_vars <- c("D_walk", "T_walk", "D_cycle", "T_cycle",
                "D_publict", "T_publict", "D_car", "T_car")
  to_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(x)))))

  d <- raw %>% mutate(across(all_of(imp_vars), to_num))

  # Winsorisation at cap_q
  # Origin-destination pairs without a public transport connection are NA and stay NA

  caps <- vapply(imp_vars, function(v) unname(quantile(d[[v]], cap_q, na.rm = TRUE)), numeric(1))

  d <- d %>%
    mutate(
      pt_none = is.na(D_publict) | is.na(T_publict),
      across(all_of(imp_vars), ~ pmin(.x, caps[cur_column()]))
    )

  if (pt_missing == "sentinel") {
    d <- d %>%
      mutate(
        D_publict = ifelse(is.na(D_publict), 1e5, D_publict),
        T_publict = ifelse(is.na(T_publict), 1e5, T_publict)
      )
  }

  d <- d %>%
    mutate(
      TM_bin = ifelse(is.na(TM_bin), 0L, as.integer(TM_bin)),
      respondent = factor(respondent),
      trip_id = factor(trip_id),
      gender = factor(gender),   										# 1 = girls / 0 = boys 
      living_place = factor(living_place, levels = 1:3,
                            labels = c("Large city (Lodz)", "Medium-sized cities",
                                       "Small towns and rural areas")),
      age_group = factor(case_when(age >= 7 & age <= 15 ~ "7–15",
                                   age >= 16 & age <= 18 ~ "16–18"),
                         levels = c("7–15", "16–18")),


      # Trip purpose: school; local recreation = (OR_N1, OR_N2); citywide recreation = (OR_C1, OR_C2)

      category1 = factor(case_when(category == "S" ~ "School",
                                   !SET$three_purposes & grepl("^OR_", category) ~ "Recreation",
                                   grepl("^OR_N", category) ~ "Local recreation",
                                   grepl("^OR_C", category) ~ "Citywide recreation"),
                         levels = PURPOSES),
      travel_mode = factor(travel_mode, levels = 1:6, labels = MODES),
      src = unname(IMP_SOURCE[as.character(travel_mode)]),


      # Mode-specific (alternative-specific) impedance

      T = case_when(src == "walk" ~ T_walk, src == "cycle" ~ T_cycle,
                    src == "publict" ~ T_publict, src == "car" ~ T_car),
      D = case_when(src == "walk" ~ D_walk, src == "cycle" ~ D_cycle,
                    src == "publict" ~ D_publict, src == "car" ~ D_car),


      # No public transport connection: indicator on the public transport row; impedance = 0

      pt_na = as.integer(pt_missing == "indicator" & travel_mode == "TM_publict" & pt_none),
      T = ifelse(pt_na == 1L, 0, T),
      D = ifelse(pt_na == 1L, 0, D),

      # Transport infrastructure at NUTS4 level, mode-specific

      TS_walk_z = as.numeric(scale(walk_lengper100km2_nuts4)),
      TS_active_z = as.numeric(scale(bike_lengper100km2_nuts4 + walk_lengper100km2_nuts4)),
      TS_PT_z = as.numeric(scale(PT_countper100km2_nuts4)),
      TS_car_z = as.numeric(scale(road_lengper100km2_nuts4)),
      transport_supply = case_when(
        travel_mode == "TM_walk" ~ TS_walk_z,
        travel_mode %in% c("TM_cycle", "TM_rollerskates_skateboard") ~ TS_active_z,
        travel_mode == "TM_publict" ~ TS_PT_z,
        TRUE ~ TS_car_z
      )
    )

  n_raw <- c(rows = nrow(d), od = n_distinct(d$trip_id), resp = n_distinct(d$respondent))

  unknown <- unique(d$category[!is.na(d$category) & is.na(d$category1)])
  if (length(unknown)) stop("Unknown trip category codes: ", paste(unknown, collapse = ", "))

  d <- d %>%
    filter(!is.na(travel_mode), !is.na(category), !is.na(respondent), !is.na(gender),
           !is.na(age), !is.na(D_walk))

  n_age_out <- d %>% filter(is.na(age_group)) %>% summarise(od = n_distinct(trip_id)) %>% pull(od)
  d <- d %>% filter(!is.na(age_group))   # the submitted model dropped these rows silently


  # Remove origin-destination pairs with no reported mode

  d <- d %>%
    group_by(trip_id) %>%
    mutate(n_modes = sum(TM_bin)) %>%
    ungroup() %>%
    filter(n_modes > 0) %>%
    mutate(multimodal = factor(ifelse(n_modes > 1, "Multiple modes", "Single mode"),
                               levels = c("Single mode", "Multiple modes"))) %>%
    droplevels()

  d$travel_mode <- factor(d$travel_mode, levels = MODES)

  # Transport infrastructure residualised on settlement type

  res_TI <- lm(transport_supply ~ living_place, data = d)
  d$transport_supply_resid_z <- as.numeric(scale(residuals(res_TI)))

  od <- d %>%
    distinct(trip_id, .keep_all = TRUE) %>%
    dplyr::select(respondent, trip_id, category, category1, gender, age_group, living_place,
                  multimodal, n_modes, pt_none, T_walk, D_walk, T_cycle, D_cycle,
                  T_publict, D_publict, T_car, D_car)

  list(
    data = d, od = od, caps = caps, cap_q = cap_q, pt_missing = pt_missing,
    info = list(n_raw = n_raw, n_age_out = n_age_out,
                n_final = c(rows = nrow(d), od = nrow(od), resp = n_distinct(d$respondent)),
                r2_TI = summary(res_TI)$r.squared)
  )
}


# --------------------------------------------------------------------------------------
# ADDITIONAL FUNCTIONS: model fitting, design matrices, coefficient draws, tests
# --------------------------------------------------------------------------------------


log_term <- function(d) {
  ok <- !(d$travel_mode == "TM_publict" & d$pt_none)
  sprintf("I(log1p(imp) - %s)", signif(mean(log1p(d$imp[ok]), na.rm = TRUE), 6))
}


# Natural cubic spline of log(1 + impedance), 3 df

spline_term <- function(d) {
  li <- log1p(d$imp[d$pt_na == 0])
  kn <- signif(quantile(li, c(1 / 3, 2 / 3), na.rm = TRUE), 6)
  bk <- signif(range(li, na.rm = TRUE), 6)
  sprintf("ns(log1p(imp), knots = c(%s, %s), Boundary.knots = c(%s, %s))", kn[1], kn[2], bk[1], bk[2])
}


# Purpose_term: impedance term interacted with trip purpose 

fixed_part <- function(imp_term, settle_slope = TRUE, purpose_slope = FALSE,
                       multimodal = FALSE, pt_indicator = (SET$pt_missing == "indicator"),
                       purpose_term = imp_term) {
  terms <- c(
    "0 + travel_mode",
    "travel_mode:gender",
    "travel_mode:age_group",
    paste0("travel_mode:", imp_term),
    if (settle_slope) paste0("travel_mode:", imp_term, ":living_place"),
    if (purpose_slope) paste0("travel_mode:", purpose_term, ":category1"),
    "travel_mode:category1",
    "travel_mode:living_place",
    "travel_mode:transport_supply_resid_z",
    if (multimodal) "travel_mode:multimodal",
    if (pt_indicator) "pt_na"
  )
  paste(terms, collapse = " + ")
}

# fast = FALSE: bobyqa, Hessian-based covariance
# fast = TRUE: same optimiser and estimates, but the Hessian is not computed, so standard errors are not usable

fit_mode_model <- function(fixed, data, fast = SET$fast_all) {
  form <- as.formula(paste("TM_bin ~", fixed, "+ (1 | trip_id) + (1 | respondent)"))
  ff <- lme4::nobars(form)
  pA <- qr(model.matrix(ff, model.frame(ff, data)))$rank
  ctrl <- lme4::glmerControl(optimizer = "bobyqa", calc.derivs = !fast,
                             optCtrl = list(maxfun = 2e5))

  # weakly informative N(0, 2.5^2) prior on all fixed effects

  eval(bquote(blme::bglmer(
    form, data = data, family = binomial(),
    fixef.prior = normal(cov = diag(.(2.5^2), .(pA))),
    control = ctrl
  )))
}

design <- function(m, newdata) {
  tt <- stats::delete.response(stats::terms(lme4::nobars(formula(m))))
  X <- model.matrix(tt, model.frame(tt, newdata, na.action = na.pass))
  b <- lme4::fixef(m)
  miss <- setdiff(names(b), colnames(X))
  if (length(miss)) stop("Design matrix lacks columns: ", paste(miss, collapse = ", "))
  X[, names(b), drop = FALSE]
}

draw_coefs <- function(m, n) {
  set.seed(SET$seed)
  MASS::mvrnorm(n, lme4::fixef(m), as.matrix(vcov(m)))
}

wald_test <- function(m, pattern, label = pattern) {
  b <- lme4::fixef(m)
  V <- as.matrix(vcov(m))
  i <- grep(pattern, names(b))
  if (!length(i)) return(NULL)
  W <- as.numeric(t(b[i]) %*% solve(V[i, i, drop = FALSE]) %*% b[i])
  data.frame(test = label, df = length(i), chi2 = W,
             p_value = pchisq(W, length(i), lower.tail = FALSE))
}

coef_table <- function(m) {
  b <- lme4::fixef(m)
  se <- sqrt(diag(as.matrix(vcov(m))))
  data.frame(term = names(b), estimate = b, std_error = se,
             OR = exp(b), OR_low = exp(b - 1.96 * se), OR_high = exp(b + 1.96 * se),
             p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE), row.names = NULL)
}

re_table <- function(m) {
  vc <- as.data.frame(lme4::VarCorr(m))
  data.frame(group = vc$grp, variance = vc$vcov)
}


# --------------------------------------------------------------------------------------
# ADDITIONAL FUNCTIONS: crossover points in predicted mode-use probabilities
# --------------------------------------------------------------------------------------

build_axis_map <- function(od, M, axis = SET$axis) {
  if (axis == "same") return(function(x, mode) x)
  xv <- od[[M$od_var]]
  fits <- list()
  for (s in setdiff(unique(IMP_SOURCE), "walk")) {
    y <- od[[paste0(M$var, "_", s)]]
    ok <- is.finite(y) & is.finite(xv) & !(s == "publict" & od$pt_none)
    fits[[s]] <- lm(log1p(y) ~ splines::ns(log1p(x), df = 3),
                    data = data.frame(y = y[ok], x = xv[ok]))
  }
  function(x, mode) {
    s <- IMP_SOURCE[[mode]]
    if (s == "walk") return(x)
    pmax(expm1(predict(fits[[s]], newdata = data.frame(x = x))), 0)
  }
}

make_profiles <- function(od, vars) {
  od %>% count(across(all_of(vars)), name = "w")
}

predict_curves <- function(m, profiles, purpose, M, map, B = NULL, levels_ref) {
  xg <- M$grid
  nd <- tidyr::expand_grid(profiles, mode_i = seq_along(MODES), xi = seq_along(xg))
  nd$travel_mode <- factor(MODES[nd$mode_i], levels = MODES)
  nd$x <- xg[nd$xi]
  nd$imp <- NA_real_
  for (k in seq_along(MODES)) {
    s <- nd$mode_i == k
    nd$imp[s] <- map(nd$x[s], MODES[k])
  }
  nd$category1 <- factor(purpose, levels = levels_ref$category1)
  nd$transport_supply_resid_z <- 0
  nd$pt_na <- 0L
  if (!"multimodal" %in% names(nd)) {
    nd$multimodal <- factor("Single mode", levels = levels_ref$multimodal)
  }

  X <- design(m, nd)
  gi <- (nd$xi - 1) * length(MODES) + nd$mode_i
  wsum <- rowsum(nd$w, gi)

  p_pt <- plogis(as.vector(X %*% lme4::fixef(m)))
  point <- t(matrix(rowsum(p_pt * nd$w, gi) / wsum, nrow = length(MODES)))
  colnames(point) <- MODES

  draws <- NULL
  if (!is.null(B)) {
    P <- plogis(X %*% t(B))
    A <- rowsum(P * nd$w, gi) / as.vector(wsum)
    draws <- aperm(array(A, dim = c(length(MODES), length(xg), nrow(B))), c(2, 1, 3))
    dimnames(draws) <- list(NULL, MODES, NULL)
  }
  list(point = point, draws = draws, x = xg)
}

cross_x <- function(d, x) {
  # first point where d moves from positive to non-positive (linear interpolation)
  i <- which(d[-length(d)] > 0 & d[-1] <= 0)[1]
  if (is.na(i)) return(NA_real_)
  x[i] + d[i] * (x[i + 1] - x[i]) / (d[i] - d[i + 1])
}

crossover_summary <- function(cv) {
  x <- cv$x
  P <- cv$point
  D <- cv$draws
  wi <- which(MODES == "TM_walk")
  oth <- setdiff(seq_along(MODES), wi)

  first_exc <- function(Pm) cross_x(Pm[, wi] - apply(Pm[, oth, drop = FALSE], 1, max), x)
  summarise_draws <- function(v) {
    c(lo = unname(quantile(v, 0.025, na.rm = TRUE)),
      hi = unname(quantile(v, 0.975, na.rm = TRUE)),
      p_exist = mean(!is.na(v)))
  }

  rows <- list()

  # Predicted probability of walking first exceeded by that of any other mode

  xp <- first_exc(P)
  exceeder <- NA_character_
  if (!is.na(xp)) {
    j <- which(x >= xp)[1]
    exceeder <- MODES[oth][which.max(P[j, oth])]
  }
  dr <- if (!is.null(D)) vapply(seq_len(dim(D)[3]), function(s) first_exc(D[, , s]), numeric(1)) else NA
  rows[[1]] <- data.frame(
    type = "first", mode_from = "TM_walk", mode_to = exceeder,
    crossover = paste0("Walking first exceeded", ifelse(is.na(exceeder), "",
                                                       paste0(" (by ", MODE_LABS[exceeder], ")"))),
    x = xp, t(summarise_draws(dr)),
    p_at = if (is.na(xp)) NA else approx(x, P[, wi], xout = xp)$y
  )


  # Pairwise crossovers

  for (pr in PAIRS) {
    xp <- cross_x(P[, pr[1]] - P[, pr[2]], x)
    dr <- if (!is.null(D)) vapply(seq_len(dim(D)[3]),
                                  function(s) cross_x(D[, pr[1], s] - D[, pr[2], s], x),
                                  numeric(1)) else NA
    rows[[length(rows) + 1]] <- data.frame(
      type = "pair", mode_from = pr[1], mode_to = pr[2],
      crossover = paste(MODE_LABS[pr[1]], "/", MODE_LABS[pr[2]]),
      x = xp, t(summarise_draws(dr)),
      p_at = if (is.na(xp)) NA else approx(x, P[, pr[1]], xout = xp)$y
    )
  }
  bind_rows(rows)
}

# Data support around each crossover

add_support <- function(cs, purpose_x, M, map, max_obs) {
  cs$n_od_within_band <- vapply(cs$x, function(z) if (is.na(z)) NA_real_ else
    sum(abs(purpose_x - z) <= M$band, na.rm = TRUE), numeric(1))
  cs$percentile_of_x <- vapply(cs$x, function(z) if (is.na(z)) NA_real_ else
    mean(purpose_x <= z, na.rm = TRUE), numeric(1))
  cs$beyond_95th <- cs$percentile_of_x > 0.95
  cs$extrapolated_mode <- mapply(function(z, a, b) {
    if (is.na(z)) return(NA_character_)
    ms <- c(a, b)
    ms <- ms[!is.na(ms)]
    ex <- ms[vapply(ms, function(mo) map(z, mo) > max_obs[[mo]], logical(1))]
    if (length(ex)) paste(MODE_LABS[ex], collapse = "; ") else ""
  }, cs$x, cs$mode_from, cs$mode_to)
  cs
}

curves_df <- function(cv, purpose, extra = list()) {
  P <- cv$point
  D <- cv$draws
  out <- expand.grid(xi = seq_along(cv$x), mode = MODES, stringsAsFactors = FALSE)
  out$x <- cv$x[out$xi]
  out$p <- P[cbind(out$xi, match(out$mode, MODES))]
  if (!is.null(D)) {
    q <- apply(D, c(1, 2), quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    out$lo <- q[1, , ][cbind(out$xi, match(out$mode, MODES))]
    out$hi <- q[2, , ][cbind(out$xi, match(out$mode, MODES))]
  }
  out$purpose <- purpose
  for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  out
}

# Full crossover analysis for fitted model: purposes, intervals, data support

run_crossovers <- function(m, prep, M, profile_vars = c("gender", "age_group", "living_place"),
                           od_filter = NULL, B = NULL, map = NULL, weighting = SET$weighting) {
  dat <- prep$data
  od <- prep$od
  if (is.null(map)) map <- build_axis_map(od, M)
  max_obs <- tapply(dat[[M$var]][dat$pt_na == 0], dat$travel_mode[dat$pt_na == 0], max, na.rm = TRUE)
  lev <- list(category1 = levels(dat$category1), multimodal = levels(dat$multimodal))
  od_sub <- if (is.null(od_filter)) od else od[od_filter(od), , drop = FALSE]

  # Population weights apply to gender x age group x settlement cells

  use_pop <- weighting == "population" && !is.null(POP) &&
    setequal(profile_vars, c("gender", "age_group", "living_place"))
  if (use_pop) {
    profiles <- POP
    if (!is.null(od_filter)) profiles <- profiles[od_filter(profiles), , drop = FALSE]
    profiles <- profiles[, c(profile_vars, "w")]
  } else {
    profiles <- make_profiles(od_sub, profile_vars)
  }

  res <- list()
  crv <- list()
  for (pp in lev$category1) {
    cv <- predict_curves(m, profiles, pp, M, map, B = B, levels_ref = lev)
    cs <- crossover_summary(cv)
    px <- od_sub[[M$od_var]][od_sub$category1 == pp]
    cs <- add_support(cs, px, M, map, max_obs)
    cs$purpose <- pp
    cs$n_resp <- n_distinct(od_sub$respondent[od_sub$category1 == pp])
    cs$n_od <- sum(od_sub$category1 == pp)
    cs$weighting <- if (use_pop) "population" else "sample"
    res[[pp]] <- cs
    crv[[pp]] <- curves_df(cv, pp)
  }
  list(cross = bind_rows(res), curves = bind_rows(crv))
}

format_cross <- function(cs, M) {
  cs %>%
    mutate(
      # interval only reported when the point estimate has a crossover in the plotted range;
      # p_exist = share of coefficient draws with a crossover in that range
      lo = ifelse(is.na(x), NA, lo),
      hi = ifelse(is.na(x), NA, hi),
      label = ifelse(is.na(x), "no crossover in range",
                     paste0(M$fmt(x), " [", M$fmt(lo), " - ", M$fmt(hi), "]")),
      flag = case_when(
        is.na(x) ~ "",
        p_exist < 0.95 ~ "crossover absent in >5% of draws",
        extrapolated_mode != "" ~ paste("extrapolated for", extrapolated_mode),
        beyond_95th ~ "beyond 95th percentile of observed pairs",
        TRUE ~ ""
      )
    )
}

# by_settlement = FALSE: one panel per trip purpose (Figures 2-3)
# by_settlement = TRUE: rows = settlement type, columns = trip purpose (Figures 4-5)


plot_crossovers <- function(curves, cross, rug, M, by_settlement = FALSE) {
  pl <- PURPOSES
  sl <- c("Large city (Lodz)", "Medium-sized cities", "Small towns and rural areas")
  curves$purpose <- factor(curves$purpose, levels = pl)
  cross$purpose <- factor(cross$purpose, levels = pl)
  rug$purpose <- factor(rug$purpose, levels = pl)
  if (by_settlement) {
    curves$living_place <- factor(curves$living_place, levels = sl)
    cross$living_place <- factor(cross$living_place, levels = sl)
    rug$living_place <- factor(rug$living_place, levels = sl)
  }
  facets <- if (by_settlement) {
    facet_grid(living_place ~ purpose, labeller = labeller(living_place = label_wrap_gen(14)))
  } else {
    facet_wrap(~ purpose)
  }


  # Crossovers flagged as uncertain (sparse data, extrapolation, absent in >5% of simulation draws)

  first <- cross %>% filter(type == "first")
  fe <- first %>% filter(!is.na(x), flag == "")
  ne <- first %>% filter(is.na(x) | flag != "") %>%
    mutate(note = ifelse(is.na(x), "Walking not exceeded\nwithin the plotted range",
                         "No reliable crossover\nwithin the observed range"))
  if (nrow(fe)) fe$lab <- paste0(M$fmt(fe$x), "\n[", M$fmt(fe$lo), " - ", M$fmt(fe$hi), "]")
  g <- ggplot(curves, aes(x = x, y = p, colour = mode, fill = mode)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_rug(data = rug, aes(x = x), inherit.aes = FALSE, sides = "b",
             alpha = 0.35, length = grid::unit(0.025, "npc")) +
    facets +
    scale_colour_discrete(name = NULL, labels = MODE_LABS) +
    scale_fill_discrete(name = NULL, labels = MODE_LABS) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = M$axis_label[[SET$axis]], y = "Predicted probability of reporting the mode") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  if (nrow(fe)) {
    g <- g +
      geom_segment(data = fe, aes(x = lo, xend = hi, y = p_at, yend = p_at),
                   inherit.aes = FALSE, linewidth = 0.6) +
      geom_point(data = fe, aes(x = x, y = p_at), inherit.aes = FALSE, size = 2) +
      geom_label(data = fe, aes(x = x, y = pmin(p_at + 0.14, 0.95), label = lab),
                 inherit.aes = FALSE, size = 2.7, label.size = 0.15, alpha = 0.85)
  }
  if (nrow(ne)) {
    g <- g + geom_text(data = ne, aes(x = max(M$grid), y = 0.97, label = note),
                       inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.7, colour = "grey30")
  }
  g
}


# --------------------------------------------------------------------------------------
# ADDITIONAL FUNCTIONS: average marginal effects (Table 1)
# --------------------------------------------------------------------------------------

ame_table <- function(m, dat, M, factors = c("gender", "age_group", "living_place", "category1")) {
  tidy <- function(x, variable) {
    x <- as.data.frame(x)
    data.frame(variable = variable, travel_mode = as.character(x$travel_mode),
               contrast = as.character(x$contrast), estimate = x$estimate,
               conf_low = x$conf.low, conf_high = x$conf.high, p_value = x$p.value)
  }
  res <- list()
  for (f in factors) {
    # trip purpose: all three pairwise contrasts (local vs school, citywide vs school,
    # citywide vs local); other factors: each level vs the reference level
    v <- if (f == "category1") list(category1 = "pairwise") else f
    res[[f]] <- tidy(avg_comparisons(m, variables = v, by = "travel_mode", type = "response"), f)
  }

  # Impedance: change in probability for +1 km / +10 min of the mode's own impedance

  d0 <- dat[dat$pt_na == 0, ]
  res$impedance <- tidy(avg_comparisons(m, variables = setNames(list(M$unit), "imp"),
                                        by = "travel_mode", newdata = d0, type = "response"),
                        paste("impedance", M$unit_label))
  res$infrastructure <- tidy(avg_slopes(m, variables = "transport_supply_resid_z",
                                        by = "travel_mode", type = "response"),
                             "transport infrastructure (residualised, per SD)")
  bind_rows(res) %>% mutate(mode = MODE_LABS[travel_mode])
}




# ======================================================================================
# Load data, prepare, describe
# ======================================================================================

raw <- read.csv(SET$data_url, stringsAsFactors = FALSE, sep = ";", dec = ".")

prep <- prepare_data(raw)
dat <- prep$data
od <- prep$od

cat("\nRows / OD pairs / respondents before filtering:", prep$info$n_raw, "\n")
cat("Final rows / OD pairs / respondents:", prep$info$n_final, "\n")

# Population shares for post-stratified curves and crossovers

POP <- read.csv(SET$population_file, colClasses = "character", sep = ";", dec = ".")

age_key <- function(x) gsub("[^0-9]", "", as.character(x))   # "7–15" and "7-15" both match

if (!is.null(POP)) {
  POP <- POP %>%
    mutate(
      gender = factor(trimws(gender), levels = levels(dat$gender)),
      age_group = factor(levels(dat$age_group)[match(age_key(age_group), age_key(levels(dat$age_group)))],
                         levels = levels(dat$age_group)),
      living_place = factor(trimws(living_place), levels = levels(dat$living_place)),
      w = as.numeric(share) / sum(as.numeric(share))
    )
  if (anyNA(POP[, c("gender", "age_group", "living_place", "w")])) {
    stop("Unmatched codes or missing shares in ", SET$population_file)
  }
  if (nrow(POP) != nlevels(dat$gender) * nlevels(dat$age_group) * nlevels(dat$living_place)) {
    stop(SET$population_file, " must have one row per gender x age group x settlement cell")
  }
}


# Multiple responses (descriptive)

combos <- dat %>%
  filter(TM_bin == 1) %>%
  group_by(trip_id, category1) %>%
  summarise(combination = paste(MODE_LABS[as.character(travel_mode)], collapse = " + "),
            n_modes = n(), .groups = "drop")
save_csv(bind_rows(
  combos %>% group_by(category1) %>%
    summarise(od_pairs = n(), share_multiple = mean(n_modes > 1), .groups = "drop") %>%
    mutate(combination = "(share of OD pairs with more than one reported mode)"),
  combos %>% filter(n_modes > 1) %>% count(category1, combination, name = "od_pairs") %>%
    group_by(category1) %>% slice_max(od_pairs, n = 8, with_ties = FALSE) %>% ungroup()
), "ESM_multiple_responses.csv")


# ======================================================================================
# Analyses: travel time and for travel distance
# ======================================================================================

results <- list()

for (M in MEASURES) {

  cat("\n\n==================== Impedance:", M$name, "====================\n")
  dat$imp <- dat[[M$var]]
  prep$data <- dat
  map <- build_axis_map(od, M)

  # ------------------------------------------------------------------------------
  # Baseline model
  # ------------------------------------------------------------------------------

  LT <- log_term(dat)       # centred log impedance (sensitivity specifications)
  ST <- spline_term(dat)    # natural spline of log impedance (primary specification)
  cat("Impedance term (primary):", ST, "\n")
  f_primary <- fixed_part(imp_term = ST, settle_slope = FALSE)
  m_primary <- fit_mode_model(f_primary, dat)
  print(summary(m_primary)$optinfo$conv$lme4$messages)

  save_csv(coef_table(m_primary), paste0("ESM_coefficients_ORs_", M$name, ".csv"))
  save_csv(re_table(m_primary), paste0("Table1_random_effects_", M$name, ".csv"))
  wald_tab <- bind_rows(          # saved after the functional-form section (adds impedance x settlement)
    wald_test(m_primary, "gender"),
    wald_test(m_primary, "age_group"),
    wald_test(m_primary, "living_place", "settlement type"),
    wald_test(m_primary, "imp", "impedance (all spline terms)"),
    wald_test(m_primary, "category1", "trip purpose"),
    wald_test(m_primary, "transport_supply", "transport infrastructure")
  )

  # Table 1: average marginal effects
  ame <- ame_table(m_primary, dat, M)
  save_csv(ame, paste0("Table1_AME_", M$name, ".csv"))
  cat("AIC primary:", AIC(m_primary), "\n")

  # ------------------------------------------------------------------------------
  # Crossovers with simulation intervals and data support
  # ------------------------------------------------------------------------------

  B <- draw_coefs(m_primary, SET$n_sim)

  # Figures 2 (time) and 3 (distance) with their crossover tables, in two variants:
  # curves averaged over the regional population and over the sample composition
  WEIGHTS <- if (is.null(POP)) "sample" else c("population", "sample")
  rug <- od %>% transmute(purpose = as.character(category1), x = .data[[M$od_var]]) %>%
    filter(x <= max(M$grid))
  fig_no <- if (M$name == "time") "Fig2" else "Fig3"
  for (wt in WEIGHTS) {
    cr_w <- run_crossovers(m_primary, prep, M, B = B, map = map, weighting = wt)
    cross_w <- format_cross(cr_w$cross, M)
    save_csv(cross_w, paste0("Fig_crossovers_table_", M$name, "_", wt, ".csv"))
    ggsave(out(paste0(fig_no, "_crossovers_", M$name, "_", wt, ".png")),
           plot_crossovers(cr_w$curves, cross_w, rug, M),
           width = 12, height = 4.8, dpi = 600, bg = "white")
    if (wt == SET$weighting) cross_main <- cross_w
  }
  print(cross_main %>% dplyr::select(purpose, crossover, label, p_exist, n_od_within_band, flag))

  # Figures 4 (time) and 5 (distance): curves by settlement type. Each panel is averaged over
  # the population of that settlement type (gender x age group shares within it); the rug
  # shows the observed origin-destination pairs of that settlement type.
  cr_set <- list()
  crv_set <- list()
  for (lp in levels(dat$living_place)) {
    r <- run_crossovers(m_primary, prep, M, B = B, map = map,
                        od_filter = function(o) o$living_place == lp)
    cr_set[[lp]] <- r$cross %>% mutate(living_place = lp)
    crv_set[[lp]] <- r$curves %>% mutate(living_place = lp)
  }
  cross_set <- format_cross(bind_rows(cr_set), M)
  rug_set <- od %>% transmute(purpose = as.character(category1), living_place = as.character(living_place),
                              x = .data[[M$od_var]]) %>% filter(x <= max(M$grid))
  fig_set <- if (M$name == "time") "Fig4" else "Fig5"
  g_set <- plot_crossovers(bind_rows(crv_set), cross_set, rug_set, M, by_settlement = TRUE)
  ggsave(out(paste0(fig_set, "_crossovers_by_settlement_", M$name, ".png")), g_set,
         width = 12, height = 10, dpi = 600, bg = "white")

  # ESM Tables 3-4: crossovers by gender, age group and settlement type, with n

  # (population and sample weights)
  for (wt in WEIGHTS) {
    sub <- list()
    for (v in c("gender", "age_group", "living_place")) {
      for (lv in levels(dat[[v]])) {
        r <- run_crossovers(m_primary, prep, M, B = B, map = map, weighting = wt,
                            od_filter = function(o) o[[v]] == lv)
        sub[[paste(v, lv)]] <- r$cross %>% mutate(variable = v, level = lv)
      }
    }
    sub_tab <- format_cross(bind_rows(sub), M) %>%
      dplyr::select(purpose, variable, level, n_resp, n_od, type, crossover, x, lo, hi, label,
                    p_exist, n_od_within_band, flag, weighting)
    save_csv(sub_tab, paste0("ESM_Table", ifelse(M$name == "time", 3, 4),
                             "_crossovers_subgroups_", M$name, "_", wt, ".csv"))
  }


  # ------------------------------------------------------------------------------
  # Trip purpose x impedance
  # ------------------------------------------------------------------------------

  # spline (as in the primary model) + purpose-specific deviation of the log-impedance slope
  m_int <- fit_mode_model(fixed_part(imp_term = ST, settle_slope = FALSE, purpose_slope = TRUE,
                                     purpose_term = LT), dat)
  ct <- coef_table(m_int)
  int_terms <- ct %>% filter(grepl("imp", term), grepl("category1", term)) %>%
    mutate(conf_low = estimate - 1.96 * std_error, conf_high = estimate + 1.96 * std_error)
  save_csv(int_terms, paste0("ESM_R1.14_interaction_coefficients_", M$name, ".csv"))

  wt <- wald_test(m_int, "imp.*category1|category1.*imp", "trip purpose x impedance (joint)")
  aic_cmp <- data.frame(model = c("primary", "primary + purpose x impedance"),
                        AIC = c(AIC(m_primary), AIC(m_int)), BIC = c(BIC(m_primary), BIC(m_int)))

  # Second differences: AME of +1 unit of impedance, one trip purpose minus another
  # (counterfactual purpose for all pairs with a PT connection; random effects at zero)
  d0 <- dat[dat$pt_na == 0, ]
  bI <- lme4::fixef(m_int)
  BI <- draw_coefs(m_int, SET$n_sim)
  eff <- list()
  for (pp in levels(dat$category1)) {
    dc <- d0
    dc$category1 <- factor(pp, levels = levels(dat$category1))
    d1 <- dc
    d1$imp <- d1$imp + M$unit
    X0 <- design(m_int, dc)
    X1 <- design(m_int, d1)
    g <- dc$travel_mode
    n_g <- as.vector(table(g))
    eff[[pp]] <- list(
      pt = as.vector(rowsum(plogis(X1 %*% bI) - plogis(X0 %*% bI), g)) / n_g,
      dr = rowsum(plogis(X1 %*% t(BI)) - plogis(X0 %*% t(BI)), g) / n_g
    )
  }
  # second differences for each pair of trip purposes
  q <- function(m, p) apply(m, 1, quantile, probs = p, na.rm = TRUE)
  purpose_pairs <- combn(rev(PURPOSES), 2, simplify = FALSE)   # later category minus earlier
  mod_tab <- bind_rows(lapply(purpose_pairs, function(pp) {
    a <- eff[[pp[1]]]
    b <- eff[[pp[2]]]
    dd <- a$dr - b$dr
    data.frame(
      comparison = paste(pp[1], "minus", pp[2]),
      mode = MODE_LABS[levels(dat$travel_mode)],
      AME_first = a$pt, AME_first_lo = q(a$dr, .025), AME_first_hi = q(a$dr, .975),
      AME_second = b$pt, AME_second_lo = q(b$dr, .025), AME_second_hi = q(b$dr, .975),
      second_difference = a$pt - b$pt,
      sd_lo95 = q(dd, .025), sd_hi95 = q(dd, .975),
      sd_lo90 = q(dd, .05), sd_hi90 = q(dd, .95)
    )
  })) %>%
    mutate(
      unit = M$unit_label,
      sesoi = SET$sesoi,
      equivalent_within_sesoi = sd_lo90 > -sesoi & sd_hi90 < sesoi,   # TOST at alpha = 0.05
      largest_difference_not_excluded = pmax(abs(sd_lo95), abs(sd_hi95))
    )
  save_csv(mod_tab, paste0("ESM_R1.14_second_differences_", M$name, ".csv"))
  save_csv(bind_rows(wt, data.frame(test = paste("AIC", aic_cmp$model), df = NA,
                                    chi2 = aic_cmp$AIC, p_value = NA)),
           paste0("ESM_R1.14_joint_test_", M$name, ".csv"))
  print(wt)
  print(mod_tab %>% dplyr::select(comparison, mode, second_difference, sd_lo95, sd_hi95,
                                  equivalent_within_sesoi))

  # ------------------------------------------------------------------------------
  # 6-8. Functional form of impedance
  # ------------------------------------------------------------------------------

  specs <- list(
    primary    = list(label = "Natural spline of log impedance, 3 df (primary)", fixed = f_primary),
    log        = list(label = "Log impedance",
                      fixed = fixed_part(imp_term = LT, settle_slope = FALSE)),
    log_settle = list(label = "Log impedance + settlement-specific slopes",
                      fixed = fixed_part(imp_term = LT)),
    linear     = list(label = paste0("Linear impedance (per ", M$unit_label, ", centred)"),
                      fixed = fixed_part(imp_term = sprintf("I((imp - %s) / %s)",
                                                            signif(mean(dat$imp[dat$pt_na == 0]), 6), M$unit),
                                         settle_slope = FALSE))
  )

  brier <- function(y, p) mean((y - p)^2)
  logloss <- function(y, p) { p <- pmin(pmax(p, 1e-6), 1 - 1e-6); -mean(y * log(p) + (1 - y) * log(1 - p)) }

  set.seed(SET$seed)
  resp <- levels(dat$respondent)
  folds <- setNames(sample(rep(seq_len(SET$k_folds), length.out = length(resp))), resp)

  ff_rows <- list()
  ff_cross <- list()
  ff_models <- list()
  for (s in names(specs)) {
    ms <- if (s == "primary") m_primary else
      tryCatch(fit_mode_model(specs[[s]]$fixed, dat),
               error = function(e) { message("Spec ", s, " failed: ", conditionMessage(e)); NULL })
    if (is.null(ms)) {
      ff_rows[[s]] <- data.frame(spec = paste(specs[[s]]$label, "(did not converge)"))
      next
    }
    p_in <- predict(ms, re.form = NA, type = "response")
    row <- data.frame(spec = specs[[s]]$label, AIC = AIC(ms), BIC = BIC(ms),
                      n_fixef = length(lme4::fixef(ms)),
                      brier_in_sample = brier(dat$TM_bin, p_in),
                      logloss_in_sample = logloss(dat$TM_bin, p_in))

    if (SET$run_cv) {
      oof <- rep(NA_real_, nrow(dat))
      for (k in seq_len(SET$k_folds)) {
        te <- folds[as.character(dat$respondent)] == k
        tr <- dat[!te, ]
        tr$respondent <- droplevels(tr$respondent)
        tr$trip_id <- droplevels(tr$trip_id)
        mk <- tryCatch(fit_mode_model(specs[[s]]$fixed, tr, fast = TRUE),
                       error = function(e) { message("CV fold ", k, " failed: ", conditionMessage(e)); NULL })
        if (!is.null(mk)) oof[te] <- predict(mk, newdata = dat[te, ], re.form = NA, type = "response")
      }
      row$cv_folds_failed <- sum(tapply(is.na(oof), folds[as.character(dat$respondent)], all))
      ok <- !is.na(oof)
      y_ok <- dat$TM_bin[ok]
      mode_ok <- dat$travel_mode[ok]
      oof <- oof[ok]
      row$brier_cv <- brier(y_ok, oof)
      row$logloss_cv <- logloss(y_ok, oof)
      for (mo in MODES) {
        i <- mode_ok == mo
        row[[paste0("brier_cv_", MODE_LABS[mo])]] <- brier(y_ok[i], oof[i])
      }
    }
    ff_rows[[s]] <- row
    ff_models[[s]] <- ms

    Bs <- if (s == "primary") B else draw_coefs(ms, SET$n_sim)
    ff_cross[[s]] <- run_crossovers(ms, prep, M, B = Bs, map = map)$cross %>%
      filter(type == "first" | (mode_from == "TM_walk" & mode_to %in% c("TM_car", "TM_publict"))) %>%
      format_cross(M) %>%
      dplyr::select(purpose, crossover, x, lo, hi, label, flag) %>%
      mutate(spec = specs[[s]]$label)
  }
  ff_tab <- bind_rows(ff_rows)
  save_csv(ff_tab, paste0("ESM_R1.8_functional_form_fit_", M$name, ".csv"))
  save_csv(bind_rows(ff_cross), paste0("ESM_R1.8_crossover_stability_", M$name, ".csv"))
  print(ff_tab[, intersect(c("spec", "AIC", "BIC", "brier_cv", "logloss_cv"), names(ff_tab))])

  if (!is.null(ff_models$log_settle)) {
    wald_tab <- bind_rows(wald_tab, wald_test(ff_models$log_settle, "imp.*living_place",
                                              "impedance x settlement (log model, settlement-specific slopes)"))
  }
  save_csv(wald_tab, paste0("ESM_wald_tests_", M$name, ".csv"))

  # ------------------------------------------------------------------------------
  # 9. Sensitivity analyses 
  # ------------------------------------------------------------------------------

  # fixed part is built from each variant's own data (its own spline knots)
  variants <- list(
    list(label = "Primary", prep = prep,
         fixed = function(d) fixed_part(imp_term = spline_term(d), settle_slope = FALSE),
         profile_vars = c("gender", "age_group", "living_place"), model = m_primary),
    list(label = "With multiple-mode indicator (submitted specification)", prep = prep,
         fixed = function(d) fixed_part(imp_term = spline_term(d), settle_slope = FALSE, multimodal = TRUE),
         profile_vars = c("gender", "age_group", "living_place", "multimodal")),
    list(label = "No PT connection coded as 100,000 (submitted coding)",
         prep = prepare_data(raw, SET$cap_q, "sentinel"),
         fixed = function(d) fixed_part(imp_term = spline_term(d), settle_slope = FALSE, pt_indicator = FALSE),
         profile_vars = c("gender", "age_group", "living_place")),
    list(label = "Winsorisation at 99th percentile",
         prep = prepare_data(raw, 0.99, SET$pt_missing),
         fixed = function(d) fixed_part(imp_term = spline_term(d), settle_slope = FALSE),
         profile_vars = c("gender", "age_group", "living_place"))
  )

  sens_ame <- list()
  sens_cross <- list()
  for (v in variants) {
    pv <- v$prep
    pv$data$imp <- pv$data[[M$var]]
    mv <- if (!is.null(v$model)) v$model else
      tryCatch(fit_mode_model(v$fixed(pv$data), pv$data),
               error = function(e) { message(v$label, " failed: ", conditionMessage(e)); NULL })
    if (is.null(mv)) next
    a <- ame_table(mv, pv$data, M, factors = "category1") %>%
      filter(variable != "transport infrastructure (residualised, per SD)")
    sens_ame[[v$label]] <- a %>% mutate(model = v$label, AIC = AIC(mv))
    Bv <- if (!is.null(v$model)) B else draw_coefs(mv, SET$n_sim)
    sens_cross[[v$label]] <- run_crossovers(mv, pv, M, profile_vars = v$profile_vars, B = Bv,
                                            map = build_axis_map(pv$od, M))$cross %>%
      format_cross(M) %>%
      dplyr::select(purpose, type, crossover, x, lo, hi, label, flag, weighting) %>%
      mutate(model = v$label)
  }
  save_csv(bind_rows(sens_ame), paste0("ESM_R1.5_sensitivity_AME_", M$name, ".csv"))
  save_csv(bind_rows(sens_cross), paste0("ESM_R1.5_sensitivity_crossovers_", M$name, ".csv"))

  results[[M$name]] <- list(m_primary = m_primary, m_int = m_int, ame = ame,
                            cross = cross_main, functional_form = ff_tab)
}

cat("\nDone. Outputs written to:", normalizePath(SET$out_dir), "\n")
print(list.files(SET$out_dir))
