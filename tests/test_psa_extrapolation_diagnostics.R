# Synthetic checks require no confidential data or production caches.
source('R/psa_extrapolation_diagnostics.R')
df_life <- expand.grid(age = 0:100, sex = c('female', 'male'), stringsAsFactors = FALSE)
df_life$qx <- ifelse(df_life$sex == 'female', .1, .2)
t <- c(0, .5, 1, 2)
m_ref <- population_reference_survival(c(30.25, 40), c('0', '1'), df_life, t)
stopifnot(max(abs(m_ref[, 1] - .9^t)) < 1e-12,
          max(abs(m_ref[, 2] - .8^t)) < 1e-12)
# Fractional baseline age must split follow-up at the attained birthday.
df_life$qx[df_life$sex == 'female' & df_life$age == 31] <- .3
m_ref <- population_reference_survival(30.5, '0', df_life, c(0, 1))
stopifnot(abs(m_ref[2, 1] - sqrt(.9 * .7)) < 1e-12)
# An unreached terminal qx=1 must not introduce NaN through zero times infinity.
df_life$qx[df_life$age == 100] <- 1
stopifnot(all(is.finite(population_reference_survival(30, '0', df_life, t))))
q <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
df_tail <- extrapolation_tail_summary(10:1, q, c(.1, .2))
stopifnot(all(df_tail$removed == c(1, 2)),
          max(abs(df_tail$contribution - c(1, 3) / 55)) < 1e-12,
          max(abs(df_tail$trimmed - c(6, 6.5))) < 1e-12)
cat('PASS: age/sex mapping, fractional attained-age integration, terminal mortality and RMST-ranked QALY accounting.\n')

# The ordering cache must be tied to both the coefficient and parameter draws.
source('R/cache_paths.R')
fixture <- tempfile(fileext = '.rds')
base <- list(time_horizon = 520, cl = 1/52, u_np = .73, u_p = .59, dr_effects = .04)
df_population <- data.frame(Age = c(60, 70), sex = c('0', '1'))
df_params <- data.frame(sim = 1:2, model_idx = c(3L, 3L))
sampling <- list(fingerprint = 'test')
inputs <- list(sampling_fingerprint = 'test', n_draws = 2L, time_horizon = 520,
  cl = 1/52, u_np = .73, u_p = .59, dr_effects = .04,
  prediction_population = cache_fingerprint(df_population),
  psa_parameters = cache_fingerprint(df_params), version = 'pfs_os_violations_v3_trapezoid')
saveRDS(list(fingerprint = cache_fingerprint(inputs)), fixture)
before <- tools::md5sum(fixture)
stopifnot(is.list(read_extrapolation_ordering(fixture, sampling, df_population, df_params, base)))
df_params$model_idx[1] <- 1L
stopifnot(inherits(try(read_extrapolation_ordering(fixture, sampling, df_population, df_params, base),
  silent = TRUE), 'try-error'), identical(before, tools::md5sum(fixture)))
unlink(fixture)
cat('PASS: stale ordering cache rejected without mutation.\n')
