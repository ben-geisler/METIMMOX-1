# EQ-5D-5L utility scoring ----------------------------------------------------
#
# Vectorised value-set functions retained from the archived QALY exploration
# notebook. Dimension inputs use EQ-5D-5L levels 1 (no problems) through 5
# (extreme problems/unable to perform).

#' Calculate EQ-5D-5L utilities using the Danish value set
#'
#' @param mo Mobility level.
#' @param sc Self-care level.
#' @param ua Usual-activities level.
#' @param pd Pain/discomfort level.
#' @param ad Anxiety/depression level.
#' @return A numeric vector of utility values.
eq5d5l_utility_dk <- function(mo, sc, ua, pd, ad) {
  1 - (
    0.041 * (mo == 2) + 0.054 * (mo == 3) + 0.157 * (mo == 4) + 0.220 * (mo == 5) +
      0.035 * (sc == 2) + 0.050 * (sc == 3) + 0.144 * (sc == 4) + 0.209 * (sc == 5) +
      0.033 * (ua == 2) + 0.040 * (ua == 3) + 0.139 * (ua == 4) + 0.174 * (ua == 5) +
      0.048 * (pd == 2) + 0.094 * (pd == 3) + 0.381 * (pd == 4) + 0.537 * (pd == 5) +
      0.072 * (ad == 2) + 0.191 * (ad == 3) + 0.430 * (ad == 4) + 0.618 * (ad == 5)
  )
}

#' Calculate EQ-5D-5L utilities using the UK value set
#'
#' @inheritParams eq5d5l_utility_dk
#' @return A numeric vector of utility values.
eq5d5l_utility_uk <- function(mo, sc, ua, pd, ad) {
  1 - (
    0.058 * (mo == 2) + 0.076 * (mo == 3) + 0.207 * (mo == 4) + 0.274 * (mo == 5) +
      0.050 * (sc == 2) + 0.080 * (sc == 3) + 0.164 * (sc == 4) + 0.203 * (sc == 5) +
      0.050 * (ua == 2) + 0.063 * (ua == 3) + 0.162 * (ua == 4) + 0.184 * (ua == 5) +
      0.063 * (pd == 2) + 0.084 * (pd == 3) + 0.276 * (pd == 4) + 0.335 * (pd == 5) +
      0.078 * (ad == 2) + 0.104 * (ad == 3) + 0.285 * (ad == 4) + 0.289 * (ad == 5)
  )
}
