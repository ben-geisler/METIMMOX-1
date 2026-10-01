# Joint OS/PFS Coefficient Sampling
Ben Geisler
2026-10-01

- [Change and scope](#change-and-scope)
- [Covariance construction](#covariance-construction)
- [Diagnostics](#diagnostics)
- [Deterministic versus probabilistic
  outcomes](#deterministic-versus-probabilistic-outcomes)
  - [Validation warnings](#validation-warnings)

# Change and scope

Issue \#159 implements the first proposed fix: estimate cross-endpoint
coefficient dependence from paired patient bootstrap fits and draw OS
and PFS coefficients together. The sampling method is mvn_joint_v2.
Every endpoint still has the multivariate-normal marginal distribution
centred at its original maximum-likelihood estimate with its
observed-information covariance. The deterministic model, parameter
means, utility distributions and common population weighting are
unchanged.

The bootstrap estimates dependence only. Its coefficient vectors are
never used as PSA draws. This preserves the normal approximation
introduced in issue \#156 while restoring an estimated relationship
between OS and PFS uncertainty. Coefficients use the unconstrained
optimisation scale (log scale for positive baseline parameters),
consistent with the [flexsurv covariance
definition](https://chjackson.github.io/flexsurv/reference/flexsurvreg.html).

# Covariance construction

For each of 1000 unstratified patient bootstrap resamples, both
endpoints are fitted to the same sampled rows. A pair is excluded if
either fit fails convergence, lacks a finite positive-definite
covariance, or changes the estimable coefficient set. There is no
rejection based on coefficient magnitude, outcome direction or survival
crossings. More than 10% failed pairs, or insufficient rank to estimate
the covariance, stops generation. Rejected pairs are never replaced with
the original estimates.

Let $B$ be the empirical covariance of the retained paired bootstrap
coefficients, with endpoint blocks $B_O$, $B_P$ and cross-block
$B_{OP}$. Let $V_O$ and $V_P$ be the original fits’ covariance matrices.
Using symmetric positive-definite matrix square roots, define

$$
C_{OP} = V_O^{1/2} B_O^{-1/2} B_{OP} B_P^{-1/2} V_P^{1/2},
\qquad
V = \begin{pmatrix}V_O & C_{OP} \\ C_{OP}^{T} & V_P\end{pmatrix}.
$$

This is a block-diagonal congruence transformation of $B$. It preserves
positive semidefiniteness and exactly retains both original marginal
covariance blocks. The cross-block is a calibrated estimate of bootstrap
dependence; it is not the unscaled empirical bootstrap cross-covariance.
Directly combining that raw cross-block with different marginal
covariance blocks need not produce a valid joint covariance matrix.

One normal vector is drawn with mean $(\hat\beta_O,\hat\beta_P)$ and
covariance $V$, then split into OS and PFS coordinates. No coefficient
or outcome centring correction is applied. Bootstrap and normal sampling
each initialise their RNG with the recorded seed, so normal variates do
not depend on the number of bootstrap attempts. All strategies share the
same endpoint vectors within a draw; the joint biomarker population
weights and economic parameters remain independent of survival
coefficients.

# Diagnostics

| Diagnostic                                 | Value         |
|:-------------------------------------------|:--------------|
| Bootstrap pairs attempted                  | 1000          |
| Successful pairs                           | 996           |
| Failed pairs                               | 4             |
| Normal PSA draws                           | 5000          |
| Minimum joint covariance eigenvalue        | 4.7652676e-07 |
| Maximum OS marginal covariance difference  | 0             |
| Maximum PFS marginal covariance difference | 0             |

Covariance-generation diagnostics

| Reason | Pairs |
|:---|---:|
| Lapack routine dgesv: system is exactly singular: U\[8,8\] = 0 | 2 |
| Lapack routine dgesv: system is exactly singular: U\[9,9\] = 0 | 1 |
| system is computationally singular: reciprocal condition number = 2.25374e-21 | 1 |

Excluded bootstrap pairs by reason

| Coefficient                  | Target | Realised |
|:-----------------------------|-------:|---------:|
| shape                        |  0.332 |    0.345 |
| rate                         |  0.682 |    0.692 |
| Age                          |  0.723 |    0.734 |
| sex1                         |  0.682 |    0.674 |
| RxExperimental arm           |  0.651 |    0.651 |
| crp1                         |  0.409 |    0.392 |
| tmb_braf1                    |  0.679 |    0.687 |
| RxExperimental arm:crp1      |  0.488 |    0.487 |
| RxExperimental arm:tmb_braf1 |  0.688 |    0.701 |

OS/PFS correlations for corresponding coefficient coordinates

The bootstrap dependence estimate has finite-sample and Monte Carlo
uncertainty that is not propagated as another PSA layer. The trial is
small and sparse treatment/biomarker cells can make bootstrap fits
unstable. Failed-pair exclusion can also affect the dependence estimate;
the counts above disclose its extent. The covariance diagnostics and
retained bootstrap coefficients are stored in the ignored sampling cache
for audit.

# Deterministic versus probabilistic outcomes

The deterministic model evaluates nonlinear survival curves at the
fitted coefficients. The PSA averages those curves over coefficient
uncertainty. In general, $E[f(\beta)] \ne f(E[\beta])$. Before the
ordering correction, discounted QALYs are a linear combination of OS and
PFS areas with independent utilities. Changing only the OS/PFS
dependence therefore cannot change their expected raw-curve QALYs while
keeping both coefficient marginals fixed. It can change crossings, the
resulting clamp correction, the distribution of net benefit and value of
information. Finite Monte Carlo results can also change when drawing
from the new covariance.

The zero-uncertainty regression sets the entire joint coefficient
covariance to zero and fixes economic parameters and population weights.
It reproduces every deterministic cost within EUR 1e-8 and QALY within
1e-12. This verifies agreement of the calculation paths, without
recentering uncertain draws or relaxing the existing diagnostic.

| Strategy | Outcome | Base   | PSA mean | MC SE  | Difference (SE) |
|:---------|:--------|:-------|:---------|:-------|:----------------|
| control  | Cost    | 25764  | 25732    | 38     | -0.83           |
| control  | QALYs   | 1.3400 | 1.3775   | 0.0042 | +8.97           |
| crp      | Cost    | 58203  | 57740    | 90     | -5.16           |
| crp      | QALYs   | 1.3775 | 1.4007   | 0.0038 | +6.11           |
| tmb_braf | Cost    | 66916  | 66302    | 97     | -6.34           |
| tmb_braf | QALYs   | 1.3389 | 1.3696   | 0.0039 | +7.89           |

PSA means versus fitted-parameter results

| Strategy | Outcome | Base    | PSA mean | MC SE  | Difference (SE) |
|:---------|:--------|:--------|:---------|:-------|:----------------|
| crp      | Cost    | 32439   | 32008    | 81     | -5.30           |
| crp      | QALYs   | 0.0375  | 0.0232   | 0.0027 | -5.39           |
| tmb_braf | Cost    | 41152   | 40570    | 88     | -6.62           |
| tmb_braf | QALYs   | -0.0011 | -0.0079  | 0.0025 | -2.76           |

Incremental PSA means versus fitted-parameter increments

At the unchanged five-Monte-Carlo-SE threshold, 5 of 6 level comparisons
exceed the criterion. The criterion tests equality of outcome means; it
is not a general consequence of correctly centred coefficient draws in a
nonlinear model. Its numerical result remains visible in the test
protocol. Expected-net-benefit decisions under uncertainty use the PSA;
the deterministic results describe the fitted-parameter scenario.

Joint normal draws do not enforce OS \>= PFS. The PSA retains its
ordering clamp and [the ordering diagnostics](pfs_os_violations.md)
quantify its effect separately. Neither nonzero correlation nor
conditioning alone guarantees ordered survival curves.

------------------------------------------------------------------------

**Report completed on:** 2026-10-01  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 1.0

## Validation warnings

No warnings recorded during rendering.
