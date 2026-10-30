# rdlocal methodology

This document explains the statistics rdlocal implements, how it was built, and the evidence that it is correct. For the source-by-source record (papers, checksums, the function-to-equation map and every deviation), see [`CLEANROOM.md`](https://github.com/mzfortech/rdlocal/blob/main/CLEANROOM.md).

## Contents

1. [Regression discontinuity in one page](#1-regression-discontinuity-in-one-page)
2. [Local polynomial estimation](#2-local-polynomial-estimation)
3. [The bandwidth trade-off](#3-the-bandwidth-trade-off)
4. [Why the naive confidence interval undercovers](#4-why-the-naive-confidence-interval-undercovers)
5. [Robust bias-corrected inference](#5-robust-bias-corrected-inference)
6. [Coverage-optimal bandwidths](#6-coverage-optimal-bandwidths)
7. [Standard errors](#7-standard-errors)
8. [Manipulation testing](#8-manipulation-testing)
9. [How rdlocal was built](#9-how-rdlocal-was-built)
10. [Validation](#10-validation)
11. [Open problem: the d̂ pilot bandwidth](#11-open-problem-the-d-pilot-bandwidth)
12. [References](#12-references)

---

## 1. Regression discontinuity in one page

Many treatments are assigned by a cutoff rule: a scholarship for students scoring at least 80, a programme for towns above a population threshold, a seat won with more than 50% of the vote. Call the score the **running variable** $X$, the cutoff $c$, and the outcome $Y$.

Units just above and just below the cutoff are comparable. The jump in the average outcome at the cutoff therefore identifies the causal effect for units at the cutoff:

$$
\tau = \lim_{x \downarrow c} \mathbb{E}[Y \mid X = x] - \lim_{x \uparrow c} \mathbb{E}[Y \mid X = x].
$$

The identifying assumption is **continuity**: absent treatment, the average outcome would pass smoothly through the cutoff.

- **Sharp RD:** treatment switches on exactly at the cutoff.
- **Fuzzy RD:** crossing the cutoff only changes the *probability* of treatment. The effect is the ratio of the jumps in the outcome and in treatment take-up, a local Wald estimator:

$$
\tau_{\text{fuzzy}} = \frac{\mu_{Y+} - \mu_{Y-}}{\mu_{T+} - \mu_{T-}}.
$$

**Running example: US Senate elections** (Cattaneo, Frandsen & Titiunik 2015). $X$ is the Democratic margin of victory, the cutoff is 0, and $Y$ is the party's vote share in the next election for the same seat. The jump measures the incumbency advantage: about 7.4 percentage points.

## 2. Local polynomial estimation

The two limits are estimated from observations near the cutoff by weighted polynomial regressions on each side, within a bandwidth $h$ and with kernel weights $K$:

$$
\hat\beta_+ = \arg\min_{\beta} \sum_{i:\,X_i \ge c} K\!\left(\frac{X_i - c}{h}\right)\Big(Y_i - \textstyle\sum_{j=0}^{p} \beta_j (X_i - c)^j\Big)^2 ,
$$

and likewise below the cutoff. Then $\hat\tau = \hat\beta_{+,0} - \hat\beta_{-,0}$.

- **Order:** local linear ($p = 1$) is the standard choice. High-order *global* polynomials behave poorly at boundaries (Gelman & Imbens 2019).
- **Kernels:** triangular (the default), Epanechnikov and uniform.

## 3. The bandwidth trade-off

A small $h$ means little bias but a noisy estimate; a large $h$ means a precise but more biased one. For order $p$, the bias is about $h^{p+1}\mathcal{B}$ and the variance about $\mathcal{V}/(nh)$. $\mathcal{B}$ depends on the $(p+1)$-th derivative of the regression function on each side, and $\mathcal{V}$ on the conditional variance and the density of $X$ at the cutoff. Minimising the mean squared error gives

$$
h_{\text{MSE}} = \left(\frac{\mathcal{V}}{2(p+1)\,\mathcal{B}^2}\right)^{1/(2p+3)} n^{-1/(2p+3)} .
$$

$\mathcal{B}$ and $\mathcal{V}$ are unknown, so a **plug-in selector** estimates them with pilot regressions and pilot bandwidths. rdlocal implements the MSE- and coverage-optimal selectors of Calonico, Cattaneo & Titiunik (2014) and Calonico, Cattaneo & Farrell (2020), following the implementation steps in the supplement to Calonico, Cattaneo, Farrell & Titiunik (2019). It also implements Imbens & Kalyanaraman (2012).

## 4. Why the naive confidence interval undercovers

At $h_{\text{MSE}}$, bias and standard deviation are of the same order by construction. An interval $\hat\tau \pm 1.96\,\widehat{\text{SE}}$ is therefore centred in the wrong place and covers the true effect less often than its nominal 95%. In rdlocal's Monte Carlo study (§10.3), conventional intervals at $h_{\text{MSE}}$ covered 85–89%.

## 5. Robust bias-corrected inference

Calonico, Cattaneo & Titiunik (2014):

1. **Estimate the bias** with a higher-order fit ($q = p + 1$) at a pilot bandwidth $b$.
2. **Subtract it:** $\hat\tau^{\text{bc}} = \hat\tau - h^{p+1}\hat{\mathcal{B}}$.
3. **Recompute the variance.** The bias estimate is itself random, so the standard error must account for it. This gives the robust standard error $\hat\sigma_{\text{rbc}}$.
4. **Report** $\hat\tau^{\text{bc}} \pm 1.96\,\hat\sigma_{\text{rbc}}$.

Standard practice reports the conventional point estimate together with the robust bias-corrected interval.

## 6. Coverage-optimal bandwidths

For inference, the right target is coverage error rather than MSE. Calonico, Cattaneo & Farrell (2020) derive

$$
h_{\text{CER}} = h_{\text{MSE}} \cdot n^{-p/((p+3)(2p+3))},
$$

which for $p = 1$ is a factor of $n^{-1/20}$. For the Senate data, $n = 1297$ and $h_{\text{MSE}} = 17.708$ give $h_{\text{CER}} = 12.374$, as published in Calonico, Cattaneo, Farrell & Titiunik (2017, p. 400).

## 7. Standard errors

- **HC0–HC3:** heteroskedasticity-robust sandwich forms. HC1 uses the per-side factor $n_s/(n_s - k)$.
- **Nearest-neighbour (NN):** each residual variance is estimated from the observations with the closest values of $X$.
- **Cluster-robust:** built from fit residuals.

## 8. Manipulation testing

If units can sort precisely around the cutoff, the running variable's density jumps there. rdlocal implements the local-polynomial density estimator of Cattaneo, Jansson & Ma (2020). It regresses the empirical distribution function on $(X - c)$ on each side, reads the density off the linear coefficient, and tests

$$
T = \frac{\hat f_+(c) - \hat f_-(c)}{\hat\sigma} \;\approx\; N(0, 1) \text{ under no manipulation.}
$$

The empirical CDF leaves each observation's own point out (CJM 2020, supplement §5.2). The standard error is the analytic estimator of CJM's Theorem 2.

## 9. How rdlocal was built

rdlocal is an **independent, clean-room implementation**:

- It was written from the published papers, their supplements and the authors' published descriptions of their algorithms, all read in full from the original PDFs. Every formula in the code cites its equation, section and page.
- The source code, help files and documentation of other RD software were not consulted.
- Other implementations were used only as **black boxes**: data in, numbers out, to produce reference values for testing.
- Wherever a paper leaves a step unspecified, the code implements the closest stated version and marks it `DEVIATION`. Every deviation is listed in [`CLEANROOM.md`](https://github.com/mzfortech/rdlocal/blob/main/CLEANROOM.md).

Two findings from implementation:

1. **The 2014 plug-in step as printed mixes scales.** In the CCT (2014) supplement's plug-in steps (§S.2.6), the regularisation term adds $3\hat V$, a variance on the derivative scale carrying a $\nu!^2$ factor, to a squared Taylor *coefficient*, which does not carry it. The implementation that reproduces published bandwidths regularises on the coefficient scale.
2. **Defaults depend on pilot details.** Default bandwidths depend heavily on pilot-step details that the main texts leave to supplements. The comparison that settles correctness is therefore at equal bandwidths, where estimators can be compared exactly (§10.2).

## 10. Validation

### 10.1 Structural identities

- With a uniform kernel and fixed $h$, the estimate equals ordinary least squares within the window, to $10^{-10}$.
- Bias correction with $h = b$ equals the order-$q$ fit.
- The fuzzy estimate equals the ratio of the two sharp estimates.
- The density-test variance matches a brute-force evaluation of its defining sum.

### 10.2 Agreement with an established implementation at equal bandwidths

On the benchmark simulation designs of CCT (2014, supplement §S.3.1) and a fuzzy design, rdlocal agrees with the reference implementation at the same $h$ and $b$:

- estimates to about $10^{-13}$;
- robust standard errors to about $10^{-14}$;
- density estimates to about $10^{-12}$.

### 10.3 Monte Carlo coverage

CCT (2014) Model 1, $n = 1000$, 500 replications:

| Bandwidth | Robust coverage | Conventional coverage |
|---|---|---|
| MSE-optimal (two seeds) | 92.6%, 91.2% | 89.0%, 84.8% |
| CER-optimal (two seeds) | 94.2%, 93.4% | 93.4%, 91.0% |

For the same design with estimated bandwidths, CCT's supplement reports 92.7% robust coverage. rdlocal matches the published behaviour of the method.

### 10.4 Published empirical results

US Senate data (Cattaneo, Frandsen & Titiunik 2015), against the results printed in Calonico, Cattaneo, Farrell & Titiunik (2017):

| Quantity | Published | rdlocal |
|---|---|---|
| Estimate at $h = 16.79369$, $b = 27.43745$ (p. 393) | 7.4253 | 7.4253043 |
| Conventional SE | 1.4954 | 1.495356 |
| Robust 95% CI | [4.06975, 10.9833] | [4.069746, 10.983319] |
| Effective N (left / right) | 343 / 310 | 343 / 310 |
| Default bandwidth $h$ / $b$ (p. 392) | 17.708 / 27.984 | 17.519 / 27.695 |
| Estimate at the default bandwidth | 7.416 | 7.4218 |

The Senate data is not redistributed with rdlocal. `tests/test_published_senate.py` runs when the file is provided and skips otherwise.

## 11. Open problem: the d̂ pilot bandwidth

Default bandwidths differ slightly from other implementations: about 1% on the Senate data, and up to 20% on one simulated design. Every remaining gap traces to one step.

The supplement to Calonico, Cattaneo, Farrell & Titiunik (2019), §12.2 (pp. 67–68), defines the pilot

$$
\hat d = \left[\frac{3 + 2q}{2}\cdot\frac{V/n}{\mathring{B}^2}\right]^{1/(5+2q)},\qquad
\mathring{B} = \hat\mu^{(1+q)}_{S,1+q}(x_\pm)\cdot O_{1+q,1+q}(v),
$$

with $x_-$ and $x_+$ the third and first quartiles of $X$ below and above the cutoff. The text does not specify how the local order-$(1+q)$ quantity is formed at a quartile. Four natural readings were tried and none reproduces the published default bandwidths:

1. the quartile as a bandwidth;
2. the quartile as the centre of a local fit;
3. a global polynomial evaluated at the quartile;
4. side-percentile fits.

rdlocal currently uses a documented substitute (`DEVIATION` in `bandwidth.py`). The test `tests/test_bandwidth.py::test_all_selectors_within_1pct` is an expected failure that will pass once the step is resolved. **Contributions are welcome, provided they are derived from the published papers.**

Every reported number is correct for the bandwidth it reports. This affects only the *default* choice. For results comparable across software, pass $h$ and $b$ explicitly.

## 12. References

- Calonico, S., Cattaneo, M. D., & Titiunik, R. (2014). Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs. *Econometrica*, 82(6), 2295–2326, and Supplemental Appendix.
- Calonico, S., Cattaneo, M. D., & Titiunik, R. (2015). rdrobust: An R Package for Robust Nonparametric Inference in Regression-Discontinuity Designs. *The R Journal*, 7(1), 38–51.
- Calonico, S., Cattaneo, M. D., Farrell, M. H., & Titiunik, R. (2017). rdrobust: Software for regression-discontinuity designs. *Stata Journal*, 17(2), 372–404.
- Calonico, S., Cattaneo, M. D., Farrell, M. H., & Titiunik, R. (2019). Regression Discontinuity Designs Using Covariates. *Review of Economics and Statistics*, 101(3), 442–451, and Supplemental Appendix.
- Calonico, S., Cattaneo, M. D., & Farrell, M. H. (2020). Optimal bandwidth choice for robust bias-corrected inference in regression discontinuity designs. *Econometrics Journal*, 23(2), 192–210, and supplement.
- Cattaneo, M. D., Frandsen, B. R., & Titiunik, R. (2015). Randomization Inference in the Regression Discontinuity Design: An Application to Party Advantages in the U.S. Senate. *Journal of Causal Inference*, 3(1), 1–24.
- Cattaneo, M. D., Jansson, M., & Ma, X. (2018). Manipulation testing based on density discontinuity. *Stata Journal*, 18(1), 234–261.
- Cattaneo, M. D., Jansson, M., & Ma, X. (2020). Simple Local Polynomial Density Estimators. *Journal of the American Statistical Association*, 115(531), 1449–1455.
- Gelman, A., & Imbens, G. (2019). Why High-Order Polynomials Should Not Be Used in Regression Discontinuity Designs. *Journal of Business & Economic Statistics*, 37(3), 447–456.
- Imbens, G., & Kalyanaraman, K. (2012). Optimal Bandwidth Choice for the Regression Discontinuity Estimator. *Review of Economic Studies*, 79(3), 933–959.
- McCrary, J. (2008). Manipulation of the running variable in the regression discontinuity design: A density test. *Journal of Econometrics*, 142(2), 698–714.
