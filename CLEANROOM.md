# Clean-room statement

`rdlocal` is an independent implementation written from published papers. It is not derived from, and
contains no code, docstrings or documentation text from, the GPL-3.0 packages `rdrobust`, `rddensity`,
`lpdensity`, `rdmulti`, `rdpower`, `rdlocrand` or `rdhte` (any language, any version).

## v0.2 scope widening (owner-approved)

Published journal articles that describe the rdpackages algorithms may be read as specifications. Still
forbidden: any source code, help files (.sthlp/.Rd), docstrings, CRAN reference manuals, vignettes, READMEs
or repository files of the rdpackages. Articles added in v0.2 are in the table below (CCFT2017, CJM2018 for
its bandwidth section, RJ2015). The Stata Journal 2014 rdrobust article (14(4):909-946) could not be
downloaded (the ageconsearch record that was tried, 340127, is a different article: Cattaneo et al. 2016 on
local randomization; its PDF is stored but unused).

**Result of the widening:** these articles state *which* bandwidth selectors exist and how they relate
(CER scaling, combination rules, the generic form h = (V/n / (2(1+p)B))^(1/(3+2p))) but defer the exact
preliminary (pilot) estimators and the regularisation scaling to the supplemental appendices of CCFT / CJM,
which are not obtainable. So the data-driven h, b of rdrobust and rddensity cannot be reproduced to 1%
from permitted text; see "Remaining deviations".

## Rules followed

1. No file of those packages was opened, grepped, viewed or copied: nothing under `site-packages/rd*` or
   `site-packages/lpdensity*`, no `inspect.getsource`, no `help()` / docstring reads, no web fetch of their
   GitHub repositories, CRAN pages or PyPI pages (one search result that was a CRAN manual of `rddensity`
   was NOT opened).
2. Allowed sources only: the papers below, general econometrics/statistics knowledge, numpy/scipy docs.
3. Black-box use: `rdrobust` 2.1.1 and `rddensity` 3.0 were executed as opaque functions on simulated data
   (`tests/make_reference.py`). Only the inputs (the simulated data in `tests/data/designs.npz` and the call
   arguments) and the returned NUMBERS (`tests/data/reference.json`) were recorded. The only introspection was
   `dir()` on a result object to learn the names of its numeric fields.
4. Design choices that were selected by comparing numbers with the black box are listed below under
   "Choices informed by black-box agreement"; each is also backed by a passage of a paper.

## Papers read (PDFs downloaded with curl and read from the PDF text / page images)

| Short name | Paper | URL | sha256 |
|---|---|---|---|
| CCT2014-S | Calonico, Cattaneo, Titiunik (2014), *Robust Nonparametric Confidence Intervals for RD Designs*, Econometrica 82(6), **Supplemental Appendix only** (the main article PDF could not be obtained from an open source) | https://mail.econometricsociety.org/publications/econometrica/2014/11/01/robust-nonparametric-confidence-intervals-regression/supp/ECTA11757SUPP.pdf | c6d9c512a72fe96debd8bec448c19ef0345d83b4fc4507c7452068edef62dab6 |
| CCF2020 | Calonico, Cattaneo, Farrell (2020), *Optimal bandwidth choice for robust bias-corrected inference in RD designs*, Econometrics Journal (arXiv version, includes supplement) | https://arxiv.org/pdf/1809.00236 | 17fc6a163ce2c3a7c87d33450452cc779f7f1014c1bdfc012fbd6d8f8e6f13ca |
| CCFT2019 | Calonico, Cattaneo, Farrell, Titiunik (2019), *RD Designs Using Covariates*, REStat (arXiv) | https://arxiv.org/pdf/1809.03904 | 6c5bb2dfbb52cc7b66c22ffee0b9e447ba59f786605d3155cf7854247a98895d |
| CJM2020 | Cattaneo, Jansson, Ma (2020), *Simple Local Polynomial Density Estimators*, JASA (arXiv version with supplement) | https://arxiv.org/pdf/1811.11512 | a7bba884a967a087e9eb1d298294fd0d1b27b931e03b8f16ba391eddd5da31a8 |
| CJM2018 | Cattaneo, Jansson, Ma (2018), *Manipulation testing based on density discontinuity*, Stata Journal 18(1) (downloaded; formulas taken from CJM2020, this paper is software-oriented) | https://ageconsearch.umn.edu/record/340234/files/Cattaneo.pdf | 9b7216e5943350224b814740d5934aefc8c12e193815764e389d61f079fafcd3 |
| IK2012 | Imbens, Kalyanaraman (2012), *Optimal bandwidth choice for the regression discontinuity estimator*, REStud (NBER WP 14726) | https://www.nber.org/system/files/working_papers/w14726/w14726.pdf | 1f86127ab69d51a07af5d36cef39415660bdfdfca63d36bfd52cbc7f7b29d60b |
| CCFT2017 | Calonico, Cattaneo, Farrell, Titiunik (2017), *rdrobust: Software for regression-discontinuity designs*, Stata Journal 17(2):372-404 (sec. 2.3 bandwidth selectors; sec. 7 published Senate output) | https://ageconsearch.umn.edu/record/340181/files/Calonico.pdf | c41e81d2b3763ea5e2ac060747f7c6a22cf6ae3f7fd27d44dfece82e2fe63c8a |
| RJ2015 | Calonico, Cattaneo, Titiunik (2015), *rdrobust: An R Package for Robust Nonparametric Inference in RD Designs*, The R Journal 7(1) | https://journal.r-project.org/archive/2015/RJ-2015-004/RJ-2015-004.pdf | b43c6d0f7188e5bcb56f92c99328a9a786ad1f2ef79f9d4ebc8ac6ba1b93b637 |
| (unused) | Cattaneo, Titiunik, Vazquez-Bare (2016), Stata Journal 16(2) - downloaded by mistake while looking for the 2014 rdrobust article | https://ageconsearch.umn.edu/record/340127/files/Cattaneo.pdf | 957f0118bfcf2559d5ced1ef31d254609086abfbd523a252ff3d38cb1f18e498 |
| CIT2019 | Cattaneo, Idrobo, Titiunik (2019), *A Practical Introduction to RD Designs: Foundations* (arXiv) - consulted for context and for a data source; no formula taken | https://arxiv.org/pdf/1911.09511 | a4f103e8518fb871503c22cceb179a4f614de43edd6d85236a489022b4c67ed3 |
| Lee2008 | Lee (2008), J. Econometrics - downloaded, only for the data search; not used | https://www.princeton.edu/~davidlee/wp/randomizerdfinal2.pdf | 49c81381c74d9560e5429ba8e279b30f9bd22b2ed3174aa02b0a3ff7f37876d9 |

Not obtained: the CCT2014 main article, McCrary (2008) (so `rd_density_test` implements only the CJM test;
the McCrary option is not provided), Lee's replication files (see "Published examples" in the README).

## Function -> source map

| Function | Source |
|---|---|
| `estimate.rd_estimate`, conventional estimator | CCT2014-S Lemma 1 / S.2.1: order-p local polynomial, intercept |
| bias-corrected linear weights (`Wb`) | CCT2014-S Thm A.1 bias expansion; CCF2020 S.5.3 equivalent-kernel display: tau_bc = tau - h^(p+1) B, B from the order-q fit at b (finite-sample form with the actual Gram matrices) |
| fuzzy ratio and linearised bias | CCT2014-S proof of Thm A.2: zeta_bc = zeta_p(h) - h^(p+1-nu) B_F, B_F = bias_Y/tau_T - tau_Y bias_T/tau_T^2; variance V_F in Lemma 2 (equals the variance of Y - tau T over tau_T) |
| `estimate._variance`, residual estimator | CCT2014-S S.2.4 (plug-in residuals; NN estimator with J neighbours J/(J+1)(Y_i - mean)^2); CCF2020 S.4 (residuals for the bc statistic: same form, bandwidth and order changed); CCF2020 S.5.2 (HC0-HC3, NN) |
| `bandwidth._cct` Step 0-2 | CCT2014-S S.2.6, "Plug-in Bandwidth Selectors", Steps 0, 1, 2 (v_n, c_n, b_n, h_n) |
| `_kernels.lead_const` | B_{nu,p} = nu! e_nu' Gamma_p^-1 vartheta_{p,p+1}, CCT2014-S S.2.6 / Lemma 1 |
| `_kernels.silverman_const` | CCT2014-S S.2.6 Step 0 (i): 2.58 (triangular); general formula printed there |
| `bandwidth.deriv_var` | V_hat_{nu,p}(h), CCT2014-S S.2.6 |
| `bandwidth._ik` | IK2012 Section 4.2, eqs (4.8)-(4.13), C_K = 3.4375 (edge kernel) |
| `bandwidth.rd_bandwidth` variants (`rd`, `two`, `sum`, `comb1`, `comb2`, `cer*`) | CCFT2017 sec. 2.3 (estimands rd/l/r/sum, h_cer = n^(-p/((3+p)(3+2p))) h_mse, comb rules; identical to CCF2020 sec. 5, G replaces n with clusters). Published Senate output on p.400 (n=1297: mserd 17.708 -> cerrd 12.374; b unchanged by CER; b side-specific for msetwo) is checked in tests |
| `density._side` | CJM2020 Sec. 2 (estimator), Thm 2 (variance estimator, evaluated as (1/n) sum_jk c_j c_k [F(x_j ^ x_k) - F_j F_k]) |
| `density.rd_density_test` statistic | CJM2020 Sec. 4 (T_p, Corollary 1, robust T_{p+1}(h_p)) |
| `density._h_mse` | CJM2020 Sec. 3.1 (h_MSE = (V/(2pB^2))^(1/(1+2p)) n^(-1/(1+2p))); see deviations |
| `_kernels.min_gamma` | CJM2020 Thm 1, matrix B(x) (integral of min{u,v}) |
| `binned.rd_binned` | generic binned means (no paper formula) |

## DEVIATIONS (where the obtainable papers do not state, or state ambiguously, what was implemented)

1. **Regularisation scale in CCT Steps 1-2.** The printed bracket adds `3 V_hat` (a variance of a
   *derivative* estimator, carrying `nu!^2`) to a squared Taylor coefficient. We use `3 V_hat / nu!^2`
   (coefficient scale) so the two terms have the same units. Used verbatim, bandwidths are far smaller
   (Lee-1 design, n=1000: h = 0.073 vs 0.110 with the coefficient-scale version; black-box value 0.114).
2. **Step 2 order.** The print writes `beta_{p+1}(b_n)` (an order p+1 fit); we use the order-q fit at b_n,
   consistent with `V_hat_{p+1,q}(b_n)` in the same bracket. Identical for the default q = p+1.
3. **`msetwo`/`msesum`.** The selectors exist per CCFT2017 sec. 2.3 and the b-per-side structure is confirmed by
   its printed table, but the formulas are ours: the same Step 1-2 expressions with the one-sided bias terms
   (two) or the sum of one-sided bias terms (sum) replacing the difference.
4. (resolved in 0.2) `cerrd`: now exactly the published rule n^(-p/((3+p)(3+2p))) h_mse with b unchanged; only
   the arg-min over Edgeworth terms of CCF2020 is not implemented (not part of rdrobust's `cerrd` per CCFT2017).
5. **Fuzzy bandwidth.** The MSE problem for the ratio is solved as the sharp problem for the pseudo-outcome
   `Y - tau0 T` (CCT2014-S Lemma 2 shows bias and variance take this form), with `tau0` from local fits at the
   Silverman pilot bandwidth. The pilot choice is ours.
6. **HC0-HC3 and cluster corrections.** CCT2014-S/CCF2020 only cite the HCk family (MacKinnon 2013) and the
   details sit in a supplement we could not obtain. We use: HC0 = 1, HC1 = n_s/(n_s-k) per side,
   HC2 = 1/(1-h_ii), HC3 = 1/(1-h_ii)^2 with the weighted hat diagonal of the fit producing the residual, and
   cluster-robust = sum of squared cluster totals x G/(G-1) x (n-1)/(n-k). Observed consequence: the
   *conventional* SE under `hc1` differs from the black box by 0.15%-1%; the other HC variants agree to 1e-14.
7. **IK**: implemented for sharp designs, p=1, triangular kernel only (the paper gives C_K for the edge kernel).
8. **Density test bandwidth (`h=None`).** CJM2020's data-driven selector for the *difference* is in an
   unobtained part of the supplement. We implement Sec. 3.1's formula with F^{(p+1)} from a global order-(p+3)
   polynomial of the empirical CDF per side (CCF2020 also uses a global fit for derivatives, order p+4 there)
   and a Silverman pilot window for the density; one common bandwidth for both sides. Our h differs from
   `rddensity`'s default (0.23 vs 0.32/0.27 on Lee-1).
9. **Empirical CDF.** F(x_i) = (1/(n-1)) sum_{j != i} 1(x_j <= x_i) (the leave-in-free form shown in CJM2020 SA
   Sec. 5.2) rather than including j = i (Sec. 2). Effect O(1/n) on f_hat.
10. **Density restricted model, jackknife and plug-in SEs** are not implemented (unrestricted model, Thm 2 SE).
11. **Pilot c_n** uses a global polynomial of degree q+2 (the print says "(m+1)th-order" for a vector
    r_m = (1,...,x^m); we read it as degree m = q+2). Trying degree q+3 made agreement with the black box worse.

## Remaining deviations after v0.2 (why 1% on h, b is not reached)

* **rdrobust h, b.** Measured gap vs the black box on the five designs: h 6%-31% low, b from 32% low to 4.5% high (see README
  table); a grid of 36 plausible readings of the printed procedure (pilot polynomial degree, c_n rule, 0/1/3
  regularisation scaling with and without the nu!^2 factor) never got below 28% max error, and feeding the
  black box's own b into our Step 2 leaves 0.2%-8% on h. CCFT2017 sec. 2.3 (pp.379-381) says the upgraded
  selector's preliminary estimates and its regularisation "are discussed in the supplemental appendix of
  CCFT"; RJ2015 pp.40-41 only gives h_MSE = C n^(-1/(2p+3)) and defers to the supplements of IK and CCT.
  That supplement is the missing specification.
* **rddensity h.** CJM2018 sec. 2.5 specifies the family (h_l, h_r, h_diff, h_sum, comb = median / min) but
  says the constants "rely on preliminary rule-of-thumb estimators" without stating them; CJM2020 sec. 3.1 only
  says F^(p+1) can come from a local polynomial or "a reference model". A normal-reference reading of that
  family gave 0.24-0.35 vs the black box's 0.27-0.33 (side-wise) and 0.59-0.96 (pooled): not within 1%.
  Our default stays the single common h of deviation 8.

## Choices informed by black-box agreement (all with a paper passage behind them)

* Robust-variance residuals taken from the order-q fit at b (CCF2020 S.4 says residual estimators have "the
  same form ... with only the bandwidth and polynomial order changed"); the alternative reading (order-p
  residuals everywhere) left a 0.15%-2% gap, the chosen one agrees to 1e-14.
* Fuzzy robust variance linearised at the conventional ratio (CCT2014-S proof of Thm A.2 linearises around
  zeta_hat_p(h)); the alternative (bias-corrected ratio) left a 0.7% gap.
* Leave-self-out empirical CDF (deviation 9); the plain form left a 0.1% gap in f_hat.
