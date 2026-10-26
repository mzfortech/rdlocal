# Clean-room statement

`rdlocal` is an independent implementation written from published papers. It is not derived from, and
contains no code, docstrings or documentation text from, the GPL-3.0 packages `rdrobust`, `rddensity`,
`lpdensity`, `rdmulti`, `rdpower`, `rdlocrand` or `rdhte` (any language, any version).

## Development iteration 2: scope widening (owner-approved)

Published journal articles that describe the rdpackages algorithms may be read as specifications. Still
forbidden: any source code, help files (.sthlp/.Rd), docstrings, CRAN reference manuals, vignettes, READMEs
or repository files of the rdpackages. Articles added in iteration 2 are in the table below (CCFT2017, CJM2018 for
its bandwidth section, RJ2015). The Stata Journal 2014 rdrobust article (14(4):909-946) could not be
downloaded (the ageconsearch record that was tried, 340127, is a different article: Cattaneo et al. 2016 on
local randomization; its PDF is stored but unused).

**Result of the iteration 2 widening:** these articles state *which* bandwidth selectors exist and how they relate
(CER scaling, combination rules, the generic form h = (V/n / (2(1+p)B))^(1/(3+2p))) but defer the exact
preliminary (pilot) estimators and the regularisation scaling to the supplemental appendices of CCFT / CJM,
which are not obtainable. So the data-driven h, b of rdrobust and rddensity cannot be reproduced to 1%
from permitted text; see "Remaining deviations".

## Development iteration 3: the missing specification

The owner supplied the journal supplements above (hosted under rdpackages.github.io/references/; only those
files were fetched, nothing else on that site). CCFT2019-S Section 12 specifies the rdrobust-family
bandwidth selector: Step 1 (pilot v), Step 2 (pilot b, with a further pilot d), Step 3 (h), using finite-sample
bias constants O(v) and variance constants at v. It does **not** print the regularisation (the factor is
identified below) and its definition of the d-pilot bias constant (B-ring) is not recoverable. The Senate
vote-share data (CCFT 2017, Stata Journal) were used **locally for verification only**; they ship inside a GPL
package and are not redistributed or committed. `tests/test_published_senate.py` skips when absent.

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
| CCFT2019-S | Calonico, Cattaneo, Farrell, Titiunik (2019), *RD Designs Using Covariates*, REStat, **Supplemental Appendix** (87 pp.; Section 12 "Sharp RD Designs: Steps 1-3" pp. 65-70 is the bandwidth implementation spec) | https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT--Supplement.pdf | 53bc60883ac3e28b0e8cefc0b1d299468d56231f9bf0b4c1febfddcc4f91398c |
| CCF2020-S | Calonico, Cattaneo, Farrell (2020), Econometrics Journal, Supplement (read for CER pilots / residuals; no new constants for the MSE selectors) | https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ--Supplement.pdf | ff2b2bed480fa2210ee1a9f430253a84fe5c2d7f1be8d82c00ba2f5fe71269ed |
| CCT2014 | Calonico, Cattaneo, Titiunik (2014), Econometrica 82(6) main text | https://rdpackages.github.io/references/cct2014_ecma.pdf | 354d9010ca1f8c3bc98e090a30a498178422b95449a0773e9366757e089ac602 |
| CFT2015 | Cattaneo, Frandsen, Titiunik (2015), J. Causal Inference (copied with the others; unused) | https://rdpackages.github.io/references/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf | 5f4ca3cc38a9ce2993aac7fcfad1e22e4967e3841d624a931d89fc125b80eabb |
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

1. (resolved in 0.3) Regularisation scale. Section 12 prints no regularisation; black-box agreement fixes it:
   B^2 -> B^2 + 3 * scaleregul * Var(B), Var(B) = sum over sides of O^2 Var(Taylor coefficient) (coefficient
   scale, no nu!^2). Given rdrobust's own b, Step 3 then reproduces its h to <=0.3% on all sharp designs with
   scaleregul = 0 and 1. The 3 and the Var(B) form are an identification from numbers, not from a page.
2. (resolved in 0.3) Step 2/3 now follow Section 12 (finite-sample O and V at v; C_K = 2.576 / 1.843; b with
   exponent 1/(3+2q) and h with 1/(3+2p)). The b-stage bias is the order-(1+q) fit coefficient at d_hat with the
   same regularisation; that structure was confirmed numerically (B(d*) equals the b implied by the black box).
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
11. **d_hat pilot (Section 12.2, p.67-68)** - THE remaining cause of selector gaps. The text defines
    d_hat = [(3+2q)/2 V/n / B_ring^2]^(1/(5+2q)) with B_ring = mu_hat^(1+q)_{1+q}(x_pm) O_{1+q,1+q}(v), x_- /
    x_+ the third / first quartile of {X<c} / {X>=c}. How mu_hat^(1+q)(x_pm) (an order-(1+q) quantity "at" a
    quartile) is formed is not recoverable: tried the quartile as a bandwidth, as a window centre, as the point of
    evaluation of a global polynomial, at v, 2v, side percentiles and with/without regularisation; none matched
    the d_hat implied by the black box (implied from rdrobust's b under scaleregul 0 and 1). We use the CCT2014
    supplement S.2.6 Step-0 rule (global polynomial of degree q+2) with the finite-sample O of Section 12.
    Consequence: 14 of 56 (design x selector x scaleregul) cells are within 1%; the Senate default is at -1.07%
    (h) and -1.03% (b); worst cells (Lee-1, 2-bandwidth selectors, p=2) are 15%-70% off.
12. **Fuzzy selector.** The pilot tau_0 (and thus the pseudo-outcome) differs from rdrobust's: fuzzy h is 12% low
    with scaleregul=1 and 75% low with scaleregul=0 (Lee-fuzzy).

## Remaining deviations after development iteration 3

* rdrobust h, b: the d_hat pilot above (CCFT2019-S Section 12.2, p.68, display for B_ring). With it fixed, the rest
  of the chain (Steps 1 and 3, regularisation, combination and CER rules) is already at <=1%.
* rddensity h: not specified by either supplement (they do not discuss the density test); the iteration 2 rule stays.
* Fuzzy selector pilot (deviation 12); `msetwo`/`msesum` inherit the d_hat error.

## Choices informed by black-box agreement (all with a paper passage behind them)

* Robust-variance residuals taken from the order-q fit at b (CCF2020 S.4 says residual estimators have "the
  same form ... with only the bandwidth and polynomial order changed"); the alternative reading (order-p
  residuals everywhere) left a 0.15%-2% gap, the chosen one agrees to 1e-14.
* Fuzzy robust variance linearised at the conventional ratio (CCT2014-S proof of Thm A.2 linearises around
  zeta_hat_p(h)); the alternative (bias-corrected ratio) left a 0.7% gap.
* Leave-self-out empirical CDF (deviation 9); the plain form left a 0.1% gap in f_hat.
