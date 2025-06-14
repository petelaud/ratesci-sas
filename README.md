# ratesci-sas

### Confidence intervals and tests for comparison of rates

ratesci-sas contains SAS macro code (%SCORECI) to compute score confidence intervals for 
rate (or risk) difference ('RD'), rate ratio ('RR', also known as relative risk), 
or odds ratio ('OR') for binomial proportions, and RD or RR for Poisson 
'exposure-adjusted' incidence  rates. Matching hypothesis tests are produced, 
including a 2-sided test for association (corresponding to the 'N-1' Chi-squared test), 
and 1-sided non-inferiority tests against any specified null value 'DELTA' (an 
'N-1' modified variant of a Farrington-Manning test). Both tests are guaranteed to be
consistent with the confidence interval.

A second macro (%PAIRBINCI) implements score-based confidence intervals and tests
for paired contrasts of binomial proportions.

Stratified calculations are catered for in %SCORECI with a range of weighting schemes, 
with direct equivalence to the Cochran-Mantel-Haenszel (CMH) test when comparing 
RD or RR with `WEIGHT=1` for MH weights, or OR with 'WEIGHT=3' for INV weights. 
The stratified analysis also produces a homogeneity test.

Note that SAS (since v9.3M2 / STAT 12.1) PROC FREQ will produce the Miettinen-Nurminen 
('MN') score interval (for binomial RD, RR or OR) for unstratified datasets only 
(and fails to produce results if there are no events or no non-events). 
The "Summary Score confidence limits" produced for a stratified analysis, e.g.
 `TABLES ... / CMH COMMONRISKDIFF(CL=SCORE TEST=SCORE);`
are *not* stratified MN intervals, and consequently can conflict with the result 
of the CMH test. 
More recently, the SAS Viya platform (version LTS 2023.10, November 2023) added 
a COMMONRISKDIFF(CL=MN) option to the TABLES statement, which gives a stratified 
Miettinen-Nurminen interval (for RD only) with Miettinen & Nurminen's iterative 
weights (or CL=MNMH for MH weighting). Those weights are not yet implemented in 
the %SCORECI macro, but are available in the ratesci package for R.

For unstratified analysis, the test and interval obtained from SAS PROC FREQ with
`TABLES ... / RISKDIFF(CL=MN EQUAL METHOD=SCORE);` 
are incoherent, because the Farrington-Manning score test omits the variance 
bias correction factor `N/(N-1)` included in the MN interval. 
This correction should be included in the test to avoid inflated type 1 error rates.

[Aside: Note that when applying one-sided tests with PROC FREQ, the confidence 
intervals in the PdiffNonInf, PdiffSup and PdiffEquiv output tables 
(with METHOD=SCORE) bear no relationship to the MN intervals in the PdiffCIs dataset, 
and actually change depending on the MARGIN value provided, which makes no sense. 
Also, only positive values of MARGIN are allowed, which means the superiority test 
result has to be inverted. Then there's the fact that the whole concept of "inferior" 
vs "superior" depends on whether the endpoint is a positive event (e.g. response rate) 
or a negative one (e.g. death). The logical solution is to use the same score statistic 
in deriving both the test and the confidence interval, output both left- and 
right-sided p-values, and let the user choose which one is relevant for their purposes.]

In addition to addressing the above issues, the %SCORECI macro incorporates 
skewness-corrected asymptotic score ('SCAS') methods, which ensure improved 
equal-tailed coverage (or central location), in other words for a nominal 95% 
confidence interval, the one-sided non-coverage probability is (on average) close 
to 2.5% on each side. 
 
Omission of the skewness correction results in the often-recommended 
'Miettinen-Nurminen' asymptotic score confidence interval, which can have inferior 
one-sided coverage, especially for the RR contrast and for unbalanced designs. 

For further details, please refer to the documentation of the ratesci package
[here](https://petelaud.github.io/ratesci)
