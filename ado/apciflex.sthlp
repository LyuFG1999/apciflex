{smcl}
{title:apciflex 1.0.0}

{p 4 4 2}
{cmd:apciflex} implements flexible-resolution Age-Period-Cohort Interaction
(APC-I) models on irregular and/or incomplete observed Lexis support. Version
1.0.0 supports Gaussian identity-link models and binary binomial models with
logit or probit links. The command estimates on the supplied raw Age x Period
support and applies Age, Period, and optional cohort grouping only after
estimation.
{p_end}

{title:Syntax}

{p 8 8 2}
{cmd:apciflex} {it:yvar} {ifin} [{it:weight}],
{cmd:age(}{it:agevar}{cmd:)}
{cmd:period(}{it:periodvar}{cmd:)}
{cmd:agegroup(}{it:agegroupvar}{cmd:)}
{cmd:periodgroup(}{it:periodgroupvar}{cmd:)}
[{cmd:cohortgroup(}{it:cohortgroupvar}{cmd:)}
 {cmd:controls(}{it:varlist}{cmd:)}
 {cmd:mincell(}{it:#}{cmd:)}
 {cmd:aggweight(}{it:equal|sample|effective|population}{cmd:)}
 {cmd:mtest(}{it:none|bonferroni|holm|bh}{cmd:)}
 {cmd:link(}{it:identity|logit|probit}{cmd:)}
 {cmd:vce(}{it:vcetype}{cmd:)}
 {cmd:detail}
 {cmd:nodisplay}
 {cmd:level(}{it:#}{cmd:)}]
{p_end}

{p 8 8 2}
Replay:
{p_end}
{p 12 12 2}{cmd:apciflex}{p_end}
{p 12 12 2}{cmd:apciflex, detail}{p_end}

{title:Core estimand}

{p 4 4 2}
Only retained observed raw Age x Period cells are estimated. Structural empty
cells are never imputed, interpolated, or treated as zeros. Raw cohort is
generated internally as period-age.
{p_end}

{p 4 4 2}
For retained support Omega, APCIFLEX first estimates the adjusted AP-cell
linear-predictor surface eta_ap and decomposes it as
{p_end}

{p 8 8 2}
eta_ap = eta_bar + alpha_a + pi_p + theta_ap.
{p_end}

{p 4 4 2}
With {cmd:link(identity)}, eta_ap is the adjusted outcome mean and the
decomposition is on the outcome scale. With {cmd:link(logit)}, eta_ap is the
log-odds. With {cmd:link(probit)}, eta_ap is the probit index. Thus the APC-I
components are always additive on the model's link scale.
{p_end}

{p 4 4 2}
Version 1.0.0 uses the equal-cell observed-support additive projection. The
degree-weighted zero-sum constraints preserve eta_bar as the equal-cell mean of
retained adjusted AP-cell linear predictors. The decomposition is an APC
surface estimand; it is not automatically a sample-composition or population
mean.
{p_end}

{title:Raw variables and grouping variables}

{p 4 4 2}
{cmd:age()} and {cmd:period()} define the raw estimation resolution. They must
be nonnegative integer-valued variables so the posted e(b) can carry stable
factor-variable metadata for Stata postestimation.
{p_end}

{p 4 4 2}
{cmd:agegroup()} and {cmd:periodgroup()} are required deterministic mappings
from raw Age and Period to theory-guided interpretation groups.
{cmd:cohortgroup()} is an optional deterministic mapping from the internally
defined raw cohort period-age to a theory-guided cohort group.
{p_end}

{title:link()}

{p 4 4 2}
{cmd:link()} selects the outcome model. Default: {cmd:link(identity)}.
{p_end}

{p 8 8 2}
{cmd:identity}: Gaussian linear model fitted by {cmd:regress}. APC-I components,
cell means, cohort deviations, and slopes are on the outcome scale. Inference
uses t statistics and F tests.
{p_end}

{p 8 8 2}
{cmd:logit}: Bernoulli/binomial GLM with logit link. The dependent variable
must be coded 0/1. APC-I components are on the log-odds scale. Pointwise
inference uses normal z statistics and joint interaction tests use Wald
chi-squared statistics.
{p_end}

{p 8 8 2}
{cmd:probit}: Bernoulli/binomial GLM with probit link. The dependent variable
must be coded 0/1. APC-I components are on the probit-index scale. Pointwise
inference uses normal z statistics and joint interaction tests use Wald
chi-squared statistics.
{p_end}

{p 4 4 2}
For logit/probit, each retained raw Age x Period cell must contain both outcome
values 0 and 1. A saturated AP-cell intercept for a cell containing only zeros
or only ones is not finite on the binomial link scale; APCIFLEX therefore
rejects such support before fitting rather than reporting an unstable
decomposition.
{p_end}

{p 4 4 2}
The inverse link is used only for prediction and margins. In particular,
invlink(eta_bar) is not the same estimand as the average of cell-level predicted
probabilities. APC-I components should therefore not be inverse-linked and
reinterpreted as additive probability components.
{p_end}

{title:Weights and mincell()}

{p 4 4 2}
The first-stage fitted model honors the supplied Stata weight. The interpretation-
stage aggregation rule is controlled separately by {cmd:aggweight()}.
{p_end}

{p 4 4 2}
For {cmd:pweight}, {cmd:mincell(#)} uses scale-invariant Kish effective sample
size within each raw AP cell:
{p_end}

{p 8 8 2}
n_eff = (sum w)^2 / sum(w^2).
{p_end}

{p 4 4 2}
For {cmd:fweight}, {cmd:mincell(#)} uses frequency-weighted N, equal to the
sum of fweights within the raw AP cell.  A compressed cell represented by
three data rows with fweights summing to 300 therefore has support 300, not 3.
For {cmd:aweight}, {cmd:iweight}, and unweighted data, {cmd:mincell(#)} uses
the number of physical data rows.
{p_end}

{p 4 4 2}
{cmd:e(cell)} stores physical row N, Kish effective N, and weight sum so the
support decision is inspectable.  {cmd:e(mincell_metric_code)} is
{cmd:kish_neff}, {cmd:fw_mass}, or {cmd:raw_n}, as appropriate.
{p_end}

{title:aggweight()}

{p 4 4 2}
{cmd:aggweight()} changes interpretation-stage averages only. It does not
change the retained AP-cell fitted model, the observed-support APC-I decomposition,
the raw interaction surface, the global interaction test, or the actual-age
cohort slope. Default: {cmd:aggweight(equal)}.
{p_end}

{p 4 4 2}
{cmd:equal}: equal temporal-unit aggregation. Raw Ages are equally weighted
within an Age group; observed raw Period occasions are equally weighted within
a Period group; AP cells are equally weighted within a raw cohort; theory
cohort groups give each retained raw cohort equal total weight after within-
cohort averaging.
{p_end}

{p 4 4 2}
{cmd:sample}: weights interpretation-stage units by retained physical data-row
counts. With fweights, use {cmd:population} when frequency mass is the desired
interpretation-stage weight.
{p_end}

{p 4 4 2}
{cmd:effective}: weights interpretation-stage units by Kish effective sample
size. With pweights, pooled effective N over disjoint AP cells is reconstructed
as
{p_end}

{p 8 8 2}
(sum_c W_c)^2 / sum_c(W_c^2 / n_eff,c),
{p_end}

{p 4 4 2}
where W_c is the AP-cell weight sum. Without weights, effective reduces to
sample weighting. Version 1.0.0 supports effective aggregation with pweights
or unweighted data.
{p_end}

{p 4 4 2}
{cmd:population}: weights by retained survey/frequency weight mass, sum(w).
With no weight it reduces to sample weighting. This mode is supported with
pweights, fweights, or unweighted data.
{p_end}

{p 4 4 2}
All grouped standard errors use the full covariance matrix through w' V w.
Standard errors are never averaged directly.
{p_end}

{title:Theory-group average and joint tests}

{p 4 4 2}
Every theory-guided Age, Period, and optional cohort group reports two
complementary tests.
{p_end}

{p 8 8 2}
The average-effect test evaluates H0: w'b_g = 0, where w is the selected
{cmd:aggweight()} aggregation vector. This is the existing group estimate,
standard error, and {cmd:Pavg}.
{p_end}

{p 8 8 2}
The joint existence test evaluates whether all underlying raw components in
the theory group are zero. For an Age group this is H0: alpha_a = 0 for every
retained raw Age a in the group; for a Period group it is H0: pi_p = 0 for
every retained raw Period p in the group.
{p_end}

{p 4 4 2}
For a cohort group, APCIFLEX first defines each retained raw cohort's average
interaction deviation Delta_c using the selected within-cohort
{cmd:aggweight()} rule, then tests H0: Delta_c = 0 for every raw cohort c in
the theory group. It deliberately does not jointly test every theta_ap cell,
because that stronger hypothesis would mix average cohort deviation with
within-cohort life-course dynamics.
{p_end}

{p 4 4 2}
Joint statistics use the full covariance matrix and rank-based degrees of
freedom. Under {cmd:link(identity)} the reported statistic is F; under
{cmd:link(logit)} or {cmd:link(probit)} it is Wald chi-squared. If
{cmd:mtest()} is not none, displayed group {cmd:Pjoint} values are adjusted
separately within Age-group, Period-group, and Cohort-group families.
{p_end}

{title:Cohort diagonal analysis}

{p 4 4 2}
For each raw cohort, APCIFLEX reports a joint nonadditivity test, an average
interaction deviation, and, when at least two retained cells exist, a linear
life-course slope. The slope is the OLS contrast of theta on actual raw Age,
not an ordinal cell-position contrast. Its unit is one link-scale interaction
deviation per one unit of Age: outcome units for identity, log-odds for logit,
and probit-index units for probit.
{p_end}

{title:mtest()}

{p 4 4 2}
{cmd:mtest()} controls multiplicity adjustment for the raw cohort diagonal
table and for the new theory-group joint tests. Default: {cmd:mtest(none)}.
{p_end}

{p 4 4 2}
Supported methods are {cmd:none}, {cmd:bonferroni}, {cmd:holm}, and
{cmd:bh} (Benjamini-Hochberg false-discovery-rate adjustment).
{p_end}

{p 4 4 2}
For the raw cohort diagonal table, adjustment remains separate within three
substantively distinct hypothesis families:
{p_end}
{p 8 8 2}1. cohort joint-deviation tests;{p_end}
{p 8 8 2}2. cohort average-deviation tests;{p_end}
{p 8 8 2}3. cohort slope tests among cohorts for which a slope is defined.{p_end}

{p 4 4 2}
Theory-group joint tests form three additional and mutually separate families:
Age-group joint tests, Period-group joint tests, and Cohort-group joint tests.
Thus Age, Period, and cohort-group hypotheses are never pooled into one
multiplicity family. Group-average p-values remain pointwise/unadjusted.
{p_end}

{p 4 4 2}
Raw p-values remain stored together with their adjusted counterparts in the
relevant e() matrices. {cmd:mtest()} does not alter the global Age x Period
interaction test.
{p_end}

{title:vce() and clustered samples}

{p 4 4 2}
{cmd:vce()} is passed to the underlying {cmd:regress} or binomial
{cmd:glm}. Version 1.0.0 explicitly supports and validates
{cmd:vce(cluster clustvar)}. Observations
with missing cluster identifiers are excluded before AP support, cell counts,
{cmd:mincell()}, and {cmd:e(sample)} are constructed. This keeps the support
reported by APCIFLEX aligned with the actual fitted-model sample.
{p_end}

{p 4 4 2}
When clustering is used, {cmd:e(clustvar)} and {cmd:e(N_clust)} are reposted.
{p_end}

{title:Replay}

{p 4 4 2}
After estimation, typing {cmd:apciflex} replays the compact results and
{cmd:apciflex, detail} replays the raw-resolution tables as well. Replay does
not refit the model. To keep replay independent of the current dataset, replay
uses the numeric group codes stored in e() rather than relying on current value
labels.
{p_end}

{title:predict}

{p 4 4 2}
{cmd:predict} after {cmd:apciflex} supports:
{p_end}
{p 8 8 2}{cmd:predict newvar, xb} -- linear predictor eta on the link scale;{p_end}
{p 8 8 2}{cmd:predict newvar, stdp} -- standard error of eta;{p_end}
{p 8 8 2}{cmd:predict newvar, pr} -- Pr(y=1), after logit/probit only;{p_end}
{p 8 8 2}{cmd:predict newvar, residuals} -- raw response residual.{p_end}

{p 4 4 2}
With {cmd:link(identity)}, the default prediction is {cmd:xb} and residuals are
y-xb. With {cmd:link(logit)} or {cmd:link(probit)}, the default prediction is
{cmd:pr} and residuals are y-Pr(y=1).
{p_end}

{p 4 4 2}
For observation i, {cmd:stdp} is always the standard error of the linear
predictor and is computed as sqrt(z_i' e(V) z_i), where z_i contains the
retained AP-cell indicator and all additive control values. Cell-control and
control-control covariance terms are included. Predictions are missing for
structural empty cells or cells removed by {cmd:mincell()}; APCIFLEX does not
extrapolate such cells.
{p_end}

{title:margins: predictive estimands versus APC-I components}

{p 4 4 2}
This distinction is essential. {cmd:margins} operates on the fitted AP-cell
prediction surface. Its results are predictive margins under Stata's chosen
standardization. They are not the APC-I decomposition components alpha_a,
pi_p, or cohort deviations.
{p_end}

{p 4 4 2}
For APC-I components use {cmd:e(age_raw)}, {cmd:e(age_group)},
{cmd:e(period_raw)}, {cmd:e(period_group)}, {cmd:e(cohort_raw)}, and optional
{cmd:e(cohort_group)}. These quantities are on the link scale selected by
{cmd:link()}.
{p_end}

{p 4 4 2}
After {cmd:link(identity)}, the default margins response is the linear
prediction. After {cmd:link(logit)} or {cmd:link(probit)}, the default margins
response is Pr(y=1). Thus, for a binary model,
{p_end}
{p 8 8 2}{cmd:margins, over(age period)}{p_end}
{p 4 4 2}
reports response-scale predictive probabilities, whereas
{p_end}
{p 8 8 2}{cmd:margins, over(age period) predict(xb)}{p_end}
{p 4 4 2}
reports predictive margins of the link-scale linear predictor.
{p_end}

{p 4 4 2}
With controls z1 and z2,
{p_end}
{p 8 8 2}{cmd:margins, over(age period) at(z1=0 z2=0) predict(xb)}{p_end}
{p 4 4 2}
reproduces the adjusted retained AP-cell linear predictors stored in
{cmd:e(cell)}. Under logit/probit,
{p_end}
{p 8 8 2}{cmd:margins, over(age period) at(z1=0 z2=0)}{p_end}
{p 4 4 2}
returns the inverse-link transformation of those cell predictors.
{p_end}

{p 4 4 2}
Commands such as {cmd:margins age}, {cmd:margins period}, and
{cmd:margins, dydx(x)} are response-scale predictive estimands under
logit/probit. They must not be reported as APC-I Age/Period/cohort components.
On incomplete support, some margins can be partly or fully nonestimable because
they require AP combinations not represented in the posted cell surface.
For example, {cmd:margins age} standardizes predictions over periods and can
request Age#Period combinations outside the retained support. Such predictive
margins are distinct from the equal-cell APC-I Age component in
{cmd:e(age_raw)}. {cmd:emptycells(reweight)} can change the predictive-margin
estimand and should not be used as a substitute for the APC-I component.
{p_end}

{p 4 4 2}
Because inverse links are nonlinear,
{p_end}
{p 8 8 2}
g^{-1}(eta_bar + alpha_a + pi_p + theta_ap)
{p_end}
{p 4 4 2}
cannot be rewritten as an additive decomposition of separate probability-scale
Age, Period, and cohort components. Probability-scale interpretation should
therefore be conducted with {cmd:predict, pr} and {cmd:margins}, while APC-I
decomposition and hypothesis tests remain on the link scale.
{p_end}

{title:Underlying model diagnostics}

{p 4 4 2}
APCIFLEX normally runs the saturated retained-cell {cmd:regress} or binomial
{cmd:glm} quietly. If the fitted model fails, version 1.0.0 reruns the same
model noisily so the original Stata diagnostic is displayed before the command
exits. Successful fits do not print the underlying model table.
{p_end}

{title:Output}

{p 4 4 2}
Default output reports support/model information, the global nonadditivity
test, and theory-guided Age, Period, and optional cohort groups. Each theory
group includes both an average-effect test and a within-group joint existence
test.
{p_end}

{p 4 4 2}
{cmd:detail} additionally reports raw Age components, raw Period components,
and raw cohort diagonal inference. If {cmd:mtest()} is not none, the displayed
cohort p-values are adjusted values; both raw and adjusted values remain stored.
{p_end}





{title:Prediction support}

{pstd}
Predictions are defined only for retained Age-by-Period cells.  If a requested
observation has nonmissing {cmd:age()} and {cmd:period()} values whose
combination is outside the retained support, {cmd:predict} exits with
{cmd:r(459)}.  This includes structural empty cells and cells excluded by
{cmd:mincell()}.

{pstd}
This rule also applies when {cmd:margins} constructs counterfactual
Age-by-Period combinations.  Thus broad commands such as {cmd:margins age} or
{cmd:margins period} may return {cmd:r(459)} on incomplete support, including
when {cmd:nose} is specified.  The explicit error prevents silent averaging
over only the subset of combinations that happen to be estimable.

{pstd}
For APC-I Age and Period components, use {cmd:e(age_raw)},
{cmd:e(age_group)}, {cmd:e(period_raw)}, and {cmd:e(period_group)}.
For predictive summaries over retained cells, use specifications such as
{cmd:margins, over(age period)}.  Manual predictions after {cmd:mincell()}
should normally be restricted with {cmd:if e(sample)} unless the requested
Age-by-Period combinations are known to be retained.


{title:Margins performance}


{pstd}
APCIFLEX uses a support-adaptive derivative policy.  Complete retained AP
support uses the faster chain-rule path.  Incomplete retained AP support uses {cmd:apcxb}/{cmd:apcpr}; raw
{cmd:xb}/{cmd:pr} are not available to {cmd:margins} there.  This routes
unsupported counterfactual combinations through the APCIFLEX prediction
support guard.


{pstd}
Version 1.0.0 no longer sets {cmd:e(marginsprop)} to {cmd:nochainrule}.
APCIFLEX predictions have the standard single-index form
{it:eta}=Xb for {cmd:link(identity)}, {cmd:link(logit)}, and
{cmd:link(probit)}.  The posted {cmd:e(b)} uses factor-variable-compatible
raw Age#Period coefficient stripes plus additive controls.  This allows
{cmd:margins} to use its chain-rule derivative path instead of numerically
perturbing every coefficient separately.

{pstd}
The prediction engine is also vectorized in Mata.  A call to
{cmd:predict} no longer executes one Stata {cmd:replace} scan for every
retained Age#Period cell.  The same engine computes {cmd:xb},
{cmd:stdp}, {cmd:pr}, and {cmd:residuals}.  In particular,
{cmd:stdp} evaluates the full quadratic form z_i' e(V) z_i in Mata,
including cell-control and control-control covariance terms.

{pstd}
These changes are computational only.  They do not change the APC-I
estimands, the posted covariance matrix, the definition of structural or
{cmd:mincell()}-excluded cells, or the distinction between predictive
margins and APC-I Age/Period/Cohort components.

{pstd}
If standard errors are not required, {cmd:margins, nose} can be faster for
predictive estimands that remain entirely inside the retained AP support, such
as observed-cell margins. On incomplete support, do not use {cmd:nose} to
bypass nonestimability diagnostics for broad commands such as
{cmd:margins age} or {cmd:margins period}; those commands may require
unobserved Age#Period combinations.



{title:Compact display}

{p 4 4 2}
Version 1.0.0 uses compact result tables designed for a 100-character Stata
line width.  Group labels are left aligned, numerical entries are displayed
with at most four decimal places and trailing zeros are suppressed, and very
small p-values are displayed as {cmd:<.0001}.  When {cmd:controls()} is
specified, their coefficients are displayed directly from the underlying
model together with their pointwise tests and a joint control-variable test.
The stored numerical results are not rounded; only the printed display is
compact.
{p_end}


{p 4 4 2}
Theory-group output is printed as two compact tables for each Age, Period, or
Cohort grouping.  The first table reports the group average, its SE, pointwise
t/z test, p-value, and confidence interval.  The second table reports the
joint F/Wald test that all raw components within the group are zero.  This
separation keeps the display within a 100-character line width and also makes
the distinction between {cmd:Pavg} and {cmd:Pjoint} explicit.
{p_end}


{title:Stored results}

{p 4 4 2}
APCIFLEX is an e-class estimation command and posts {cmd:e(b)}, {cmd:e(V)}, and
{cmd:e(sample)}. The AP-cell coefficients use factor-variable-compatible raw
Age#Period stripes followed by additive control coefficients.
{p_end}



{p 4 4 2}
Sample-size scalars distinguish Stata estimation N from physical data rows.
{cmd:e(N)} is the N posted by the underlying estimation command; with
{cmd:fweight} it equals the sum of frequency weights.  {cmd:e(N_rows)} is the
number of retained physical rows marked by {cmd:e(sample)}.
{cmd:e(dropped_rows)} is the number of physical rows removed by
{cmd:mincell()}; {cmd:e(dropped_observations)} is retained as a backward-
compatible alias for the same physical-row count.
{p_end}
{p 4 4 2}
Primary matrices are {cmd:e(age_raw)}, {cmd:e(period_raw)}, {cmd:e(cell)},
{cmd:e(cohort_raw)}, {cmd:e(age_group)}, {cmd:e(period_group)}, optional
{cmd:e(cohort_group)}, {cmd:e(V_cellmean)}, {cmd:e(V_age_raw)},
{cmd:e(V_period_raw)}, {cmd:e(V_interaction)}, {cmd:e(global_test)},
{cmd:e(control_table)}, and {cmd:e(control_test)} when additive controls are
specified.
{p_end}

{p 4 4 2}
{cmd:e(control_table)} copies the additive-control coefficient results from the
underlying {cmd:regress} or {cmd:glm} fit before APCIFLEX reparameterizes the
AP-cell surface.  Its columns are estimate, SE, test statistic, p-value, lower
confidence limit, and upper confidence limit.  {cmd:e(control_test)} stores the
joint test of all additive controls.  {cmd:e(n_controls)} gives the number of
controls.
{p_end}

{p 4 4 2}
{cmd:e(cell)} columns are age, period, cohort, raw N, Kish effective N, weight
sum, adjusted AP-cell linear predictor, its SE, interaction deviation,
interaction SE, e(b) coefficient index, Age group, Period group, and optional
cohort group. Under identity the linear predictor is the adjusted outcome mean;
under logit/probit it is on the link scale.
{p_end}

{p 4 4 2}
{cmd:e(cohort_raw)} columns are cohort, cells, joint statistic, joint df, raw
joint p, average deviation, average SE, average test statistic, raw average p,
slope/age, slope SE, slope test statistic, raw slope p, adjusted joint p,
adjusted average p, and adjusted slope p. The joint statistic is F under
identity and Wald chi-squared under logit/probit; pointwise statistics are t
under identity and z under logit/probit.
{p_end}

{p 4 4 2}
{cmd:e(age_group)}, {cmd:e(period_group)}, and optional
{cmd:e(cohort_group)} retain their original first nine columns:
group code, number of raw levels, cells, average estimate, SE, pointwise
statistic, raw average p, lower CI, and upper CI. Version 1.0.0 appends four
columns: joint statistic, joint df, raw joint p, and multiplicity-adjusted
joint p. The first nine columns are unchanged for backward compatibility.
{p_end}

{p 4 4 2}
Important locals include {cmd:e(aggweight)}, {cmd:e(mtest)},
{cmd:e(age_group_weighting)}, {cmd:e(period_group_weighting)},
{cmd:e(raw_cohort_weighting)}, {cmd:e(cohort_group_weighting)},
{cmd:e(grand_mean_weighting)}, {cmd:e(group_joint_test)},
{cmd:e(cohort_group_joint_test)}, and {cmd:e(group_joint_mtest_scope)}.
{p_end}

{title:Examples}

{p 8 8 2}
{cmd:apciflex y [pw=w], age(age) period(year) agegroup(age5) periodgroup(periodg) cohortgroup(cohort5) controls(x1 x2) mincell(30) aggweight(equal) mtest(holm) vce(cluster id)}
{p_end}

{p 8 8 2}{cmd:apciflex}{p_end}
{p 8 8 2}{cmd:apciflex, detail}{p_end}
{p 8 8 2}{cmd:predict double xb, xb}{p_end}
{p 8 8 2}{cmd:predict double se_xb, stdp}{p_end}
{p 8 8 2}{cmd:predict double resid, residuals}{p_end}

{p 8 8 2}
{cmd:apciflex ybin, age(age) period(year) agegroup(age5) periodgroup(periodg) controls(x1 x2) link(logit) vce(robust)}
{p_end}
{p 8 8 2}{cmd:predict double p, pr}{p_end}
{p 8 8 2}{cmd:margins, over(age period)}{p_end}
{p 8 8 2}{cmd:margins, dydx(x1)}{p_end}
{p 8 8 2}{cmd:margins, over(age period) predict(xb)}{p_end}

{title:Scope and limitations}

{p 4 4 2}
Version 1.0.0 supports identity, logit, and probit links. The nonlinear
extension is deliberately link-scale for APC-I decomposition and response-scale
for probability prediction/margins. APCIFLEX does not create additive
probability-scale Age, Period, or cohort components because the inverse link is
nonlinear.
{p_end}

{p 4 4 2}
Controls are currently numeric additive variables. Adjusted AP-cell linear predictors and
the APC-I grand link mean are evaluated at controls=0. Center controls before
estimation when a mean-profile interpretation is desired.
{p_end}


{title:Nonlinear stored results}

{p 4 4 2}
{cmd:e(link)} stores identity, logit, or probit; {cmd:e(family)} stores gaussian
or binomial; {cmd:e(inference)} stores t or z; and {cmd:e(component_scale)}
identifies the scale of the APC-I components.
{p_end}

{p 4 4 2}
{cmd:e(global_stat)} stores the global interaction statistic. For identity,
{cmd:e(global_F)} is populated and {cmd:e(global_chi2)} is missing. For
logit/probit, {cmd:e(global_chi2)} is populated, {cmd:e(global_F)} is missing,
and {cmd:e(global_df2)} is missing because large-sample chi-squared inference
is used.
{p_end}

{p 4 4 2}
For logit/probit, {cmd:e(ll)} and {cmd:e(deviance)} record the fitted binomial
GLM diagnostics when available. {cmd:e(r2)} and {cmd:e(rmse)} are retained as
missing for nonlinear links.
{p_end}


{title:Margins prediction interface on incomplete support}

{pstd}
When retained Age-by-Period support is incomplete, APCIFLEX uses
margins-specific prediction statistics.  The default response for
{cmd:margins} is {cmd:apcxb} under {cmd:link(identity)} and {cmd:apcpr}
under {cmd:link(logit)} or {cmd:link(probit)}.  These statistics are
numerically identical to {cmd:xb} and {cmd:pr} on retained AP cells, but they
are routed through {cmd:apciflex_p} so that unsupported counterfactual
Age-by-Period combinations are detected.

{pstd}
On incomplete support, {cmd:margins, predict(xb)} and
{cmd:margins, predict(pr)} are intentionally disallowed.  Use the default
margins response, {cmd:predict(apcxb)}, or {cmd:predict(apcpr)} as
appropriate.  Broad requests such as {cmd:margins age} or
{cmd:margins period} return {cmd:r(459)} when they require unretained AP
combinations, including with {cmd:nose}.  Retained-cell summaries such as
{cmd:margins, over(age period)} remain supported.


{title:Finite-cell inference note}

{pstd}
Under {cmd:link(identity)}, APCIFLEX currently uses the model residual degrees
of freedom {cmd:e(df_r)} for pointwise t-based confidence intervals of raw
Age, Period, interaction, and grouped linear-contrast estimates.  The command
does not assign a separate cell-specific degrees of freedom to each
{it:theta}_{ap}.  When retained AP cells are small, cell-level interaction
intervals should therefore be interpreted as asymptotic/large-cell
approximations.  {cmd:mincell()} and the reported cell support diagnostics
are especially important for fine-resolution interaction inference.


{title:Author}

{p 4 4 2}
Fengguang Lyu, Department of Sociology, Xi'an Jiaotong University.
Email: {browse "mailto:lvfg1999@126.com":lvfg1999@126.com}.
{p_end}
