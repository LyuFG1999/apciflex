APCIFLEX 1.0.0
==============

Author: Fengguang Lyu
Department of Sociology, Xi'an Jiaotong University, Xi'an, China
Email: lvfg1999@126.com

Purpose
-------
apciflex fits age-period-cohort interaction models on a connected observed
age-by-period grid and aggregates the components into separately defined age,
period, and cohort groups. Identity, logit, and probit links are supported.

Requirements and installation
-----------------------------
Stata 16 or later. Extract the archive and set the extracted directory as
Stata's working directory. To install the command, type:

    net install apciflex, from("ado") replace
    which apciflex
    help apciflex

Reproducing the article
----------------------
From the extracted directory, type:

    do apciflex_article_example.do

The do-file loads the bundled ado directory, uses the mt64 random-number
generator with seed 12345, creates apciflex_example.dta, and reproduces the
article's age-by-period counts, model fit, and grouped estimates and tests.
It writes apciflex_article_example.log. A generated log is supplied for
comparison. The simulated data are generated entirely by the do-file; no
external dataset or download is required. Running the do-file replaces its
previously generated dataset and log.

Program contents
----------------

    ado/apciflex.ado                  Estimation command
    ado/apciflex_p.ado                Prediction program
    ado/apciflex.sthlp                Detailed help
    ado/apciflex.pkg                  Installation manifest
    ado/stata.toc                     Installation catalog
    ado/README_APCIFLEX.txt           This guide
    ado/LICENSE.txt                  MIT license

Weights and interpretation
--------------------------
Estimation weights follow the underlying regress/glm rules. mincell() uses
Kish effective sample size for pweights, frequency totals for fweights, and
physical row counts otherwise. The additive decomposition gives equal weight
to retained cells. aggweight() controls aggregation after decomposition:
equal, sample, effective, or population. See help apciflex for definitions,
permitted weight combinations, statistical tests, and stored results.

All components are on the selected link scale. predict supports xb, stdp,
and residuals, and pr for binary links. Prediction requests requiring cells
outside the retained grid return r(459). See help apciflex for margins.

Citation and license
--------------------
Lyu, Fengguang. 2026. apciflex: Age-period-cohort interaction analysis with
flexible grouping on irregular observed Lexis grids. Version 1.0.0.
Stata software. Cite the accompanying article when published.
MIT License; see LICENSE.txt.
