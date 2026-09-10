# apciflex

**Age–period–cohort interaction analysis with flexible grouping on irregular observed Lexis grids**

This repository contains the data accompanying the paper introducing `apciflex`, a Stata command that combines estimation on an observed age–period grid with flexible grouping of age, period, and cohort.

📄 **[Read the paper on SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7440381)**

## Overview

In applied research, surveys may be conducted at unevenly spaced intervals, age coverage may differ across waves, and some age–period combinations may be unobserved. Researchers may also need different groupings for different questions—for example, ages grouped by life stage, periods by historical events, and cohorts by birth decade.

`apciflex` brings these choices into one estimation and grouping procedure. It fits the outcome surface at the original age and period values, decomposes the fitted surface, and then constructs summaries for separately defined age, period, and cohort groups.

## Method

### 1. Estimate and decompose the observed surface

The method builds on the age–period–cohort interaction (APC-I) model developed by [Luo and Hodges (2022)](https://doi.org/10.1177/0049124119882451). For each retained age–period cell, the fitted outcome surface is decomposed as

$$
\hat{\eta}_{ap}=\hat{\mu}+\hat{\alpha}_a+\hat{\pi}_p+\hat{\theta}_{ap}
$$

Here, $\hat{\mu}$ is the overall level, $\hat{\alpha}_a$ and $\hat{\pi}_p$ are the age and period components, and $\hat{\theta}_{ap}$ is the cell's deviation from their additive structure. The surface can incorporate control variables and is expressed on the scale of the selected link function.

An equal-cell least-squares projection determines the additive age–period structure over the retained cells. On a complete grid, this gives the usual effect-coded APC-I decomposition for the same fitted surface. On an incomplete grid, estimation and decomposition use the observed combinations. The retained cells must connect all included ages and periods, and the grid must retain interaction degrees of freedom.

Cohort membership follows actual time values:

$$
c=p-a
$$

Interaction components belonging to the same cohort are then summarized to describe that cohort's average deviation from the additive age–period baseline. Actual survey years preserve the temporal relationship between observations when survey intervals differ.

### 2. Define reporting groups after estimation

Age, period, and cohort groups can be specified separately. Each original time value maps to one reporting group in its dimension. For example, a study can combine five-year age groups, survey periods defined by policy changes, and ten-year birth cohorts.

Age and period group means aggregate their corresponding original components. Cohort group means first aggregate interaction components within each original cohort and then aggregate those cohort means within the reporting group. Aggregation can use equal weights, sample sizes, effective sample sizes, or population weights, subject to the selected estimation weights.

Holding the estimation sample and model specification fixed, changing the reporting groups preserves the underlying fitted surface and decomposition. This allows researchers to compare summaries that reflect different theoretical definitions of the time dimensions.

### 3. Obtain estimates and statistical tests

Grouped estimates are linear contrasts of the fitted cell surface. Their standard errors use the full covariance matrix, including covariances between components.

The analysis provides three types of tests:

| Test                      | Question                                                                                          |
| ------------------------- | ------------------------------------------------------------------------------------------------- |
| Global interaction test   | Do the retained cells contain interaction variation beyond the additive age and period structure? |
| Group mean test           | Does a group's average component differ from zero?                                                |
| Joint test within a group | Are the original components included in the group jointly zero?                                   |

For cohort groups, the joint test concerns the mean deviations of the original cohorts within the group. A group can have an average close to zero while still containing detectable deviations among its constituent cohorts.

Bonferroni, Holm, and Benjamini–Hochberg adjustments are available for joint tests, with separate test families for age, period, and cohort groups. Group mean tests and confidence intervals use pointwise inference.

## Relationship to existing APC-I software

The project builds on the APC-I framework and its existing software implementations. [Xu and Luo (2022)](https://journal.r-project.org/articles/RJ-2022-026/) introduce the R package `APCI` and Stata command; their R implementation also accommodates unequal age and period intervals. The focus of `apciflex` is the combination of decomposition on an irregular observed grid and separately defined reporting groups, using a shared fitted surface and covariance matrix.

## Example data

The paper uses simulated data generated in Stata with the `mt64` random-number generator and seed **12345**. The example includes:

- Eight observed ages: 25, 30, 35, 40, 45, 50, 55, and 60.
- Five survey years: 2010, 2012, 2015, 2018, and 2022.
- Three unobserved age–period cells and one cell below the minimum sample-size threshold.
- Separately defined age, period, and cohort reporting groups.

The generated dataset contains **4,650 observations**. Applying the example's minimum cell size of 30 retains **4,640 observations in 36 cells**. The data illustrate the estimation and grouping procedure; they are simulated and do not represent an empirical population.

The paper provides the model specification, command syntax, and interpretation of the results.

## Paper and references

**Accompanying paper**

Lyu, Fengguang. *apciflex: Age-period-cohort interaction analysis with flexible grouping on irregular observed Lexis grids.* [SSRN working paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7440381).

**Methodological background**

- Luo, L., and J. S. Hodges. 2022. The age-period-cohort-interaction model for describing and investigating inter-cohort deviations and intra-cohort life-course dynamics. *Sociological Methods & Research* 51: 1164–1210. [DOI](https://doi.org/10.1177/0049124119882451).
- Xu, J., and L. Luo. 2022. APCI: An R and Stata package for visualizing and analyzing age-period-cohort data. *The R Journal* 14: 77–95. [Article](https://journal.r-project.org/articles/RJ-2022-026/).

## Contact

Fengguang Lyu · Xi'an Jiaotong University
Email: [lvfg1999@126.com](mailto:lvfg1999@126.com)
