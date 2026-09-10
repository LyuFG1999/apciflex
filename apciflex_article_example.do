********************************************************************************
* APCIFLEX: data generation and examples for the accompanying manuscript
* Run from the extracted submission folder:
*     do apciflex_article_example.do
* Generates apciflex_example.dta and apciflex_article_example.log.
********************************************************************************
version 16.0
clear all
set more off
set linesize 200
capture log close _all
log using "apciflex_article_example.log", text replace
set rng mt64
set seed 12345
adopath ++ "ado"
which apciflex

* 4.1 Generate the simulated data and inspect the observed grid.
clear
input byte age
25
30
35
40
45
50
55
60
end
tempfile ages
save "`ages'", replace

clear
input int period
2010
2012
2015
2018
2022
end

cross using "`ages'"
sort age period

* Unequal cell sizes make the example closer to survey data.
gen int ncell = 110 + 10*mod(age/5 + period, 5)
expand ncell
bysort age period: gen int within=_n
drop ncell

* Three unobserved age-by-period cells.
drop if age==25 & period==2012
drop if age==40 & period==2015
drop if age==55 & period==2018

* One observed but weak cell; mincell(30) will remove it.
drop if age==60 & period==2022 & within>10
drop within

* Theory-guided interpretation groups.
gen byte agegrp=cond(age<=30,1,cond(age<=40,2,cond(age<=50,3,4)))
label define agegrp_lbl 1 "25-30" 2 "35-40" 3 "45-50" 4 "55-60"
label values agegrp agegrp_lbl

gen byte periodgrp=cond(period<=2012,1,cond(period<=2018,2,3))
label define periodgrp_lbl 1 "2010-2012" 2 "2015-2018" 3 "2022"
label values periodgrp periodgrp_lbl

* Define raw cohorts for the reporting-group variables.
gen int cohort=period-age
gen byte cohortgrp = cond(cohort<=1969,1, ///
    cond(cohort<=1979,2,cond(cohort<=1989,3,4)))
label define cohortgrp_lbl 1 "<=1969" 2 "1970-79" 3 "1980-89" 4 ">=1990"
label values cohortgrp cohortgrp_lbl

* Centered additive controls.
gen byte female=runiform()<.5
gen double female_c=female-.5
gen double ses_c=rnormal()

* Continuous outcome with a nonadditive cohort-related component.
gen double wellbeing = 5 ///
    + .012*(age-42.5) ///
    - .12*(period==2010) ///
    - .05*(period==2012) ///
    + .04*(period==2018) ///
    + .10*(period==2022) ///
    + .20*sin((cohort-1975)/5) ///
    + .06*((age-42.5)/17.5)*((period-2015)/7) ///
    + .30*female_c ///
    + .20*ses_c ///
    + rnormal(0,.55)

label variable wellbeing "Well-being score"
label variable female_c  "Female, centered"
label variable ses_c     "Socioeconomic covariate"

save "apciflex_example.dta", replace

* 4.1 Load the generated data and inspect the grid.

use "apciflex_example.dta", clear
count
tabulate age period

* 4.2 Estimate the model and report grouped results.
apciflex wellbeing, age(age) period(period) ///
    agegroup(agegrp) periodgroup(periodgrp) ///
    cohortgroup(cohortgrp) controls(female_c ses_c) ///
    mincell(30) mtest(holm) vce(robust)
	
display "ARTICLE EXAMPLES COMPLETED; SEED = 12345"
log close
