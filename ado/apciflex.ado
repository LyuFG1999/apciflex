*! apciflex 1.0.0 08sep2026
*! Flexible-resolution APC-I on irregular/incomplete Lexis support
*! identity/logit/probit links + theory-guided post-estimation aggregation

program define apciflex, eclass sortpreserve
    version 16.0

    * Standard estimation-command replay.
    if replay() {
        if "`e(cmd)'"!="apciflex" error 301
        syntax [, DETAIL]
        if "`detail'"!="" _apciflex_display, detail replay
        else _apciflex_display, replay
        exit
    }

    local cmdline `"apciflex `0'"'

    syntax varname(numeric) [if] [in] [fw aw pw iw], ///
        AGE(varname numeric) ///
        PERIOD(varname numeric) ///
        AGEGroup(varname numeric) ///
        PERIODGroup(varname numeric) ///
        [ COHORTGroup(varname numeric) ///
          CONTROLS(varlist numeric) ///
          MINCELL(integer 1) ///
          AGGWeight(string) ///
          MTEST(string) ///
          LINK(string) ///
          VCE(string asis) ///
          DETAIL ///
          NODISPLAY ///
          Level(cilevel) ]

    local y `varlist'
    marksample touse
    markout `touse' `y' `age' `period'
    if "`controls'"!="" markout `touse' `controls'

    if `mincell'<1 {
        di as err "mincell() must be an integer >= 1."
        exit 198
    }

    if "`aggweight'"=="" local aggweight "equal"
    local aggweight=lower("`aggweight'")

    if !inlist("`aggweight'","equal","sample","effective","population") {
        di as err "aggweight() must be equal, sample, effective, or population."
        exit 198
    }

    * effective aggregation has a clear Kish interpretation for survey
    * pweights; without weights it reduces to raw-N weighting.
    if "`aggweight'"=="effective" & !inlist("`weight'","","pweight") {
        di as err "aggweight(effective) is supported with pweights or unweighted data."
        di as err "For aweights/iweights/fweights, choose equal, sample, or an appropriate population specification."
        exit 198
    }

    * population aggregation uses the supplied survey/frequency mass.
    * With no weights it reduces to sample weighting.
    if "`aggweight'"=="population" & !inlist("`weight'","","pweight","fweight") {
        di as err "aggweight(population) is supported with pweights, fweights, or unweighted data."
        exit 198
    }

    if "`mtest'"=="" local mtest "none"
    local mtest=lower("`mtest'")
    if !inlist("`mtest'","none","bonferroni","holm","bh") {
        di as err "mtest() must be none, bonferroni, holm, or bh."
        exit 198
    }

    if "`link'"=="" local link "identity"
    local link=lower("`link'")
    if !inlist("`link'","identity","logit","probit") {
        di as err "link() must be identity, logit, or probit."
        exit 198
    }

    local family "gaussian"
    local inference "t"
    local component_scale "outcome"
    if inlist("`link'","logit","probit") {
        local family "binomial"
        local inference "z"
        if "`link'"=="logit" local component_scale "log_odds"
        else local component_scale "probit_index"
    }

    * Parse cluster VCE early enough that missing cluster identifiers are
    * excluded before AP support, mincell(), and e(sample) are constructed.
    local requested_clustvar ""
    if `"`vce'"'!="" {
        local __vcecopy `"`vce'"'
        gettoken __vcehead __vcetail : __vcecopy
        local __vcehead=lower("`__vcehead'")
        if "`__vcehead'"=="cluster" {
            gettoken requested_clustvar __vceextra : __vcetail
            if "`requested_clustvar'"=="" | `"`__vceextra'"'!="" {
                di as err "vce(cluster ...) requires exactly one cluster variable."
                exit 198
            }
            capture confirm variable `requested_clustvar'
            if _rc {
                di as err "cluster variable `requested_clustvar' not found."
                exit 111
            }
            markout `touse' `requested_clustvar'
        }
    }

    * Weight validity is part of the analytical-sample definition.
    tempvar __w
    if "`weight'"!="" {
        quietly gen double `__w' `exp' if `touse'
        quietly replace `touse'=0 if missing(`__w') | `__w'<=0
    }
    else {
        quietly gen double `__w'=1 if `touse'
    }

    quietly count if `touse'
    if r(N)==0 {
        di as err "no observations remain in the analytical sample."
        exit 2000
    }

    if inlist("`link'","logit","probit") {
        quietly count if `touse' & !inlist(`y',0,1)
        if r(N)>0 {
            di as err "link(`link') requires a binary outcome coded 0/1."
            di as err r(N) " analytical observation(s) have values outside {0,1}."
            exit 459
        }
        quietly summarize `y' if `touse', meanonly
        if r(min)==r(max) {
            di as err "link(`link') requires both outcome values 0 and 1 in the analytical sample."
            exit 459
        }
    }

    * Standard Stata factor-variable metadata used by margins requires
    * nonnegative integer-valued categorical levels.
    quietly count if `touse' & (`age'<0 | `age'!=floor(`age'))
    if r(N)>0 {
        di as err "age() must contain nonnegative integer-valued raw levels."
        di as err "This is required for factor-variable-compatible postestimation."
        exit 459
    }
    quietly count if `touse' & (`period'<0 | `period'!=floor(`period'))
    if r(N)>0 {
        di as err "period() must contain nonnegative integer-valued raw levels."
        di as err "This is required for factor-variable-compatible postestimation."
        exit 459
    }

    * Group variables must be defined on the analytical sample; do not silently
    * drop observations because a reporting group is missing.
    foreach g in `agegroup' `periodgroup' {
        quietly count if `touse' & missing(`g')
        if r(N)>0 {
            di as err "`g' is missing for " r(N) " analytical observation(s)."
            exit 459
        }
    }
    if "`cohortgroup'"!="" {
        quietly count if `touse' & missing(`cohortgroup')
        if r(N)>0 {
            di as err "`cohortgroup' is missing for " r(N) " analytical observation(s)."
            exit 459
        }
    }

    * Raw cohort is always defined internally.
    tempvar __cohort __cohortgroup
    quietly gen double `__cohort'=`period'-`age' if `touse'

    if "`cohortgroup'"!="" {
        quietly gen double `__cohortgroup'=`cohortgroup' if `touse'
        local hascg=1
    }
    else {
        quietly gen double `__cohortgroup'=. if `touse'
        local hascg=0
    }

    * Deterministic mapping from raw temporal values to theory groups.
    quietly _apciflex_checkmap if `touse', ///
        base(`age') group(`agegroup') name("agegroup()")
    quietly _apciflex_checkmap if `touse', ///
        base(`period') group(`periodgroup') name("periodgroup()")
    if `hascg' {
        quietly _apciflex_checkmap if `touse', ///
            base(`__cohort') group(`cohortgroup') name("cohortgroup()")
    }

    * Source support before mincell().
    quietly levelsof `age' if `touse', local(alevels0)
    quietly levelsof `period' if `touse', local(plevels0)
    local A0 : word count `alevels0'
    local P0 : word count `plevels0'
    local potential_cells=`A0'*`P0'

    tempvar __nraw __tag __w2 __wsum __w2sum __neff __support __small
    quietly bysort `age' `period': egen long `__nraw'=total(`touse')
    quietly egen byte `__tag'=tag(`age' `period') if `touse'

    quietly gen double `__w2'=`__w'^2 if `touse'
    quietly bysort `age' `period': egen double `__wsum'=total(`__w') if `touse'
    quietly bysort `age' `period': egen double `__w2sum'=total(`__w2') if `touse'
    quietly gen double `__neff'=(`__wsum'^2)/`__w2sum' if `touse'

    quietly count if `__tag'
    local observed_cells=r(N)
    local empty_cells=`potential_cells'-`observed_cells'

    * Support metric for mincell().  pweights are arbitrary up to multiplicative
    * scale, so use scale-invariant Kish effective N.  fweights represent
    * repeated observations, so use within-cell frequency mass (sum fweight),
    * not the number of physical data rows.  Other weights use physical row N.
    if "`weight'"=="pweight" {
        quietly gen double `__support'=`__neff' if `touse'
        local support_metric "Kish effective N"
        local support_code "kish_neff"
    }
    else if "`weight'"=="fweight" {
        quietly gen double `__support'=`__wsum' if `touse'
        local support_metric "frequency-weighted N"
        local support_code "fw_mass"
    }
    else {
        quietly gen double `__support'=`__nraw' if `touse'
        local support_metric "physical row N"
        local support_code "raw_n"
    }

    quietly gen byte `__small'=(`touse' & `__support'<`mincell')
    quietly count if `__tag' & `__support'<`mincell'
    local dropped_cells=r(N)
    quietly count if `__small'
    local dropped_rows=r(N)

    if `dropped_rows'>0 quietly replace `touse'=0 if `__small'

    quietly count if `touse'
    if r(N)==0 {
        di as err "mincell(`mincell') removes all observed Age x Period cells."
        exit 2000
    }

    * Retained raw support.
    quietly levelsof `age' if `touse', local(alevels)
    quietly levelsof `period' if `touse', local(plevels)
    local A : word count `alevels'
    local P : word count `plevels'

    if `A'<2 | `P'<2 {
        di as err "at least two retained Age levels and two retained Period levels are required."
        exit 459
    }

    tempvar __apcell
    quietly egen long `__apcell'=group(`age' `period') if `touse'
    quietly summarize `__apcell' if `touse', meanonly
    local K=r(max)

    if inlist("`link'","logit","probit") {
        tempvar __ymin __ymax __ytag
        quietly bysort `__apcell': egen double `__ymin'=min(`y') if `touse'
        quietly bysort `__apcell': egen double `__ymax'=max(`y') if `touse'
        quietly egen byte `__ytag'=tag(`__apcell') if `touse'
        quietly count if `__ytag' & `__ymin'==`__ymax'
        if r(N)>0 {
            di as err "link(`link') cannot estimate the saturated AP-cell surface:"
            di as err r(N) " retained Age x Period cell(s) contain only one binary outcome value."
            di as err "Increase cell support, combine substantively defensible raw levels, or use link(identity) for a linear probability specification."
            exit 459
        }
    }

    if `K'<=`A'+`P'-1 {
        di as err "the retained support has no residual Age x Period interaction degrees of freedom."
        exit 459
    }

    * A common observed-support additive baseline requires connected support.
    tempname CONNECTED
    mata: st_numscalar("`CONNECTED'", ///
        _apciflex_connected("`age'","`period'","`touse'"))
    if scalar(`CONNECTED')!=1 {
        di as err "the retained Age x Period support is disconnected."
        di as err "A common additive Age+Period baseline is not identified across components."
        di as err "Reduce mincell(), add support, or analyze connected components separately."
        exit 459
    }

    * Saturated observed AP-cell model with optional additive controls.
    * identity: Gaussian linear model fitted by regress.
    * logit/probit: binomial GLM; APC-I decomposition is on the link scale.
    local wgt ""
    if "`weight'"!="" local wgt "[`weight'`exp']"

    local comma ""
    if `"`vce'"'!="" local comma `", vce(`vce')"'

    if "`link'"=="identity" {
        capture quietly regress `y' i.`__apcell' `controls' if `touse' `wgt' `comma'
        local fit_rc=_rc
        if `fit_rc' {
            di as err "Underlying regress failed; original Stata diagnostics follow:"
            capture noisily regress `y' i.`__apcell' `controls' if `touse' `wgt' `comma'
            local fit_rc2=_rc
            if `fit_rc2' exit `fit_rc2'
            exit `fit_rc'
        }
    }
    else {
        local glmopts "family(binomial) link(`link')"
        if `"`vce'"'!="" local glmopts "`glmopts' vce(`vce')"
        capture quietly glm `y' i.`__apcell' `controls' if `touse' `wgt', `glmopts'
        local fit_rc=_rc
        if `fit_rc' {
            di as err "Underlying binomial GLM failed; original Stata diagnostics follow:"
            capture noisily glm `y' i.`__apcell' `controls' if `touse' `wgt', `glmopts'
            local fit_rc2=_rc
            if `fit_rc2' exit `fit_rc2'
            exit `fit_rc'
        }
    }

    tempvar __fitsample
    quietly gen byte `__fitsample'=e(sample)

    quietly count if `touse'
    local NROWS=r(N)
    quietly count if `__fitsample'
    if r(N)!=`NROWS' {
        di as err "the fitted model changed the retained analytical sample."
        exit 498
    }
    quietly count if `touse' & !`__fitsample'
    if r(N)>0 {
        di as err "the fitted model changed the retained analytical sample."
        exit 498
    }

    tempname BASEB BASEV DFRES NOBS RSQ RMSE LL DEVIANCE NCLUST
    matrix `BASEB'=e(b)
    matrix `BASEV'=e(V)
    scalar `NOBS'=e(N)

    scalar `DFRES'=.
    scalar `RSQ'=.
    scalar `RMSE'=.
    scalar `LL'=.
    scalar `DEVIANCE'=.

    if "`link'"=="identity" {
        scalar `DFRES'=e(df_r)
        scalar `RSQ'=e(r2)
        scalar `RMSE'=e(rmse)
    }
    else {
        capture scalar `LL'=e(ll)
        capture scalar `DEVIANCE'=e(deviance)
    }

    local base_vce "`e(vce)'"
    local base_vcetype "`e(vcetype)'"
    local base_clustvar "`e(clustvar)'"
    local has_nclust=0
    capture scalar `NCLUST'=e(N_clust)
    if !_rc local has_nclust=1

    * Additive-control coefficients and tests are copied directly from the
    * underlying regress/glm fit, before APCIFLEX reparameterizes e(b)/e(V).
    local nctrl : word count `controls'
    tempname CTRLTAB CTRLTEST CTRLCRIT CB CS CT CP
    if `nctrl'>0 {
        matrix `CTRLTAB'=J(`nctrl',6,.)
        matrix colnames `CTRLTAB'=estimate se stat pvalue lb ub
        matrix rownames `CTRLTAB'=`controls'

        if "`inference'"=="t" {
            scalar `CTRLCRIT'=invttail(scalar(`DFRES'),(100-`level')/200)
        }
        else {
            scalar `CTRLCRIT'=invnormal(1-(100-`level')/200)
        }

        local rr=0
        foreach z of local controls {
            local ++rr
            local jj=colnumb(`BASEB',"`z'")
            if `jj'>=. {
                di as err "control coefficient `z' could not be mapped to the underlying model."
                exit 498
            }

            scalar `CB'=`BASEB'[1,`jj']
            scalar `CS'=sqrt(`BASEV'[`jj',`jj'])
            scalar `CT'=.
            scalar `CP'=.
            if scalar(`CS')>0 & scalar(`CS')<. {
                scalar `CT'=scalar(`CB')/scalar(`CS')
                if "`inference'"=="t" scalar `CP'=2*ttail(scalar(`DFRES'),abs(scalar(`CT')))
                else scalar `CP'=2*normal(-abs(scalar(`CT')))
            }

            matrix `CTRLTAB'[`rr',1]=scalar(`CB')
            matrix `CTRLTAB'[`rr',2]=scalar(`CS')
            matrix `CTRLTAB'[`rr',3]=scalar(`CT')
            matrix `CTRLTAB'[`rr',4]=scalar(`CP')
            matrix `CTRLTAB'[`rr',5]=scalar(`CB')-scalar(`CTRLCRIT')*scalar(`CS')
            matrix `CTRLTAB'[`rr',6]=scalar(`CB')+scalar(`CTRLCRIT')*scalar(`CS')
        }

        matrix `CTRLTEST'=J(1,4,.)
        matrix colnames `CTRLTEST'=stat df1 df2 pvalue
        capture quietly testparm `controls'
        if !_rc {
            if "`inference'"=="t" {
                matrix `CTRLTEST'[1,1]=r(F)
                matrix `CTRLTEST'[1,2]=r(df)
                matrix `CTRLTEST'[1,3]=r(df_r)
                matrix `CTRLTEST'[1,4]=r(p)
            }
            else {
                matrix `CTRLTEST'[1,1]=r(chi2)
                matrix `CTRLTEST'[1,2]=r(df)
                matrix `CTRLTEST'[1,3]=.
                matrix `CTRLTEST'[1,4]=r(p)
            }
        }
    }

    * Recover retained adjusted AP-cell linear predictors at controls=0.
    tempname LMU MU VMU
    matrix `LMU'=J(`K',colsof(`BASEB'),0)

    local consj=colnumb(`BASEB',"_cons")
    if `consj'>=. {
        di as err "internal parameter mapping failed: fitted-model intercept not found."
        exit 498
    }

    forvalues k=1/`K' {
        matrix `LMU'[`k',`consj']=1
        if `k'>1 {
            local jj=colnumb(`BASEB',"`k'.`__apcell'")
            if `jj'>=. {
                di as err "internal AP-cell parameter mapping failed for retained cell `k'."
                di as err "This is not an ordinary control-variable collinearity message."
                exit 498
            }
            matrix `LMU'[`k',`jj']=1
        }
    }

    mata: _apciflex_surface( ///
        "`BASEB'","`BASEV'","`LMU'","`MU'","`VMU'")

    * Cell metadata in exactly the same order as MU.
    tempname CELLINFO
    mata: _apciflex_cellinfo( ///
        "`__apcell'","`touse'","`age'","`period'","`__cohort'", ///
        "`agegroup'","`periodgroup'","`__cohortgroup'", ///
        "`__nraw'","`__neff'","`__wsum'","`CELLINFO'")

    * Observed-support APC-I decomposition.
    tempname AGERAW PERIODRAW VAGE VPER THETA VTH CELL GLOBAL COHORTRAW
    tempname GRAND GRANDSE GRANDP

    mata: _apciflex_decompose( ///
        "`CELLINFO'","`MU'","`VMU'","`DFRES'","`inference'",`level', ///
        "`AGERAW'","`VAGE'","`PERIODRAW'","`VPER'", ///
        "`THETA'","`VTH'","`CELL'", ///
        "`GRAND'","`GRANDSE'","`GRANDP'")

    local interaction_df=`K'-`A'-`P'+1
    mata: _apciflex_global_test( ///
        "`THETA'","`VTH'","`DFRES'","`inference'","`GLOBAL'")

    mata: _apciflex_cohort_table( ///
        "`CELLINFO'","`THETA'","`VTH'","`DFRES'","`inference'", ///
        "`aggweight'","`mtest'","`COHORTRAW'")

    * Theory-guided post-estimation aggregation.
    tempname AGEGROUP PERIODGROUP COHORTGROUP
    mata: _apciflex_group_ageperiod( ///
        1,"`CELLINFO'","`AGERAW'","`VAGE'","`DFRES'","`inference'",`level',"`aggweight'","`mtest'","`AGEGROUP'")
    mata: _apciflex_group_ageperiod( ///
        2,"`CELLINFO'","`PERIODRAW'","`VPER'","`DFRES'","`inference'",`level',"`aggweight'","`mtest'","`PERIODGROUP'")
    if `hascg' {
        mata: _apciflex_group_cohort( ///
            "`CELLINFO'","`THETA'","`VTH'","`DFRES'","`inference'",`level',"`aggweight'","`mtest'","`COHORTGROUP'")
    }

    matrix colnames `AGERAW'=value estimate se stat pvalue lb ub
    matrix colnames `PERIODRAW'=value estimate se stat pvalue lb ub
    matrix colnames `CELL'=age period cohort n_raw n_eff weight_sum mean mean_se interaction interaction_se coef_index age_group period_group cohort_group
    matrix colnames `COHORTRAW'=cohort cells joint_stat joint_df1 joint_p average average_se average_stat average_p slope_age slope_se slope_stat slope_p joint_p_adj average_p_adj slope_p_adj
    matrix colnames `GLOBAL'=stat df1 df2 pvalue
    matrix colnames `AGEGROUP'=group_code raw_levels cells estimate se stat pvalue lb ub joint_stat joint_df1 joint_p joint_p_adj
    matrix colnames `PERIODGROUP'=group_code raw_levels cells estimate se stat pvalue lb ub joint_stat joint_df1 joint_p joint_p_adj
    if `hascg' matrix colnames `COHORTGROUP'=group_code raw_levels cells estimate se stat pvalue lb ub joint_stat joint_df1 joint_p joint_p_adj

    * --------------------------------------------------------------------------
    * Proper e-class posting.
    *
    * e(b) is an invertible reparameterization of the saturated fitted model:
    *   [raw Age#Period cell means] = adjusted retained AP-cell linear predictors at controls=0
    *   [control:z]                 = additive control coefficients
    * The joint e(V) contains their full cross-covariance.
    * --------------------------------------------------------------------------
    tempname LPOST CSEL BPOST VPOST
    matrix `LPOST'=`LMU'
    if `nctrl'>0 {
        matrix `CSEL'=J(`nctrl',colsof(`BASEB'),0)
        local rr=0
        foreach z of local controls {
            local ++rr
            local jj=colnumb(`BASEB',"`z'")
            if `jj'>=. {
                di as err "control coefficient `z' could not be mapped to the fitted model."
                exit 498
            }
            matrix `CSEL'[`rr',`jj']=1
        }
        matrix `LPOST'=`LPOST' \ `CSEL'
    }

    matrix `BPOST'=(`LPOST'*`BASEB'')'
    matrix `VPOST'=`LPOST'*`BASEV'*`LPOST''

    * Factor-variable-compatible coefficient stripes.
    *
    * The K cell-mean parameters are represented as a no-base cell-means
    * interaction of the ORIGINAL raw variables:
    *
    *     ibn.age # ibn.period
    *
    * rather than synthetic c1...cK names.  This allows margins to map the
    * parameter vector back to the data variables age() and period().
    local post_names ""
    local age_base : word 1 of `alevels'
    local period_base : word 1 of `plevels'

    forvalues k=1/`K' {
        local aa=`CELLINFO'[`k',1]
        local pp=`CELLINFO'[`k',2]

        local aterm "`aa'.`age'"
        local pterm "`pp'.`period'"

        if `aa'==`age_base' local aterm "`aa'bn.`age'"
        if `pp'==`period_base' local pterm "`pp'bn.`period'"

        local post_names "`post_names' `aterm'#`pterm'"
    }

    foreach z of local controls {
        local post_names "`post_names' `z'"
    }

    matrix colnames `BPOST'=`post_names'
    matrix rownames `VPOST'=`post_names'
    matrix colnames `VPOST'=`post_names'

    local Npost=scalar(`NOBS')
    local DFpost=scalar(`DFRES')

    * buildfvinfo is required so margins/contrast can identify factor-variable
    * structure and estimability from e(b).
    tempname POSTRANK
    mata: st_numscalar("`POSTRANK'", rank(st_matrix("`VPOST'")))

    ereturn clear

    if "`inference'"=="t" {
        ereturn post `BPOST' `VPOST', ///
            esample(`touse') depname(`y') obs(`Npost') dof(`DFpost') buildfvinfo
    }
    else {
        ereturn post `BPOST' `VPOST', ///
            esample(`touse') depname(`y') obs(`Npost') buildfvinfo
    }

    ereturn scalar rank=scalar(`POSTRANK')
    ereturn scalar age_base=`age_base'
    ereturn scalar period_base=`period_base'

    * Result matrices.
    ereturn matrix age_raw=`AGERAW', copy
    ereturn matrix period_raw=`PERIODRAW', copy
    ereturn matrix cell=`CELL', copy
    ereturn matrix cohort_raw=`COHORTRAW', copy
    ereturn matrix age_group=`AGEGROUP', copy
    ereturn matrix period_group=`PERIODGROUP', copy
    if `hascg' ereturn matrix cohort_group=`COHORTGROUP', copy

    ereturn matrix V_cellmean=`VMU', copy
    ereturn matrix V_age_raw=`VAGE', copy
    ereturn matrix V_period_raw=`VPER', copy
    ereturn matrix V_interaction=`VTH', copy
    ereturn matrix global_test=`GLOBAL', copy
    if `nctrl'>0 {
        ereturn matrix control_table=`CTRLTAB', copy
        ereturn matrix control_test=`CTRLTEST', copy
    }

    * Scalars and metadata.
    ereturn scalar grand_mean=scalar(`GRAND')
    ereturn scalar grand_mean_se=scalar(`GRANDSE')
    ereturn scalar grand_mean_p=scalar(`GRANDP')

    ereturn scalar N_rows=`NROWS'
    ereturn scalar A_source=`A0'
    ereturn scalar P_source=`P0'
    ereturn scalar potential_cells=`potential_cells'
    ereturn scalar observed_cells=`observed_cells'
    ereturn scalar empty_cells=`empty_cells'
    ereturn scalar A=`A'
    ereturn scalar P=`P'
    ereturn scalar K=`K'
    ereturn scalar interaction_df=`interaction_df'
    ereturn scalar mincell=`mincell'
    ereturn scalar dropped_cells=`dropped_cells'
    ereturn scalar dropped_rows=`dropped_rows'
    * Backward-compatible alias: physical data rows removed by mincell().
    ereturn scalar dropped_observations=`dropped_rows'
    ereturn scalar has_cohortgroup=`hascg'
    ereturn scalar n_controls=`nctrl'
    ereturn scalar r2=scalar(`RSQ')
    ereturn scalar rmse=scalar(`RMSE')
    ereturn scalar ll=scalar(`LL')
    ereturn scalar deviance=scalar(`DEVIANCE')
    ereturn scalar level=`level'

    ereturn scalar global_stat=`GLOBAL'[1,1]
    ereturn scalar global_df1=`GLOBAL'[1,2]
    ereturn scalar global_df2=`GLOBAL'[1,3]
    ereturn scalar global_p=`GLOBAL'[1,4]
    if "`inference'"=="t" {
        ereturn scalar global_F=`GLOBAL'[1,1]
        ereturn scalar global_chi2=.
    }
    else {
        ereturn scalar global_F=.
        ereturn scalar global_chi2=`GLOBAL'[1,1]
    }

    if "`link'"=="identity" {
        ereturn local title "Flexible-resolution linear APC-I"
    }
    else {
        ereturn local title "Flexible-resolution binomial-`link' APC-I"
    }
    ereturn local cmd "apciflex"
    ereturn local cmdline `"`cmdline'"'
    ereturn local version "1.0.0"
    ereturn local predict "apciflex_p"
    ereturn local predict_support_policy "error_outside_retained_support"
    * margins policy:
    * Complete retained AP support keeps the standard fast xb/pr interface.
    * Incomplete retained AP support uses margins-specific prediction aliases
    * (apcxb/apcpr).  These aliases must pass through apciflex_p, where the
    * retained-support guard rejects unsupported counterfactual AP combinations.
    if `K' < (`A'*`P') {
        if "`link'"=="identity" {
            ereturn local marginsok "APCXB default"
            ereturn local marginsdefault "predict(apcxb)"
        }
        else {
            ereturn local marginsok "APCXB APCPR default"
            ereturn local marginsdefault "predict(apcpr)"
        }
        ereturn local marginsnotok "XB PR"
        ereturn local margins_derivative_policy "custom_predict_chainrule_candidate"
    }
    else {
        if "`link'"=="identity" ereturn local marginsok "XB APCXB default"
        else ereturn local marginsok "XB PR APCXB APCPR default"
        ereturn local margins_derivative_policy "chainrule_complete_support"
    }

    ereturn local fvops "true"

    ereturn local family "`family'"
    ereturn local link "`link'"
    ereturn local inference "`inference'"
    ereturn local component_scale "`component_scale'"
    ereturn local cell_scale "`component_scale'"
    ereturn local grand_scale "`component_scale'"

    ereturn local outcome_var "`y'"
    ereturn local age_var "`age'"
    ereturn local period_var "`period'"
    ereturn local agegroup_var "`agegroup'"
    ereturn local periodgroup_var "`periodgroup'"
    ereturn local cohortgroup_var "`cohortgroup'"
    ereturn local cohort_definition "period-age"
    ereturn local controls "`controls'"
    ereturn local mincell_metric "`support_metric'"
    ereturn local mincell_metric_code "`support_code'"

    * Interpretation-stage aggregation convention.
    ereturn local aggweight "`aggweight'"
    ereturn local mtest "`mtest'"
    ereturn local age_group_weighting "`aggweight'"
    ereturn local period_group_weighting "`aggweight'"
    ereturn local raw_cohort_weighting "`aggweight'"
    ereturn local cohort_group_weighting "`aggweight'"
    ereturn local grand_mean_weighting "equal_retained_ap_cells"
    ereturn local group_joint_test "all_raw_components_within_group_zero"
    ereturn local cohort_group_joint_test "all_raw_cohort_average_deviations_within_group_zero"
    ereturn local group_joint_mtest_scope "separate_age_period_cohort_families"

    if "`weight'"!="" {
        ereturn local wtype "`weight'"
        ereturn local wexp "`exp'"
    }
    if "`base_vce'"!="" ereturn local vce "`base_vce'"
    if "`base_vcetype'"!="" ereturn local vcetype "`base_vcetype'"
    if "`base_clustvar'"!="" ereturn local clustvar "`base_clustvar'"
    if `has_nclust' ereturn scalar N_clust=scalar(`NCLUST')

    * Compact output. detail controls ALL raw tables.
    if "`nodisplay'"=="" {
        if "`detail'"!="" _apciflex_display, detail
        else _apciflex_display
    }

end



program define _apciflex_display
    version 16.0
    syntax [, DETAIL REPLAY]

    if "`e(cmd)'"!="apciflex" error 301

    tempname AGEGROUP PERIODGROUP COHORTGROUP AGERAW PERIODRAW COHORTRAW
    tempname CTRLTAB CTRLTEST
    matrix `AGEGROUP'=e(age_group)
    matrix `PERIODGROUP'=e(period_group)
    matrix `AGERAW'=e(age_raw)
    matrix `PERIODRAW'=e(period_raw)
    matrix `COHORTRAW'=e(cohort_raw)
    if e(has_cohortgroup) matrix `COHORTGROUP'=e(cohort_group)
    if e(n_controls)>0 {
        matrix `CTRLTAB'=e(control_table)
        matrix `CTRLTEST'=e(control_test)
    }

    local y "`e(outcome_var)'"
    local age "`e(age_var)'"
    local period "`e(period_var)'"
    local agegroup "`e(agegroup_var)'"
    local periodgroup "`e(periodgroup_var)'"
    local cohortgroup "`e(cohortgroup_var)'"
    local aggweight "`e(aggweight)'"
    local mtest "`e(mtest)'"
    local link "`e(link)'"
    if "`link'"=="" local link "identity"
    local family "`e(family)'"
    if "`family'"=="" local family "gaussian"
    local inference "`e(inference)'"
    if "`inference'"=="" local inference "t"
    local level=e(level)
    local hascg=e(has_cohortgroup)

    local nolabel ""
    if "`replay'"!="" local nolabel "nolabel"

    di as txt _newline "=============================================================================="
    di as txt "APCIFLEX 1.0.0: FLEXIBLE-RESOLUTION APC-I"
    di as txt "=============================================================================="
    di as txt "Outcome                         : " as res "`y'"
    di as txt "Family / link                   : " as res "`family' / `link'"
    if "`link'"=="identity" di as txt "APC-I component scale           : " as res "outcome (identity-link)"
    else if "`link'"=="logit" di as txt "APC-I component scale           : " as res "log-odds (linear predictor)"
    else di as txt "APC-I component scale           : " as res "probit index (linear predictor)"
    di as txt "Raw Age / Period variables      : " as res "`age' / `period'"
    di as txt "Theory Age / Period groups      : " as res "`agegroup' / `periodgroup'"
    if `hascg' di as txt "Theory Cohort group             : " as res "`cohortgroup'"
    else di as txt "Theory Cohort group             : " as res "(not specified)"
    di as txt "Cohort definition               : " as res "period - age"
    di as txt "Analytical rows retained        : " as res %12.0fc e(N_rows)
    di as txt "Source Age / Period levels      : " as res e(A_source) " / " e(P_source)
    di as txt "Potential / observed AP cells   : " as res e(potential_cells) " / " e(observed_cells)
    di as txt "Structural empty AP cells       : " as res e(empty_cells)
    di as txt "mincell() support metric        : " as res "`e(mincell_metric)'"
    di as txt "mincell() threshold             : " as res e(mincell)
    di as txt "Cells / rows removed            : " as res e(dropped_cells) " / " e(dropped_rows)
    di as txt "Retained Age / Period / AP cells: " as res e(A) " / " e(P) " / " e(K)
    di as txt "Interaction df (support)        : " as res e(interaction_df)
    if "`e(vce)'"!="" {
        local vceshow "`e(vcetype)'"
        if "`vceshow'"=="" local vceshow "`e(vce)'"
        di as txt "VCE                             : " as res "`vceshow'"
        if "`e(clustvar)'"!="" {
            di as txt "Cluster variable / clusters     : " as res "`e(clustvar)'" ///
                as txt " / " as res %9.0f e(N_clust)
        }
    }
    tempname __HEAD
    matrix `__HEAD'=J(1,4,.)
    matrix `__HEAD'[1,1]=e(grand_mean)
    matrix `__HEAD'[1,2]=e(grand_mean_se)
    matrix `__HEAD'[1,3]=e(grand_mean_p)
    matrix `__HEAD'[1,4]=e(global_p)

    _apciflex_fmt4, matrix(`__HEAD') row(1) col(1)
    local __gm "`r(s)'"
    _apciflex_fmt4, matrix(`__HEAD') row(1) col(2)
    local __gse "`r(s)'"
    _apciflex_fmt4, matrix(`__HEAD') row(1) col(3) pvalue
    local __gp "`r(s)'"
    _apciflex_fmt4, matrix(`__HEAD') row(1) col(4) pvalue
    local __globalp "`r(s)'"

    if "`link'"=="identity" {
        di as txt "R-squared / Root MSE            : " as res %8.4f e(r2) " / " %8.4f e(rmse)
        di as txt "Grand mean (equal-cell; ctrls=0): " ///
            as res "`__gm'" as txt "   SE=" as res "`__gse'" ///
            as txt "   p=" as res "`__gp'"
    }
    else {
        di as txt "Log likelihood / Deviance       : " as res %11.4f e(ll) " / " %11.4f e(deviance)
        local glabel "log-odds"
        if "`link'"=="probit" local glabel "probit index"
        di as txt "Grand link mean (equal-cell; controls=0; `glabel'): " ///
            as res "`__gm'" as txt "   SE=" as res "`__gse'" ///
            as txt "   p=" as res "`__gp'"
    }

    di as txt _newline "GLOBAL NONADDITIVE AGE x PERIOD TEST"
    if "`inference'"=="t" {
        di as txt "F(" as res %6.0f e(global_df1) as txt ", " ///
            as res %10.0f e(global_df2) as txt ") = " ///
            as res %10.4f e(global_F) ///
            as txt "   p = " as res "`__globalp'"
    }
    else {
        di as txt "Wald chi2(" as res %6.0f e(global_df1) as txt ") = " ///
            as res %10.4f e(global_chi2) ///
            as txt "   p = " as res "`__globalp'"
    }

    if e(n_controls)>0 {
        di as txt _newline "ADDITIVE CONTROLS (UNDERLYING MODEL)"
        _apciflex_print_controls, matrix(`CTRLTAB') testmatrix(`CTRLTEST') ///
            level(`level') inference("`inference'")
    }

    di as txt _newline "THEORY-GUIDED AGE GROUPS"
    _apciflex_print_group, matrix(`AGEGROUP') gvar("`agegroup'") ///
        second("Raw") level(`level') inference("`inference'") mtest("`mtest'") `nolabel'

    di as txt _newline "THEORY-GUIDED PERIOD GROUPS"
    _apciflex_print_group, matrix(`PERIODGROUP') gvar("`periodgroup'") ///
        second("Raw") level(`level') inference("`inference'") mtest("`mtest'") `nolabel'

    if `hascg' {
        di as txt _newline "THEORY-GUIDED COHORT GROUPS"
        _apciflex_print_group, matrix(`COHORTGROUP') gvar("`cohortgroup'") ///
            second("Raw") level(`level') inference("`inference'") mtest("`mtest'") `nolabel'
    }

    if "`mtest'"!="none" {
        di as txt "Note: Pjoint* is adjusted by mtest(`mtest') within each displayed dimension."
        di as txt "      Pavg is unadjusted."
    }

    if "`detail'"!="" {
        di as txt _newline "RAW AGE COMPONENTS"
        _apciflex_print_effect, matrix(`AGERAW') gvar("`age'") level(`level') inference("`inference'") `nolabel'

        di as txt _newline "RAW PERIOD COMPONENTS"
        _apciflex_print_effect, matrix(`PERIODRAW') gvar("`period'") level(`level') inference("`inference'") `nolabel'

        di as txt _newline "RAW COHORT DIAGONAL ANALYSIS"
        _apciflex_print_cohort, matrix(`COHORTRAW') mtest("`mtest'") inference("`inference'")
        di as txt "Note: one-cell cohorts have no slope; their joint test is a one-cell interaction test."
        if e(K)<e(A)*e(P) {
            di as txt "Note: retained Age x Period support is incomplete; broad margins requests may return r(459)"
            di as txt "      when they require unretained Age x Period combinations."
        }
        if inlist("`link'","logit","probit") {
            di as txt "Note: the grand mean and APC-I components are on the model link scale; use predict/margins"
            di as txt "      for response-probability estimands."
        }
        if "`mtest'"=="none" {
            di as txt "      Raw cohort p-values are unadjusted; specify mtest() to adjust three hypothesis families separately."
        }
        else {
            di as txt "      Displayed cohort p-values use mtest(`mtest'); raw and adjusted p-values are both stored in e(cohort_raw)."
        }
    }
    else {
        di as txt _newline "Raw Age, Period, and cohort diagnostics are stored in e()."
        di as txt "Specify detail to display all raw-resolution tables."
    }

    di as txt _newline "Postestimation: predict xb/stdp; see help apciflex for margins on incomplete support."
    di as txt "Aggregation    : link(`link'), aggweight(`aggweight'), mtest(`mtest')."
    if "`detail'"!="" {
        di as txt "  APC-I components: e(age_raw), e(period_raw), e(cell), e(cohort_raw)."
        di as txt "  Unsupported Age x Period prediction requests return r(459)."
    }

    if "`replay'"!="" {
        di as txt "  Replay uses stored numeric group codes so it remains valid after the original dataset is cleared."
    }
end


program define _apciflex_checkmap
    version 16.0
    syntax [if], BASE(varname numeric) GROUP(varname numeric) NAME(string)

    marksample touse
    tempvar __min __max

    quietly bysort `base': egen double `__min'=min(`group') if `touse'
    quietly bysort `base': egen double `__max'=max(`group') if `touse'

    quietly count if `touse' & `__min'!=`__max'
    if r(N)>0 {
        di as err "`name' is not a valid deterministic grouping variable."
        di as err "At least one raw value of `base' is assigned to more than one group."
        exit 459
    }
end


program define _apciflex_print_effect
    version 16.0
    syntax, MATRIX(name) GVAR(string) Level(cilevel) [ INFERENCE(string) NOLABEL ]

    if "`inference'"=="" local inference "t"
    local statlabel "t"
    if "`inference'"!="t" local statlabel "z"

    local vl ""
    if "`nolabel'"=="" & "`gvar'"!="" {
        capture confirm variable `gvar'
        if !_rc local vl : value label `gvar'
    }

    di as txt "--------------------------------------------------------------------------"
    di as txt %-16s "Level" %10s "Estimate" %9s "SE" %8s "`statlabel'" %9s "P" ///
        %10s "LB" %10s "UB"

    forvalues i=1/`=rowsof(`matrix')' {
        local code=`matrix'[`i',1]
        local lab ""
        if "`vl'"!="" local lab : label (`gvar') `code'
        if `"`lab'"'=="" local lab = string(`code',"%16.0g")

        _apciflex_fmt4, matrix(`matrix') row(`i') col(2)
        local b "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(3)
        local se "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(4)
        local st "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(5) pvalue
        local p "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(6)
        local lb "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(7)
        local ub "`r(s)'"

        di as res %-16s `"`lab'"' %10s "`b'" %9s "`se'" %8s "`st'" ///
            %9s "`p'" %10s "`lb'" %10s "`ub'"
    }
end


program define _apciflex_fmt4, rclass
    version 16.0
    syntax, MATRIX(name) ROW(integer) COL(integer) [ PVALUE ]

    tempname __x
    scalar `__x'=`matrix'[`row',`col']

    if missing(scalar(`__x')) {
        return local s "."
        exit
    }

    if "`pvalue'"!="" & scalar(`__x')>=0 & scalar(`__x')<.0001 {
        return local s "<.0001"
        exit
    }

    local s = strtrim(string(scalar(`__x'),"%14.4f"))
    if inlist("`s'","-0.0000","+0.0000") local s "0.0000"

    return local s "`s'"
end


program define _apciflex_print_controls
    version 16.0
    syntax, MATRIX(name) TESTMATRIX(name) Level(cilevel) [ INFERENCE(string) ]

    if "`inference'"=="" local inference "t"
    local statlabel "t"
    local jointlabel "F"
    if "`inference'"!="t" {
        local statlabel "z"
        local jointlabel "Wald X2"
    }

    di as txt "----------------------------------------------------------------------------"
    di as txt %-18s "Variable" %10s "Coef." %9s "SE" %10s "`statlabel'" ///
        %9s "P>|`statlabel'|" %10s "LB" %10s "UB"

    local rn : rownames `matrix'
    local i=0
    foreach z of local rn {
        local ++i
        _apciflex_fmt4, matrix(`matrix') row(`i') col(1)
        local b "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(2)
        local se "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(3)
        local st "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(4) pvalue
        local p "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(5)
        local lb "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(6)
        local ub "`r(s)'"

        di as res %-18s "`z'" %10s "`b'" %9s "`se'" %10s "`st'" ///
            %9s "`p'" %10s "`lb'" %10s "`ub'"
    }

    if `testmatrix'[1,1]<. {
        _apciflex_fmt4, matrix(`testmatrix') row(1) col(1)
        local jt "`r(s)'"
        _apciflex_fmt4, matrix(`testmatrix') row(1) col(4) pvalue
        local jp "`r(s)'"

        local df1 = string(`testmatrix'[1,2],"%9.0f")
        local df1 = strtrim("`df1'")
        local df2 = string(`testmatrix'[1,3],"%12.0f")
        local df2 = strtrim("`df2'")

        if "`inference'"=="t" {
            local __pshow = cond("`jp'"=="<.0001","`jp'","= `jp'")
            di as txt "Joint controls: F(`df1', `df2') = " ///
                as res "`jt'" as txt ", p " as res "`__pshow'"
        }
        else {
            local __pshow = cond("`jp'"=="<.0001","`jp'","= `jp'")
            di as txt "Joint controls: Wald X2(`df1') = " ///
                as res "`jt'" as txt ", p " as res "`__pshow'"
        }
    }
end


program define _apciflex_print_group
    version 16.0
    syntax, MATRIX(name) GVAR(string) [ SECOND(string) Level(cilevel) INFERENCE(string) MTEST(string) NOLABEL ]

    if "`second'"=="" local second "Raw"
    if "`inference'"=="" local inference "t"
    if "`mtest'"=="" local mtest "none"

    local statlabel "t"
    local jointlabel "Joint F"
    if "`inference'"!="t" {
        local statlabel "z"
        local jointlabel "Joint X2"
    }

    local pjcol=12
    local ptag ""
    if "`mtest'"!="none" {
        local pjcol=13
        local ptag "*"
    }

    local vl ""
    if "`nolabel'"=="" & "`gvar'"!="" {
        capture confirm variable `gvar'
        if !_rc local vl : value label `gvar'
    }

    * Table A: average component and pointwise t/z test.
    di as txt "AVERAGE COMPONENT"
    di as txt "--------------------------------------------------------------------------------------"
    di as txt %-14s "Group" %5s "`second'" %6s "Cells" ///
        %9s "Average" %8s "SE" %10s "`statlabel'" %9s "Pavg" ///
        %9s "LB" %9s "UB"

    forvalues i=1/`=rowsof(`matrix')' {
        local code=`matrix'[`i',1]
        local lab ""
        if "`vl'"!="" local lab : label (`gvar') `code'
        if `"`lab'"'=="" local lab = string(`code',"%14.0g")

        _apciflex_fmt4, matrix(`matrix') row(`i') col(4)
        local av "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(5)
        local se "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(6)
        local st "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(7) pvalue
        local pa "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(8)
        local lb "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(9)
        local ub "`r(s)'"

        di as res %-14s `"`lab'"' ///
            %5.0f `matrix'[`i',2] ///
            %6.0f `matrix'[`i',3] ///
            %9s "`av'" ///
            %8s "`se'" ///
            %10s "`st'" ///
            %9s "`pa'" ///
            %9s "`lb'" ///
            %9s "`ub'"
    }

    * Table B: joint test of all raw components in the group.
    di as txt _newline "JOINT TEST"
    di as txt "-------------------------------------------------"
    di as txt %-14s "Group" %11s "`jointlabel'" %6s "df" %11s "Pjoint`ptag'"

    forvalues i=1/`=rowsof(`matrix')' {
        local code=`matrix'[`i',1]
        local lab ""
        if "`vl'"!="" local lab : label (`gvar') `code'
        if `"`lab'"'=="" local lab = string(`code',"%14.0g")

        _apciflex_fmt4, matrix(`matrix') row(`i') col(10)
        local js "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(`pjcol') pvalue
        local pj "`r(s)'"

        di as res %-14s `"`lab'"' ///
            %11s "`js'" ///
            %6.0f `matrix'[`i',11] ///
            %11s "`pj'"
    }
end


program define _apciflex_print_cohort
    version 16.0
    syntax, MATRIX(name) [ MTEST(string) INFERENCE(string) ]

    if "`mtest'"=="" local mtest "none"
    if "`inference'"=="" local inference "t"

    local pjcol=5
    local pacol=9
    local pscol=13
    local ptag ""
    if "`mtest'"!="none" {
        local pjcol=14
        local pacol=15
        local pscol=16
        local ptag "*"
    }

    local jointlabel "Joint F"
    if "`inference'"!="t" local jointlabel "Joint X2"

    di as txt "------------------------------------------------------------------------------------------------"
    di as txt %-10s "Cohort" %6s "Cells" %9s "`jointlabel'" %4s "df" %9s "Pjoint`ptag'" ///
        %9s "Average" %8s "SEavg" %9s "Pavg`ptag'" ///
        %9s "Slope" %8s "SE" %9s "Pslope`ptag'"

    forvalues i=1/`=rowsof(`matrix')' {
        _apciflex_fmt4, matrix(`matrix') row(`i') col(3)
        local js "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(`pjcol') pvalue
        local jp "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(6)
        local av "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(7)
        local ase "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(`pacol') pvalue
        local ap "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(10)
        local sl "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(11)
        local sse "`r(s)'"
        _apciflex_fmt4, matrix(`matrix') row(`i') col(`pscol') pvalue
        local sp "`r(s)'"

        local clab=string(`matrix'[`i',1],"%10.0g")
        di as res %-10s "`clab'" ///
            %6.0f `matrix'[`i',2] ///
            %9s "`js'" ///
            %4.0f `matrix'[`i',4] ///
            %9s "`jp'" ///
            %9s "`av'" ///
            %8s "`ase'" ///
            %9s "`ap'" ///
            %9s "`sl'" ///
            %8s "`sse'" ///
            %9s "`sp'"
    }

    if "`mtest'"!="none" {
        di as txt "* adjusted separately within joint, average, and slope families by mtest(`mtest')."
    }
end


mata:

real scalar _apciflex_connected(
    string scalar agevar,
    string scalar periodvar,
    string scalar touse)
{
    real colvector use, idx, a, p, alev, plev, ai, pi
    real matrix z, G
    real scalar K, A, P, k, i, j

    use=st_data(.,touse)
    idx=selectindex(use:==1)
    a=st_data(idx,agevar)
    p=st_data(idx,periodvar)
    z=uniqrows(sort((a,p),(1,2)))

    K=rows(z)
    if (K<=1) return(1)

    alev=uniqrows(sort(z[.,1],1))
    plev=uniqrows(sort(z[.,2],1))
    A=rows(alev)
    P=rows(plev)

    ai=J(K,1,.)
    pi=J(K,1,.)
    for (k=1;k<=K;k++) {
        for (i=1;i<=A;i++) if (z[k,1]==alev[i]) ai[k]=i
        for (j=1;j<=P;j++) if (z[k,2]==plev[j]) pi[k]=j
    }

    // A bipartite Age-Period support graph is connected iff the
    // incidence-style additive design has rank A+P-1.
    // This avoids the previous repeated K-by-K flood-fill scans.
    G=J(K,A+P,0)
    for (k=1;k<=K;k++) {
        G[k,ai[k]]=1
        G[k,A+pi[k]]=1
    }

    return(rank(G)==A+P-1)
}




void _apciflex_surface(
    string scalar B0name,
    string scalar V0name,
    string scalar Lname,
    string scalar MUname,
    string scalar Vname)
{
    real rowvector b
    real matrix V0, L, V
    real colvector mu

    b=st_matrix(B0name)
    V0=st_matrix(V0name)
    L=st_matrix(Lname)

    mu=L*b'
    V=L*V0*L'
    V=(V+V')/2

    st_matrix(MUname,mu)
    st_matrix(Vname,V)
}


void _apciflex_cellinfo(
    string scalar cellvar,
    string scalar touse,
    string scalar agevar,
    string scalar periodvar,
    string scalar cohortvar,
    string scalar agegroupvar,
    string scalar periodgroupvar,
    string scalar cohortgroupvar,
    string scalar nrawvar,
    string scalar neffvar,
    string scalar wsumvar,
    string scalar OUTname)
{
    real colvector c, use, age, period, cohort, ag, pg, cg
    real colvector nr, ne, ws, idx
    real matrix out
    real scalar K, k

    c=st_data(.,cellvar)
    use=st_data(.,touse)
    age=st_data(.,agevar)
    period=st_data(.,periodvar)
    cohort=st_data(.,cohortvar)
    ag=st_data(.,agegroupvar)
    pg=st_data(.,periodgroupvar)
    cg=st_data(.,cohortgroupvar)
    nr=st_data(.,nrawvar)
    ne=st_data(.,neffvar)
    ws=st_data(.,wsumvar)

    K=max(select(c,use:==1))
    out=J(K,9,.)

    for (k=1;k<=K;k++) {
        idx=selectindex((use:==1):&(c:==k))
        out[k,.]=(age[idx[1]],period[idx[1]],cohort[idx[1]],
            nr[idx[1]],ne[idx[1]],ws[idx[1]],
            ag[idx[1]],pg[idx[1]],cg[idx[1]])
    }

    st_matrix(OUTname,out)
}



real scalar _apciflex_crit(
    real scalar level,
    string scalar inference,
    real scalar df)
{
    real scalar alpha2
    alpha2=(100-level)/200

    if (inference=="t") return(invttail(df,alpha2))
    return(invnormal(1-alpha2))
}


real scalar _apciflex_twosided_p(
    real scalar stat,
    string scalar inference,
    real scalar df)
{
    if (inference=="t") return(2*ttail(df,abs(stat)))
    return(2*normal(-abs(stat)))
}


void _apciflex_decompose(
    string scalar CELLINFOname,
    string scalar MUname,
    string scalar Vname,
    string scalar DFname,
    string scalar inference,
    real scalar level,
    string scalar AGEname,
    string scalar VAGEname,
    string scalar PERname,
    string scalar VPERname,
    string scalar THname,
    string scalar VTHname,
    string scalar CELLname,
    string scalar GRANDname,
    string scalar GRANDSEname,
    string scalar GRANDPname)
{
    real matrix ci, V, G, Cc, Q, RHS, Tall, T, Vd, R, Vth
    real matrix ageout, perout, cell
    real colvector mu, alev, plev, ai, pi, da, dp, d, th
    real colvector vdiag, thdiag, coefindex
    real scalar K, A, P, q, k, i, j, df, crit
    real scalar grand, gse, gp, vv, se, st, pv

    ci=st_matrix(CELLINFOname)
    mu=st_matrix(MUname)
    V=st_matrix(Vname)
    V=(V+V')/2
    df=st_numscalar(DFname)

    K=rows(ci)
    alev=uniqrows(sort(ci[.,1],1))
    plev=uniqrows(sort(ci[.,2],1))
    A=rows(alev)
    P=rows(plev)

    ai=J(K,1,.)
    pi=J(K,1,.)
    da=J(A,1,0)
    dp=J(P,1,0)

    for (k=1;k<=K;k++) {
        for (i=1;i<=A;i++) {
            if (ci[k,1]==alev[i]) ai[k]=i
        }
        for (j=1;j<=P;j++) {
            if (ci[k,2]==plev[j]) pi[k]=j
        }
        da[ai[k]]=da[ai[k]]+1
        dp[pi[k]]=dp[pi[k]]+1
    }

    q=1+A+P
    G=J(K,q,0)
    G[.,1]=J(K,1,1)

    for (k=1;k<=K;k++) {
        G[k,1+ai[k]]=1
        G[k,1+A+pi[k]]=1
    }

    // Current APCIFLEX estimand: equal-cell observed-support projection.
    // Degree-weighted zero-sum constraints preserve the intercept as the
    // equal-cell mean of retained AP-cell adjusted means.
    Cc=J(2,q,0)
    Cc[1,2..(1+A)]=da'
    Cc[2,(2+A)..q]=dp'

    Q=J(q+2,q+2,0)
    Q[1..q,1..q]=G'*G
    Q[1..q,(q+1)..(q+2)]=Cc'
    Q[(q+1)..(q+2),1..q]=Cc

    RHS=J(q+2,K,0)
    RHS[1..q,.]=G'

    Tall=luinv(Q)*RHS
    T=Tall[1..q,.]

    d=T*mu
    Vd=T*V*T'
    Vd=(Vd+Vd')/2

    R=I(K)-G*T
    th=R*mu
    Vth=R*V*R'
    Vth=(Vth+Vth')/2

    crit=_apciflex_crit(level,inference,df)

    ageout=J(A,7,.)
    for (i=1;i<=A;i++) {
        vv=Vd[1+i,1+i]
        if (vv<0 & abs(vv)<1e-12) vv=0
        se=.
        if (vv>=0) se=sqrt(vv)
        st=.; pv=.
        if (se>0) {
            st=d[1+i]/se
            pv=_apciflex_twosided_p(st,inference,df)
        }
        ageout[i,.]=(alev[i],d[1+i],se,st,pv,
            d[1+i]-crit*se,d[1+i]+crit*se)
    }

    perout=J(P,7,.)
    for (j=1;j<=P;j++) {
        vv=Vd[1+A+j,1+A+j]
        if (vv<0 & abs(vv)<1e-12) vv=0
        se=.
        if (vv>=0) se=sqrt(vv)
        st=.; pv=.
        if (se>0) {
            st=d[1+A+j]/se
            pv=_apciflex_twosided_p(st,inference,df)
        }
        perout[j,.]=(plev[j],d[1+A+j],se,st,pv,
            d[1+A+j]-crit*se,d[1+A+j]+crit*se)
    }

    grand=d[1]
    vv=Vd[1,1]
    if (vv<0 & abs(vv)<1e-12) vv=0
    gse=.
    if (vv>=0) gse=sqrt(vv)
    gp=.
    if (gse>0) gp=_apciflex_twosided_p(grand/gse,inference,df)

    vdiag=diagonal(V)
    thdiag=diagonal(Vth)
    for (k=1;k<=K;k++) {
        if (vdiag[k]<0 & abs(vdiag[k])<1e-12) vdiag[k]=0
        if (thdiag[k]<0 & abs(thdiag[k])<1e-12) thdiag[k]=0
    }

    coefindex=1::K
    cell=(ci[.,1..6],mu,sqrt(vdiag),th,sqrt(thdiag),
        coefindex,ci[.,7..9])

    st_matrix(AGEname,ageout)
    st_matrix(VAGEname,Vd[2..(1+A),2..(1+A)])
    st_matrix(PERname,perout)
    st_matrix(VPERname,Vd[(2+A)..q,(2+A)..q])
    st_matrix(THname,th)
    st_matrix(VTHname,Vth)
    st_matrix(CELLname,cell)

    st_numscalar(GRANDname,grand)
    st_numscalar(GRANDSEname,gse)
    st_numscalar(GRANDPname,gp)
}


void _apciflex_global_test(
    string scalar THname,
    string scalar VTHname,
    string scalar DFRESname,
    string scalar inference,
    string scalar OUTname)
{
    real colvector th
    real matrix V, out
    real scalar df2, rk, wald, stat, p

    th=st_matrix(THname)
    V=st_matrix(VTHname)
    df2=st_numscalar(DFRESname)

    rk=rank(V)
    if (rk<=0) {
        out=(.,0,.,.)
        if (inference=="t") out[1,3]=df2
        st_matrix(OUTname,out)
        return
    }

    wald=th'*pinv(V)*th
    if (wald<0 & abs(wald)<1e-10) wald=0

    if (inference=="t") {
        stat=wald/rk
        p=Ftail(rk,df2,stat)
        out=(stat,rk,df2,p)
    }
    else {
        stat=wald
        p=chi2tail(rk,stat)
        out=(stat,rk,.,p)
    }

    st_matrix(OUTname,out)
}


real scalar _apciflex_unit_support(
    real matrix ci,
    real colvector idx,
    string scalar mode)
{
    real scalar W, Q

    if (mode=="equal") return(1)

    if (mode=="sample") {
        return(sum(ci[idx,4]))
    }

    if (mode=="population") {
        return(sum(ci[idx,6]))
    }

    // effective: pooled Kish effective N across disjoint retained cells.
    W=sum(ci[idx,6])
    Q=sum((ci[idx,6]:^2):/ci[idx,5])

    if (Q<=0) return(.)
    return((W^2)/Q)
}


real colvector _apciflex_cell_support(
    real matrix ci,
    real colvector idx,
    string scalar mode)
{
    if (mode=="equal") return(J(rows(idx),1,1))
    if (mode=="sample") return(ci[idx,4])
    if (mode=="effective") return(ci[idx,5])
    if (mode=="population") return(ci[idx,6])

    return(J(rows(idx),1,.))
}


real colvector _apciflex_adjust_p(
    real colvector p,
    string scalar method)
{
    real colvector out, idx, pv, ord, ps, adj
    real scalar m, i, v

    out=J(rows(p),1,.)
    idx=selectindex(p:<.)
    m=rows(idx)
    if (m==0) return(out)

    pv=p[idx]

    if (method=="none") {
        out[idx]=pv
        return(out)
    }

    if (method=="bonferroni") {
        for (i=1;i<=m;i++) out[idx[i]]=min((1,m*pv[i]))
        return(out)
    }

    ord=order(pv,1)
    ps=pv[ord]
    adj=J(m,1,.)

    if (method=="holm") {
        v=0
        for (i=1;i<=m;i++) {
            v=max((v,(m-i+1)*ps[i]))
            adj[i]=min((1,v))
        }
    }
    else if (method=="bh") {
        v=1
        for (i=m;i>=1;i--) {
            v=min((v,(m/i)*ps[i]))
            adj[i]=min((1,v))
        }
    }
    else {
        return(out)
    }

    for (i=1;i<=m;i++) {
        out[idx[ord[i]]]=adj[i]
    }

    return(out)
}


void _apciflex_cohort_table(
    string scalar CELLINFOname,
    string scalar THname,
    string scalar VTHname,
    string scalar DFname,
    string scalar inference,
    string scalar mode,
    string scalar mtest,
    string scalar OUTname)
{
    real matrix ci, V, out, Vk
    real colvector th, cohort, clev, idx, ord, est, x, cw
    real rowvector w, ws
    real scalar df, C, i, m, rk, wald, F, pj
    real scalar avg, vv, se, tt, pa
    real scalar slope, sv, sse, st, ps, den, sw

    ci=st_matrix(CELLINFOname)
    th=st_matrix(THname)
    V=st_matrix(VTHname)
    df=st_numscalar(DFname)

    cohort=ci[.,3]
    clev=uniqrows(sort(cohort,1))
    C=rows(clev)
    out=J(C,16,.)

    for (i=1;i<=C;i++) {
        idx=selectindex(cohort:==clev[i])

        // Order by actual age, not ordinal cell position.
        ord=order(ci[idx,1],1)
        idx=idx[ord]
        m=rows(idx)

        est=th[idx]
        Vk=V[idx,idx]
        rk=rank(Vk)

        F=.; pj=.
        if (rk>0) {
            wald=est'*pinv(Vk)*est
            if (wald<0 & abs(wald)<1e-10) wald=0
            if (inference=="t") {
                F=wald/rk
                pj=Ftail(rk,df,F)
            }
            else {
                F=wald
                pj=chi2tail(rk,F)
            }
        }

        // Raw cohort average uses the selected aggregation support.
        cw=_apciflex_cell_support(ci,idx,mode)
        sw=sum(cw)

        w=J(1,rows(th),0)
        if (sw>0) w[1,idx']=(cw:/sw)'

        avg=w*th
        vv=w*V*w'
        if (vv<0 & abs(vv)<1e-12) vv=0
        se=.
        if (vv>=0) se=sqrt(vv)
        tt=.; pa=.
        if (se>0) {
            tt=avg/se
            pa=_apciflex_twosided_p(tt,inference,df)
        }

        // Life-course slope remains an actual-age OLS contrast.
        // aggweight() changes averages, not the slope estimand.
        slope=.; sse=.; st=.; ps=.
        if (m>=2) {
            x=ci[idx,1]
            x=x:-mean(x)
            den=quadcross(x,x)

            if (den>0) {
                ws=J(1,rows(th),0)
                ws[1,idx']=x'/den

                slope=ws*th
                sv=ws*V*ws'
                if (sv<0 & abs(sv)<1e-12) sv=0
                if (sv>=0) sse=sqrt(sv)

                if (sse>0) {
                    st=slope/sse
                    ps=_apciflex_twosided_p(st,inference,df)
                }
            }
        }

        out[i,1..13]=(clev[i],m,F,rk,pj,avg,se,tt,pa,slope,sse,st,ps)
    }

    // Multiplicity adjustment is defined separately for the three distinct
    // hypothesis families: joint deviation, average deviation, and slope.
    out[.,14]=_apciflex_adjust_p(out[.,5],mtest)
    out[.,15]=_apciflex_adjust_p(out[.,9],mtest)
    out[.,16]=_apciflex_adjust_p(out[.,13],mtest)

    st_matrix(OUTname,out)
}


void _apciflex_group_ageperiod(
    real scalar type,
    string scalar CELLINFOname,
    string scalar RAWname,
    string scalar VRAWname,
    string scalar DFname,
    string scalar inference,
    real scalar level,
    string scalar mode,
    string scalar mtest,
    string scalar OUTname)
{
    real matrix ci, raw, V, out, Vg
    real colvector levels, effects, cellbase, cellgroup, rawgroup, groups
    real colvector idx, gidx, bg
    real rowvector w
    real scalar df, i, j, G, nlev, ncells, code, support, sw
    real scalar est, vv, se, st, p, crit
    real scalar rk, wald, jstat, jp

    ci=st_matrix(CELLINFOname)
    raw=st_matrix(RAWname)
    V=st_matrix(VRAWname)
    df=st_numscalar(DFname)

    levels=raw[.,1]
    effects=raw[.,2]

    if (type==1) {
        cellbase=ci[.,1]
        cellgroup=ci[.,7]
    }
    else {
        cellbase=ci[.,2]
        cellgroup=ci[.,8]
    }

    rawgroup=J(rows(levels),1,.)
    for (i=1;i<=rows(levels);i++) {
        for (j=1;j<=rows(ci);j++) {
            if (cellbase[j]==levels[i]) {
                rawgroup[i]=cellgroup[j]
                break
            }
        }
    }

    groups=uniqrows(sort(rawgroup,1))
    G=rows(groups)
    out=J(G,13,.)
    crit=_apciflex_crit(level,inference,df)

    for (i=1;i<=G;i++) {
        code=groups[i]
        w=J(1,rows(levels),0)
        nlev=0
        sw=0

        for (j=1;j<=rows(levels);j++) {
            if (rawgroup[j]==code) {
                nlev=nlev+1
                idx=selectindex(cellbase:==levels[j])
                support=_apciflex_unit_support(ci,idx,mode)
                w[j]=support
                sw=sw+support
            }
        }

        if (sw>0) w=w:/sw
        ncells=sum(cellgroup:==code)

        // Average group component: H0 w'b = 0.
        est=w*effects
        vv=w*V*w'
        if (vv<0 & abs(vv)<1e-12) vv=0
        se=.
        if (vv>=0) se=sqrt(vv)
        st=.; p=.
        if (se>0) {
            st=est/se
            p=_apciflex_twosided_p(st,inference,df)
        }

        // Joint existence test:
        // H0: every retained raw Age/Period component in this theory group = 0.
        gidx=selectindex(rawgroup:==code)
        bg=effects[gidx]
        Vg=V[gidx,gidx]
        rk=rank(Vg)
        jstat=.; jp=.
        if (rk>0) {
            wald=bg'*pinv(Vg)*bg
            if (wald<0 & abs(wald)<1e-10) wald=0

            if (inference=="t") {
                jstat=wald/rk
                jp=Ftail(rk,df,jstat)
            }
            else {
                jstat=wald
                jp=chi2tail(rk,jstat)
            }
        }

        out[i,1..12]=(code,nlev,ncells,est,se,st,p,
            est-crit*se,est+crit*se,jstat,rk,jp)
    }

    // mtest() for the group-joint family only. Age and Period call this
    // function separately, so their multiplicity families remain separate.
    out[.,13]=_apciflex_adjust_p(out[.,12],mtest)

    st_matrix(OUTname,out)
}


void _apciflex_group_cohort(
    string scalar CELLINFOname,
    string scalar THname,
    string scalar VTHname,
    string scalar DFname,
    string scalar inference,
    real scalar level,
    string scalar mode,
    string scalar mtest,
    string scalar OUTname)
{
    real matrix ci, V, out, L, Vd
    real colvector th, cohort, cg, groups, cidx, gcohorts, cellw, d
    real rowvector w
    real scalar df, G, i, j, code, nlev, ncells, cval
    real scalar est, vv, se, st, p, crit
    real scalar cohort_support, total_cohort_support, swcell
    real scalar rk, wald, jstat, jp

    ci=st_matrix(CELLINFOname)
    th=st_matrix(THname)
    V=st_matrix(VTHname)
    df=st_numscalar(DFname)

    cohort=ci[.,3]
    cg=ci[.,9]
    groups=uniqrows(sort(cg,1))
    G=rows(groups)
    out=J(G,13,.)
    crit=_apciflex_crit(level,inference,df)

    for (i=1;i<=G;i++) {
        code=groups[i]

        gcohorts=uniqrows(sort(select(cohort,cg:==code),1))
        nlev=rows(gcohorts)
        ncells=sum(cg:==code)

        w=J(1,rows(th),0)
        L=J(nlev,rows(th),0)
        total_cohort_support=0

        // Between-cohort support used only for the group-average estimand.
        for (j=1;j<=nlev;j++) {
            cval=gcohorts[j]
            cidx=selectindex(cohort:==cval)
            total_cohort_support=total_cohort_support +
                _apciflex_unit_support(ci,cidx,mode)
        }

        // First construct each raw cohort's average deviation:
        // Delta_c = sum_j w_{cj} theta_{cj}.
        // The same within-cohort weighting rule is used by e(cohort_raw).
        for (j=1;j<=nlev;j++) {
            cval=gcohorts[j]
            cidx=selectindex(cohort:==cval)

            cohort_support=_apciflex_unit_support(ci,cidx,mode)
            cellw=_apciflex_cell_support(ci,cidx,mode)
            swcell=sum(cellw)

            if (swcell>0) {
                L[j,cidx']=(cellw:/swcell)'
            }

            // Group average is a second-stage weighted average of raw cohorts.
            if (total_cohort_support>0 & swcell>0) {
                w[1,cidx']=
                    (cohort_support/total_cohort_support) :*
                    (cellw:/swcell)'
            }
        }

        // Average cohort-group deviation: H0 w'theta = 0.
        est=w*th
        vv=w*V*w'
        if (vv<0 & abs(vv)<1e-12) vv=0
        se=.
        if (vv>=0) se=sqrt(vv)
        st=.; p=.
        if (se>0) {
            st=est/se
            p=_apciflex_twosided_p(st,inference,df)
        }

        // Joint cohort-group existence test:
        // H0: Delta_c = 0 for every retained raw cohort c in the theory group.
        // This deliberately tests raw-cohort average deviations, not every
        // individual theta_ap cell on the covered Lexis diagonals.
        d=L*th
        Vd=L*V*L'
        rk=rank(Vd)
        jstat=.; jp=.
        if (rk>0) {
            wald=d'*pinv(Vd)*d
            if (wald<0 & abs(wald)<1e-10) wald=0

            if (inference=="t") {
                jstat=wald/rk
                jp=Ftail(rk,df,jstat)
            }
            else {
                jstat=wald
                jp=chi2tail(rk,jstat)
            }
        }

        out[i,1..12]=(code,nlev,ncells,est,se,st,p,
            est-crit*se,est+crit*se,jstat,rk,jp)
    }

    // Cohort groups form their own joint-test multiplicity family.
    out[.,13]=_apciflex_adjust_p(out[.,12],mtest)

    st_matrix(OUTname,out)
}

end
