*! apciflex_p 1.0.0 08sep2026
*! optimized predict program for apciflex
*! associative-array AP lookup; lazy e(V); explicit retained-support guard

program define apciflex_p
    version 16.0

    if "`e(cmd)'"!="apciflex" {
        error 301
    }

    syntax newvarname [if] [in] [, XB STDP PR RESiduals APCXB APCPR ]

    local nopt=0
    if "`xb'"!="" local ++nopt
    if "`stdp'"!="" local ++nopt
    if "`pr'"!="" local ++nopt
    if "`residuals'"!="" local ++nopt
    if "`apcxb'"!="" local ++nopt
    if "`apcpr'"!="" local ++nopt

    if `nopt'>1 {
        di as err "only one of xb, stdp, pr, residuals, apcxb, or apcpr may be specified."
        exit 198
    }

    local link "`e(link)'"
    if "`link'"=="" local link "identity"

    if ("`pr'"!="" | "`apcpr'"!="") & "`link'"=="identity" {
        di as err "option pr/apcpr is available only after link(logit) or link(probit)."
        exit 198
    }

    local mode "xb"
    if "`stdp'"!="" local mode "stdp"
    else if "`residuals'"!="" local mode "residuals"
    else if "`pr'"!="" local mode "pr"
    else if "`apcpr'"!="" local mode "pr"
    else if "`apcxb'"!="" local mode "xb"
    else if `nopt'==0 & inlist("`link'","logit","probit") local mode "pr"

    marksample touse, novarlist

    local avar "`e(age_var)'"
    local pvar "`e(period_var)'"
    local controls "`e(controls)'"
    local depvar "`e(depvar)'"

    * Compute into a temporary variable so a support error does not leave
    * a partially generated user variable behind.
    tempvar __apcpred
    quietly gen double `__apcpred'=.

    mata: _apciflex_p_core( ///
        "`__apcpred'","`touse'","`avar'","`pvar'", ///
        "`controls'","`depvar'","`mode'","`link'")

    rename `__apcpred' `varlist'
    label variable `varlist' "APCIFLEX prediction: `mode'"
end

mata:

void _apciflex_p_core(
    string scalar outvar,
    string scalar tousev,
    string scalar avar,
    string scalar pvar,
    string scalar ctrls,
    string scalar depvar,
    string scalar mode,
    string scalar link)
{
    real matrix C, V, X, Xh, Vcc
    real colvector b, a, p, eta, out, mu, y, idx, qv, bad, kidx, ks
    real rowvector Vkc
    string rowvector zvars
    real scalar K, nc, n, k, m, vkk, i, j
    transmorphic APmap

    C=st_matrix("e(cell)")
    b=st_matrix("e(b)")'

    K=st_numscalar("e(K)")
    zvars=tokens(ctrls)
    nc=cols(zvars)

    if (rows(C)!=K | rows(b)<K+nc) {
        errprintf("internal APCIFLEX prediction mapping is inconsistent with e(b).\n")
        _error(498)
    }

    a=st_data(.,avar,tousev)
    p=st_data(.,pvar,tousev)
    n=rows(a)

    /*
      Map (age, period) -> retained AP-cell coefficient index once.
      A two-dimensional real key avoids arbitrary scalar-key multipliers.
    */
    APmap=asarray_create("real",2)
    asarray_notfound(APmap,0)

    for (k=1; k<=K; k++) {
        asarray(APmap,(C[k,1],C[k,2]),k)
    }

    kidx=J(n,1,0)
    for (i=1; i<=n; i++) {
        if (a[i]<. & p[i]<.) {
            kidx[i]=asarray(APmap,(a[i],p[i]))
        }
    }

    bad=selectindex((a:<.) :& (p:<.) :& (kidx:==0))
    if (rows(bad)>0) {
        errprintf("APCIFLEX prediction requested %g observation(s) outside the retained Age x Period support.\n", rows(bad))
        errprintf("The fitted AP-cell surface is undefined for these Age x Period combinations.\n")
        errprintf("Use predict ... if e(sample), retained-cell margins such as margins, over(age period),\n")
        errprintf("or stored APC-I component matrices for Age/Period/cohort inference.\n")
        _error(459)
    }

    if (nc>0) {
        X=st_data(.,zvars,tousev)
    }
    else {
        X=J(n,0,.)
    }

    if (mode=="stdp") {
        /*
          e(V) is needed only for stdp.  xb/pr/residuals no longer copy the
          full covariance matrix on every prediction call.
        */
        V=st_matrix("e(V)")
        if (rows(V)<K+nc | cols(V)<K+nc) {
            errprintf("internal APCIFLEX prediction mapping is inconsistent with e(V).\n")
            _error(498)
        }

        out=J(n,1,.)
        if (nc>0) Vcc=V[(K+1)..(K+nc),(K+1)..(K+nc)]

        idx=selectindex(kidx:>0)
        if (rows(idx)>0) {
            ks=uniqrows(sort(kidx[idx],1))

            for (j=1; j<=rows(ks); j++) {
                k=ks[j]
                idx=selectindex(kidx:==k)
                m=rows(idx)
                if (m==0) continue

                vkk=V[k,k]

                if (nc==0) {
                    qv=J(m,1,vkk)
                }
                else {
                    Xh=X[idx,.]
                    Vkc=V[k,(K+1)..(K+nc)]
                    qv=J(m,1,vkk) +
                        2:*Xh*Vkc' +
                        rowsum((Xh*Vcc):*Xh)
                }

                out[idx]=qv
            }
        }

        idx=selectindex(out:<0 :& abs(out):<1e-12)
        if (rows(idx)>0) out[idx]=J(rows(idx),1,0)

        idx=selectindex(out:<0 :& out:<.)
        if (rows(idx)>0) out[idx]=J(rows(idx),1,.)

        out=sqrt(out)
        st_store(.,outvar,tousev,out)
        return
    }

    eta=J(n,1,.)
    idx=selectindex(kidx:>0)
    if (rows(idx)>0) {
        eta[idx]=b[kidx[idx]]
    }

    if (nc>0) {
        eta=eta + X*b[(K+1)..(K+nc)]
    }

    if (link=="logit") {
        mu=1:/(1:+exp(-eta))
    }
    else if (link=="probit") {
        mu=normal(eta)
    }
    else {
        mu=eta
    }

    if (mode=="residuals") {
        y=st_data(.,depvar,tousev)
        out=y-mu
    }
    else if (mode=="pr") {
        out=mu
    }
    else {
        out=eta
    }

    st_store(.,outvar,tousev,out)
}

end
