clear all
set more off
program drop _all
mata: mata clear
version 16
args root core
if "`root'" == "" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"
quietly do "`root'/ado/sdea.ado"

program define _make_hetero_data
    version 16.0
    syntax, Scenario(string)

    clear
    set obs 480
    gen id = _n
    gen period = ceil(_n/40)
    gen group = 1

    if "`scenario'" == "composition" {
        bysort period: gen within = _n
        replace group = 2 if period < 7 & within > 30
        replace group = 2 if period >= 7 & within > 10
        drop within
    }
    else {
        bysort period: gen within = _n
        replace group = 2 if within > 20
        drop within
    }

    gen double x1 = runiform(10, 100)
    gen double x2 = runiform(5, 50)
    gen double tech = 1
    replace tech = 1.45 if group == 2

    if "`scenario'" == "within_drift" {
        replace tech = tech * 1.35 if group == 2 & period >= 7
    }

    gen double y = tech * x1^.40 * x2^.50 * exp(-abs(rnormal(0,.06)))
end

program define _run_case, rclass
    version 16.0
    syntax, Scenario(string) [BYGROUP]

    _make_hetero_data, scenario(`scenario')
    local byopt ""
    if "`bygroup'" != "" local byopt "by(group)"

    sdea y, inputs(x1 x2) time(period) dmu(id) method(buildhull) ///
        starters(40) drifttest driftalgo(frontier) driftwindow(3) ///
        thrvol(1.5) `byopt' generate(alpha frame drift) replace

    local n_drifts = r(N_drifts)
    quietly levelsof period if drift > 0 & !missing(drift), local(drift_periods) clean

    di as text "SCENARIO `scenario' `byopt': drifts=`n_drifts' periods=`drift_periods'"

    return scalar N_drifts = `n_drifts'
    return local drift_periods "`drift_periods'"
end

set seed 20260528
_run_case, scenario(composition)
local comp_pooled = r(N_drifts)

set seed 20260528
_run_case, scenario(composition) bygroup
local comp_by = r(N_drifts)

set seed 20260528
_run_case, scenario(within_drift) bygroup
local within_by = r(N_drifts)
local within_periods "`r(drift_periods)'"

di as text _n "Concept drift heterogeneity validation"
di as text "  composition shift, pooled:     drifts=`comp_pooled'"
di as text "  composition shift, by(group):  drifts=`comp_by'"
di as text "  within-group drift, by(group): drifts=`within_by' periods=`within_periods'"

assert `comp_by' == 0
assert `within_by' > 0
