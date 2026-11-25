clear all
set more off
program drop _all
mata: mata clear
version 16
args root core
if "`root'" == "" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"
quietly do "`root'/ado/sdea.ado"

program define _run_sdea_drift_case, rclass
    version 16.0
    syntax, Case(string) [ALGorithm(string)]
    if "`algorithm'" == "" local algorithm "frontier"

    clear
    set obs 360
    gen id = _n
    gen period = ceil(_n/30)
    gen double x1 = runiform(10, 100)
    gen double x2 = runiform(5, 50)
    gen double tech = 1 + .01 * period

    if "`case'" == "support" {
        replace x1 = runiform(60, 160) if period >= 7
        replace x2 = runiform(30, 80) if period >= 7
    }
    if "`case'" == "drift" {
        replace tech = tech * 1.35 if period >= 7
    }

    gen double y = tech * x1^.40 * x2^.50 * exp(-abs(rnormal(0,.08)))

    sdea y, inputs(x1 x2) time(period) dmu(id) method(buildhull) ///
        starters(30) drifttest driftalgo(`algorithm') driftwindow(3) ///
        thrvol(1.5) threxp(.5) thrturn(.7) ///
        generate(alpha frame drift) replace

    local n_drifts = r(N_drifts)
    quietly levelsof period if drift > 0 & !missing(drift), local(drift_periods) clean
    quietly count if frame == 1
    local n_frame = r(N)
    quietly summarize alpha, meanonly

    di as text "CASE `case': period_drifts=`n_drifts' frame=`n_frame' mean_alpha=" ///
        as result %9.6f r(mean) as text " drift_periods=`drift_periods'"

    if "`algorithm'" == "multi" {
        preserve
        keep period _sdea_drift_V _sdea_drift_D _sdea_drift_T _sdea_drift_triggers drift
        collapse (mean) V=_sdea_drift_V D=_sdea_drift_D T=_sdea_drift_T ///
            triggers=_sdea_drift_triggers drift=drift, by(period)
        di as text "Diagnostics for `case' (`algorithm')"
        list period V D T triggers drift if period >= 4, noobs abbrev(12)
        restore
    }
    if "`algorithm'" == "frontier" {
        preserve
        keep period _sdea_drift_frontier drift
        collapse (mean) frontier=_sdea_drift_frontier drift=drift, by(period)
        di as text "Diagnostics for `case' (`algorithm')"
        list period frontier drift if period >= 4, noobs abbrev(12)
        restore
    }

    return scalar N_drifts = `n_drifts'
    return scalar N_frame = `n_frame'
    return scalar mean_alpha = r(mean)
    return local drift_periods "`drift_periods'"
end

set seed 20260528
_run_sdea_drift_case, case(stationary) algorithm(frontier)
local stationary_drifts = r(N_drifts)
local stationary_periods "`r(drift_periods)'"

set seed 20260528
_run_sdea_drift_case, case(support) algorithm(frontier)
local support_drifts = r(N_drifts)
local support_periods "`r(drift_periods)'"

set seed 20260528
_run_sdea_drift_case, case(drift) algorithm(frontier)
local drift_drifts = r(N_drifts)
local drift_periods "`r(drift_periods)'"

di as text _n "Concept drift scenario validation"
di as text "  stationary DGP:       drifts=`stationary_drifts' periods=`stationary_periods'"
di as text "  support expansion:    drifts=`support_drifts' periods=`support_periods'"
di as text "  technology drift:     drifts=`drift_drifts' periods=`drift_periods'"

assert `stationary_drifts' == 0
assert `support_drifts' == 0
assert `drift_drifts' > 0
