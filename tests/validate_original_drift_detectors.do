clear all
set more off
program drop _all
mata: mata clear
version 16
args root core
if "`root'" == "" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"
quietly do "`root'/ado/sdea.ado"

program define _run_original_detector, rclass
    version 16.0
    syntax, ALGorithm(string)

    clear
    set obs 800
    gen id = _n
    gen period = ceil(_n/40)
    gen double x1 = runiform(10, 100)
    gen double x2 = runiform(5, 50)
    gen double y = x1^.45 * x2^.45 * exp(-abs(rnormal(0,.05)))

    gen byte err = (runiform() < .02)
    replace err = (runiform() < .35) if period >= 11

    sdea y, inputs(x1 x2) time(period) dmu(id) method(buildhull) ///
        starters(40) drifttest driftalgo(`algorithm') ///
        driftmode(original) drifterror(err) ///
        generate(alpha frame drift) replace

    quietly count if drift > 0 & !missing(drift)
    local n_drift_obs = r(N)
    quietly levelsof period if drift > 0 & !missing(drift), local(drift_periods) clean
    quietly summarize err if period < 11, meanonly
    local p0 = r(mean)
    quietly summarize err if period >= 11, meanonly
    local p1 = r(mean)

    di as text "ORIGINAL `algorithm': error_before=" as result %6.3f `p0' ///
        as text " error_after=" as result %6.3f `p1' ///
        as text " drift_obs=" as result `n_drift_obs' ///
        as text " drift_periods=`drift_periods'"

    return scalar N_drift_obs = `n_drift_obs'
    return scalar p_before = `p0'
    return scalar p_after = `p1'
    return local drift_periods "`drift_periods'"
end

set seed 20260528
_run_original_detector, algorithm(ddm)
local ddm_drifts = r(N_drift_obs)
local ddm_periods "`r(drift_periods)'"

set seed 20260528
_run_original_detector, algorithm(eddm)
local eddm_drifts = r(N_drift_obs)
local eddm_periods "`r(drift_periods)'"

di as text _n "Original DDM/EDDM validation"
di as text "  DDM:  drift_obs=`ddm_drifts' periods=`ddm_periods'"
di as text "  EDDM: drift_obs=`eddm_drifts' periods=`eddm_periods'"

assert `ddm_drifts' > 0
assert `eddm_drifts' > 0
