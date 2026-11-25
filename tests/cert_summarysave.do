*! cert_summarysave.do
*! Certification for sdea summarysave() monitoring output

clear all
set more off
program drop _all
mata: mata clear

version 16
args root core
if "`root'" == "" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"
quietly do "`root'/ado/sdea.ado"

set seed 12062026
set obs 60
gen long dmu = _n
gen int period = ceil(_n/10)
gen double x1 = runiform(20, 100)
gen double x2 = runiform(10, 60)
gen double y1 = 0.7*x1 + 0.2*x2 + runiform(5, 30)
gen double y2 = 0.3*x1 + 0.6*x2 + runiform(5, 25)

tempfile monitor
sdea y1 y2, inputs(x1 x2) time(period) method(buildhull) ///
    starters(12) gate(update_gate) summarysave("`monitor'") replace

preserve
use "`monitor'", clear

assert _N == 6
confirm variable period
confirm variable N_obs
confirm variable mean_score
confirm variable max_score
confirm variable frame_members
confirm variable frame_entries
confirm variable frame_pruned
confirm variable mean_drift
confirm variable max_drift
confirm variable drift_points
confirm variable max_frame_size
confirm variable gate_enveloped
confirm variable gate_updates
confirm variable gate_rejected
confirm variable gate_starters

summarize N_obs, meanonly
assert r(sum) == 60

egen long gate_total = rowtotal(gate_enveloped gate_updates gate_rejected gate_starters)
assert gate_total == N_obs

summarize max_frame_size, meanonly
assert r(max) >= 12

restore

di as result "cert_summarysave.do passed"
