version 16.0
clear all
set more off
program drop _all
mata: mata clear

version 16
args root core
if "`root'" == "" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"
quietly do "`root'/ado/sdea.ado"

set obs 80
set seed 12062026
gen long id = _n
gen int t = _n
gen double x1 = 1 + runiform()*4
gen double x2 = 1 + runiform()*3
gen double y1 = 1.5*x1 + 0.8*x2 + rnormal(0, .05)
replace y1 = y1 + 1.5 if inlist(id, 45, 60, 75)
gen double frontier_score = cond(inlist(id, 45, 60, 75), .95, .20)

sdea y1, inputs(x1 x2) time(t) method(buildhull) starters(20) ///
    generate(alpha0 frame0 drift0) replace
assert !missing(alpha0)
assert !missing(frame0)
assert !missing(drift0)

sdea y1, inputs(x1 x2) time(t) method(buildhull) starters(20) ///
    gate(update_gate) mlscore(frontier_score) mlthreshold(.50) ///
    generate(alpha1 frame1 drift1) replace
assert !missing(alpha1)
assert !missing(frame1)
assert !missing(drift1)
assert inlist(update_gate, 0, 1, 2, 3)
assert r(N_gate_starters) == 20
assert r(N_gate_updates) >= 1
assert r(N_gate_rejected) >= 0

di as result "cert_frontier_gate.do passed"
