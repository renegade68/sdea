version 16
args root core output
if "`root'" == "" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"
clear
set rng mt64
set seed 20260528
set obs 360
generate long arrival = _n
generate int period = ceil(arrival/30)
generate double x1 = runiform(10,100)
generate double x2 = runiform(5,50)
generate double y = (1+.01*period) * cond(period>=7,1.35,1) * ///
    x1^.40 * x2^.50 * exp(-abs(rnormal(0,.08)))
sort period arrival
sdea y, inputs(x1 x2) time(period) dmu(arrival) method(buildhull) ///
    starters(30) drifttest driftalgo(frontier) driftwindow(3) ///
    thrvol(1.5) threxp(.5) thrturn(.7) generate(alpha frame drift)
summarize alpha frame drift
levelsof period if drift>0 & !missing(drift)
if "`output'" != "" export delimited arrival period alpha frame drift using "`output'/streaming.csv"
