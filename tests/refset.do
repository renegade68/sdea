version 16
args root core
if "`root'" == "" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"

clear
set obs 72
set seed 56124
gen long id = _n
gen str12 labelid = string(mod(id,7))
gen double x = runiform() + .1
gen double x2 = 1e9*(runiform()+.1)
gen double y = 1e-9*(runiform()+.1)
gen double y2 = runiform()+.1
gen double zero = 0
gen byte period = mod(id,3)
gen int _sdea_frame_size = id
gen byte sample = period!=1 & inrange(id,5,69)
local shape "inputs(x x2 zero) method(standard) model(te) dmu(labelid)"
local out "generate(te member drift)"
foreach rts in vrs crs {
    foreach ort in input output {
        quietly dea x x2 zero = y y2 if sample, rts(`rts') orientation(`ort') ///
            stage(1) engine(simplex) dmu(id) generate(reference)
        foreach variant in modified dula {
            foreach solver in simplex ipm auto {
                local extra
                if "`variant'"=="modified" local extra prune(5)
                quietly sdea y y2 if period!=1 in 5/69, `shape' refset(buildhull) ///
                    orientation(`ort') rts(`rts') bhvariant(`variant') solver(`solver') `extra' `out'
                assert "`r(refset)'"=="buildhull" & "`r(score_solver)'"=="simplex"
                assert "`r(hull_solver)'"=="`solver'"
                if "`solver'"=="simplex" assert r(N_hull_simplex)>0 & r(N_hull_ipm)==0
                else assert r(N_hull_ipm)>0
                assert r(N_hull_failed)>=0
                quietly count if sample & abs(te-reference)>1e-7
                assert r(N)==0
                assert inlist(member,0,1) & drift==0 if sample
                assert missing(te) & missing(member) & missing(drift) if !sample
                assert id==_n & _sdea_frame_size==id
                drop te member drift
            }
        }
        drop reference
        if "`rts'"=="vrs" drop CRS_TE VRS_TE SCALE RTS_NUM RTS
    }
}

* Ties, zero-valued coordinates and two different input proportions.
clear
input double(x x2 y)
1 3 5
1 3 5
2 4 4
3 5 2
end
local opts "inputs(x x2) method(standard) refset(buildhull) rts(vrs) generate(te member drift)"
quietly sdea y, `opts' orientation(input)
assert r(N_frame)==1 & r(N_candidates)==1 & r(N_efficient)==2
assert te==1 in 1/2
assert abs(te-.75)<1e-7 in 3
assert member==(_n==1)
quietly sdea y, `opts' orientation(output) replace
assert abs(te-y/5)<1e-7
quietly sdea y in 3, `opts' replace
assert r(N)==1 & r(N_frame)==1
assert te==1 & member==1 in 3
assert missing(te) if _n!=3
capture sdea y, `opts'
assert _rc==110
capture sdea y, inputs(x x2) method(standard) refset(buildhull) gen(x dummy dummy2) replace
assert _rc!=0 & x[1]==1
capture sdea y, inputs(x x2) method(standard) refset(buildhull) gen(dup dup drift) replace
assert _rc!=0 & te[3]==1
capture sdea y, inputs(x x2) method(standard) refset(buildhull) time(x)
assert _rc==198
capture sdea y, inputs(x x2) method(standard) refset(buildhull) model(sbm)
assert _rc==198
capture sdea y, inputs(x x2) method(buildhull) time(x) refset(buildhull)
assert _rc==198
gen long firm = _n
quietly sdea y, inputs(x x2) method(standard) refset(buildhull) generate(f m d) replace
assert firm==_n

* Standard full-sample and period-specific behaviour with and without the new option.
gen byte year = 1+(_n>2)
quietly sdea y, inputs(x x2) method(standard) time(year) generate(orig oframe odrift)
quietly sdea y, inputs(x x2) method(standard) time(year) refset(full) generate(copy cf cd)
assert "`r(refset)'"=="full"
assert orig==copy & oframe==cf & odrift==cd

* Individual zero coordinates require feasibility in each input and output.
clear
input double(x1 x2 y1 y2)
0 1 1 0
1 0 0 1
0 2 1 0
2 0 0 1
end
quietly sdea y1 y2, inputs(x1 x2) method(standard) refset(buildhull) ///
    orientation(input) rts(vrs) generate(te m d)
assert r(N_frame)==2
assert te==1 in 1/2
assert abs(te-.5)<1e-7 in 3/4
gen double old_te = te
capture sdea y1 y2 if x1<0, inputs(x1 x2) method(standard) refset(buildhull) ///
    generate(te m d) replace
assert _rc==2000 & old_te==te

* Dense VRS frontier, and score invariance to reversal of the observation order.
clear
set obs 48
gen long id = _n
gen double x = id
gen double y = sqrt(x)
gen double xzero = 0
local dense "inputs(x xzero) method(standard) refset(buildhull) rts(vrs) orientation(input)"
quietly sdea y, `dense' generate(te nf dr)
assert te>=1-1e-7 & te<=1 & nf==1
gsort -id
quietly sdea y, `dense' generate(rev rf rd)
assert abs(te-rev)<=1e-7 & id==_N+1-_n
replace xzero=-1 in 2
capture sdea y, `dense' generate(te nf dr) replace
assert _rc==459 & nf[2]==1 & te[2]>=1-1e-7

* Exercise conservative recovery using deliberately incomplete proposed frames.
clear
input double(x y) byte(candidate member)
2 1 1 1
1 2 1 0
end
gen double score = .
program drop _all
mata: mata clear
quietly do "`root'/ado/sdea.ado"
mata: st_matrix("checks",_sdea_static_evaluate("x","y","x","y","candidate","member","score","input","vrs"))
assert checks[1,1]==1 & checks[1,2]==0
assert member==1 & abs(score-cond(_n==1,.5,1))<1e-7
replace member=0 in 2
replace candidate=0 in 2
mata: st_matrix("checks",_sdea_static_evaluate("x","y","x","y","candidate","member","score","input","vrs"))
assert checks[1,1]==0 & checks[1,2]==1 & score[2]==1
matrix drop checks
display "PASS refset regression tests"
