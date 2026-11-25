version 16
args root core repeats
if "`root'" == "" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"
if "`repeats'" == "" local repeats 1
import delimited using "`root'/tests/fixtures/full_dea_9040.csv", clear asdouble
isid firm
tempfile reference
save "`reference'"
do "`root'/examples/generate_9040.do"
merge 1:1 firm using "`reference'", assert(match) nogen
assert _N==9040 & firm==_n
forvalues run=1/`repeats' {
    foreach rts in crs vrs {
        foreach ort in input output {
            local side = cond("`ort'"=="input","in","out")
            timer clear 1
            timer on 1
            quietly sdea y1 y2, inputs(x1 x2 x3 x4 x5) method(standard) ///
                refset(buildhull) dmu(firm) rts(`rts') orientation(`ort') generate(te isref drift)
            timer off 1
            local candidates = r(N_candidates)
            local refs = r(N_frame)
            local kept = r(N_uncertain)
            local retries = r(N_full_retry)
            assert r(N)==_N & r(N_efficient)>0
            quietly timer list 1
            local seconds = r(t1)
            assert te>0 & te<=1 & inlist(isref,0,1) & drift==0
            quietly count if isref
            assert r(N)==`refs'
            generate double difference = abs(te - `rts'_`side')
            assert difference<=1e-7
            quietly summarize difference, meanonly
            di "RESULT `run' `rts' `ort' `seconds' `candidates' `refs' `kept' `retries' " %12.5g r(max)
            drop te isref drift difference
        }
    }
}
