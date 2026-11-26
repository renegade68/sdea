version 16
args root core
if "`root'"=="" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"
clear
set obs 71
set seed 8197
generate long id = _N-_n+1
foreach var in x1 x2 y1 y2 {
    generate double `var' = 1+runiform()
}
generate byte sample = mod(id,3)!=0
local spec inputs(x1 x2) method(standard) refset(buildhull) dmu(id)
foreach rts in crs vrs {
    foreach ort in input output {
        quietly sdea y1 y2 if sample, `spec' rts(`rts') orientation(`ort') generate(s m d)
        assert r(workers)==1
        local refs = r(N_frame)
        local retries = r(N_full_retry)
        quietly sdea y1 y2 if sample, `spec' rts(`rts') orientation(`ort') workers(2) generate(s2 m2 d2)
        assert r(workers)==2 & r(N_frame)==`refs' & r(N_full_retry)==`retries'
        assert id==_N-_n+1 & s==s2 & m==m2 & d==d2
        assert !missing(s,m,d) if sample
        assert missing(s2) & missing(m2) & missing(d2) if !sample
        drop s m d s2 m2 d2
    }
}
quietly sdea y1 y2 in 1, `spec' workers(2) generate(s m d)
assert r(workers)==1 & r(N)==1 & s[1]==1
capture sdea y1 y2, `spec' workers(0) generate(s m d) replace
assert _rc==198 & s[1]==1
capture sdea y1 y2, inputs(x1 x2) method(standard) refset(full) workers(2)
assert _rc==198
