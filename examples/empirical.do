version 16
args root core data
if "`root'" == "" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"
foreach example in pharmacy banks {
    if "`example'"=="pharmacy" {
        use "`data'/pharmacy_yearly_1980_1989.dta" if year==1989, clear
        local X x1 x2 x3 x4
        local Y y1 y2 y3 y4
    }
    else {
        use "`data'/us_banks_quarterly_2005_2012.dta" if year==2005 & quarter==1, clear
        keep if !missing(deposits,employees,assets,loans) & min(deposits,employees,assets,loans)>0
        sort id
        keep in 1/200
        local X deposits employees
        local Y assets loans
    }
    generate long row_id = _n
    foreach rts in crs vrs {
        foreach orient in input output {
            quietly dea `X' = `Y', dmu(row_id) stage(1) engine(simplex) ///
                rts(`rts') orientation(`orient') generate(baseline)
            quietly sdea `Y', inputs(`X') dmu(row_id) method(standard) refset(buildhull) ///
                rts(`rts') orientation(`orient') generate(te member drift)
            local nf = r(N_frame)
            gen double error = abs(te-baseline)
            assert error<=1e-7
            quietly summarize error, meanonly
            display "EXAMPLE `example' `rts' `orient' " _N " references `nf' max " r(max)
            drop baseline te member drift error
            if "`rts'"=="vrs" drop CRS_TE VRS_TE SCALE RTS_NUM RTS
        }
    }
}
