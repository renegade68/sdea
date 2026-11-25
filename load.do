version 16
args root core
if "`root'" == "" local root "`c(pwd)'"
mata: st_local("root", pathresolve(pwd(), st_local("root")))
local revision "ad77287466b22c3edd1ae3cc074981225c95d424"
if "`core'" == "" {
    local core "`root'/DEA-Stata-`revision'/ado"
}
mata: st_local("core", pathresolve(pwd(), st_local("core")))
confirm file "`root'/ado/sdea.ado"
foreach f in dea.ado _dea_simplex.ado _dea_ipm.ado _dea_lpcheck.ado _dea_nonte_lp.ado _dea_teaudit.ado _dea_recovery.ado {
    capture confirm file "`core'/`f'"
    if _rc {
        display as error "DEA-Stata 2.0.73 requires `f' in `core'."
        display as error "See Installation in README.md to extract the complete pinned archive."
        exit 601
    }
}
tempname header
file open `header' using "`core'/dea.ado", read text
file read `header' firstline
file close `header'
if !regexm("`firstline'", "^[*]! version 2[.]0[.]73( |$)") {
    display as error "This sdea release requires DEA-Stata 2.0.73. See Installation in README.md."
    exit 459
}
adopath ++ "`core'"
adopath ++ "`root'/ado"
which dea
which sdea
