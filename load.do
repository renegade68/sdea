version 16
args root core
if "`root'" == "" local root "`c(pwd)'"
if "`core'" == "" {
    display as error "Supply the sdea repository root and the pinned DEA-Stata ado directory."
    exit 198
}
foreach f in dea.ado _dea_simplex.ado _dea_ipm.ado {
    confirm file "`core'/`f'"
}
confirm file "`root'/ado/sdea.ado"
adopath ++ "`core'"
adopath ++ "`root'/ado"
which dea
which sdea
