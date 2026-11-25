version 16
clear
set rng mt64
set seed 119
set obs 9040
generate long firm = _n
foreach v in x1 x2 x3 x4 x5 y1 y2 {
    generate double `v' = 1 + 10*runiform()
}
