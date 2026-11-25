version 16
args root core output
if "`root'" == "" local root "`c(pwd)'"
do "`root'/load.do" "`root'" "`core'"
clear
input long firm double(x1 x2 y)
1 1 3 5
2 1 3 5
3 2 4 4
4 3 5 2
end
sdea y, inputs(x1 x2) dmu(firm) method(standard) refset(buildhull) ///
    rts(vrs) orientation(input) generate(te isref drift)
list firm te isref drift, noobs
if "`output'" != "" export delimited firm te isref drift using "`output'/cross_section.csv"
