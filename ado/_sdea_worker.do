*! Batch entry point for sdea workers(). The parent owns the job directory.
version 16
args worker workers
program define _sdea_worker_main
    args worker workers
    set processors 1
    mata: fh = fopen("../sample.bin","r")
    mata: paths = fgetmatrix(fh)
    mata: fclose(fh)
    mata: st_local("source",paths[1]); st_local("core",paths[2])
    adopath ++ "`macval(core)'"
    quietly run "`macval(source)'"
    mata: _sdea_worker_run(`worker',`workers')
end
capture noisily _sdea_worker_main `worker' `workers'
local rc = _rc
tempname status
file open `status' using status, write text
file write `status' "`rc'" _n
file close `status'
exit `rc'
