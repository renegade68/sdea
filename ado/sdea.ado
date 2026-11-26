*! version 5.9.0  26nov2026
*! Static radial TE with refset(buildhull), using DEA-Stata 2.0.73.
*! sdea - Streaming DEA Package
*!
*! FIX: Orientation bug in super-efficiency efficient DMU identification
*!      (output-oriented supereff was using input criterion)
*! FIX: Super-SBM replaced with Tone (2002, EJOR) original formulation
*!      using Charnes-Cooper linearization of ratio form
*! FIX: Infeasible super-efficiency cases now return missing (not 999/0.001)
*! FIX: tf() option now accepts real values (for decimal time like 2019.75)
*! ADD: TDW (Time-Decay Weighted) option for TFDEA aggregate ROC
*! FIX: TFDEA ROC calculation using t_eff (Inman et al. 2006 Eq. 10)
*! FIX: TFDEA t_eff calculation using lambda weights (Inman et al. 2006 Eq.9)
*! FIX: TFDEA tf_exp formula using t_eff + ln(1/φ)/ln(γ) (Inman et al. 2006 Eq.13)
*! ADD: _sdea_solve_dea_with_lambda function for ROC t_eff calculation
*! ADD: Frontier update gate for BuildHull with optional ML/rule score
*! ADD: summarysave() period-level monitoring output for Stata workflows
*!
*! FEATURES:
*!   Streaming Methods: buildhull, fullstream, standard, time, window, malmquist, tfdea
*!   TFDEA: Technology forecasting with Rate of Change (ROC)
*!   TDW-TFDEA: Time-decay weighted ROC aggregation (decay option)
*!   DEA Models: te, sbm, supersbm, additive, cost, revenue, allocative, ddf, supereff, multiplier
*!   Malmquist: Adjacent period and Global MPI with TFPCH, EFFCH, TECHCH, PECH, SECH (CRS/VRS selectable)
*!   Drift Diagnostics: frontier relationship, multi-indicator, and univariate detectors
*!   Visualization: Evolution, drift, and frame graphs
*!   Time-based: Sequential evaluation with frontier tracking
*!   Window DEA: Moving window analysis for panel data
*!   NEW: Dual Fallback - Automatic primal-to-dual switching for robustness
*!   NEW: Frontier update gate - event-triggered updates for streaming frontiers
*!
*! TFDEA algorithm based on: Lim et al. (2014), Inman et al. (2006)
*! TDW-TFDEA: Lee (2026) - Time-decay weighted technology forecasting
*!
*! DEPENDENCY NOTE:
*!   sdea is organized as an extension of the current dea core. The
*!   development build uses the full ado directory from DEA-Stata 2.0.73.
*!   Streaming-specific state remains here: BuildHull frame maintenance,
*!   update gates, drift diagnostics, and time-indexed frontier management.
*!   Long-term maintenance should factor shared DEA model/solver routines
*!   into the common dea layer where feasible.
*!
*! SYNTAX: sdea outputs [if] [in], inputs(varlist) [options]
*!
*! Author: Choonjoo Lee (email: bloom.rampike@gmail.com)
*! Copyright (c) 2026 Choonjoo Lee.
*! Licensed under the Academic Free License version 3.0.
*! SPDX-License-Identifier: AFL-3.0

program define sdea, rclass
    version 16.0
    syntax varlist(min=1 numeric) [if] [in], ///
        INputs(varlist numeric min=1) ///
        [Time(varname numeric)] ///
        [DMU(varname)] ///
        [STARTers(integer 0)] ///
        [METhod(string)] ///
        [MODel(string)] ///
        [ORIENTation(string)] ///
        [RTS(string)] ///
        [WINdow(integer 3)] ///
        [GLOBal] ///
        [TF(real 0)] ///
        [DECAY(real 0)] ///
        [DECAYType(string)] ///
        [GENerate(namelist max=3)] ///
        [DRIFTtest] ///
        [DRIFTAlgo(string)] ///
        [DRIFTMode(string)] ///
        [DRIFTError(varname numeric)] ///
        [BY(varname)] ///
        [DRIFTWINdow(integer 3)] ///
        [THRESHold(real 5)] ///
        [THRVol(real 0.1)] ///
        [THRExp(real 0.2)] ///
        [THRTurn(real 0.3)] ///
        [GRaph] ///
        [SAVing(string)] ///
        [SAVERef(string)] ///
        [SAVEHistory(string)] ///
        [SUMMARYSave(string)] ///
        [REPlace] ///
        [VERBose] ///
        [GATE(name)] ///
        [GATEEps(real 1e-8)] ///
        [MLSCore(varname numeric)] ///
        [MLTHReshold(real 0.5)] ///
        [SOLver(string)] ///
        [W(varlist numeric)] ///
        [P(varlist numeric)] ///
        [GX(varlist numeric)] ///
        [GY(varlist numeric)] ///
        [Bounded(varlist numeric)] ///
        [Upper(varlist numeric)] ///
        [SEModel(string)] ///
        [REFSet(string)] ///
        [BHVariant(string)] ///
        [PRUNE(integer 0)] ///
        [WORKers(integer 1)]
    
    marksample touse
    local outputs `varlist'
    
    // Initialize LP solver based on Stata version
    
    quietly count if `touse'
    if r(N) == 0 error 2000
    local n_total = r(N)
    
    // Decay type validation for TDW-TFDEA
    if "`decaytype'" == "" local decaytype "exponential"
    if !inlist("`decaytype'", "exponential", "hyperbolic", "linear", "gaussian", "logistic") {
        di as error "decaytype() must be one of: exponential hyperbolic linear gaussian logistic"
        exit 198
    }
    
    // Method validation
    if "`method'" == "" local method "buildhull"
    if !inlist("`method'", "buildhull", "fullstream", "cumulative", "standard", "time") & ///
       !inlist("`method'", "window", "moving", "malmquist", "tfdea") {
        di as error "method() must be one of: buildhull fullstream cumulative standard time window moving malmquist tfdea"
        exit 198
    }
    if "`bhvariant'" == "" local bhvariant "modified"
    if !inlist("`bhvariant'", "modified", "dula") {
        di as error "bhvariant() must be one of: modified dula"
        exit 198
    }
    local refset = lower("`refset'")
    if `workers'<1 | missing(`workers') {
        di as err "workers() must be a positive integer"
        exit 198
    }
    if `workers'>1 {
        if "`method'"!="standard" | "`refset'"!="buildhull" {
            di as err "workers() requires method(standard) refset(buildhull)"
            exit 198
        }
        if !inlist("`c(os)'","MacOSX","Unix") | !c(MP) {
            di as err "parallel scoring currently requires Stata/MP on macOS or Unix"
            exit 198
        }
        if `workers'>c(processors_mach) {
            di as err "workers() exceeds the number of processors on this machine"
            exit 198
        }
        local workers = min(`workers',`n_total')
    }
    if "`refset'"!="" & (!inlist("`refset'","full","buildhull") | "`method'"!="standard") {
        di as err "refset(full|buildhull) requires method(standard)"
        exit 198
    }
    if "`refset'"=="buildhull" & "`solver'"=="" local solver "simplex"
    if "`solver'" == "" local solver "ipm"
    local solver = lower("`solver'")
    if !inlist("`solver'", "ipm", "simplex", "auto") {
        di as error "solver() must be one of: ipm simplex auto"
        exit 198
    }
    if "`summarysave'" != "" & "`time'" == "" {
        di as error "summarysave() requires time()"
        exit 198
    }
    
    // Time variable requirement
    if inlist("`method'", "buildhull", "fullstream", "cumulative", "time", "window", "moving", "malmquist", "tfdea") & "`time'" == "" {
        di as error "time() required for method(`method')"
        exit 198
    }
    
    // DMU variable for malmquist and tfdea
    if inlist("`method'", "malmquist", "tfdea") {
        if "`dmu'" == "" {
            capture confirm variable dmu
            if _rc == 0 {
                local dmu "dmu"
                di as text "Note: Using 'dmu' as DMU identifier"
            }
            else {
                tempvar dmu_temp
                quietly gen `dmu_temp' = _n
                local dmu "`dmu_temp'"
                di as text "Note: No dmu() specified, using observation number"
            }
        }
    }
    
    // Window size validation
    if "`method'" == "window" {
        if `window' < 1 {
            di as error "window() must be at least 1"
            exit 198
        }
    }
    
    // TFDEA specific validation
    if "`method'" == "tfdea" {
        if `tf' == 0 {
            quietly summarize `time' if `touse'
            local tf = r(max)
            di as text "Note: tf() not specified, using maximum time = `tf'"
        }
    }
    
    // Orientation validation
    if "`orientation'" == "" local orientation "output"
    if !inlist("`orientation'", "input", "output", "in", "out") {
        di as error "orientation() must be 'input' or 'output'"
        exit 198
    }
    if "`orientation'" == "in" local orientation "input"
    if "`orientation'" == "out" local orientation "output"
    
    // RTS validation
    if "`rts'" == "" local rts "vrs"
    if !inlist("`rts'", "crs", "vrs", "CRS", "VRS") {
        di as error "rts() must be 'crs' or 'vrs'"
        exit 198
    }
    local rts = lower("`rts'")
    
    // Model validation
    if "`model'" == "" local model "te"
    local model = lower("`model'")
    if !inlist("`model'", "te", "sbm", "supersbm", "additive", "cost") & ///
       !inlist("`model'", "revenue", "allocative", "ddf", "supereff", "multiplier") & ///
       !inlist("`model'", "fdh", "imprecise") {
        di as error "model() must be one of: te sbm supersbm additive cost revenue allocative ddf supereff multiplier fdh imprecise"
        exit 198
    }
    
    if "`refset'"!="" & "`model'"!="te" {
        di as err "refset() requires model(te)"
        exit 198
    }
    if "`refset'"=="buildhull" {
        if "`time'`drifttest'`graph'`saveref'`savehistory'`saving'`by'"!="" {
            di as err "refset(buildhull) supports static, ungrouped radial scores"
            exit 198
        }
        if `prune'<0 | `starters'<0 | ("`bhvariant'"=="dula" & (`prune'>0 | `starters'>0)) {
            di as err "prune() and starters() require nonnegative counts and apply to bhvariant(modified)"
            exit 198
        }
    }

    // Model-specific option validation
    if "`model'" == "cost" | ("`model'" == "allocative" & "`w'" != "") {
        if "`w'" == "" {
            di as error "model(`model') requires w(varlist) for input prices"
            exit 198
        }
        local nw : word count `w'
        local ninput_check : word count `inputs'
        if `nw' != `ninput_check' {
            di as error "w() must have same number of variables as inputs()"
            exit 198
        }
    }
    if "`model'" == "revenue" | ("`model'" == "allocative" & "`p'" != "") {
        if "`p'" == "" & "`model'" == "revenue" {
            di as error "model(revenue) requires p(varlist) for output prices"
            exit 198
        }
        if "`p'" != "" {
            local np : word count `p'
            local noutput_check : word count `outputs'
            if `np' != `noutput_check' {
                di as error "p() must have same number of variables as outputs"
                exit 198
            }
        }
    }
    if "`model'" == "allocative" & "`w'" == "" & "`p'" == "" {
        di as error "model(allocative) requires w(varlist) or p(varlist)"
        exit 198
    }
    if "`model'" == "ddf" {
        if "`gx'" == "" & "`gy'" == "" {
            di as error "model(ddf) requires gx() and/or gy() direction vectors"
            exit 198
        }
    }
    if "`model'" == "imprecise" {
        if "`bounded'" == "" & "`upper'" == "" {
            di as error "model(imprecise) requires bounded() or upper() option"
            exit 198
        }
        if "`upper'" == "" local upper "`bounded'"
    }
    
    // Non-oriented models ignore orientation
    if inlist("`model'", "sbm", "supersbm", "additive") {
        local orientation_display "non-oriented"
    }
    else {
        local orientation_display "`orientation'"
    }
    
    // Drift algorithm
    if "`driftalgo'" == "" local driftalgo "multi"
    if !inlist("`driftalgo'", "multi", "cusum", "adwin", "ph", "ddm", "eddm", "frontier") {
        di as error "driftalgo() must be one of: multi cusum adwin ph ddm eddm frontier"
        exit 198
    }
    if "`driftmode'" == "" local driftmode "dea"
    local driftmode = lower("`driftmode'")
    if !inlist("`driftmode'", "dea", "original") {
        di as error "driftmode() must be one of: dea original"
        exit 198
    }
    if "`driftmode'" == "original" {
        if !inlist("`driftalgo'", "ddm", "eddm") {
            di as error "driftmode(original) is currently supported only with driftalgo(ddm) or driftalgo(eddm)"
            exit 198
        }
        if "`drifterror'" == "" {
            di as error "driftmode(original) requires drifterror(varname), a 0/1 error indicator"
            exit 198
        }
        quietly count if missing(`drifterror') & `touse'
        if r(N) > 0 {
            di as error "drifterror() has missing values in the estimation sample"
            exit 459
        }
        quietly count if !inlist(`drifterror', 0, 1) & `touse'
        if r(N) > 0 {
            di as error "drifterror() must be a 0/1 binary error indicator"
            exit 198
        }
        quietly summarize `drifterror' if `touse', meanonly
        if r(N) == 0 | r(min) < 0 | r(max) > 1 {
            di as error "drifterror() must be nonmissing and bounded in [0,1] for the estimation sample"
            exit 198
        }
    }
    if "`by'" != "" {
        capture confirm numeric variable `by'
        if _rc {
            di as error "by() must specify a numeric group variable"
            exit 198
        }
    }
    local _N_drifts .
    local _drift_algo "`driftalgo'"
    
    local ninput : word count `inputs'
    local noutput : word count `outputs'

    // Frontier update gate validation
    if ("`gate'" != "" | "`mlscore'" != "") & "`method'" != "" & "`method'" != "buildhull" {
        di as error "gate() and mlscore() are currently supported only with method(buildhull)"
        exit 198
    }
    if ("`gate'" != "" | "`mlscore'" != "") & "`bhvariant'" == "dula" {
        di as error "gate() and mlscore() are currently supported only with bhvariant(modified)"
        exit 198
    }
    if `gateeps' < 0 {
        di as error "gateeps() must be nonnegative"
        exit 198
    }
    if `mlthreshold' < 0 | `mlthreshold' > 1 {
        di as error "mlthreshold() must be between 0 and 1"
        exit 198
    }
    
    *===========================================================================
    * Data Validation
    *===========================================================================
    
    // Check for missing values
    foreach v of varlist `inputs' `outputs' {
        quietly count if missing(`v') & `touse'
        if r(N) > 0 {
            di as error "Variable `v' has `r(N)' missing values"
            exit 459
        }
    }
    
    // Check outputs: each observation must have at least one positive output
    tempvar has_output has_input
    egen `has_output' = rowmax(`outputs') if `touse'
    quietly count if `has_output' <= 0 & `touse'
    if r(N) > 0 {
        di as error "Error: `r(N)' observations have all outputs <= 0"
        exit 459
    }
    
    // Check inputs: each observation must have at least one positive input
    egen `has_input' = rowmax(`inputs') if `touse'
    quietly count if `has_input' <= 0 & `touse'
    if r(N) > 0 {
        di as error "Error: `r(N)' observations have all inputs <= 0"
        exit 459
    }
    drop `has_output' `has_input'
    
    *===========================================================================
    * Header Display
    *===========================================================================
    
    di as text _n "{hline 70}"
    di as text "Streaming DEA Analysis (sdea v5.9.0)"
    di as text "{hline 70}"
    di as text "Method: " as result "`method'"
    di as text "Model: " as result "`model'"
    di as text "Observations: " as result `n_total'
    di as text "Inputs: " as result `ninput' " (" "`inputs'" ")"
    di as text "Outputs: " as result `noutput' " (" "`outputs'" ")"
    if "`time'" != "" di as text "Time variable: " as result "`time'"
    if inlist("`method'", "standard", "time") | "`model'" != "te" {
        di as text "Orientation: " as result "`orientation_display'"
        di as text "Returns to scale: " as result upper("`rts'")
    }
    if "`method'" == "buildhull" {
        di as text "LP solver: " as result "`solver'"
    }
    
    *===========================================================================
    * Variable Names Setup
    *===========================================================================
    
    if "`generate'" == "" {
        local gen_alpha "_sdea_alpha"
        local gen_frame "_sdea_frame"
        local gen_drift "_sdea_drift"
    }
    else {
        tokenize `generate'
        local gen_alpha "`1'"
        if "`2'" != "" local gen_frame "`2'"
        else local gen_frame "_sdea_frame"
        if "`3'" != "" local gen_drift "`3'"
        else local gen_drift "_sdea_drift"
    }
    
    if "`refset'"=="buildhull" {
        _sdea_refset `outputs', inputs(`inputs') touse(`touse') ///
            generate(`gen_alpha' `gen_frame' `gen_drift') dmu(`dmu') ///
            orientation(`orientation') rts(`rts') bhvariant(`bhvariant') ///
            solver(`solver') prune(`prune') starters(`starters') workers(`workers') `replace'
        return add
        return local method "standard"
        return local model "te"
        return local refset "buildhull"
        return scalar workers = `workers'
        return local orientation "`orientation'"
        return local rts "`rts'"
        return local bhvariant "`bhvariant'"
        return local hull_solver "`solver'"
        return local score_solver "simplex"
        exit
    }

    // Check if variables exist
    foreach v in `gen_alpha' `gen_frame' `gen_drift' {
        capture confirm variable `v'
        if _rc == 0 {
            if "`replace'" != "" {
                drop `v'
            }
            else {
                di as error "Variable `v' already exists. Use replace option."
                exit 110
            }
        }
    }
    if "`gate'" != "" {
        capture confirm variable `gate'
        if _rc == 0 {
            if "`replace'" != "" {
                drop `gate'
            }
            else {
                di as error "Variable `gate' already exists. Use replace option."
                exit 110
            }
        }
    }
    local gate_internal "`gate'"
    if "`gate_internal'" == "" & "`mlscore'" != "" {
        tempvar gate_internal
    }
    
    *===========================================================================
    * Method Dispatch
    *===========================================================================
    
    // Pre-create output variables for methods that need them
    if inlist("`method'", "fullstream", "cumulative", "moving") {
        quietly gen double `gen_alpha' = . if `touse'
        quietly gen byte `gen_frame' = 0 if `touse'
        quietly gen double `gen_drift' = 0 if `touse'
    }
    
    if inlist("`method'", "buildhull", "fullstream") & "`method'" == "buildhull" {
        // BuildHull method
        if `starters' == 0 {
            local min1 = 30
            local min2 = 3 * (`ninput' + `noutput')
            local min3 = floor(sqrt(`n_total'))
            local starters = max(`min1', `min2', `min3')
            if `starters' > 0.8 * `n_total' local starters = floor(0.8 * `n_total')
        }
        di as text "Starters: " as result `starters' " (" %4.1f 100*`starters'/`n_total' "%)"
        if "`gate'" != "" | "`mlscore'" != "" {
            if "`gate'" != "" di as text "Frontier update gate: " as result "`gate'"
            else di as text "Frontier update gate: " as result "enabled"
            di as text "Gate tolerance: " as result %9.2e `gateeps'
            if "`mlscore'" != "" {
                di as text "ML/rule score: " as result "`mlscore'" ///
                    as text " (threshold " as result %5.3f `mlthreshold' as text ")"
            }
        }
        
        if `starters' < 10 {
            di as error "Too few starters (minimum 10)"
            exit 198
        }
        if `starters' >= `n_total' {
            di as error "Starters must be less than total observations"
            exit 198
        }
        
        di as text "{hline 70}"
        local buildhull_sortvars "`dmu' `inputs' `outputs'"
        local gate_opts ""
        if "`gate_internal'" != "" {
            local gate_opts "`gate_opts' gate(`gate_internal') gateeps(`gateeps')"
        }
        if "`mlscore'" != "" {
            local gate_opts "`gate_opts' mlscore(`mlscore') mlthreshold(`mlthreshold')"
        }
        _sdea_buildhull `inputs', outputs(`outputs') time(`time') ///
            starters(`starters') alpha(`gen_alpha') frame(`gen_frame') ///
            drift(`gen_drift') touse(`touse') prune(`prune') ///
            sortvars(`buildhull_sortvars') bhvariant(`bhvariant') solver(`solver') ///
            `gate_opts'
        
        local _eff_type "buildhull"
        
        // Drift test
        if "`drifttest'" != "" {
            di as text _n "{hline 70}"
            di as text "Drift Detection"
            di as text "{hline 70}"
            _sdea_drift `gen_alpha', time(`time') drift(`gen_drift') ///
                frame(`gen_frame') touse(`touse') inputs(`inputs') outputs(`outputs') ///
                by(`by') driftmode(`driftmode') drifterror(`drifterror') ///
                window(`driftwindow') algorithm(`driftalgo') ///
                thrvol(`thrvol') threxp(`threxp') thrturn(`thrturn') ///
                threshold(`threshold') `replace'
            local _N_drifts = r(N_drifts)
            local _drift_algo "`r(algo_name)'"
        }
        
        // Graph
        if "`graph'" != "" {
            di as text _n "{hline 70}"
            di as text "Creating visualizations..."
            di as text "{hline 70}"
            _sdea_graph `gen_alpha', time(`time') drift(`gen_drift') ///
                frame(`gen_frame') touse(`touse') saving(`saving')
        }
        
        // Return values
        quietly count if `gen_frame' == 1 & `touse'
        return scalar N_frame = r(N)
        quietly count if abs(`gen_alpha') < 1e-6 & `touse'
        return scalar N_efficient = r(N)
        if "`drifttest'" != "" {
            quietly count if `gen_drift' > 0 & !missing(`gen_drift') & `touse'
            return scalar N_drifts = r(N)
        }
        else {
            return scalar N_drifts = 0
        }
        if "`gate_internal'" != "" {
            quietly count if `gate_internal' == 0 & `touse'
            return scalar N_gate_enveloped = r(N)
            quietly count if `gate_internal' == 1 & `touse'
            return scalar N_gate_updates = r(N)
            quietly count if `gate_internal' == 2 & `touse'
            return scalar N_gate_rejected = r(N)
            quietly count if `gate_internal' == 3 & `touse'
            return scalar N_gate_starters = r(N)
        }
    }
    else if inlist("`method'", "fullstream", "cumulative", "moving") {
        // Unified streaming methods
        di as text "{hline 70}"
        if "`method'" == "fullstream" {
            di as text "FullStream Streaming DEA"
        }
        else if "`method'" == "cumulative" {
            di as text "Cumulative Streaming DEA (Expanding Frontier)"
        }
        else if "`method'" == "moving" {
            di as text "Moving Window Streaming DEA"
            di as text "Window size: `window'"
        }
        di as text "{hline 70}"
        
        // DEBUG: Show actual values being passed
        
        // Sort by time before calling Mata
        tempvar sort_order
        quietly gen long `sort_order' = _n
        sort `time' `sort_order'
        
        // Call unified Mata function
        mata: _sdea_streaming_unified("`inputs'", "`outputs'", "`time'", "`method'", `window', "`gen_alpha'", "`gen_frame'", "", "`orientation'", "`rts'", "`verbose'")
        
        // Restore original order
        sort `sort_order'
        
        // Label variables
        label variable `gen_alpha' "Efficiency score (`method')"
        label variable `gen_frame' "Frontier member"
        
        local _eff_type "radial"
        
        // Drift test (common for streaming methods)
        if "`drifttest'" != "" {
            di as text _n "{hline 70}"
            di as text "Drift Detection"
            di as text "{hline 70}"
            _sdea_drift `gen_alpha', time(`time') drift(`gen_drift') ///
                frame(`gen_frame') touse(`touse') inputs(`inputs') outputs(`outputs') ///
                by(`by') driftmode(`driftmode') drifterror(`drifterror') ///
                window(`driftwindow') algorithm(`driftalgo') ///
                thrvol(`thrvol') threxp(`threxp') thrturn(`thrturn') ///
                threshold(`threshold') `replace'
            local _N_drifts = r(N_drifts)
            local _drift_algo "`r(algo_name)'"
        }
        
        // Return values (common)
        quietly count if `gen_frame' == 1 & `touse'
        return scalar N_frame = r(N)
        quietly count if abs(`gen_alpha' - 1) < 1e-6 & `touse'
        return scalar N_efficient = r(N)
        if "`drifttest'" != "" {
            quietly count if `gen_drift' > 0 & !missing(`gen_drift') & `touse'
            return scalar N_drifts = r(N)
        }
        else {
            return scalar N_drifts = 0
        }
    }
    else if "`method'" == "standard" {
        // Standard DEA (period by period)
        di as text "{hline 70}"
        
        // Use dea_lp for extended models, internal function for te
        if "`model'" == "te" {
            _sdea_standard `outputs', inputs(`inputs') time(`time') ///
                orientation(`orientation') rts(`rts') efficiency(`gen_alpha') touse(`touse')
        }
        else {
            // Use dea_lp with selected model
            _sdea_standard_model `outputs', inputs(`inputs') time(`time') ///
                model(`model') orientation(`orientation') rts(`rts') ///
                efficiency(`gen_alpha') touse(`touse') ///
                w(`w') p(`p') gx(`gx') gy(`gy') ///
                bounded(`bounded') upper(`upper')
        }
        
        quietly gen byte `gen_frame' = 0 if `touse'
        quietly gen double `gen_drift' = 0 if `touse'
        
        // Count efficient DMUs based on model
        if inlist("`model'", "additive", "ddf") {
            quietly count if abs(`gen_alpha') < 1e-6 & `touse'
        }
        else if "`model'" == "supereff" {
            // Orientation-dependent: input SE>=1, output SE<=1
            if "`orientation'" == "input" {
                quietly count if `gen_alpha' >= 1 - 1e-6 & `touse' & !missing(`gen_alpha')
            }
            else {
                quietly count if `gen_alpha' <= 1 + 1e-6 & `touse' & !missing(`gen_alpha')
            }
        }
        else if "`model'" == "supersbm" {
            // Super-SBM (Tone 2002): SE >= 1 always means super-efficient
            quietly count if `gen_alpha' >= 1 - 1e-6 & `touse' & !missing(`gen_alpha')
        }
        else {
            quietly count if abs(`gen_alpha' - 1) < 1e-6 & `touse'
        }
        return scalar N_efficient = r(N)
        return scalar N_frame = 0
        return scalar N_drifts = 0
    }
    else if "`method'" == "time" {
        // Time-based streaming DEA
        di as text "{hline 70}"
        
        // Use dea_lp for extended models
        if "`model'" == "te" {
            _sdea_time `outputs', inputs(`inputs') timevar(`time') ///
                orientation(`orientation') rts(`rts') ///
                saveref(`saveref') savehistory(`savehistory') `verbose'
        }
        else {
            _sdea_time_model `outputs', inputs(`inputs') timevar(`time') ///
                model(`model') orientation(`orientation') rts(`rts') ///
                saveref(`saveref') savehistory(`savehistory') `verbose' ///
                w(`w') p(`p') gx(`gx') gy(`gy')
        }
        
        // Rename standard output variables
        capture rename theta_t `gen_alpha'
        capture rename is_efficient_t `gen_frame'
        quietly gen double `gen_drift' = 0 if `touse'
        
        quietly count if `gen_frame' == 1 & `touse'
        return scalar N_efficient = r(N)
        return scalar N_frame = r(N)
        return scalar N_drifts = 0
    }
    else if "`method'" == "window" {
        // Window DEA
        di as text "{hline 70}"
        di as text "Window DEA Analysis"
        di as text "Window size: `window'"
        di as text "{hline 70}"
        
        // Get unique time periods
        tempvar time_id
        egen `time_id' = group(`time') if `touse'
        quietly summarize `time_id' if `touse'
        local T = r(max)
        local n_windows = `T' - `window' + 1
        
        if `n_windows' < 1 {
            di as error "Window size (`window') is larger than number of periods (`T')"
            exit 198
        }
        
        di as text "Number of time periods: `T'"
        di as text "Number of windows: `n_windows'"
        di as text "{hline 70}"
        
        // Drop existing _eff_w* variables if replace option specified
        if "`replace'" != "" {
            capture drop _eff_w*
        }
        
        // Save price options to avoid conflict with loop variable
        local w_prices "`w'"
        local p_prices "`p'"
        
        // Create efficiency variables for each window
        forvalues wnum = 1/`n_windows' {
            quietly gen double _eff_w`wnum' = . if `touse'
        }
        
        // Also create average efficiency variable
        quietly gen double `gen_alpha' = . if `touse'
        
        // Process each window
        forvalues wnum = 1/`n_windows' {
            local t_start = `wnum'
            local t_end = `wnum' + `window' - 1
            
            if "`verbose'" != "" {
                di as text "Processing window `wnum': periods `t_start' to `t_end'"
            }
            
            // Create window indicator
            tempvar in_window
            quietly gen byte `in_window' = (`time_id' >= `t_start' & `time_id' <= `t_end') & `touse'
            
            // Run DEA on pooled window data
            tempvar eff_temp
            quietly gen double `eff_temp' = . if `in_window'
            
            // Call appropriate DEA function based on model
            if "`model'" == "te" {
                mata: _sdea_run_dea("`inputs'", "`outputs'", "`eff_temp'", "`in_window'", "`orientation'", "`rts'")
            }
            else {
                // For other models, use internal batch solver
                mata: _sdea_solve_model_batch("`inputs'", "`outputs'", "`eff_temp'", ///
                    "`in_window'", "`model'", "`orientation'", "`rts'", "`w_prices'", "`p_prices'", "`gx'", "`gy'")
            }
            
            // Store window efficiency
            quietly replace _eff_w`wnum' = `eff_temp' if `in_window'
            
            drop `eff_temp' `in_window'
        }
        
        // Calculate average efficiency across windows for each observation
        // Each DMU-time may appear in multiple windows
        quietly {
            egen __n_windows__ = rownonmiss(_eff_w*) if `touse'
            egen __sum_eff__ = rowtotal(_eff_w*) if `touse', missing
            replace `gen_alpha' = __sum_eff__ / __n_windows__ if __n_windows__ > 0 & `touse'
            drop __n_windows__ __sum_eff__
        }
        
        // Create frame indicator (efficient = avg efficiency >= threshold)
        quietly gen byte `gen_frame' = (`gen_alpha' >= 1 - 1e-6) if `touse' & !missing(`gen_alpha')
        quietly gen double `gen_drift' = 0 if `touse'
        
        // Label variables
        label var `gen_alpha' "Window DEA average efficiency"
        forvalues wnum = 1/`n_windows' {
            local t_start = `wnum'
            local t_end = `wnum' + `window' - 1
            label var _eff_w`wnum' "Window `wnum' efficiency (t=`t_start'-`t_end')"
        }
        
        quietly count if `gen_frame' == 1 & `touse'
        return scalar N_efficient = r(N)
        return scalar N_frame = r(N)
        return scalar N_drifts = 0
        return scalar N_windows = `n_windows'
    }
    else if "`method'" == "malmquist" {
        // Malmquist Productivity Index
        di as text "{hline 70}"
        di as text "Malmquist Productivity Index Analysis"
        di as text "{hline 70}"
        
        // Drop existing Malmquist variables if replace option specified
        if "`replace'" != "" {
            capture drop _tfpch _effch _techch _pech _sech _eff_t _eff_t1 _eff_global
        }
        
        // Get unique time periods
        tempvar time_id
        egen `time_id' = group(`time') if `touse'
        quietly summarize `time_id' if `touse'
        local T = r(max)
        
        if `T' < 2 {
            di as error "At least 2 time periods required for Malmquist index"
            exit 198
        }
        
        di as text "Time periods: `T'"
        if "`global'" != "" {
            di as text "Method: Global Malmquist"
        }
        else {
            di as text "Method: Adjacent Period Malmquist"
        }
        di as text "Returns to scale: " as result upper("`rts'")
        di as text "{hline 70}"
        
        // Get DMU identifier
        if "`dmu'" == "" {
            capture confirm variable dmu
            if _rc == 0 {
                local dmu "dmu"
            }
            else {
                tempvar dmu_temp
                bysort `time' (`inputs' `outputs'): gen `dmu_temp' = _n if `touse'
                local dmu "`dmu_temp'"
            }
        }
        
        // Run Malmquist analysis
        if "`global'" != "" {
            mata: _sdea_malmquist_global("`inputs'", "`outputs'", "`time'", "`dmu'", ///
                "`touse'", "`orientation'", "`rts'", "`verbose'")
        }
        else {
            mata: _sdea_malmquist_adjacent("`inputs'", "`outputs'", "`time'", "`dmu'", ///
                "`touse'", "`orientation'", "`rts'", "`verbose'")
        }
        
        // Create standard sdea output variables from malmquist results
        quietly gen double `gen_alpha' = _tfpch if `touse'
        quietly gen byte `gen_frame' = (_tfpch >= 1) if `touse' & !missing(_tfpch)
        quietly gen double `gen_drift' = 0 if `touse'
        
        label var `gen_alpha' "Malmquist TFP Change"
        
        quietly count if `gen_frame' == 1 & `touse' & !missing(`gen_frame')
        return scalar N_efficient = r(N)
        return scalar N_frame = r(N)
        return scalar N_drifts = 0
        return scalar T = `T'
    }
    else if "`method'" == "tfdea" {
        // Technology Forecasting DEA
        di as text "{hline 70}"
        if `decay' > 0 {
            di as text "TDW-TFDEA: Time-Decay Weighted Technology Forecasting DEA"
        }
        else {
            di as text "Technology Forecasting using DEA (TFDEA)"
        }
        di as text "{hline 70}"
        
        // Set default semodel if not specified
        if "`semodel'" == "" local semodel "supereff"
        if !inlist("`semodel'", "supereff", "supersbm") {
            di as error "semodel() must be supereff or supersbm"
            exit 198
        }
        
        di as text "Forecast year (tf): `tf'"
        di as text "Super-efficiency model: `semodel'"
        if `decay' > 0 {
            di as text "TDW decay parameter (δ): `decay'"
            di as text "TDW decay function: `decaytype'"
        }
        di as text "{hline 70}"
        
        // Create TFDEA result variables
        foreach v in _tfdea_theta_tk _tfdea_theta_tf _tfdea_tk_eff _tfdea_roc ///
                     _tfdea_se _tfdea_tf_eff _tfdea_tfexp _tfdea_t_eff_roc ///
                     _tfdea_mean_roc_std _tfdea_mean_roc_tdw {
            capture drop `v'
            quietly gen double `v' = .
        }
        
        // Create main output variable first
        quietly gen double `gen_alpha' = . if `touse'
        
        // Call TFDEA Mata function with decay parameter for TDW
        mata: _sdea_tfdea("`inputs'", "`outputs'", "`time'", "`dmu'", ///
                         "`touse'", `tf', "`orientation'", "`rts'", ///
                         "`semodel'", "`verbose'" != "", `decay', "`decaytype'")
        
        // Set main efficiency variable from theta_tf
        quietly replace `gen_alpha' = _tfdea_theta_tf if `touse' & _tfdea_theta_tf != .
        quietly gen byte `gen_frame' = (abs(`gen_alpha' - 1) < 1e-6) if `touse' & !missing(`gen_alpha')
        quietly gen double `gen_drift' = 0 if `touse'
        
        // Return scalars
        quietly summarize _tfdea_roc if _tfdea_roc != . & _tfdea_roc > 0 & _tfdea_roc < 2
        if r(N) > 0 {
            return scalar mean_roc = r(mean)
            return scalar N_roc_dmu = r(N)
        }
        
        // Return TDW ROC if available
        quietly summarize _tfdea_mean_roc_std if _n == 1
        if r(N) > 0 {
            return scalar mean_roc_std = r(mean)
        }
        quietly summarize _tfdea_mean_roc_tdw if _n == 1
        if r(N) > 0 {
            return scalar mean_roc_tdw = r(mean)
        }
        
        return scalar tf = `tf'
        return scalar decay = `decay'
        return local decaytype "`decaytype'"
        return scalar N_efficient = 0
        return scalar N_frame = 0
        return scalar N_drifts = 0
    }
    
    *===========================================================================
    * Save Results
    *===========================================================================
    
    if "`saving'" != "" & "`graph'" == "" {
        if "`replace'" == "" {
            local timestamp = string(year(date(c(current_date), "DMY"))) + ///
                string(month(date(c(current_date), "DMY")), "%02.0f") + ///
                string(day(date(c(current_date), "DMY")), "%02.0f")
            local saving_file = "`saving'_`timestamp'.dta"
            save "`saving_file'"
            di as text _n "Results saved: " as result "`saving_file'"
        }
        else {
            save "`saving'", replace
            di as text _n "Results saved: " as result "`saving'"
        }
    }

    if "`summarysave'" != "" {
        local summary_gate_opt ""
        if "`gate_internal'" != "" {
            local summary_gate_opt "gate(`gate_internal')"
        }
        local summary_replace_opt ""
        if "`replace'" != "" local summary_replace_opt "replace"
        _sdea_save_monitoring_summary, time(`time') alpha(`gen_alpha') ///
            frame(`gen_frame') drift(`gen_drift') touse(`touse') ///
            saving("`summarysave'") `summary_gate_opt' `summary_replace_opt'
        return local summarysave "`summarysave'"
    }
    
    *===========================================================================
    * Summary
    *===========================================================================
    
    di as text _n "{hline 70}"
    di as text "Summary"
    di as text "{hline 70}"
    
    if inlist("`method'", "buildhull", "fullstream") {
        // BuildHull: alpha is a non-oriented hull distance/exteriority measure
        quietly summarize `gen_frame' if `touse'
        di as text "Frame members: " as result %4.0f r(sum) " (" %5.2f 100*r(mean) "%)"
        quietly count if `gen_frame' == 0 & `touse'
        di as text "Interior (dominated): " as result %4.0f r(N) " (" %5.2f 100*r(N)/`n_total' "%)"
        quietly summarize `gen_alpha' if `touse' & `gen_alpha' > 0, detail
        if r(N) > 0 {
            if "`method'" == "buildhull" {
                di as text "BuildHull alpha (mean, exterior only): " as result %8.4f r(mean)
                di as text "BuildHull alpha (max): " as result %8.4f r(max)
            }
            else {
                di as text "Streaming distance (mean, exterior only): " as result %8.4f r(mean)
                di as text "Streaming distance (max): " as result %8.4f r(max)
            }
        }
    }
    else if "`method'" == "time" {
        // Time-based: efficiency score 0-1
        quietly summarize `gen_frame' if `touse'
        di as text "Efficient at arrival: " as result %4.0f r(sum) " (" %5.2f 100*r(mean) "%)"
        quietly summarize `gen_alpha' if `touse', detail
        di as text "Efficiency score (mean): " as result %6.4f r(mean)
        di as text "Efficiency score (median): " as result %6.4f r(p50)
        di as text "Efficiency score (min): " as result %6.4f r(min)
    }
    else if "`method'" == "standard" {
        // Standard DEA: model-dependent interpretation
        quietly summarize `gen_alpha' if `touse', detail
        if inlist("`model'", "additive", "ddf") {
            di as text "Inefficiency score (mean): " as result %6.4f r(mean)
            di as text "Inefficiency score (max): " as result %6.4f r(max)
            quietly count if abs(`gen_alpha') < 1e-6 & `touse'
            di as text "Fully efficient DMUs (score=0): " as result %4.0f r(N) " (" %5.2f 100*r(N)/`n_total' "%)"
        }
        else if inlist("`model'", "supereff", "supersbm") {
            di as text "Super-efficiency score (mean): " as result %6.4f r(mean)
            di as text "Super-efficiency score (max): " as result %6.4f r(max)
            // Orientation-dependent super-efficiency identification
            if "`model'" == "supereff" & "`orientation'" == "output" {
                quietly count if `gen_alpha' <= 1 + 1e-6 & `touse' & !missing(`gen_alpha')
                di as text "Super-efficient DMUs (score<=1): " as result %4.0f r(N) " (" %5.2f 100*r(N)/`n_total' "%)"
            }
            else {
                quietly count if `gen_alpha' >= 1 - 1e-6 & `touse' & !missing(`gen_alpha')
                di as text "Super-efficient DMUs (score>=1): " as result %4.0f r(N) " (" %5.2f 100*r(N)/`n_total' "%)"
            }
        }
        else {
            di as text "Efficiency score (mean): " as result %6.4f r(mean)
            di as text "Efficiency score (median): " as result %6.4f r(p50)
            di as text "Efficiency score (min): " as result %6.4f r(min)
            quietly count if abs(`gen_alpha' - 1) < 1e-6 & `touse'
            di as text "Fully efficient DMUs: " as result %4.0f r(N) " (" %5.2f 100*r(N)/`n_total' "%)"
        }
    }
    else if "`method'" == "window" {
        // Window DEA: average efficiency across windows
        quietly summarize `gen_alpha' if `touse', detail
        di as text "Average efficiency (mean): " as result %6.4f r(mean)
        di as text "Average efficiency (median): " as result %6.4f r(p50)
        di as text "Average efficiency (min): " as result %6.4f r(min)
        quietly count if abs(`gen_alpha' - 1) < 1e-6 & `touse' & !missing(`gen_alpha')
        di as text "Fully efficient (avg=1): " as result %4.0f r(N) " (" %5.2f 100*r(N)/`n_total' "%)"
        
        // Window-specific stats
        di as text _n "Window-specific statistics:"
        local n_windows = `T' - `window' + 1
        forvalues w = 1/`n_windows' {
            quietly summarize _eff_w`w' if `touse' & !missing(_eff_w`w')
            local t_start = `w'
            local t_end = `w' + `window' - 1
            di as text "  Window `w' (t=`t_start'-`t_end'): mean=" as result %6.4f r(mean) as text ", N=" as result r(N)
        }
    }
    else if "`method'" == "malmquist" {
        // Malmquist Productivity Index
        di as text "Malmquist Productivity Index Results:"
        di as text ""
        
        quietly summarize _tfpch if `touse' & !missing(_tfpch)
        local n_malm = r(N)
        di as text "TFPCH (Total Factor Productivity Change)"
        di as text "  Mean: " as result %9.4f r(mean)
        quietly count if _tfpch > 1 & !missing(_tfpch) & `touse'
        di as text "  Improved (>1): " as result %4.0f r(N) " (" %5.2f 100*r(N)/`n_malm' "%)"
        
        quietly summarize _effch if `touse' & !missing(_effch)
        di as text "EFFCH (Efficiency Change - catch-up)"
        di as text "  Mean: " as result %9.4f r(mean)
        
        quietly summarize _techch if `touse' & !missing(_techch)
        di as text "TECHCH (Technical Change - frontier shift)"
        di as text "  Mean: " as result %9.4f r(mean)
        
        capture confirm variable _pech
        if _rc == 0 {
            quietly summarize _pech if `touse' & !missing(_pech)
            di as text "PECH (Pure Efficiency Change)"
            di as text "  Mean: " as result %9.4f r(mean)
            
            quietly summarize _sech if `touse' & !missing(_sech)
            di as text "SECH (Scale Efficiency Change)"
            di as text "  Mean: " as result %9.4f r(mean)
        }
        
        di as text ""
        di as text "Note: Values > 1 indicate improvement, < 1 indicate decline"
        if "`rts'" == "vrs" {
            di as text "Note: Under VRS, PECH = EFFCH and SECH = 1 (no scale decomposition)"
        }
    }
    else if "`method'" == "tfdea" {
        // TFDEA summary
        di as text "TFDEA Results:"
        di as text ""
        
        quietly summarize _tfdea_roc if _tfdea_roc != . & _tfdea_roc > 0 & _tfdea_roc < 2 & `touse'
        if r(N) > 0 {
            di as text "Rate of Change (ROC) Calculation:"
            di as text "  Annual ROC (mean): " as result %9.6f r(mean)
            di as text "  Valid DMUs for ROC: " as result %4.0f r(N)
            
            // Show TDW results if decay > 0
            if `decay' > 0 {
                quietly summarize _tfdea_mean_roc_std if _n == 1
                local roc_std = r(mean)
                quietly summarize _tfdea_mean_roc_tdw if _n == 1
                local roc_tdw = r(mean)
                di as text ""
                di as text "TDW-TFDEA (δ=`decay', `decaytype'):"
                di as text "  Aggregate ROC (Std): " as result %9.6f `roc_std' " (" %5.2f (`roc_std'-1)*100 "% annual)"
                di as text "  Aggregate ROC (TDW): " as result %9.6f `roc_tdw' " (" %5.2f (`roc_tdw'-1)*100 "% annual)"
            }
        }
        
        quietly count if _tfdea_se != . & `touse'
        if r(N) > 0 {
            di as text ""
            di as text "Technology Forecasting:"
            quietly summarize _tfdea_se if _tfdea_se != . & `touse'
            di as text "  Super-efficiency (mean): " as result %9.4f r(mean)
            di as text "  DMUs forecasted: " as result %4.0f r(N)
            
            quietly summarize _tfdea_tfexp if _tfdea_tfexp != . & `touse'
            if r(N) > 0 {
                di as text "  Forecasted year (mean): " as result %9.2f r(mean)
                di as text "  Forecasted year (min): " as result %9.2f r(min)
                di as text "  Forecasted year (max): " as result %9.2f r(max)
            }
        }
        
        di as text ""
        di as text "Variables created:"
        di as text "  _tfdea_theta_tk: Initial efficiency (when DMU first appeared)"
        di as text "  _tfdea_theta_tf: Final efficiency (at forecast year tf)"
        di as text "  _tfdea_roc: Rate of change"
        di as text "  _tfdea_t_eff_roc: t_eff for ROC (λ-weighted time)"
        di as text "  _tfdea_se: Super-efficiency (for forecasted DMUs)"
        di as text "  _tfdea_tf_eff: t_eff for forecasting"
        di as text "  _tfdea_tfexp: Forecasted release year"
        if `decay' > 0 {
            di as text "  _tfdea_mean_roc_std: Aggregate ROC (simple mean)"
            di as text "  _tfdea_mean_roc_tdw: Aggregate ROC (TDW weighted)"
        }
    }
    
    if "`drifttest'" != "" {
        quietly count if `gen_drift' > 0 & !missing(`gen_drift') & `touse'
        di as text "Drift points detected: " as result r(N)
    }
    di as text "{hline 70}"
    
    return scalar N = `n_total'
    if "`drifttest'" != "" {
        return scalar N_drifts = `_N_drifts'
        return local algo_name "`_drift_algo'"
    }
    if "`method'"=="standard" & "`model'"=="te" return local refset "full"
    return local method "`method'"
    return local bhvariant "`bhvariant'"
    return local model "`model'"
    return local rts "`rts'"
    return local orientation "`orientation_display'"
    
    di as text _n "Analysis complete!"
end

*===============================================================================
* INTERNAL PROGRAMS
*===============================================================================

*-------------------------------------------------------------------------------
* Static stage-1 scores and reference selection
*-------------------------------------------------------------------------------
program define _sdea_refset, rclass sortpreserve
    version 16
    syntax varlist(numeric), Inputs(varlist numeric) TOUSE(varname) ///
        GENerate(namelist min=3 max=3) [DMU(varname)] ORIentation(string) ///
        RTS(string) BHVariant(string) SOLver(string) PRUNE(integer) STARTers(integer) WORKers(integer) [REPlace]

    local raw : list inputs | varlist
    local protected `raw' `dmu'
    local unique : list uniq generate
    if "`unique'" != "`generate'" {
        di as err "generate() requires three distinct variable names"
        exit 198
    }
    local existing
    foreach name of local generate {
        if `: list name in protected' {
            di as err "output `name' also identifies a source variable"
            exit 198
        }
        capture confirm variable `name', exact
        if !_rc {
            if "`replace'"=="" {
                di as err "variable `name' already exists"
                exit 110
            }
            local existing `existing' `name'
        }
        else confirm new variable `name'
    }
    foreach v of local raw {
        quietly count if `touse' & `v'<0
        if r(N) {
            di as err "refset(buildhull) requires nonnegative inputs and outputs"
            exit 459
        }
    }

    tempname work
    tempvar obs link score member drift
    quietly generate long `obs' = _n
    frame put `raw' `obs' if `touse', into(`work')
    frame `work': quietly _sdea_refset_fit `varlist', inputs(`inputs') ///
        score(`score') member(`member') orientation(`orientation') rts(`rts') ///
        bhvariant(`bhvariant') solver(`solver') prune(`prune') starters(`starters') workers(`workers')
    local hull_ipm = r(hull_ipm)
    local hull_simplex = r(hull_simplex)
    local hull_failed = r(hull_failed)
    local candidates = r(candidates)
    local retained = r(retained)
    local retries = r(retries)

    quietly frlink 1:1 `obs', frame(`work') generate(`link')
    quietly generate double `score' = frval(`link',`score')
    quietly generate byte `member' = frval(`link',`member')
    quietly generate byte `drift' = 0 if `touse'
    quietly count if `touse'
    local n = r(N)
    quietly count if `member'==1
    local nf = r(N)
    quietly count if abs(`score'-1)<1e-6 & `touse'
    local ne = r(N)
    label variable `score' "Radial technical efficiency (0-1)"
    label variable `member' "Final static reference set"
    label variable `drift' "Static analysis (zero drift)"
    * Commit all outputs together. A failed rename restores replaced variables.
    preserve
    if "`existing'"!="" drop `existing'
    rename (`score' `member' `drift') (`generate')
    restore, not
    di as txt "Static radial TE: " as res `n' as txt " firms, " as res `nf' ///
        as txt " references, " as res `retries' as txt " full-reference retries"
    return scalar N = `n'
    return scalar N_hull_ipm = `hull_ipm'
    return scalar N_hull_simplex = `hull_simplex'
    return scalar N_hull_failed = `hull_failed'
    return scalar N_candidates = `candidates'
    return scalar N_frame = `nf'
    return scalar N_efficient = `ne'
    return scalar N_drifts = 0
    return scalar N_uncertain = `retained'
    return scalar N_full_retry = `retries'
end

program define _sdea_refset_fit, rclass
    version 16
    syntax varlist(numeric), Inputs(varlist numeric) SCORE(name) MEMBER(name) ///
        ORIentation(string) RTS(string) BHVariant(string) SOLver(string) ///
        PRUNE(integer) STARTers(integer) WORKers(integer)
    * Keep source columns out of the legacy builder's fixed output namespace.
    local xs
    local ys
    local nx
    local ny
    foreach side in x y {
        if "`side'"=="x" local vars `inputs'
        else local vars `varlist'
        foreach v of local vars {
            tempvar raw norm
            clonevar `raw' = `v'
            quietly summarize `raw', meanonly
            generate double `norm' = `raw' / cond(r(max)>0,r(max),1)
            local `side's ``side's' `raw'
            local n`side' `n`side'' `norm'
        }
    }
    tempvar candidates alpha streamdrift period
    generate byte `candidates' = 0
    mata: _sdea_static_pareto("`xs'", "`ys'", "`candidates'")
    quietly count if `candidates'
    local nc = r(N)
    tempname hullcounts
    matrix `hullcounts' = J(1,3,0)
    if `nc'==1 {
        generate byte `member' = `candidates'
    }
    else {
        if `starters'==0 local starters = max(30,3*wordcount("`nx' `ny'"),floor(sqrt(`nc')))
        local starters = min(`starters',`nc')
        generate byte `period' = 1
        _sdea_buildhull `nx', outputs(`ny') time(`period') touse(`candidates') ///
            starters(`starters') alpha(`alpha') frame(`member') drift(`streamdrift') ///
            bhvariant(`bhvariant') solver(`solver') prune(`prune') audit(`hullcounts')
        replace `member' = 0 if !`candidates'
    }
    generate double `score' = .
    tempname stats
    if `workers'>1 {
        tempfile stub
        local job "`stub'.sdea"
        quietly findfile sdea.ado
        mata: st_local("source",pathresolve(pwd(),st_global("r(fn)")))
        quietly findfile _dea_simplex.ado
        mata: st_local("core",pathgetparent(pathresolve(pwd(),st_global("r(fn)"))))
        local executable "`c(sysdir_stata)'StataMP.app/Contents/MacOS/stata-mp"
        capture confirm file "`macval(executable)'"
        if _rc local executable "`c(sysdir_stata)'stata-mp"
        confirm file "`macval(executable)'"
        mkdir "`macval(job)'"
    }
    capture noisily mata: st_matrix("`stats'", _sdea_static_evaluate("`xs'","`ys'","`nx'","`ny'", ///
        "`candidates'","`member'","`score'","`orientation'","`rts'",`workers'))
    local rc = _rc
    if `workers'>1 {
        mata: _sdea_worker_cleanup(st_local("job"),`workers')
    }
    if `rc' exit `rc'
    return scalar hull_ipm = `hullcounts'[1,1]
    return scalar hull_simplex = `hullcounts'[1,2]
    return scalar hull_failed = `hullcounts'[1,3]
    return scalar candidates = `nc'
    return scalar retained = `stats'[1,1]
    return scalar retries = `stats'[1,2]
end


* Period-level Stata monitoring output
program define _sdea_save_monitoring_summary
    version 16.0
    syntax, Time(varname numeric) Alpha(varname numeric) FRame(varname numeric) ///
        Drift(varname numeric) Touse(varname) SAVing(string) ///
        [GATE(varname numeric) REPlace]
    
    preserve
    quietly {
        keep if `touse'
        
        tempvar obs score driftval frame_member frame_entry frame_pruned frame_size drift_point
        gen long `obs' = 1
        gen double `score' = `alpha'
        gen double `driftval' = `drift'
        gen byte `frame_member' = (`frame' == 1) if !missing(`frame')
        gen byte `drift_point' = (`drift' > 0) if !missing(`drift')
        
        capture confirm variable _sdea_frame_entry
        if _rc == 0 gen double `frame_entry' = _sdea_frame_entry
        else gen double `frame_entry' = .
        
        capture confirm variable _sdea_frame_pruned
        if _rc == 0 gen double `frame_pruned' = _sdea_frame_pruned
        else gen double `frame_pruned' = .
        
        capture confirm variable _sdea_frame_size
        if _rc == 0 gen double `frame_size' = _sdea_frame_size
        else gen double `frame_size' = .
        
        if "`gate'" != "" {
            tempvar gate_enveloped gate_updates gate_rejected gate_starters
            gen byte `gate_enveloped' = (`gate' == 0) if !missing(`gate')
            gen byte `gate_updates' = (`gate' == 1) if !missing(`gate')
            gen byte `gate_rejected' = (`gate' == 2) if !missing(`gate')
            gen byte `gate_starters' = (`gate' == 3) if !missing(`gate')
            collapse (sum) N_obs=`obs' frame_members=`frame_member' ///
                    frame_entries=`frame_entry' frame_pruned=`frame_pruned' ///
                    drift_points=`drift_point' gate_enveloped=`gate_enveloped' ///
                    gate_updates=`gate_updates' gate_rejected=`gate_rejected' ///
                    gate_starters=`gate_starters' ///
                    (mean) mean_score=`score' mean_drift=`driftval' ///
                    (max) max_score=`score' max_drift=`driftval' ///
                    max_frame_size=`frame_size', by(`time')
        }
        else {
            collapse (sum) N_obs=`obs' frame_members=`frame_member' ///
                    frame_entries=`frame_entry' frame_pruned=`frame_pruned' ///
                    drift_points=`drift_point' ///
                    (mean) mean_score=`score' mean_drift=`driftval' ///
                    (max) max_score=`score' max_drift=`driftval' ///
                    max_frame_size=`frame_size', by(`time')
        }
        
        label variable `time' "Time period"
        label variable N_obs "Number of stream observations"
        label variable mean_score "Mean sdea score or hull distance"
        label variable max_score "Maximum sdea score or hull distance"
        label variable frame_members "Final frame members observed in period"
        label variable frame_entries "Accepted frame entries in period"
        label variable frame_pruned "Frame members pruned in period"
        label variable mean_drift "Mean drift magnitude"
        label variable max_drift "Maximum drift magnitude"
        label variable drift_points "Observations with positive drift"
        label variable max_frame_size "Maximum maintained frame size in period"
        capture label variable gate_enveloped "Gate: enveloped/no update"
        capture label variable gate_updates "Gate: accepted frontier updates"
        capture label variable gate_rejected "Gate: rejected candidate updates"
        capture label variable gate_starters "Gate: starter observations"
    }
    
    if "`replace'" != "" {
        save "`saving'", replace
    }
    else {
        save "`saving'"
    }
    restore
    di as text "Monitoring summary saved: " as result "`saving'"
end

*-------------------------------------------------------------------------------
* _sdea_fullstream: FullStream Algorithm (O(n) single-pass)
* Each observation is evaluated only once against the current frontier
*-------------------------------------------------------------------------------
*-------------------------------------------------------------------------------
* _sdea_buildhull: BuildHull streaming algorithm
*-------------------------------------------------------------------------------
program define _sdea_buildhull, sortpreserve
    version 16.0
    syntax varlist(numeric), OUtputs(varlist numeric) Time(varname numeric) ///
        STARTers(integer) Alpha(name) FRame(name) DRift(name) Touse(varname) ///
        [PRUNE(integer 0)] [SORTvars(varlist)] [BHVariant(string)] [SOLver(string)] [AUDIT(name)] ///
        [GATE(name) GATEEps(real 1e-8) MLSCore(varname numeric) MLTHReshold(real 0.5)]
    
    local _sdea_static_meter `audit'
    local inputs `varlist'
    local frame_size "_sdea_frame_size"
    local frame_entry "_sdea_frame_entry"
    local frame_pruned "_sdea_frame_pruned"
    
    tempvar sort_order obs_id stream_id is_starter
    quietly {
        gen long `sort_order' = _n if `touse'
        sort `time' `sortvars' `sort_order'
        gen long `obs_id' = _n if `touse'
        gen long `stream_id' = sum(`touse') if `touse'
        gen byte `is_starter' = (`stream_id' <= `starters') if `touse'
    }
    
    if "`bhvariant'" == "" local bhvariant "modified"
    if "`solver'" == "" local solver "ipm"
    if "`bhvariant'" == "dula" {
        mata: _sdea_buildhull_dula_main("`inputs'", "`outputs'", "`time'", ///
            "`obs_id'", "`alpha'", "`frame'", "`drift'", "`touse'", "`solver'")
    }
    else {
        mata: _sdea_buildhull_main("`inputs'", "`outputs'", "`time'", ///
            "`is_starter'", "`obs_id'", "`alpha'", "`frame'", "`drift'", "`touse'", ///
            `prune', "`solver'", "`gate'", `gateeps', "`mlscore'", `mlthreshold')
    }
    
    label variable `alpha' "BuildHull alpha (hull distance)"
    label variable `frame' "Frame member (1=yes, 0=no)"
    label variable `drift' "Drift magnitude (log scale)"
    label variable `frame_size' "BuildHull frame size after arrival"
    label variable `frame_entry' "BuildHull frame entry at arrival"
    label variable `frame_pruned' "BuildHull frame members pruned at arrival"
    if "`gate'" != "" {
        label variable `gate' "Frontier update gate"
        label define _sdea_gate 0 "enveloped/no update" 1 "accepted update" ///
            2 "candidate rejected" 3 "starter", replace
        label values `gate' _sdea_gate
    }
end

*-------------------------------------------------------------------------------
* _sdea_standard: Standard DEA (period by period)
*-------------------------------------------------------------------------------
program define _sdea_standard
    version 16.0
    syntax varlist(min=1), Inputs(varlist numeric) [Time(varname numeric)] ///
        ORIENTation(string) RTS(string) EFFiciency(name) Touse(varname)
    
    local outputs `varlist'
    quietly gen double `efficiency' = . if `touse'
    
    if "`time'" == "" {
        mata: _sdea_run_dea("`inputs'", "`outputs'", "`efficiency'", ///
            "`touse'", "`orientation'", "`rts'")
    }
    else {
        quietly levelsof `time' if `touse', local(time_levels)
        tempvar touse_t
        foreach t of local time_levels {
            quietly gen byte `touse_t' = (`touse' & `time'==`t')
            mata: _sdea_run_dea("`inputs'", "`outputs'", "`efficiency'", ///
                "`touse_t'", "`orientation'", "`rts'")
            drop `touse_t'
        }
    }
    
    // Convert output orientation to 0-1 scale
    if "`orientation'" == "output" {
        quietly replace `efficiency' = 1/`efficiency' if `touse' & `efficiency' > 0
    }
    
    label variable `efficiency' "DEA efficiency score"
    quietly summ `efficiency' if `touse', detail
    di as text "Efficiency statistics:"
    di as text "  Mean: " as result %6.4f r(mean)
    di as text "  Min: " as result %6.4f r(min)
    di as text "  Max: " as result %6.4f r(max)
end

*-------------------------------------------------------------------------------
* _sdea_standard_model: Standard DEA with extended models (internal solver)
*-------------------------------------------------------------------------------
program define _sdea_standard_model
    version 16.0
    syntax varlist(min=1), Inputs(varlist numeric) [Time(varname numeric)] ///
        MODel(string) ORIENTation(string) RTS(string) EFFiciency(name) Touse(varname) ///
        [W(varlist numeric)] [P(varlist numeric)] [GX(varlist numeric)] [GY(varlist numeric)] ///
        [Bounded(varlist numeric)] [Upper(varlist numeric)]
    
    local outputs `varlist'
    
    quietly gen double `efficiency' = . if `touse'
    
    if "`time'" == "" {
        // Single period analysis - use Mata solver directly
        mata: _sdea_solve_model_batch("`inputs'", "`outputs'", "`efficiency'", ///
            "`touse'", "`model'", "`orientation'", "`rts'", "`w'", "`p'", "`gx'", "`gy'")
    }
    else {
        // Period-by-period analysis
        quietly levelsof `time' if `touse', local(time_levels)
        foreach t of local time_levels {
            tempvar touse_t
            quietly gen byte `touse_t' = `touse' & `time' == `t'
            mata: _sdea_solve_model_batch("`inputs'", "`outputs'", "`efficiency'", ///
                "`touse_t'", "`model'", "`orientation'", "`rts'", "`w'", "`p'", "`gx'", "`gy'")
            drop `touse_t'
        }
    }
    
    // Label based on model
    if inlist("`model'", "additive", "ddf") {
        label variable `efficiency' "DEA inefficiency score (`model')"
    }
    else if inlist("`model'", "supereff", "supersbm") {
        label variable `efficiency' "DEA super-efficiency score (`model')"
    }
    else {
        label variable `efficiency' "DEA efficiency score (`model')"
    }
    
    quietly summ `efficiency' if `touse', detail
    di as text "Efficiency statistics (`model'):"
    di as text "  Mean: " as result %6.4f r(mean)
    di as text "  Min: " as result %6.4f r(min)
    di as text "  Max: " as result %6.4f r(max)
end

*-------------------------------------------------------------------------------
* _sdea_time_model: Time-based streaming DEA with extended models
*-------------------------------------------------------------------------------
program define _sdea_time_model, rclass sortpreserve
    version 16.0
    syntax varlist(min=1 numeric), ///
        Inputs(varlist numeric) ///
        TIMEvar(varname numeric) ///
        MODel(string) ///
        [ORIENTation(string)] ///
        [RTS(string)] ///
        [SAVEref(string)] ///
        [SAVEhistory(string)] ///
        [VERBose] ///
        [W(varlist numeric)] [P(varlist numeric)] [GX(varlist numeric)] [GY(varlist numeric)]
    
    local outputs "`varlist'"
    
    if "`rts'" == "" local rts "vrs"
    if "`orientation'" == "" local orientation "output"
    
    // Check for existing variables
    foreach v in theta_t is_efficient_t ref_size_t time_period {
        capture confirm variable `v'
        if _rc == 0 {
            drop `v'
        }
    }
    
    // Create result variables
    quietly gen double theta_t = .
    quietly gen byte is_efficient_t = .
    quietly gen int ref_size_t = .
    quietly gen int time_period = .
    
    // Get time periods
    quietly levelsof `timevar', local(time_levels)
    local n_periods : word count `time_levels'
    
    di as text _n "Time-based Streaming DEA (model: `model')"
    di as text "{hline 50}"
    di as text "Periods: " as result `n_periods'
    
    local period_num = 0
    local cumulative_n = 0
    
    foreach t of local time_levels {
        local ++period_num
        
        // Count observations up to and including this period
        quietly count if `timevar' <= `t'
        local cumulative_n = r(N)
        
        // Evaluate all observations up to this period using internal solver
        tempvar _eff_temp touse_t
        quietly gen byte `touse_t' = (`timevar' <= `t')
        quietly gen double `_eff_temp' = . if `touse_t'
        mata: _sdea_solve_model_batch("`inputs'", "`outputs'", "`_eff_temp'", ///
            "`touse_t'", "`model'", "`orientation'", "`rts'", "`w'", "`p'", "`gx'", "`gy'")
        
        // Store results for current period observations
        quietly replace theta_t = `_eff_temp' if `timevar' == `t'
        quietly replace time_period = `period_num' if `timevar' == `t'
        quietly replace ref_size_t = `cumulative_n' if `timevar' == `t'
        
        // Determine efficiency based on model
        if inlist("`model'", "additive", "ddf") {
            quietly replace is_efficient_t = (abs(`_eff_temp') < 1e-6) if `timevar' == `t'
        }
        else if inlist("`model'", "supereff", "supersbm") {
            // Super-efficiency: orientation-dependent identification
            if "`model'" == "supereff" & "`orientation'" == "output" {
                quietly replace is_efficient_t = (`_eff_temp' <= 1 + 1e-6) if `timevar' == `t' & !missing(`_eff_temp')
            }
            else {
                quietly replace is_efficient_t = (`_eff_temp' >= 1 - 1e-6) if `timevar' == `t' & !missing(`_eff_temp')
            }
        }
        else {
            quietly replace is_efficient_t = (abs(`_eff_temp' - 1) < 1e-6) if `timevar' == `t'
        }
        
        drop `_eff_temp' `touse_t'
        
        if "`verbose'" != "" {
            quietly count if is_efficient_t == 1 & `timevar' == `t'
            local n_eff = r(N)
            quietly count if `timevar' == `t'
            local n_t = r(N)
            di as text "  Period `period_num' (t=`t'): " as result "`n_t' obs, `n_eff' efficient, ref_size=`cumulative_n'"
        }
    }
    
    label variable theta_t "Efficiency score at arrival"
    label variable is_efficient_t "Efficient at time of arrival"
    label variable ref_size_t "Reference set size"
    label variable time_period "Time period number"
    
    di as text "{hline 50}"
    quietly summ theta_t, detail
    di as text "Final efficiency (mean): " as result %6.4f r(mean)
    quietly count if is_efficient_t == 1
    di as text "Total efficient DMUs: " as result r(N)
end

*-------------------------------------------------------------------------------
* _sdea_time: Time-based streaming DEA
*-------------------------------------------------------------------------------
program define _sdea_time, rclass sortpreserve
    version 16.0
    syntax varlist(min=1 numeric), ///
        Inputs(varlist numeric) ///
        TIMEvar(varname numeric) ///
        [ORIENTation(string)] ///
        [RTS(string)] ///
        [SAVEref(string)] ///
        [SAVEhistory(string)] ///
        [VERBose]
    
    local outputs "`varlist'"
    
    if "`rts'" == "" local rts "vrs"
    if "`orientation'" == "" local orientation "output"
    
    local rtsU = upper("`rts'")
    if "`orientation'" == "input" local ortU "IN"
    else local ortU "OUT"
    
    // Check for existing variables
    foreach v in theta_t is_efficient_t ref_size_t time_period {
        capture confirm variable `v'
        if _rc == 0 {
            drop `v'
        }
    }
    
    // Create result variables
    quietly gen double theta_t = .
    quietly gen byte is_efficient_t = 0
    quietly gen long ref_size_t = .
    quietly gen long time_period = .
    
    label var theta_t "Time-based efficiency score"
    label var is_efficient_t "Efficient at arrival time (1=yes)"
    label var ref_size_t "Reference set size at evaluation"
    label var time_period "Time period of evaluation"
    
    // Sort by time
    tempvar orig_order touse_temp
    gen long `orig_order' = _n
    gen byte `touse_temp' = 1
    sort `timevar'
    
    if "`verbose'" != "" {
        di _n as text "{hline 78}"
        di as text "Time-based Streaming DEA"
        di as text "{hline 78}"
        di as text "RTS:         `rtsU'"
        di as text "Orientation: `ortU'"
        di as text "{hline 78}" _n
    }
    
    // Run streaming DEA
    mata: _sdea_time_streaming("`inputs'", "`outputs'", "`touse_temp'", ///
        "`timevar'", "`rtsU'", "`ortU'", ///
        "theta_t", "is_efficient_t", "ref_size_t", "time_period", ///
        "`saveref'", "`savehistory'", "`verbose'")
    
    // Restore original order
    sort `orig_order'
    
    // Summary
    quietly count if is_efficient_t == 1
    local n_eff = r(N)
    quietly count
    local N = r(N)
    quietly sum theta_t, detail
    
    if "`verbose'" != "" {
        di _n as text "{hline 78}"
        di as text "Results Summary"
        di as text "{hline 78}"
        di as text "Total observations:   " as result `N'
        di as text "Efficient at arrival: " as result `n_eff' " (" %5.2f 100*`n_eff'/`N' "%)"
        di as text "Mean efficiency:      " as result %6.4f r(mean)
        di as text "{hline 78}" _n
    }
    
    return scalar N = `N'
    return scalar n_efficient = `n_eff'
    return scalar mean_eff = r(mean)
end

*-------------------------------------------------------------------------------
* _sdea_drift: Drift diagnostics
* 
* Implements production-relationship and DEA-signal diagnostics, including:
*   1. Efficiency Volatility (V_t): Large deviations in mean efficiency
*   2. Frontier Expansion Rate (D_t): Abnormal growth in efficient units
*   3. Reference Set Turnover (T_t): Instability in benchmark composition
*   4. Rolling input-output relationship checks and univariate detectors
*
* A drift signal is raised when at least 2 of 3 indicators exceed thresholds
*-------------------------------------------------------------------------------
program define _sdea_drift, rclass
    version 16.0
    syntax varname(numeric), Time(varname numeric) DRift(varname numeric) ///
        FRame(varname numeric) Touse(varname) ///
        [INputs(varlist numeric)] ///
        [OUtputs(varlist numeric)] ///
        [BY(varname)] ///
        [DRIFTMode(string)] ///
        [DRIFTError(varname numeric)] ///
        [WINdow(integer 3)] ///
        [THRvol(real 0.1)] ///
        [THRexp(real 0.2)] ///
        [THRturn(real 0.3)] ///
        [THRESHold(real 5)] ///
        [ALGOrithm(string)] ///
        [REPlace]
    
    local alpha `varlist'
    
    // Default algorithm
    if "`algorithm'" == "" local algorithm "multi"
    if "`driftmode'" == "" local driftmode "dea"
    local driftmode = lower("`driftmode'")
    if !inlist("`driftmode'", "dea", "original") {
        di as error "driftmode() must be one of: dea original"
        exit 198
    }
    if "`driftmode'" == "original" {
        if !inlist("`algorithm'", "ddm", "eddm") {
            di as error "driftmode(original) is supported only with algorithm(ddm) or algorithm(eddm)"
            exit 198
        }
        if "`drifterror'" == "" {
            di as error "driftmode(original) requires drifterror(varname), a 0/1 error indicator"
            exit 198
        }
        quietly count if missing(`drifterror') & `touse'
        if r(N) > 0 {
            di as error "drifterror() has missing values in the estimation sample"
            exit 459
        }
        quietly count if !inlist(`drifterror', 0, 1) & `touse'
        if r(N) > 0 {
            di as error "drifterror() must be a 0/1 binary error indicator"
            exit 198
        }
        quietly summarize `drifterror' if `touse', meanonly
        if r(N) == 0 | r(min) < 0 | r(max) > 1 {
            di as error "drifterror() must be nonmissing and bounded in [0,1] for the estimation sample"
            exit 198
        }
    }
    
    // Ensure time ordering
    tempvar __ord
    quietly gen long `__ord' = _n
    quietly sort `time' `__ord'
    
    // Get distinct time periods
    quietly levelsof `time' if `touse', local(time_levels)
    local T : word count `time_levels'
    
    if `T' < `window' + 2 & "`algorithm'" == "multi" {
        di as text "Warning: Insufficient time periods for multi-indicator drift detection"
        di as text "         (need at least window + 2 = `=`window'+2' periods, have `T')"
        di as text "         Falling back to CUSUM algorithm"
        local algorithm "cusum"
    }
    
    // Create period-level statistics
    capture confirm variable _sdea_frame_size
    local has_frame_history = (_rc == 0)
    if `has_frame_history' {
        capture confirm variable _sdea_frame_entry
        if _rc local has_frame_history = 0
    }
    if `has_frame_history' {
        capture confirm variable _sdea_frame_pruned
        if _rc local has_frame_history = 0
    }
    
    tempvar period_tag mean_eff sd_eff n_eff n_obs n_frame_t
    tempvar hist_frame_size hist_frame_entry hist_frame_pruned
    quietly {
        egen `period_tag' = tag(`time') if `touse'
        bysort `time': egen `mean_eff' = mean(`alpha') if `touse'
        bysort `time': egen `sd_eff' = sd(`alpha') if `touse'
        bysort `time': egen `n_eff' = total(`frame') if `touse'
        bysort `time': egen `n_obs' = count(`alpha') if `touse'
        bysort `time': egen `n_frame_t' = total(`frame') if `touse'
        if `has_frame_history' {
            bysort `time': egen `hist_frame_size' = max(_sdea_frame_size) if `touse'
            bysort `time': egen `hist_frame_entry' = total(_sdea_frame_entry) if `touse'
            bysort `time': egen `hist_frame_pruned' = total(_sdea_frame_pruned) if `touse'
        }
    }
    
    // =========================================================================
    // CUSUM Algorithm
    // =========================================================================
    
    if "`algorithm'" == "cusum" {
        di as text "Algorithm: CUSUM"
        di as text "{hline 50}"
        
        // Check for variation
        quietly summarize `alpha' if `touse', detail
        local mu = r(mean)
        local sd = r(sd)
        
        if `sd' <= 1e-10 | missing(`sd') {
            di as text "Warning: Alpha has near-zero variance. No drift detected."
            quietly sort `__ord'
            return scalar N_drifts = 0
            return scalar threshold = `thrvol'
            return local algo_name "cusum"
            exit 0
        }
        
        // CUSUM threshold
        local cusum_thr = `threshold'
        
        // Create CUSUM statistics
        tempvar cusum_z cusum_stat cusum_alarm
        quietly {
            gen double `cusum_z' = (`alpha' - `mu')/`sd' if `touse'
            gen double `cusum_stat' = sum(`cusum_z') if `touse'
            gen byte `cusum_alarm' = (abs(`cusum_stat') > `cusum_thr') if `touse'
            replace `drift' = `cusum_alarm' if `touse'
        }
        
        quietly sort `__ord'
        
        quietly count if `cusum_alarm' == 1 & `period_tag' == 1 & `touse'
        local N_drifts = r(N)
        
        di as text "CUSUM threshold: " as result `cusum_thr'
        di as text "Drift events detected: " as result `N_drifts'
        
        if `N_drifts' > 0 & `N_drifts' <= 20 {
            di as text _n "Drift time points:"
            quietly levelsof `time' if `drift' > 0 & !missing(`drift') & `period_tag' == 1 & `touse', local(drift_times) clean
            di as result "  `drift_times'"
        }
        
        return scalar N_drifts = `N_drifts'
        return scalar threshold = `cusum_thr'
        return local algo_name "cusum"
        exit 0
    }
    
    // =========================================================================
    // ADWIN-style Algorithm (Adaptive Windowing)
    // Bifet and Gavalda (2007), adapted to a DEA-derived univariate signal
    // =========================================================================
    
    if "`algorithm'" == "adwin" {
        di as text "Algorithm: ADWIN (Adaptive Windowing)"
        di as text "{hline 50}"
        
        // ADWIN parameters
        local delta = 0.002  // confidence parameter
        local min_window = max(10, `window')
        
        tempvar adwin_alarm adwin_stat
        quietly {
            gen byte `adwin_alarm' = 0 if `touse'
            gen double `adwin_stat' = . if `touse'
        }
        
        // Get time-sorted efficiency values
        tempvar time_order
        quietly gen long `time_order' = _n
        quietly sort `time' `time_order'
        
        quietly count if `touse'
        local n_obs = r(N)
        
        // Process each observation
        mata: _sdea_adwin("`alpha'", "`touse'", "`adwin_alarm'", "`adwin_stat'", `delta', `min_window')
        
        quietly replace `drift' = `adwin_alarm' if `touse'
        
        quietly sort `time_order'
        
        quietly count if `adwin_alarm' == 1 & `period_tag' == 1 & `touse'
        local N_drifts = r(N)
        
        di as text "ADWIN delta: " as result `delta'
        di as text "Min window: " as result `min_window'
        di as text "Drift events detected: " as result `N_drifts'
        
        if `N_drifts' > 0 & `N_drifts' <= 20 {
            di as text _n "Drift time points:"
            quietly levelsof `time' if `drift' > 0 & !missing(`drift') & `period_tag' == 1 & `touse', local(drift_times) clean
            di as result "  `drift_times'"
        }
        
        return scalar N_drifts = `N_drifts'
        return scalar delta = `delta'
        return local algo_name "adwin"
        exit 0
    }
    
    // =========================================================================
    // Page-Hinkley Test
    // Page (1954) - sequential change detection
    // =========================================================================
    
    if "`algorithm'" == "ph" {
        di as text "Algorithm: Page-Hinkley Test"
        di as text "{hline 50}"
        
        // Page-Hinkley parameters
        local ph_delta = 0.005  // magnitude of allowed change
        local ph_lambda = 50    // detection threshold
        
        quietly summarize `alpha' if `touse', detail
        local mu = r(mean)
        
        tempvar ph_m ph_M ph_alarm
        quietly {
            gen double `ph_m' = 0 if `touse'
            gen double `ph_M' = 0 if `touse'
            gen byte `ph_alarm' = 0 if `touse'
        }
        
        // Get time-sorted order
        tempvar time_order obs_num
        quietly gen long `time_order' = _n
        quietly sort `time' `time_order'
        quietly gen long `obs_num' = _n if `touse'
        
        // Compute Page-Hinkley statistics
        quietly {
            // m_t = Σ(x_i - μ - δ)
            gen double __ph_cum__ = sum(`alpha' - `mu' - `ph_delta') if `touse'
            replace `ph_m' = __ph_cum__ if `touse'
            
            // M_t = max(m_1, ..., m_t)
            gen double __ph_max__ = `ph_m'[1] if `touse' & `obs_num' == 1
            replace __ph_max__ = max(__ph_max__[_n-1], `ph_m') if `touse' & `obs_num' > 1
            replace `ph_M' = __ph_max__ if `touse'
            
            // PH_t = M_t - m_t, alarm if PH_t > λ
            gen double __ph_stat__ = `ph_M' - `ph_m' if `touse'
            replace `ph_alarm' = (__ph_stat__ > `ph_lambda') if `touse'
            
            drop __ph_cum__ __ph_max__ __ph_stat__
        }
        
        quietly replace `drift' = `ph_alarm' if `touse'
        
        quietly sort `time_order'
        
        quietly count if `ph_alarm' == 1 & `period_tag' == 1 & `touse'
        local N_drifts = r(N)
        
        di as text "PH delta (tolerance): " as result `ph_delta'
        di as text "PH lambda (threshold): " as result `ph_lambda'
        di as text "Drift events detected: " as result `N_drifts'
        
        if `N_drifts' > 0 & `N_drifts' <= 20 {
            di as text _n "Drift time points:"
            quietly levelsof `time' if `drift' > 0 & !missing(`drift') & `period_tag' == 1 & `touse', local(drift_times) clean
            di as result "  `drift_times'"
        }
        
        return scalar N_drifts = `N_drifts'
        return scalar ph_delta = `ph_delta'
        return scalar ph_lambda = `ph_lambda'
        return local algo_name "ph"
        exit 0
    }
    
    // =========================================================================
    // DDM Drift Detection Method
    // Gama et al. (2004). Default mode adapts the detector from classification
    // errors to DEA signal deviations; original mode uses a supplied 0/1
    // error stream.
    // =========================================================================
    
    if "`algorithm'" == "ddm" {
        di as text "Algorithm: DDM (Drift Detection Method)"
        di as text "{hline 50}"
        
        // DDM thresholds
        local ddm_warning = 2.0   // warning level (std devs)
        local ddm_drift = 3.0     // drift level (std devs)
        local min_samples = 30
        
        tempvar ddm_error ddm_p ddm_s ddm_pmin ddm_smin ddm_alarm ddm_warning_flag
        
        quietly {
            if "`driftmode'" == "original" {
                gen double `ddm_error' = `drifterror' if `touse'
            }
            else {
                gen double `ddm_error' = abs(1 - `alpha') if `touse'
            }
            
            // Running mean and std dev of error
            gen double `ddm_p' = . if `touse'
            gen double `ddm_s' = . if `touse'
            gen double `ddm_pmin' = . if `touse'
            gen double `ddm_smin' = . if `touse'
            gen byte `ddm_alarm' = 0 if `touse'
            gen byte `ddm_warning_flag' = 0 if `touse'
        }
        
        tempvar time_order obs_num
        quietly gen long `time_order' = _n
        quietly sort `time' `time_order'
        quietly gen long `obs_num' = _n if `touse'
        
        // Compute DDM statistics using Mata for efficiency
        mata: _sdea_ddm("`ddm_error'", "`touse'", "`ddm_alarm'", "`ddm_warning_flag'", ///
                        `ddm_warning', `ddm_drift', `min_samples')
        
        quietly replace `drift' = `ddm_alarm' if `touse'
        
        quietly sort `time_order'
        
        tempvar ddm_period_alarm
        bysort `time': egen `ddm_period_alarm' = max(`ddm_alarm') if `touse'
        quietly count if `ddm_period_alarm' == 1 & `period_tag' == 1 & `touse'
        local N_drifts = r(N)
        
        quietly count if `ddm_warning_flag' == 1 & `period_tag' == 1 & `touse'
        local N_warnings = r(N)
        
        di as text "DDM warning level: " as result `ddm_warning' " std devs"
        di as text "DDM drift level: " as result `ddm_drift' " std devs"
        if "`driftmode'" == "original" {
            di as text "Input stream: " as result "`drifterror'" as text " (0/1 original DDM error indicator)"
        }
        else {
            di as text "Input stream: DEA signal deviation |1 - alpha|"
        }
        di as text "Min samples: " as result `min_samples'
        di as text "Warning events: " as result `N_warnings'
        di as text "Drift events detected: " as result `N_drifts'
        
        if `N_drifts' > 0 & `N_drifts' <= 20 {
            di as text _n "Drift time points:"
            quietly levelsof `time' if `ddm_period_alarm' == 1 & `period_tag' == 1 & `touse', local(drift_times) clean
            di as result "  `drift_times'"
        }
        
        return scalar N_drifts = `N_drifts'
        return scalar N_warnings = `N_warnings'
        return scalar ddm_warning = `ddm_warning'
        return scalar ddm_drift = `ddm_drift'
        return local algo_name "ddm"
        exit 0
    }
    
    // =========================================================================
    // EDDM Early Drift Detection Method
    // Baena-Garcia et al. (2006). Default mode adapts distances between errors
    // to DEA signal deviations; original mode uses a supplied 0/1 error stream.
    // =========================================================================
    
    if "`algorithm'" == "eddm" {
        di as text "Algorithm: EDDM (Early Drift Detection Method)"
        di as text "{hline 50}"
        
        // EDDM thresholds (based on ratio p' + 2s' / p'_max + 2s'_max)
        local eddm_alpha = 0.95   // warning threshold
        local eddm_beta = 0.90    // drift threshold
        local min_samples = 30
        
        tempvar eddm_dist eddm_alarm eddm_warning_flag
        
        quietly {
            // Distance between "errors" (deviations from 1)
            gen double `eddm_dist' = . if `touse'
            gen byte `eddm_alarm' = 0 if `touse'
            gen byte `eddm_warning_flag' = 0 if `touse'
        }
        
        tempvar time_order
        quietly gen long `time_order' = _n
        quietly sort `time' `time_order'
        
        if "`driftmode'" == "original" {
            local eddm_signal "`drifterror'"
            local eddm_binary = 1
        }
        else {
            local eddm_signal "`alpha'"
            local eddm_binary = 0
        }
        
        // Compute EDDM statistics using Mata
        mata: _sdea_eddm("`eddm_signal'", "`touse'", "`eddm_alarm'", "`eddm_warning_flag'", ///
                         `eddm_alpha', `eddm_beta', `min_samples', `eddm_binary')
        
        quietly replace `drift' = `eddm_alarm' if `touse'
        
        quietly sort `time_order'
        
        tempvar eddm_period_alarm
        bysort `time': egen `eddm_period_alarm' = max(`eddm_alarm') if `touse'
        quietly count if `eddm_period_alarm' == 1 & `period_tag' == 1 & `touse'
        local N_drifts = r(N)
        
        quietly count if `eddm_warning_flag' == 1 & `period_tag' == 1 & `touse'
        local N_warnings = r(N)
        
        di as text "EDDM alpha (warning): " as result `eddm_alpha'
        di as text "EDDM beta (drift): " as result `eddm_beta'
        if "`driftmode'" == "original" {
            di as text "Input stream: " as result "`drifterror'" as text " (0/1 original EDDM error indicator)"
        }
        else {
            di as text "Input stream: DEA signal deviations from alpha = 1"
        }
        di as text "Min samples: " as result `min_samples'
        di as text "Warning events: " as result `N_warnings'
        di as text "Drift events detected: " as result `N_drifts'
        
        if `N_drifts' > 0 & `N_drifts' <= 20 {
            di as text _n "Drift time points:"
            quietly levelsof `time' if `eddm_period_alarm' == 1 & `period_tag' == 1 & `touse', local(drift_times) clean
            di as result "  `drift_times'"
        }
        
        return scalar N_drifts = `N_drifts'
        return scalar N_warnings = `N_warnings'
        return scalar eddm_alpha = `eddm_alpha'
        return scalar eddm_beta = `eddm_beta'
        return local algo_name "eddm"
        exit 0
    }

    // =========================================================================
    // Frontier Relationship Drift
    // Distinguishes frontier expansion/covariate support changes from concept
    // drift by testing whether the conditional input-output relation shifts.
    // =========================================================================

    if "`algorithm'" == "frontier" {
        if "`inputs'" == "" | "`outputs'" == "" {
            di as error "driftalgo(frontier) requires inputs() and outputs()"
            exit 198
        }

        di as text "Algorithm: Frontier relationship drift"
        di as text "{hline 50}"
        local frontier_thr = max(1.5, `thrvol')
        di as text "Window size:              " as result `window'
        di as text "Frontier z threshold:     " as result %5.2f `frontier_thr'
        if "`by'" != "" di as text "Heterogeneity groups:     " as result "`by'"
        di as text "{hline 50}"

        tempvar frontier_stat frontier_alarm
        quietly {
            gen double `frontier_stat' = . if `touse'
            gen byte `frontier_alarm' = 0 if `touse'
        }

        local lninputs
        foreach x of varlist `inputs' {
            tempvar lnx
            quietly gen double `lnx' = ln(`x') if `touse' & `x' > 0
            local lninputs `lninputs' `lnx'
        }
        local k_inputs : word count `lninputs'

        local group_levels "__pooled__"
        if "`by'" != "" {
            quietly levelsof `by' if `touse', local(group_levels) clean
        }

        foreach g of local group_levels {
            if "`by'" == "" local group_cond "1"
            else local group_cond "(`by' == `g')"

            local t_num = 0
            foreach t of local time_levels {
                local t_num = `t_num' + 1
                if `t_num' > `window' {
                    local w_start = `t_num' - `window'
                    local train_cond ""
                    forvalues w = `w_start'/`=`t_num'-1' {
                        local tw : word `w' of `time_levels'
                        if "`train_cond'" == "" local train_cond "(`time' == `tw')"
                        else local train_cond "`train_cond' | (`time' == `tw')"
                    }

                    local max_z = 0
                    foreach y of varlist `outputs' {
                        tempvar lny xb res
                        quietly gen double `lny' = ln(`y') if `touse' & `y' > 0
                        capture quietly regress `lny' `lninputs' if `touse' & (`train_cond') & `group_cond'
                        if _rc == 0 & e(N) > `k_inputs' + 5 & e(rmse) > 1e-10 {
                            quietly predict double `xb' if `touse' & `time' == `t' & `group_cond', xb
                            quietly gen double `res' = `lny' - `xb' if `touse' & `time' == `t' & `group_cond'
                            quietly summarize `res' if `touse' & `time' == `t' & `group_cond', meanonly
                            if r(N) > 0 & !missing(r(mean)) {
                                local z = abs(r(mean)) / e(rmse)
                                local max_z = max(`max_z', `z')
                            }
                        }
                        capture drop `lny' `xb' `res'
                    }

                    quietly replace `frontier_stat' = max(cond(missing(`frontier_stat'), 0, `frontier_stat'), `max_z') ///
                        if `time' == `t' & `touse' & `group_cond'
                    quietly replace `frontier_alarm' = (`frontier_stat' > `frontier_thr') ///
                        if `time' == `t' & `touse' & `group_cond'
                }
            }
        }

        quietly replace `drift' = `frontier_stat' * `frontier_alarm' if `touse'
        quietly sort `__ord'

        quietly levelsof `time' if `frontier_alarm' == 1 & `touse', local(frontier_drift_times) clean
        local N_drifts : word count `frontier_drift_times'

        capture drop _sdea_drift_frontier
        quietly gen double _sdea_drift_frontier = `frontier_stat' if `touse'
        label variable _sdea_drift_frontier "Drift diagnostic: conditional frontier shift"

        quietly summarize `frontier_stat' if `period_tag' == 1 & `touse'
        if r(N) > 0 {
            di as text "Frontier shift: mean=" as result %6.3f r(mean) ///
                as text "  max=" as result %6.3f r(max)
        }
        di as text "Drift events detected: " as result `N_drifts'
        if `N_drifts' > 0 & `N_drifts' <= 20 {
            di as text _n "Drift time points:"
            di as result "  `frontier_drift_times'"
        }

        return scalar N_drifts = `N_drifts'
        return scalar window = `window'
        return scalar frontier_threshold = `frontier_thr'
        return local algo_name "frontier"
        exit 0
    }
    
    // =========================================================================
    // Multi-indicator Framework (default)
    // =========================================================================
    
    di as text "Algorithm: Multi-indicator"
    di as text "{hline 50}"
    di as text "Window size:              " as result `window'
    di as text "Volatility threshold:     " as result %5.2f `thrvol'
    di as text "Expansion threshold:      " as result %5.2f `threxp'
    di as text "Turnover threshold:       " as result %5.2f `thrturn'
    di as text "{hline 50}"
    
    // Create temporary variables for indicators
    tempvar V_t D_t T_t drift_signal n_triggered
    tempvar lag_n_eff lag_n_obs exp_rate
    
    quietly {
        // Initialize
        gen double `V_t' = . if `touse'
        gen double `D_t' = . if `touse'
        gen double `T_t' = . if `touse'
        gen byte `drift_signal' = 0 if `touse'
        
        // =====================================================================
        // Indicator 1: Efficiency Volatility (V_t)
        // V_t = |μ_t - MA^w_t| / σ(μ_{t-w}, ..., μ_t)
        // =====================================================================
        
        local t_num = 0
        foreach t of local time_levels {
            local t_num = `t_num' + 1
            
            if `t_num' > `window' {
                // Get window of past mean efficiencies
                local w_start = `t_num' - `window'
                local sum_mu = 0
                local sum_sq = 0
                local w_count = 0
                
                forvalues w = `w_start'/`t_num' {
                    local tw : word `w' of `time_levels'
                    summ `mean_eff' if `time' == `tw' & `period_tag' == 1, meanonly
                    if r(N) > 0 {
                        local sum_mu = `sum_mu' + r(mean)
                        local sum_sq = `sum_sq' + r(mean)^2
                        local w_count = `w_count' + 1
                    }
                }
                
                if `w_count' > 1 {
                    local ma_val = `sum_mu' / `w_count'
                    local var_val = (`sum_sq' / `w_count') - (`ma_val')^2
                    local sd_val = sqrt(max(`var_val', 1e-10))
                    
                    // Current mean efficiency
                    summ `mean_eff' if `time' == `t' & `period_tag' == 1, meanonly
                    local mu_t = r(mean)
                    
                    // Volatility indicator
                    local v_val = abs(`mu_t' - `ma_val') / `sd_val'
                    replace `V_t' = `v_val' if `time' == `t' & `touse'
                }
            }
        }
        
        // =====================================================================
        // Indicator 2: Frontier Expansion Rate (D_t)
        // With BuildHull history: entry share among arrivals in period t.
        // Fallback: change in final frame counts, retained for non-history methods.
        // =====================================================================
        
        // Compute expansion rate for each period
        local prev_n_eff = .
        local prev_n_obs = .
        local t_num = 0
        
        gen double `exp_rate' = . if `touse'
        
        foreach t of local time_levels {
            local t_num = `t_num' + 1
            
            summ `n_obs' if `time' == `t' & `period_tag' == 1, meanonly
            local curr_n_obs = r(mean)
            
            if `has_frame_history' {
                summ `hist_frame_entry' if `time' == `t' & `period_tag' == 1, meanonly
                local curr_entries = r(mean)
                if `curr_n_obs' > 0 {
                    local r_t = `curr_entries' / `curr_n_obs'
                    replace `exp_rate' = `r_t' if `time' == `t' & `touse'
                }
            }
            else {
                summ `n_eff' if `time' == `t' & `period_tag' == 1, meanonly
                local curr_n_eff = r(mean)
                
                if `t_num' > 1 & !missing(`prev_n_eff') {
                    local denom = `curr_n_obs' - `prev_n_obs'
                    if abs(`denom') > 0.5 {
                        local r_t = (`curr_n_eff' - `prev_n_eff') / `denom'
                        replace `exp_rate' = `r_t' if `time' == `t' & `touse'
                    }
                }
                local prev_n_eff = `curr_n_eff'
            }
            
            local prev_n_obs = `curr_n_obs'
        }
        
        // Compute D_t based on deviation from expected expansion rate
        local t_num = 0
        foreach t of local time_levels {
            local t_num = `t_num' + 1
            
            if `t_num' > `window' + 1 {
                local w_start = `t_num' - `window'
                local sum_r = 0
                local w_count = 0
                
                forvalues w = `w_start'/`=`t_num'-1' {
                    local tw : word `w' of `time_levels'
                    summ `exp_rate' if `time' == `tw' & `period_tag' == 1, meanonly
                    if r(N) > 0 & !missing(r(mean)) {
                        local sum_r = `sum_r' + r(mean)
                        local w_count = `w_count' + 1
                    }
                }
                
                if `w_count' > 0 {
                    local r_expected = `sum_r' / `w_count'
                    
                    summ `exp_rate' if `time' == `t' & `period_tag' == 1, meanonly
                    if r(N) > 0 & !missing(r(mean)) & abs(`r_expected') > 1e-6 {
                        local r_actual = r(mean)
                        local d_val = abs(`r_actual' - `r_expected') / abs(`r_expected')
                        replace `D_t' = min(`d_val', 10) if `time' == `t' & `touse'
                    }
                }
            }
        }
        
        // =====================================================================
        // Indicator 3: Reference Set Turnover (T_t)
        // With BuildHull history: (entries + exits) / union size.
        // Fallback: approximate by final frame-size change ratio.
        // =====================================================================
        
        local prev_n_frame = .
        local t_num = 0
        
        foreach t of local time_levels {
            local t_num = `t_num' + 1
            
            if `has_frame_history' {
                summ `hist_frame_size' if `time' == `t' & `period_tag' == 1, meanonly
                local curr_n_frame = r(mean)
                summ `hist_frame_entry' if `time' == `t' & `period_tag' == 1, meanonly
                local curr_entries = r(mean)
                summ `hist_frame_pruned' if `time' == `t' & `period_tag' == 1, meanonly
                local curr_pruned = r(mean)
                
                if `t_num' > 1 & !missing(`prev_n_frame') & `prev_n_frame' > 0 {
                    local union_size = `prev_n_frame' + `curr_entries'
                    if `union_size' > 0 {
                        local turnover = (`curr_entries' + `curr_pruned') / `union_size'
                        replace `T_t' = min(`turnover', 1) if `time' == `t' & `touse'
                    }
                }
            }
            else {
                summ `n_frame_t' if `time' == `t' & `period_tag' == 1, meanonly
                local curr_n_frame = r(mean)
            
                if `t_num' > 1 & !missing(`prev_n_frame') & `prev_n_frame' > 0 {
                    local min_size = min(`curr_n_frame', `prev_n_frame')
                    local max_size = max(`curr_n_frame', `prev_n_frame')
                    
                    if `max_size' > 0 {
                        local overlap_est = `min_size' / `max_size'
                        local turnover = 1 - `overlap_est'
                        replace `T_t' = `turnover' if `time' == `t' & `touse'
                    }
                }
            }
            
            local prev_n_frame = `curr_n_frame'
        }
        
        // =====================================================================
        // Combined drift signal: at least 2 of 3 indicators exceed thresholds
        // =====================================================================
        
        gen int `n_triggered' = 0 if `touse'
        replace `n_triggered' = `n_triggered' + (`V_t' > `thrvol' & !missing(`V_t')) if `touse'
        replace `n_triggered' = `n_triggered' + (`D_t' > `threxp' & !missing(`D_t')) if `touse'
        replace `n_triggered' = `n_triggered' + (`T_t' > `thrturn' & !missing(`T_t')) if `touse'
        
        // Drift if >= 2 indicators triggered
        replace `drift_signal' = (`n_triggered' >= 2) if `touse'
        
        // Store drift magnitude (max normalized indicator)
        tempvar drift_magnitude
        gen double `drift_magnitude' = 0 if `touse'
        replace `drift_magnitude' = max(`V_t'/`thrvol', 0) if !missing(`V_t') & `touse'
        replace `drift_magnitude' = max(`drift_magnitude', `D_t'/`threxp') if !missing(`D_t') & `touse'
        replace `drift_magnitude' = max(`drift_magnitude', `T_t'/`thrturn') if !missing(`T_t') & `touse'
        
        replace `drift' = `drift_magnitude' * `drift_signal' if `touse'

        capture drop _sdea_drift_V _sdea_drift_D _sdea_drift_T _sdea_drift_triggers
        gen double _sdea_drift_V = `V_t' if `touse'
        gen double _sdea_drift_D = `D_t' if `touse'
        gen double _sdea_drift_T = `T_t' if `touse'
        gen int _sdea_drift_triggers = `n_triggered' if `touse'
        label variable _sdea_drift_V "Drift diagnostic: efficiency volatility"
        label variable _sdea_drift_D "Drift diagnostic: frontier expansion"
        label variable _sdea_drift_T "Drift diagnostic: frame turnover"
        label variable _sdea_drift_triggers "Drift diagnostic: triggered indicators"
    }
    
    // Sort back
    quietly sort `__ord'
    
    // =========================================================================
    // Report results
    // =========================================================================
    
    // Count drifts at period level
    quietly count if `drift' > 0 & !missing(`drift') & `period_tag' == 1 & `touse'
    local N_drifts = r(N)
    
    // Count triggered indicators
    quietly count if `V_t' > `thrvol' & `period_tag' == 1 & `touse'
    local n_vol_triggered = r(N)
    quietly count if `D_t' > `threxp' & `period_tag' == 1 & `touse'
    local n_exp_triggered = r(N)
    quietly count if `T_t' > `thrturn' & `period_tag' == 1 & `touse'
    local n_turn_triggered = r(N)
    
    di as text _n "Indicator Summary"
    di as text "{hline 50}"
    
    quietly summ `V_t' if `period_tag' == 1 & `touse'
    if r(N) > 0 {
        di as text "Volatility (V):   mean=" as result %6.3f r(mean) ///
            as text "  max=" as result %6.3f r(max) ///
            as text "  triggered=" as result `n_vol_triggered'
    }
    
    quietly summ `D_t' if `period_tag' == 1 & `touse'
    if r(N) > 0 {
        di as text "Expansion (D):    mean=" as result %6.3f r(mean) ///
            as text "  max=" as result %6.3f r(max) ///
            as text "  triggered=" as result `n_exp_triggered'
    }
    
    quietly summ `T_t' if `period_tag' == 1 & `touse'
    if r(N) > 0 {
        di as text "Turnover (T):     mean=" as result %6.3f r(mean) ///
            as text "  max=" as result %6.3f r(max) ///
            as text "  triggered=" as result `n_turn_triggered'
    }
    
    di as text "{hline 50}"
    di as text "Drift events (>=2 indicators): " as result `N_drifts'
    
    if `N_drifts' > 0 & `N_drifts' <= 20 {
        di as text _n "Drift time points:"
        quietly levelsof `time' if `drift' > 0 & !missing(`drift') & `period_tag' == 1 & `touse', local(drift_times) clean
        di as result "  `drift_times'"
    }
    
    // Return values
    return scalar N_drifts = `N_drifts'
    return scalar window = `window'
    return scalar thr_volatility = `thrvol'
    return scalar thr_expansion = `threxp'
    return scalar thr_turnover = `thrturn'
    return scalar n_vol_triggered = `n_vol_triggered'
    return scalar n_exp_triggered = `n_exp_triggered'
    return scalar n_turn_triggered = `n_turn_triggered'
    return local algo_name "multi"
end

*-------------------------------------------------------------------------------
* _sdea_graph: Visualization
*-------------------------------------------------------------------------------
program define _sdea_graph
    version 16.0
    syntax varname(numeric), Time(varname numeric) DRift(varname numeric) ///
        FRame(varname) Touse(varname) [SAVing(string)]
    
    local alpha `varlist'
    
    // 시간별 통계 계산
    tempvar mean_drift_by_time max_drift_by_time time_tag drift_detected
    quietly {
        // 시간별 평균/최대 drift 계산
        bysort `time': egen `mean_drift_by_time' = mean(`drift') if `touse'
        bysort `time': egen `max_drift_by_time' = max(`drift') if `touse'
        egen `time_tag' = tag(`time') if `touse'
        // drift가 감지된 시점 표시 (해당 시간의 max drift > 0)
        gen byte `drift_detected' = (`time_tag' == 1 & `max_drift_by_time' > 0) if `touse'
    }
    
    // Panel 1: Evolution (increased height ratio for y-axis labels)
    twoway (scatter `alpha' `time' if `frame'==0 & `touse', ///
                msize(vsmall) mcolor(navy%30)) ///
           (scatter `alpha' `time' if `frame'==1 & `touse', ///
                msize(small) mcolor(blue) msymbol(D)) ///
           (scatter `alpha' `time' if `drift_detected'==1 & `touse', ///
                msize(large) mcolor(red) msymbol(Oh)), ///
           legend(order(1 "Interior" 2 "Frame" 3 "Drift") ///
                  pos(1) ring(0) col(1) size(vsmall)) ///
           ytitle("α (radial distance)", size(small)) ///
           xtitle("Time", size(small)) ///
           ylabel(, labsize(vsmall) angle(0) nogrid) ///
           title("Streaming DEA Evolution", size(medsmall) pos(11)) ///
           name(sdea_evolution, replace)
    
    // Panel 2: Drift magnitude
    quietly summ `drift' if `touse' & `drift' > 0 & !missing(`drift')
    if r(N) > 0 {
        twoway (line `mean_drift_by_time' `time' if `time_tag'==1 & `touse', lcolor(blue) lwidth(medium)) ///
               (scatter `max_drift_by_time' `time' if `drift_detected'==1 & `touse', ///
                    msize(medium) mcolor(red) msymbol(Oh)), ///
               legend(order(1 "Drift magnitude" 2 "Major drifts") ///
                      pos(1) ring(0) col(1) size(vsmall)) ///
               ytitle("Drift magnitude", size(small)) ///
               xtitle("Time", size(small)) ///
               ylabel(, labsize(vsmall) angle(0) nogrid) ///
               title("Drift Diagnostics", size(medsmall) pos(11)) ///
               name(sdea_drift, replace)
    }
    
    // Panel 3: Frame size evolution
    quietly {
        tempvar frame_count tag
        bysort `time': egen `frame_count' = total(`frame')
        egen `tag' = tag(`time')
    }
    twoway (line `frame_count' `time' if `tag'==1, ///
                lcolor(blue) lwidth(medium)), ///
           ytitle("Frame size", size(small)) ///
           xtitle("Time", size(small)) ///
           ylabel(, labsize(vsmall) angle(0) nogrid) ///
           title("Frame Evolution", size(medsmall) pos(11)) ///
           name(sdea_frame, replace)
    
    // Combine (increased ysize for better y-axis label spacing)
    quietly count if `drift' > 0 & !missing(`drift') & `touse'
    if r(N) > 0 {
        graph combine sdea_evolution sdea_drift sdea_frame, ///
            cols(1) ysize(10) xsize(6) ///
            title("Streaming DEA Analysis", size(medium)) ///
            name(sdea_combined, replace)
    }
    else {
        graph combine sdea_evolution sdea_frame, ///
            cols(1) ysize(7) xsize(6) ///
            title("Streaming DEA Analysis", size(medium)) ///
            name(sdea_combined, replace)
    }
    
    // Save
    if "`saving'" != "" {
        graph export "`saving'.png", replace width(1200) name(sdea_combined)
        di as text "(Graph saved as `saving'.png)"
    }
    else {
        di as text "(Graphs: sdea_evolution, sdea_drift, sdea_frame, sdea_combined)"
    }
end
   

*===============================================================================
* MATA FUNCTIONS
*===============================================================================

// Drop existing functions to prevent conflicts
capture mata: mata drop _sdea_solve_lp()
capture mata: mata drop _sdea_prune_frame()
capture mata: mata drop _sdea_buildhull_dual_plane()
capture mata: mata drop _sdea_buildhull_dula_main()
capture mata: mata drop _sdea_buildhull_eval()
capture mata: mata drop _sdea_buildhull_main()
capture mata: mata drop _sdea_standard_eval()
capture mata: mata drop _sdea_run_dea()
capture mata: mata drop _sdea_time_streaming()
capture mata: mata drop _sdea_save_reference()
capture mata: mata drop _sdea_save_history()

capture mata: mata which _rsm_solve_lp()
if _rc {
    capture findfile _dea_simplex.ado
    if _rc {
        di as error "_dea_simplex.ado not found. It is required for solver(simplex) and solver(auto)."
        exit 601
    }
    quietly do "`r(fn)'"
}

mata:
mata set matastrict on

//==============================================================================
// LP Solver - IPM only (requires Stata 16+)
//==============================================================================

struct _sdea_lp_result {
    real scalar converged
    real scalar value
    real rowvector solution
}

struct _sdea_lp_result scalar _sdea_solve_lp(
    real rowvector c,
    real matrix A_ineq,
    real colvector b_ineq,
    real matrix A_eq,
    real colvector b_eq,
    real rowvector lb,
    real rowvector ub,
    string scalar minmax,
    string scalar solver)
{
    struct _sdea_lp_result scalar result
    struct RSM_Result scalar simplex_result
    class LinearProgram scalar lp
    string scalar meter
    real scalar used_simplex
    real rowvector simplex_c

    meter = st_local("_sdea_static_meter")
    used_simplex = 0
    result.converged = 0
    result.value = .
    result.solution = J(1, cols(c), .)

    if (solver == "") solver = "ipm"

    if (solver == "simplex") {
        if (min(lb) >= 0 & hasmissing(ub)) {
            used_simplex = 1
            simplex_c = c
            if (minmax == "max") simplex_c = -c
            simplex_result = _rsm_solve_lp(A_ineq, b_ineq, A_eq, b_eq,
                simplex_c, 1e-8)
            if (simplex_result.status == 0) {
                result.converged = 1
                result.value = simplex_result.objval
                if (minmax == "max") result.value = -simplex_result.objval
                result.solution = simplex_result.x'
            }
        }
        if (meter!="") st_matrix(meter,st_matrix(meter)+(0,used_simplex,!result.converged))
        return(result)
    }
    
    lp = LinearProgram()
    lp.setCoefficients(c)
    lp.setMaxOrMin(minmax)
    if (rows(A_ineq) > 0) lp.setInequality(A_ineq, b_ineq)
    if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
    lp.setBounds(lb, ub)
    (void) lp.optimize()
    
    result.converged = (lp.converged() == 1 & lp.errorcode() == 0)
    if (result.converged) {
        result.value = lp.value()
        result.solution = lp.parameters()
    }
    else {
        result.value = .
        result.solution = J(1, cols(c), .)
    }

    if (!result.converged & solver == "auto" & min(lb) >= 0 & hasmissing(ub)) {
        used_simplex = 1
        simplex_c = c
        if (minmax == "max") simplex_c = -c
        simplex_result = _rsm_solve_lp(A_ineq, b_ineq, A_eq, b_eq,
            simplex_c, 1e-8)
        if (simplex_result.status == 0) {
            result.converged = 1
            result.value = simplex_result.objval
            if (minmax == "max") result.value = -simplex_result.objval
            result.solution = simplex_result.x'
        }
    }
    
    if (meter!="") st_matrix(meter,st_matrix(meter)+(1,used_simplex,!result.converged))
    return(result)
}

//==============================================================================
// Prune Frame: Remove dominated DMUs from reference set
// Returns number of pruned DMUs; modifies Xf, Yf, frame_idx in place
//==============================================================================

real scalar _sdea_prune_frame(
    real matrix Xf,
    real matrix Yf,
    real colvector frame_idx,
    string scalar solver
) {
    real scalar p, k, n_pruned, a_val
    real matrix Xtemp, Ytemp
    real colvector keep, idx_temp
    
    p = rows(Xf)
    if (p <= 1) return(0)
    
    keep = J(p, 1, 1)
    
    // Check each DMU in frame: is it dominated by others?
    for (k=1; k<=p; k++) {
        // Build reference set excluding DMU k
        Xtemp = Xf[selectindex((1::p) :!= k), ]
        Ytemp = Yf[selectindex((1::p) :!= k), ]
        
        // Evaluate DMU k against remaining frame
        a_val = _sdea_buildhull_eval(Xtemp, Ytemp, Xf[k,]', Yf[k,]', solver)
        
        // If alpha <= 0 (or very small), DMU k is dominated
        if (!missing(a_val) & a_val < 1e-8) {
            keep[k] = 0
        }
    }
    
    n_pruned = p - sum(keep)
    
    if (n_pruned > 0) {
        idx_temp = selectindex(keep)
        Xf = Xf[idx_temp, ]
        Yf = Yf[idx_temp, ]
        frame_idx = frame_idx[idx_temp]
    }
    
    return(n_pruned)
}

//==============================================================================
// BuildHull Evaluation (Dulá et al., 2011)
// min α s.t. Σx_kλ_k - α ≤ x_0, Σy_kλ_k + α ≥ y_0, Σλ_k = 1
//==============================================================================

real scalar _sdea_buildhull_eval(
    real matrix Xframe,
    real matrix Yframe,
    real colvector x0,
    real colvector y0,
    string scalar solver, | real rowvector lambda
) {
    struct _sdea_lp_result scalar result
    real scalar p, m, s, i, r
    real matrix A
    real colvector b
    real rowvector c
    real matrix A_eq
    real colvector b_eq
    
    p = rows(Xframe)
    m = cols(Xframe)
    s = cols(Yframe)
    
    if (p == 0) return(0)
    
    // min α
    c = (J(1, p, 0), 1)
    
    // Σ x_k λ_k - α ≤ x_0
    A = J(m, p+1, 0)
    b = J(m, 1, 0)
    for (i=1; i<=m; i++) {
        A[i, 1..p] = Xframe[,i]'
        A[i, p+1] = -1
        b[i] = x0[i]
    }
    
    // Σ y_k λ_k + α ≥ y_0 → -Σ y_k λ_k - α ≤ -y_0
    A = A \ J(s, p+1, 0)
    b = b \ (-y0)
    for (r=1; r<=s; r++) {
        A[m+r, 1..p] = -Yframe[,r]'
        A[m+r, p+1] = -1
    }
    
    // Equality: Σλ_k = 1
    A_eq = (J(1,p,1), 0)
    b_eq = 1
    
    // Solve LP
    result = _sdea_solve_lp(c, A, b, A_eq, b_eq, J(1,p+1,0), J(1,p+1,.), "min", solver)
    
    if (result.converged) {
        if (args()==6) lambda = result.solution[1..p]
        return(result.solution[p+1])
    }
    return(.)
}

real rowvector _sdea_buildhull_dual_plane(
    real matrix Aframe,
    real rowvector b0,
    string scalar solver
) {
    real scalar p, d, i
    real matrix Aineq
    real colvector bineq
    real rowvector c, lb, ub, sol, plane
    struct _sdea_lp_result scalar result

    p = rows(Aframe)
    d = cols(Aframe)
    if (p == 0) return(J(1, d, .))

    // Dual of min delta: A_F lambda + e delta >= b, sum(lambda)=1.
    // Variables are u >= 0 and w = w_plus - w_minus.
    // max b'u + w subject to A_F u + w <= 0 and e'u <= 1.
    c = (b0, 1, -1)
    Aineq = J(p + 1, d + 2, 0)
    bineq = J(p + 1, 1, 0)
    for (i=1; i<=p; i++) {
        Aineq[i, 1..d] = Aframe[i,]
        Aineq[i, d+1] = 1
        Aineq[i, d+2] = -1
    }
    Aineq[p+1, 1..d] = J(1, d, 1)
    bineq[p+1] = 1

    lb = J(1, d+2, 0)
    ub = J(1, d+2, .)
    result = _sdea_solve_lp(c, Aineq, bineq, J(0, d+2, .), J(0, 1, .), lb, ub, "max", solver)

    if (result.converged) {
        sol = result.solution
        plane = sol[1..d]
        return(plane)
    }
    return(J(1, d, .))
}

//==============================================================================
// Dulá-style BuildHull Phase 1
// Identifies a VRS frame by hyperplane translation. Inputs are converted to
// undesirable attributes, so the all-desirable vector is a=(-x,y).
//==============================================================================

void _sdea_buildhull_dula_main(string scalar inputvars, string scalar outputvars,
    string scalar timevar, string scalar obsvar, string scalar alphavar,
    string scalar framevar, string scalar driftvar, string scalar touse,
    string scalar solver) {

    real matrix X, Y, A, Af
    real colvector obs_id, frame_idx, active, all_alpha, all_drift
    real colvector all_frame_size, all_frame_entry, all_frame_pruned
    real scalar n, d, h, k, bpos, bidx, a_val, n_interior, n_exterior, nframe
    real scalar best_pos, best_idx, best_score, score, newvar
    real rowvector plane
    real colvector candidate

    X = st_data(., tokens(inputvars), touse)
    Y = st_data(., tokens(outputvars), touse)
    obs_id = st_data(., obsvar, touse)
    A = (-X, Y)
    n = rows(A)
    d = cols(A)

    printf("\nDula BuildHull Phase 1\n")
    printf("------------------------------------------------------------\n")
    printf("Total: %g, Attributes: %g\n", n, d)

    all_alpha = J(n, 1, .)
    all_drift = J(n, 1, 0)
    all_frame_size = J(n, 1, 0)
    all_frame_entry = J(n, 1, 0)
    all_frame_pruned = J(n, 1, 0)
    frame_idx = J(0, 1, .)

    // Initialize with coordinate-maximizing observations. These are natural
    // extreme candidates for the all-desirable representation.
    for (h=1; h<=d; h++) {
        candidate = selectindex(A[,h] :== max(A[,h]))
        if (rows(candidate) > 0) {
            k = candidate[1]
            if (sum(frame_idx :== k) == 0) {
                frame_idx = frame_idx \ k
                all_frame_entry[k] = 1
                all_frame_size[k] = rows(frame_idx)
            }
        }
    }
    if (rows(frame_idx) == 0) {
        frame_idx = 1
        all_frame_entry[1] = 1
        all_frame_size[1] = 1
    }

    active = J(n, 1, 1)
    active[frame_idx] = J(rows(frame_idx), 1, 0)
    n_interior = 0
    n_exterior = 0

    while (sum(active) > 0) {
        bpos = selectindex(active)[1]
        bidx = bpos
        Af = A[frame_idx,]
        a_val = _sdea_buildhull_eval(-Af[,1..cols(X)], Af[,(cols(X)+1)..d], -A[bidx,1..cols(X)]', A[bidx,(cols(X)+1)..d]', solver)
        all_alpha[bidx] = a_val

        if (!missing(a_val) & a_val < 1e-8) {
            active[bidx] = 0
            n_interior++
            all_frame_size[bidx] = rows(frame_idx)
        }
        else if (!missing(a_val)) {
            plane = _sdea_buildhull_dual_plane(Af, A[bidx,], solver)
            if (missing(plane[1])) {
                best_idx = bidx
            }
            else {
                best_score = .
                best_idx = bidx
                candidate = selectindex(active)
                for (k=1; k<=rows(candidate); k++) {
                    score = A[candidate[k],] * plane'
                    if (missing(best_score) | score > best_score + 1e-10) {
                        best_score = score
                        best_idx = candidate[k]
                    }
                }
            }

            frame_idx = frame_idx \ best_idx
            active[best_idx] = 0
            all_frame_entry[best_idx] = 1
            all_frame_size[best_idx] = rows(frame_idx)
            if (best_idx != bidx) all_frame_size[bidx] = rows(frame_idx)
            n_exterior++
        }
        else {
            active[bidx] = 0
            all_frame_size[bidx] = rows(frame_idx)
        }
    }

    nframe = rows(frame_idx)
    printf("------------------------------------------------------------\n")
    printf("Completed: Frame=%g, Interior tests=%g, Exterior tests=%g\n", nframe, n_interior, n_exterior)
    printf("------------------------------------------------------------\n")

    stata("capture drop " + alphavar)
    newvar = st_addvar("double", alphavar)
    st_store(., alphavar, touse, all_alpha)

    stata("capture drop " + framevar)
    newvar = st_addvar("byte", framevar)
    st_store(., framevar, touse, J(n, 1, 0))
    if (rows(frame_idx) > 0) {
        st_store(obs_id[frame_idx], framevar, J(rows(frame_idx), 1, 1))
    }

    stata("capture drop " + driftvar)
    newvar = st_addvar("double", driftvar)
    st_store(., driftvar, touse, all_drift)

    stata("capture drop _sdea_frame_size")
    newvar = st_addvar("long", "_sdea_frame_size")
    st_store(., "_sdea_frame_size", touse, all_frame_size)

    stata("capture drop _sdea_frame_entry")
    newvar = st_addvar("byte", "_sdea_frame_entry")
    st_store(., "_sdea_frame_entry", touse, all_frame_entry)

    stata("capture drop _sdea_frame_pruned")
    newvar = st_addvar("long", "_sdea_frame_pruned")
    st_store(., "_sdea_frame_pruned", touse, all_frame_pruned)
}

//==============================================================================
// Unified Streaming DEA Function
// Combines: fullstream, cumulative, moving
//==============================================================================

void _sdea_streaming_unified(string scalar inputvars, string scalar outputvars,
    string scalar timevar, string scalar method, real scalar window,
    string scalar alphavar, string scalar framevar, string scalar touse,
    string scalar orientation, string scalar rts, string scalar verbose)
{
    real matrix X, Y, X_ref, Y_ref
    real colvector time_vec, alpha, frame_flag, periods
    real scalar n, i, t, p, theta, n_efficient, n_frame
    real scalar current_period, win_start, win_end
    real rowvector x0, y0
    
    // Load data
    X = st_data(., tokens(inputvars), touse)
    Y = st_data(., tokens(outputvars), touse)
    time_vec = st_data(., timevar, touse)
    
    n = rows(X)
    periods = uniqrows(time_vec)
    
    alpha = J(n, 1, .)
    frame_flag = J(n, 1, 0)
    
    if (verbose != "") {
        printf("{txt}%s: %g observations, %g periods\n", strupper(method), n, rows(periods))
        if (method == "moving") printf("{txt}  Window size: %g\n", window)
    }
    
    // Calculate minimum DMU requirement: max(m*s, 3*(m+s))
    real scalar m_inp, s_out, min_dmu
    m_inp = cols(X)
    s_out = cols(Y)
    min_dmu = max((m_inp * s_out, 3 * (m_inp + s_out)))
    
    if (verbose != "") {
        printf("{txt}  Minimum DMU rule of thumb: %g (inputs=%g, outputs=%g)\n", min_dmu, m_inp, s_out)
    }
    
    //==========================================================================
    // FULLSTREAM: Sequential single-pass, expanding frontier
    //==========================================================================
    if (method == "fullstream") {
        // Initialize frontier with first min_dmu observations (or first period)
        // This ensures sufficient DMUs for DEA estimation
        real scalar init_size
        init_size = min((min_dmu, n))
        
        // Initialize with first init_size observations
        X_ref = X[1::init_size, .]
        Y_ref = Y[1::init_size, .]
        
        // Mark initial observations as frontier members with efficiency = 1
        for (i = 1; i <= init_size; i++) {
            alpha[i] = 1
            frame_flag[i] = 1
        }
        n_frame = init_size
        
        if (verbose != "") {
            printf("{txt}  Initial frontier size: %g (min DMU rule)\n", init_size)
        }
        
        // Process remaining observations
        for (i = init_size + 1; i <= n; i++) {
            x0 = X[i, .]
            y0 = Y[i, .]
            
            // Add new DMU to frontier FIRST (self-evaluation)
            X_ref = X_ref \ x0
            Y_ref = Y_ref \ y0
            
            // Evaluate against frontier INCLUDING itself
            theta = _sdea_solve_dea_basic(X_ref, Y_ref, x0, y0, orientation, rts)
            
            if (theta == .) {
                alpha[i] = .
                // Remove from frontier if infeasible
                X_ref = X_ref[1::rows(X_ref)-1, .]
                Y_ref = Y_ref[1::rows(Y_ref)-1, .]
            }
            else {
                alpha[i] = theta
                // Check if efficient (keep in frontier) or not (remove)
                // Input: theta=1 is efficient (minimize), theta<1 means can reduce inputs
                // Output: theta=1 is efficient (maximize), theta>1 means can expand outputs
                if ((orientation == "input" & theta >= 1 - 1e-6) |
                    (orientation != "input" & theta <= 1 + 1e-6)) {
                    // Efficient - keep in frontier
                    frame_flag[i] = 1
                    n_frame++
                }
                else {
                    // Not efficient - remove from frontier
                    X_ref = X_ref[1::rows(X_ref)-1, .]
                    Y_ref = Y_ref[1::rows(Y_ref)-1, .]
                }
            }
            
            if (verbose != "" & mod(i, 500) == 0) {
                printf("{txt}  Processed %g/%g, frame size: %g\n", i, n, n_frame)
            }
        }
        
        printf("{txt}FullStream Results:\n")
        printf("{txt}  Observations: %g\n", n)
        printf("{txt}  Frame members: %g (%5.2f%%)\n", n_frame, 100*n_frame/n)
    }
    //==========================================================================
    // CUMULATIVE: Expanding frontier (all observations up to current period)
    //==========================================================================
    else if (method == "cumulative") {
        X_ref = J(0, cols(X), .)
        Y_ref = J(0, cols(Y), .)
        
        for (t = 1; t <= rows(periods); t++) {
            current_period = periods[t]
            
            // Add all observations in current period to cumulative frontier
            for (i = 1; i <= n; i++) {
                if (time_vec[i] == current_period) {
                    X_ref = X_ref \ X[i, .]
                    Y_ref = Y_ref \ Y[i, .]
                }
            }
            
            // Only evaluate if frontier has minimum required DMUs
            if (rows(X_ref) >= min_dmu) {
                // Evaluate all observations in current period
                for (i = 1; i <= n; i++) {
                    if (time_vec[i] == current_period) {
                        theta = _sdea_solve_dea_basic(X_ref, Y_ref, X[i, .], Y[i, .], orientation, rts)
                        alpha[i] = theta
                        if (theta != . & abs(theta - 1) < 1e-6) frame_flag[i] = 1
                    }
                }
            }
            else {
                // Not enough DMUs yet - mark as efficient (initial period)
                for (i = 1; i <= n; i++) {
                    if (time_vec[i] == current_period) {
                        alpha[i] = 1
                        frame_flag[i] = 1
                    }
                }
            }
            
            if (verbose != "") {
                printf("{txt}  Period %g: frontier size = %g\n", current_period, rows(X_ref))
            }
        }
        
        printf("{txt}Cumulative DEA Results:\n")
        printf("{txt}  Observations: %g, Periods: %g\n", n, rows(periods))
    }
    //==========================================================================
    // MOVING: Sliding window frontier
    //==========================================================================
    else if (method == "moving") {
        for (t = 1; t <= rows(periods); t++) {
            current_period = periods[t]
            win_start = max((1, t - window + 1))
            win_end = t
            
            // Build window frontier
            X_ref = J(0, cols(X), .)
            Y_ref = J(0, cols(Y), .)
            
            for (i = 1; i <= n; i++) {
                for (p = 1; p <= rows(periods); p++) {
                    if (time_vec[i] == periods[p] & p >= win_start & p <= win_end) {
                        X_ref = X_ref \ X[i, .]
                        Y_ref = Y_ref \ Y[i, .]
                        break
                    }
                }
            }
            
            // Only evaluate if frontier has minimum required DMUs
            if (rows(X_ref) >= min_dmu) {
                // Evaluate observations in current period
                for (i = 1; i <= n; i++) {
                    if (time_vec[i] == current_period) {
                        theta = _sdea_solve_dea_basic(X_ref, Y_ref, X[i, .], Y[i, .], orientation, rts)
                        alpha[i] = theta
                        if (theta != . & abs(theta - 1) < 1e-6) frame_flag[i] = 1
                    }
                }
            }
            else {
                // Not enough DMUs - mark as efficient
                for (i = 1; i <= n; i++) {
                    if (time_vec[i] == current_period) {
                        alpha[i] = 1
                        frame_flag[i] = 1
                    }
                }
                if (verbose != "") {
                    printf("{txt}  Period %g: insufficient DMUs (%g < %g), skipped\n", 
                           current_period, rows(X_ref), min_dmu)
                }
            }
            
            if (verbose != "") {
                printf("{txt}  Period %g (window %g-%g): frontier = %g\n", 
                       current_period, periods[win_start], periods[win_end], rows(X_ref))
            }
        }
        
        printf("{txt}Moving Window DEA Results:\n")
        printf("{txt}  Observations: %g, Periods: %g, Window: %g\n", n, rows(periods), window)
    }
    
    // Store results
    st_store(., st_varindex(alphavar), touse, alpha)
    st_store(., st_varindex(framevar), touse, frame_flag)
    
    // Common summary
    n_efficient = sum(frame_flag)
    printf("{txt}  Efficient DMUs: %g (%5.2f%%)\n", n_efficient, 100*n_efficient/n)
    printf("{txt}  Mean efficiency: %9.6f\n", mean(select(alpha, alpha :< .)))
}

//==============================================================================
// ADWIN-style Algorithm (Adaptive Windowing)
// Bifet and Gavalda (2007), adapted to a DEA-derived univariate signal
//==============================================================================

void _sdea_adwin(string scalar alphavar, string scalar touse,
                  string scalar alarmvar, string scalar statvar,
                  real scalar delta, real scalar min_window)
{
    real colvector alpha, alarm, stat
    real scalar n, i, j, n0, n1, cut_found
    real scalar mu0, mu1, mu_total, m, epsilon_cut
    real scalar var0, var1, n_harm
    
    alpha = st_data(., alphavar, touse)
    n = rows(alpha)
    
    alarm = J(n, 1, 0)
    stat = J(n, 1, .)
    
    // ADWIN maintains a window W and looks for a cut point
    // where the averages of two sub-windows differ significantly
    
    for (i = min_window; i <= n; i++) {
        cut_found = 0
        
        // Try different cut points in the window [max(1, i-W)..i]
        for (j = max((1, i - 100)); j < i - min_window/2; j++) {
            // Sub-window 0: [j, cut), Sub-window 1: [cut, i]
            n0 = j
            n1 = i - j
            
            if (n0 < min_window/2 | n1 < min_window/2) continue
            
            mu0 = mean(alpha[1..j])
            mu1 = mean(alpha[(j+1)..i])
            mu_total = mean(alpha[1..i])
            
            // Compute Hoeffding bound
            m = 1 / (1/n0 + 1/n1)
            n_harm = 2 * n0 * n1 / (n0 + n1)
            
            // Epsilon cut threshold
            epsilon_cut = sqrt(ln(4/delta) / (2 * m))
            
            // Check if means differ significantly
            if (abs(mu0 - mu1) > epsilon_cut) {
                cut_found = 1
                alarm[i] = 1
                stat[i] = abs(mu0 - mu1) / epsilon_cut
                break
            }
        }
    }
    
    st_store(., st_varindex(alarmvar), touse, alarm)
    st_store(., st_varindex(statvar), touse, stat)
}

//==============================================================================
// DDM-style Drift Detection Method
// Gama et al. (2004), adapted to DEA signal deviations
//==============================================================================

void _sdea_ddm(string scalar errorvar, string scalar touse,
               string scalar alarmvar, string scalar warningvar,
               real scalar warning_level, real scalar drift_level,
               real scalar min_samples)
{
    real colvector err, alarm, warning
    real scalar n, i, p_i, s_i, p_min, s_min
    real scalar sum_err, sum_sq
    
    err = st_data(., errorvar, touse)
    n = rows(err)
    
    alarm = J(n, 1, 0)
    warning = J(n, 1, 0)
    
    // Initialize
    p_min = 1
    s_min = 0
    sum_err = 0
    sum_sq = 0
    
    for (i = 1; i <= n; i++) {
        // Update running statistics
        sum_err = sum_err + err[i]
        sum_sq = sum_sq + err[i]^2
        
        // Running mean (p_i) and std dev (s_i)
        p_i = sum_err / i
        if (i > 1) {
            s_i = sqrt((sum_sq/i - p_i^2) / i)
        }
        else {
            s_i = 0
        }
        
        // Update minimum
        if (i >= min_samples) {
            if (p_min + s_min == 0 | p_i + s_i < p_min + s_min) {
                p_min = p_i
                s_min = s_i
            }
            else if (p_i + s_i >= p_min + drift_level * s_min) {
                alarm[i] = 1
                // Reset after drift
                sum_err = 0
                sum_sq = 0
                p_min = 1
                s_min = 0
            }
            else if (p_i + s_i >= p_min + warning_level * s_min) {
                warning[i] = 1
            }
        }
    }
    
    st_store(., st_varindex(alarmvar), touse, alarm)
    st_store(., st_varindex(warningvar), touse, warning)
}

//==============================================================================
// EDDM Early Drift Detection Method
// Baena-Garcia et al. (2006). binary_error=1 uses the original 0/1 error stream;
// binary_error=0 adapts the error event to large DEA signal deviations.
//==============================================================================

void _sdea_eddm(string scalar signalvar, string scalar touse,
                string scalar alarmvar, string scalar warningvar,
                real scalar alpha_thr, real scalar beta_thr,
                real scalar min_samples, real scalar binary_error)
{
    real colvector signal_vec, alarm, warning
    real scalar n, i, dist, p_prime, s_prime
    real scalar p_prime_max, s_prime_max, ratio
    real scalar sum_dist, sum_sq_dist, n_errors
    real scalar last_error_pos, is_error
    
    signal_vec = st_data(., signalvar, touse)
    n = rows(signal_vec)
    
    alarm = J(n, 1, 0)
    warning = J(n, 1, 0)
    
    // Initialize
    p_prime_max = 0
    s_prime_max = 0
    sum_dist = 0
    sum_sq_dist = 0
    n_errors = 0
    last_error_pos = 0
    
    for (i = 1; i <= n; i++) {
        if (binary_error) {
            is_error = (signal_vec[i] != 0)
        }
        else {
            is_error = (abs(signal_vec[i] - 1) > 0.01)
        }
        
        if (is_error) {
            n_errors = n_errors + 1
            
            // Distance between consecutive errors
            if (last_error_pos > 0) {
                dist = i - last_error_pos
                sum_dist = sum_dist + dist
                sum_sq_dist = sum_sq_dist + dist^2
            }
            last_error_pos = i
            
            // Compute p' and s' (mean and std dev of distances)
            if (n_errors >= 2) {
                p_prime = sum_dist / (n_errors - 1)
                if (n_errors > 2) {
                    s_prime = sqrt(sum_sq_dist/(n_errors-1) - p_prime^2)
                }
                else {
                    s_prime = 0
                }
                
                // Update maximum
                if (p_prime + 2*s_prime > p_prime_max + 2*s_prime_max) {
                    p_prime_max = p_prime
                    s_prime_max = s_prime
                }
                
                // Check for warning/drift using ratio
                if (n_errors >= min_samples & p_prime_max + 2*s_prime_max > 0) {
                    ratio = (p_prime + 2*s_prime) / (p_prime_max + 2*s_prime_max)
                    
                    if (ratio < beta_thr) {
                        alarm[i] = 1
                        // Reset
                        p_prime_max = 0
                        s_prime_max = 0
                        sum_dist = 0
                        sum_sq_dist = 0
                        n_errors = 0
                        last_error_pos = 0
                    }
                    else if (ratio < alpha_thr) {
                        warning[i] = 1
                    }
                }
            }
        }
    }
    
    st_store(., st_varindex(alarmvar), touse, alarm)
    st_store(., st_varindex(warningvar), touse, warning)
}

//==============================================================================
// BuildHull Main Processing Loop
//==============================================================================

void _sdea_buildhull_main(string scalar inputvars, string scalar outputvars,
    string scalar timevar, string scalar startervar, string scalar obsvar, string scalar alphavar,
    string scalar framevar, string scalar driftvar, string scalar touse,
    real scalar prune_interval, string scalar solver,
    string scalar gatevar, real scalar gate_eps,
    string scalar mlscorevar, real scalar mlthreshold) {
    
    real matrix X, Y, Xf, Yf
    real colvector time, is_start, obs_id, all_alpha, all_drift, frame_idx
    real colvector all_frame_size, all_frame_entry, all_frame_pruned
    real colvector gate_status, mlscore
    real scalar n, nstart, j, a_val, last_exterior, nframe, n_interior, n_exterior
    real scalar n_pruned, total_pruned, prune_count, n_rejected, accept_update
    real scalar newvar
    
    X = st_data(., tokens(inputvars), touse)
    Y = st_data(., tokens(outputvars), touse)
    time = st_data(., timevar, touse)
    is_start = st_data(., startervar, touse)
    obs_id = st_data(., obsvar, touse)
    if (mlscorevar != "") {
        mlscore = st_data(., mlscorevar, touse)
    }
    else {
        mlscore = J(rows(X), 1, .)
    }
    n = rows(X)
    nstart = sum(is_start)
    
    printf("\nFrameStream with BuildHull\n")
    printf("------------------------------------------------------------\n")
    printf("Total: %g, Starters: %g, Streamers: %g\n", n, nstart, n-nstart)
    if (prune_interval > 0) {
        printf("Pruning interval: every %g DMUs\n", prune_interval)
    }
    if (gatevar != "" | mlscorevar != "") {
        printf("Frontier update gate enabled (epsilon=%9.2e)\n", gate_eps)
        if (mlscorevar != "") {
            printf("Candidate updates require %s >= %5.3f\n", mlscorevar, mlthreshold)
        }
    }
    printf("\nBuilding initial frame...\n")
    
    frame_idx = J(0, 1, .)
    all_alpha = J(n, 1, 0)
    all_drift = J(n, 1, 0)
    all_frame_size = J(n, 1, 0)
    all_frame_entry = J(n, 1, 0)
    all_frame_pruned = J(n, 1, 0)
    gate_status = J(n, 1, .)
    
    for (j=1; j<=nstart; j++) {
        if (rows(frame_idx) == 0) {
            Xf = X[j,]; Yf = Y[j,]; frame_idx = j
            all_frame_entry[j] = 1
        }
        else {
            a_val = _sdea_buildhull_eval(Xf, Yf, X[j,]', Y[j,]', solver)
            if (!missing(a_val) & a_val > 1e-8) {
                Xf = Xf \ X[j,]; Yf = Yf \ Y[j,]; frame_idx = frame_idx \ j
                all_frame_entry[j] = 1
            }
        }
        all_frame_size[j] = rows(frame_idx)
    }
    
    // Prune initial frame if prune_interval > 0
    total_pruned = 0
    if (prune_interval > 0 & rows(frame_idx) > 1) {
        n_pruned = _sdea_prune_frame(Xf, Yf, frame_idx, solver)
        total_pruned = total_pruned + n_pruned
        if (n_pruned > 0) {
            printf("Initial frame pruned: %g -> %g DMUs\n", rows(frame_idx)+n_pruned, rows(frame_idx))
            all_frame_pruned[nstart] = all_frame_pruned[nstart] + n_pruned
            all_frame_size[nstart] = rows(frame_idx)
        }
    }
    
    printf("Initial frame: %g DMUs\n", rows(frame_idx))
    printf("\nProcessing streamers...\n")
    
    // Starters have drift = 0
    for (j=1; j<=nstart; j++) {
        all_alpha[j] = 0
        all_drift[j] = 0
        gate_status[j] = 3
    }
    
    last_exterior = 0
    n_interior = 0
    n_exterior = 0
    n_rejected = 0
    prune_count = 0
    
    for (j=nstart+1; j<=n; j++) {
        a_val = _sdea_buildhull_eval(Xf, Yf, X[j,]', Y[j,]', solver)
        all_alpha[j] = a_val
        
        if (!missing(a_val)) {
            if (a_val <= gate_eps) {
                n_interior++
                all_drift[j] = 0
                gate_status[j] = 0
            }
            else {
                accept_update = 1
                if (mlscorevar != "") {
                    accept_update = (!missing(mlscore[j]) & mlscore[j] >= mlthreshold)
                }
                
                if (accept_update) {
                    n_exterior++
                    gate_status[j] = 1
                    if (last_exterior > gate_eps) {
                        all_drift[j] = abs(log(a_val + 1) - log(last_exterior + 1))
                    }
                    else {
                        all_drift[j] = abs(log(a_val + 1))
                    }
                    last_exterior = a_val
                    Xf = Xf \ X[j,]; Yf = Yf \ Y[j,]; frame_idx = frame_idx \ j
                    all_frame_entry[j] = 1
                    prune_count++
                    
                    // Prune frame periodically
                    if (prune_interval > 0 & prune_count >= prune_interval) {
                        n_pruned = _sdea_prune_frame(Xf, Yf, frame_idx, solver)
                        total_pruned = total_pruned + n_pruned
                        all_frame_pruned[j] = all_frame_pruned[j] + n_pruned
                        prune_count = 0
                    }
                }
                else {
                    n_rejected++
                    gate_status[j] = 2
                    all_drift[j] = 0
                }
            }
        }
        all_frame_size[j] = rows(frame_idx)
        
        if (mod(j-nstart, 100) == 0) {
            printf("  %g/%g (Frame: %g)\n", j-nstart, n-nstart, rows(frame_idx))
        }
    }
    
    // Final prune
    if (prune_interval > 0 & rows(frame_idx) > 1) {
        n_pruned = _sdea_prune_frame(Xf, Yf, frame_idx, solver)
        total_pruned = total_pruned + n_pruned
        all_frame_pruned[n] = all_frame_pruned[n] + n_pruned
        all_frame_size[n] = rows(frame_idx)
    }
    
    nframe = rows(frame_idx)
    printf("------------------------------------------------------------\n")
    printf("Completed: Frame=%g, Interior=%g, Exterior=%g", nframe, n_interior, n_exterior)
    if (gatevar != "" | mlscorevar != "") {
        printf(", Rejected=%g", n_rejected)
    }
    printf("\n")
    if (prune_interval > 0) {
        printf("Total pruned from frame: %g DMUs\n", total_pruned)
    }
    printf("------------------------------------------------------------\n")
    
    // Store results
    stata("capture drop " + alphavar)
    newvar = st_addvar("double", alphavar)
    st_store(., alphavar, touse, all_alpha)
    
    stata("capture drop " + framevar)
    newvar = st_addvar("byte", framevar)
    st_store(., framevar, touse, J(n, 1, 0))
    if (rows(frame_idx) > 0) {
        st_store(obs_id[frame_idx], framevar, J(rows(frame_idx), 1, 1))
    }
    
    stata("capture drop " + driftvar)
    newvar = st_addvar("double", driftvar)
    st_store(., driftvar, touse, all_drift)
    
    stata("capture drop _sdea_frame_size")
    newvar = st_addvar("long", "_sdea_frame_size")
    st_store(., "_sdea_frame_size", touse, all_frame_size)
    
    stata("capture drop _sdea_frame_entry")
    newvar = st_addvar("byte", "_sdea_frame_entry")
    st_store(., "_sdea_frame_entry", touse, all_frame_entry)
    
    stata("capture drop _sdea_frame_pruned")
    newvar = st_addvar("long", "_sdea_frame_pruned")
    st_store(., "_sdea_frame_pruned", touse, all_frame_pruned)
    
    if (gatevar != "") {
        stata("capture drop " + gatevar)
        newvar = st_addvar("byte", gatevar)
        st_store(., gatevar, touse, gate_status)
    }
}

//==============================================================================
// Standard DEA Evaluation (Output-oriented VRS)
//==============================================================================

real scalar _sdea_standard_eval(
    real matrix Xref,
    real matrix Yref,
    real colvector x0,
    real colvector y0
) {
    struct _sdea_lp_result scalar result
    real scalar p, m, s, i, r
    real matrix A
    real colvector b
    real rowvector c
    real matrix A_eq
    real colvector b_eq
    
    // Add self to reference set
    Xref = Xref \ x0'
    Yref = Yref \ y0'
    
    p = rows(Xref)
    m = cols(Xref)
    s = cols(Yref)
    
    // max φ
    c = (J(1, p, 0), 1)
    
    // Σ x_k λ_k ≤ x_0
    A = J(m, p+1, 0)
    b = J(m, 1, 0)
    for (i=1; i<=m; i++) {
        A[i, 1..p] = Xref[,i]'
        b[i] = x0[i]
    }
    
    // Σ y_k λ_k ≥ φ y_0 → -Σ y_k λ_k + φ y_0 ≤ 0
    A = A \ J(s, p+1, 0)
    b = b \ J(s, 1, 0)
    for (r=1; r<=s; r++) {
        A[m+r, 1..p] = -Yref[,r]'
        A[m+r, p+1] = y0[r]
    }
    
    // Equality: Σλ_k = 1
    A_eq = (J(1,p,1), 0)
    b_eq = 1
    
    // Solve LP
    result = _sdea_solve_lp(c, A, b, A_eq, b_eq, J(1,p+1,0), J(1,p+1,.), "max", "ipm")
    
    if (result.converged) {
        return(result.solution[p+1])
    }
    return(.)
}

//==============================================================================
// Run DEA for all DMUs
//==============================================================================

void _sdea_run_dea(string scalar inputs, string scalar outputs,
    string scalar efficiency, string scalar touse,
    string scalar orientation, string scalar rts) {
    
    struct _sdea_lp_result scalar result
    real matrix X, Y
    real colvector x0, y0, obs_indices
    real scalar n, m_in, m_out, k, obs_k, eff_score
    real matrix A
    real colvector b
    real rowvector c
    real matrix A_eq
    real colvector b_eq
    
    X = st_data(., inputs, touse)
    Y = st_data(., outputs, touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m_in = cols(X)
    m_out = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k,]'
        y0 = Y[k,]'
        
        if (orientation == "output") {
            // max φ: Xλ ≤ x0, Yλ ≥ φy0, Σλ=1 (VRS)
            c = (J(1, n, 0), -1)
            A = (X', J(m_in, 1, 0)) \ (-Y', y0)
            b = x0 \ J(m_out, 1, 0)
            
            if (rts == "vrs") {
                A_eq = (J(1,n,1), 0)
                b_eq = 1
            }
            else {
                A_eq = J(0, n+1, 0)
                b_eq = J(0, 1, 0)
            }
            
            result = _sdea_solve_lp(c, A, b, A_eq, b_eq, J(1,n+1,0), J(1,n+1,.), "min", "ipm")
            
            if (result.converged) {
                eff_score = -result.value
                st_store(obs_k, efficiency, eff_score)
            }
            else {
                st_store(obs_k, efficiency, .)
            }
        }
        else {
            // min θ: Xλ ≤ θx0, Yλ ≥ y0, Σλ=1 (VRS)
            c = (J(1, n, 0), 1)
            A = (X', -x0) \ (-Y', J(m_out, 1, 0))
            b = J(m_in, 1, 0) \ (-y0)
            
            if (rts == "vrs") {
                A_eq = (J(1,n,1), 0)
                b_eq = 1
            }
            else {
                A_eq = J(0, n+1, 0)
                b_eq = J(0, 1, 0)
            }
            
            result = _sdea_solve_lp(c, A, b, A_eq, b_eq, J(1,n+1,0), J(1,n+1,.), "min", "ipm")
            
            if (result.converged) {
                eff_score = result.value
                st_store(obs_k, efficiency, eff_score)
            }
            else {
                st_store(obs_k, efficiency, .)
            }
        }
    }
}

//==============================================================================
// Time-based Streaming DEA
//==============================================================================

void _sdea_time_streaming(string scalar inputs, string scalar outputs,
    string scalar touse, string scalar timevar,
    string scalar rts, string scalar orient,
    string scalar theta_var, string scalar eff_var,
    string scalar refsize_var, string scalar period_var,
    string scalar saveref, string scalar savehistory, string scalar verbose)
{
    real matrix DI, DO, REF_DI, REF_DO
    real colvector time_vals, unique_times
    real scalar i, t, n, n_times, n_ref, eff, prev_ref_size
    real scalar start_idx, end_idx
    real matrix history
    real colvector theta_vec, is_eff_vec, ref_size_vec, period_vec
    string scalar msg
    real matrix A
    real colvector b
    real rowvector c
    real scalar m_in, m_out, p
    
    // Load data
    DI = st_data(., tokens(inputs), touse)
    DO = st_data(., tokens(outputs), touse)
    time_vals = st_data(., timevar, touse)
    n = rows(DI)
    m_in = cols(DI)
    m_out = cols(DO)
    
    // Get unique time periods
    unique_times = uniqrows(time_vals)
    n_times = rows(unique_times)
    
    // Initialize results
    theta_vec = J(n, 1, .)
    is_eff_vec = J(n, 1, 0)
    ref_size_vec = J(n, 1, .)
    period_vec = J(n, 1, .)
    
    // Initialize reference set (transposed: inputs × DMUs)
    REF_DI = J(m_in, 0, .)
    REF_DO = J(m_out, 0, .)
    
    // Track history
    history = J(n_times, 5, .)
    prev_ref_size = 0
    
    // Process each time period
    for (t = 1; t <= n_times; t++) {
        // Find observations at this time
        start_idx = 0
        end_idx = 0
        for (i = 1; i <= n; i++) {
            if (time_vals[i] == unique_times[t]) {
                if (start_idx == 0) start_idx = i
                end_idx = i
            }
        }
        
        if (start_idx == 0) continue
        
        n_ref = cols(REF_DI)
        
        if (verbose != "") {
            msg = sprintf("Time %g: %g new obs, ref size = %g", 
                         unique_times[t], end_idx - start_idx + 1, n_ref)
            printf("{txt}%s\n", msg)
            displayflush()
        }
        
        // Evaluate each new observation
        for (i = start_idx; i <= end_idx; i++) {
            if (n_ref == 0) {
                // First observation is efficient
                eff = 1
                is_eff_vec[i] = 1
                REF_DI = DI[i, .]'
                REF_DO = DO[i, .]'
            }
            else {
                // Evaluate against current reference set
                struct _sdea_lp_result scalar lp_result
                real matrix X_all, Y_all, A_in, A_out, A_eq_t
                real colvector x0, y0, b_eq_t
                
                // Combine reference set with current DMU
                X_all = REF_DI' \ DI[i, .]
                Y_all = REF_DO' \ DO[i, .]
                p = rows(X_all)
                
                x0 = DI[i, .]'
                y0 = DO[i, .]'
                
                if (orient == "IN") {
                    // Input-oriented: min θ
                    c = (J(1, p, 0), 1)
                    A_in = (X_all', -x0)
                    A_out = (-Y_all', J(m_out, 1, 0))
                    A = A_in \ A_out
                    b = J(m_in, 1, 0) \ (-y0)
                    
                    if (rts == "VRS") {
                        A_eq_t = (J(1,p,1), 0)
                        b_eq_t = 1
                    }
                    else {
                        A_eq_t = J(0, p+1, 0)
                        b_eq_t = J(0, 1, 0)
                    }
                    
                    lp_result = _sdea_solve_lp(c, A, b, A_eq_t, b_eq_t, 
                                               J(1,p+1,0), J(1,p+1,.), "min", "ipm")
                    
                    if (lp_result.converged) {
                        eff = lp_result.value
                    }
                    else {
                        eff = .
                    }
                }
                else {
                    // Output-oriented: max φ
                    c = (J(1, p, 0), -1)
                    A_in = (X_all', J(m_in, 1, 0))
                    A_out = (-Y_all', y0)
                    A = A_in \ A_out
                    b = x0 \ J(m_out, 1, 0)
                    
                    if (rts == "VRS") {
                        A_eq_t = (J(1,p,1), 0)
                        b_eq_t = 1
                    }
                    else {
                        A_eq_t = J(0, p+1, 0)
                        b_eq_t = J(0, 1, 0)
                    }
                    
                    lp_result = _sdea_solve_lp(c, A, b, A_eq_t, b_eq_t,
                                               J(1,p+1,0), J(1,p+1,.), "min", "ipm")
                    
                    if (lp_result.converged) {
                        // φ = -value (≥1), convert to 0-1 scale: 1/φ
                        eff = 1 / (-lp_result.value)
                    }
                    else {
                        eff = .
                    }
                }
                
                // Check if efficient
                if (!missing(eff) & abs(eff - 1) < 1e-6) {
                    is_eff_vec[i] = 1
                    REF_DI = (REF_DI, DI[i, .]')
                    REF_DO = (REF_DO, DO[i, .]')
                }
            }
            
            theta_vec[i] = eff
            ref_size_vec[i] = cols(REF_DI)
            period_vec[i] = t
        }
        
        // Record history
        history[t, 1] = unique_times[t]
        history[t, 2] = end_idx - start_idx + 1
        history[t, 3] = sum(is_eff_vec[start_idx::end_idx])
        history[t, 4] = cols(REF_DI)
        if (prev_ref_size > 0) {
            history[t, 5] = (cols(REF_DI) - prev_ref_size) / prev_ref_size
        }
        prev_ref_size = cols(REF_DI)
    }
    
    // Store results
    st_store(., theta_var, touse, theta_vec)
    st_store(., eff_var, touse, is_eff_vec)
    st_store(., refsize_var, touse, ref_size_vec)
    st_store(., period_var, touse, period_vec)
    
    // Save reference set if requested
    if (saveref != "") {
        _sdea_save_reference(saveref, REF_DI', REF_DO', tokens(inputs), tokens(outputs))
    }
    
    // Save history if requested
    if (savehistory != "") {
        _sdea_save_history(savehistory, history)
    }
}

//==============================================================================
// Save Reference Set to File
//==============================================================================

void _sdea_save_reference(string scalar filename, real matrix DI, real matrix DO,
    string rowvector input_names, string rowvector output_names)
{
    real scalar fh, i, j, n_ref, n_inputs, n_outputs
    string scalar line
    
    n_ref = rows(DI)
    n_inputs = cols(DI)
    n_outputs = cols(DO)
    
    // Remove existing file if it exists
    unlink(filename)
    
    fh = fopen(filename, "w")
    
    // Header
    fput(fh, "Reference Set - " + strofreal(n_ref) + " efficient DMUs")
    
    line = "Inputs: "
    for (i = 1; i <= n_inputs; i++) {
        line = line + input_names[i] + " "
    }
    fput(fh, line)
    
    line = "Outputs: "
    for (i = 1; i <= n_outputs; i++) {
        line = line + output_names[i] + " "
    }
    fput(fh, line)
    fput(fh, "")
    
    // Data
    for (i = 1; i <= n_ref; i++) {
        line = ""
        for (j = 1; j <= n_inputs; j++) {
            line = line + strofreal(DI[i, j]) + " "
        }
        for (j = 1; j <= n_outputs; j++) {
            line = line + strofreal(DO[i, j]) + " "
        }
        fput(fh, line)
    }
    
    fclose(fh)
    printf("{txt}Reference set saved to: %s\n", filename)
}

//==============================================================================
// Save History to CSV
//==============================================================================

void _sdea_save_history(string scalar filename, real matrix history)
{
    real scalar fh, i, n_times
    string scalar line
    
    n_times = rows(history)
    
    // Remove existing file if it exists
    unlink(filename)
    
    fh = fopen(filename, "w")
    
    // Header
    fput(fh, "Time,N_New,N_Efficient,Ref_Size,Expansion_Rate")
    
    // Data
    for (i = 1; i <= n_times; i++) {
        line = strofreal(history[i,1]) + "," + 
               strofreal(history[i,2]) + "," + 
               strofreal(history[i,3]) + "," + 
               strofreal(history[i,4]) + "," + 
               strofreal(history[i,5])
        fput(fh, line)
    }
    
    fclose(fh)
    printf("{txt}History saved to: %s\n", filename)
}

//==============================================================================
// Malmquist Productivity Index - Adjacent Period
//==============================================================================
void _sdea_malmquist_adjacent(string scalar inputs,
                              string scalar outputs,
                              string scalar period_var,
                              string scalar dmu_var,
                              string scalar touse,
                              string scalar orientation,
                              string scalar rts,
                              string scalar verbose)
{
    real matrix X, Y
    real colvector period_vec, dmu_vec
    real scalar T, n, m, s, t
    real colvector unique_periods
    real matrix X_t, Y_t, X_t1, Y_t1
    real colvector idx_t, idx_t1
    real scalar n_dmu, i
    real scalar d_tt, d_t1t1, d_t1t, d_tt1, vrs_t, vrs_t1
    real scalar mi_tfpch, mi_effch, mi_techch, mi_pech, mi_sech
    real scalar period_t, period_t1, obs_idx
    
    // Load data
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    period_vec = st_data(., period_var, touse)
    dmu_vec = st_data(., dmu_var, touse)
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    unique_periods = uniqrows(period_vec)
    T = rows(unique_periods)
    
    // Create output variables - safe drop using _st_varindex (returns . if not found)
    if (_st_varindex("_tfpch") != .) st_dropvar("_tfpch")
    if (_st_varindex("_effch") != .) st_dropvar("_effch")
    if (_st_varindex("_techch") != .) st_dropvar("_techch")
    if (_st_varindex("_pech") != .) st_dropvar("_pech")
    if (_st_varindex("_sech") != .) st_dropvar("_sech")
    if (_st_varindex("_eff_t") != .) st_dropvar("_eff_t")
    if (_st_varindex("_eff_t1") != .) st_dropvar("_eff_t1")
    
    st_addvar("double", "_tfpch")
    st_addvar("double", "_effch")
    st_addvar("double", "_techch")
    st_addvar("double", "_pech")
    st_addvar("double", "_sech")
    st_addvar("double", "_eff_t")
    st_addvar("double", "_eff_t1")
    
    // Process each adjacent period pair
    for (t = 1; t <= T - 1; t++) {
        period_t = unique_periods[t]
        period_t1 = unique_periods[t + 1]
        
        if (verbose != "") {
            printf("Processing period %g -> %g\n", period_t, period_t1)
        }
        
        // Get data for period t and t+1
        idx_t = selectindex(period_vec :== period_t)
        idx_t1 = selectindex(period_vec :== period_t1)
        
        X_t = X[idx_t, .]
        Y_t = Y[idx_t, .]
        
        X_t1 = X[idx_t1, .]
        Y_t1 = Y[idx_t1, .]
        
        n_dmu = rows(X_t)
        
        // Calculate for each DMU
        for (i = 1; i <= n_dmu; i++) {
            // D^t(x^t, y^t)
            d_tt = _sdea_calc_distance(X_t, Y_t, X_t[i,.], Y_t[i,.], orientation, rts)
            
            // D^{t+1}(x^{t+1}, y^{t+1})
            d_t1t1 = _sdea_calc_distance(X_t1, Y_t1, X_t1[i,.], Y_t1[i,.], orientation, rts)
            
            // D^t(x^{t+1}, y^{t+1}) - cross period
            d_t1t = _sdea_calc_distance(X_t, Y_t, X_t1[i,.], Y_t1[i,.], orientation, rts)
            
            // D^{t+1}(x^t, y^t) - cross period
            d_tt1 = _sdea_calc_distance(X_t1, Y_t1, X_t[i,.], Y_t[i,.], orientation, rts)
            
            // VRS for PECH/SECH decomposition (only meaningful when base rts is CRS)
            if (rts == "crs") {
                vrs_t = _sdea_calc_distance(X_t, Y_t, X_t[i,.], Y_t[i,.], orientation, "vrs")
                vrs_t1 = _sdea_calc_distance(X_t1, Y_t1, X_t1[i,.], Y_t1[i,.], orientation, "vrs")
            }
            else {
                // When rts is VRS, PECH = EFFCH, SECH = 1
                vrs_t = d_tt
                vrs_t1 = d_t1t1
            }
            
            // Calculate indices
            if (d_tt > 0 & d_t1t1 > 0 & d_t1t > 0 & d_tt1 > 0) {
                mi_tfpch = sqrt((d_t1t / d_tt) * (d_t1t1 / d_tt1))
                mi_effch = d_t1t1 / d_tt
                mi_techch = sqrt((d_t1t / d_t1t1) * (d_tt / d_tt1))
                
                if (vrs_t > 0 & vrs_t1 > 0) {
                    mi_pech = vrs_t1 / vrs_t
                    mi_sech = mi_effch / mi_pech
                }
                else {
                    mi_pech = .
                    mi_sech = .
                }
            }
            else {
                mi_tfpch = .
                mi_effch = .
                mi_techch = .
                mi_pech = .
                mi_sech = .
            }
            
            // Store results
            obs_idx = idx_t1[i]
            st_store(obs_idx, st_varindex("_tfpch"), mi_tfpch)
            st_store(obs_idx, st_varindex("_effch"), mi_effch)
            st_store(obs_idx, st_varindex("_techch"), mi_techch)
            st_store(obs_idx, st_varindex("_pech"), mi_pech)
            st_store(obs_idx, st_varindex("_sech"), mi_sech)
            st_store(obs_idx, st_varindex("_eff_t"), d_tt)
            st_store(obs_idx, st_varindex("_eff_t1"), d_t1t1)
        }
    }
    
    stata(`"label var _tfpch "Total Factor Productivity Change""')
    stata(`"label var _effch "Efficiency Change (catch-up)""')
    stata(`"label var _techch "Technical Change (frontier shift)""')
    stata(`"label var _pech "Pure Efficiency Change""')
    stata(`"label var _sech "Scale Efficiency Change""')
    stata(`"label var _eff_t "Efficiency at period t""')
    stata(`"label var _eff_t1 "Efficiency at period t+1""')
}

//==============================================================================
// Malmquist Productivity Index - Global (Pooled Frontier)
//==============================================================================
void _sdea_malmquist_global(string scalar inputs,
                            string scalar outputs,
                            string scalar period_var,
                            string scalar dmu_var,
                            string scalar touse,
                            string scalar orientation,
                            string scalar rts,
                            string scalar verbose)
{
    real matrix X, Y, X_all, Y_all
    real colvector period_vec, dmu_vec
    real scalar T, n, t
    real colvector unique_periods
    real matrix X_t, Y_t, X_t1, Y_t1
    real colvector idx_t, idx_t1
    real scalar n_dmu, i
    real scalar d_t_global, d_t1_global, d_t_t, d_t1_t1
    real scalar mi_tfpch, mi_effch, mi_techch
    real scalar period_t, period_t1, obs_idx
    
    // Load data
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    period_vec = st_data(., period_var, touse)
    dmu_vec = st_data(., dmu_var, touse)
    
    // Global frontier (all data pooled)
    X_all = X
    Y_all = Y
    
    n = rows(X)
    unique_periods = uniqrows(period_vec)
    T = rows(unique_periods)
    
    // Create output variables - safe drop using _st_varindex (returns . if not found)
    if (_st_varindex("_tfpch") != .) st_dropvar("_tfpch")
    if (_st_varindex("_effch") != .) st_dropvar("_effch")
    if (_st_varindex("_techch") != .) st_dropvar("_techch")
    if (_st_varindex("_eff_t") != .) st_dropvar("_eff_t")
    if (_st_varindex("_eff_t1") != .) st_dropvar("_eff_t1")
    if (_st_varindex("_eff_global") != .) st_dropvar("_eff_global")
    
    st_addvar("double", "_tfpch")
    st_addvar("double", "_effch")
    st_addvar("double", "_techch")
    st_addvar("double", "_eff_t")
    st_addvar("double", "_eff_t1")
    st_addvar("double", "_eff_global")
    
    for (t = 1; t <= T - 1; t++) {
        period_t = unique_periods[t]
        period_t1 = unique_periods[t + 1]
        
        if (verbose != "") {
            printf("Processing period %g -> %g (Global)\n", period_t, period_t1)
        }
        
        idx_t = selectindex(period_vec :== period_t)
        idx_t1 = selectindex(period_vec :== period_t1)
        
        X_t = X[idx_t, .]
        Y_t = Y[idx_t, .]
        X_t1 = X[idx_t1, .]
        Y_t1 = Y[idx_t1, .]
        
        n_dmu = rows(X_t)
        
        for (i = 1; i <= n_dmu; i++) {
            // Distance to global frontier
            d_t_global = _sdea_calc_distance(X_all, Y_all, X_t[i,.], Y_t[i,.], orientation, rts)
            d_t1_global = _sdea_calc_distance(X_all, Y_all, X_t1[i,.], Y_t1[i,.], orientation, rts)
            
            // Distance to own-period frontier
            d_t_t = _sdea_calc_distance(X_t, Y_t, X_t[i,.], Y_t[i,.], orientation, rts)
            d_t1_t1 = _sdea_calc_distance(X_t1, Y_t1, X_t1[i,.], Y_t1[i,.], orientation, rts)
            
            if (d_t_global > 0 & d_t1_global > 0 & d_t_t > 0 & d_t1_t1 > 0) {
                mi_tfpch = d_t1_global / d_t_global
                mi_effch = d_t1_t1 / d_t_t
                mi_techch = mi_tfpch / mi_effch
            }
            else {
                mi_tfpch = .
                mi_effch = .
                mi_techch = .
            }
            
            obs_idx = idx_t1[i]
            st_store(obs_idx, st_varindex("_tfpch"), mi_tfpch)
            st_store(obs_idx, st_varindex("_effch"), mi_effch)
            st_store(obs_idx, st_varindex("_techch"), mi_techch)
            st_store(obs_idx, st_varindex("_eff_t"), d_t_t)
            st_store(obs_idx, st_varindex("_eff_t1"), d_t1_t1)
            st_store(obs_idx, st_varindex("_eff_global"), d_t1_global)
        }
    }
    
    stata(`"label var _tfpch "Global Total Factor Productivity Change""')
    stata(`"label var _effch "Efficiency Change""')
    stata(`"label var _techch "Technical Change (Best Practice Gap)""')
    stata(`"label var _eff_global "Efficiency relative to global frontier""')
}

//==============================================================================
// Calculate Distance Function using LinearProgram()
//==============================================================================
real scalar _sdea_calc_distance(real matrix X_ref,
                                real matrix Y_ref,
                                real rowvector x0,
                                real rowvector y0,
                                string scalar orientation,
                                string scalar rts)
{
    real scalar n, m, s, eff
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub
    class LinearProgram scalar lp
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (orientation == "input") {
        // min θ: X'λ <= θx0, Y'λ >= y0
        c = (J(1, n, 0), 1)
        
        A_ineq = (X_ref', -x0')
        b_ineq = J(m, 1, 0)
        
        A_ineq = A_ineq \ (-Y_ref', J(s, 1, 0))
        b_ineq = b_ineq \ (-y0')
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
        }
        else {
            A_eq = J(0, n + 1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lb = J(1, n + 1, 0)
        ub = J(1, n + 1, .)
        
        lp = LinearProgram()
        lp.setCoefficients(c)
        lp.setMaxOrMin("min")
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = lp.value()
            if (eff < 0) eff = 0
        }
        else {
            eff = .
        }
    }
    else {
        // max φ: X'λ <= x0, Y'λ >= φy0
        c = (J(1, n, 0), -1)
        
        A_ineq = (X_ref', J(m, 1, 0))
        b_ineq = x0'
        
        A_ineq = A_ineq \ (-Y_ref', y0')
        b_ineq = b_ineq \ J(s, 1, 0)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
        }
        else {
            A_eq = J(0, n + 1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lb = J(1, n + 1, 0)
        ub = J(1, n + 1, .)
        
        lp = LinearProgram()
        lp.setCoefficients(c)
        lp.setMaxOrMin("min")
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = -lp.value()
            if (eff > 0) eff = 1 / eff
        }
        else {
            eff = .
        }
    }
    
    return(eff)
}

//==============================================================================
// TFDEA: Technology Forecasting using DEA
// With Time-Decay Weighted (TDW) option for aggregate ROC
//==============================================================================

void _sdea_tfdea(string scalar inputs_str, string scalar outputs_str,
                  string scalar time_str, string scalar dmu_str,
                  string scalar touse_str, real scalar tf,
                  string scalar orientation, string scalar rts,
                  string scalar model, real scalar verbose,
                  real scalar decay_delta, string scalar decay_type)
{
    real matrix X, Y, X_t, Y_t
    real colvector time_var, dmu_var, touse
    real colvector unique_times, theta_all
    real colvector theta_tk, theta_tf, tk_eff, roc_vec, t_eff_roc_vec
    real scalar n, i, t, n_times, theta, theta_init
    real scalar tk_val, elapsed
    string rowvector input_vars, output_vars
    
    // Get variable names
    input_vars = tokens(inputs_str)
    output_vars = tokens(outputs_str)
    
    // Load data
    time_var = st_data(., time_str)
    dmu_var = st_data(., dmu_str)
    touse = st_data(., touse_str)
    X = st_data(., input_vars)
    Y = st_data(., output_vars)
    
    n = rows(X)
    
    // Get unique times up to tf
    unique_times = uniqrows(select(time_var, (time_var :<= tf) :& touse))
    n_times = rows(unique_times)
    
    printf("{txt}Step 1: Computing Rate of Change (ROC)...\n")
    printf("{txt}  Time periods for ROC: %g (up to tf=%g)\n", n_times, tf)
    if (decay_delta > 0) {
        printf("{txt}  TDW enabled: δ=%g, type=%s\n", decay_delta, decay_type)
    }
    
    // Initialize result vectors
    theta_tk = J(n, 1, .)
    theta_tf = J(n, 1, .)
    tk_eff = J(n, 1, .)
    roc_vec = J(n, 1, .)
    t_eff_roc_vec = J(n, 1, .)  // Store t_eff for each ROC observation
    
    // Declare loop variables for ROC calculation
    real colvector time_t
    real rowvector lambda_roc
    real scalar t_eff_roc, was_efficient_tk, is_inefficient_tf, lambda_sum_roc
    
    // For each DMU, compute efficiency at each time period
    // and track when it became inefficient
    
    for (i = 1; i <= n; i++) {
        if (touse[i] == 0 || time_var[i] > tf) continue
        
        tk_val = time_var[i]  // When this DMU first appeared
        
        // Compute efficiency at initial time (when DMU appeared)
        // Reference set: all DMUs at or before tk_val
        X_t = select(X, (time_var :<= tk_val) :& touse)
        Y_t = select(Y, (time_var :<= tk_val) :& touse)
        
        theta_init = _sdea_solve_dea_basic(X_t, Y_t, X[i,.], Y[i,.], orientation, rts)
        if (theta_init != .) {
            theta_tk[i] = theta_init
        }
        
        // Compute efficiency at final time (tf) with lambda extraction for t_eff
        // Reference set: all DMUs at or before tf
        X_t = select(X, (time_var :<= tf) :& touse)
        Y_t = select(Y, (time_var :<= tf) :& touse)
        time_t = select(time_var, (time_var :<= tf) :& touse)
        
        // Get efficiency and lambda vector
        theta = _sdea_solve_dea_with_lambda(X_t, Y_t, X[i,.], Y[i,.], 
                                             orientation, rts, lambda_roc)
        if (theta != .) {
            theta_tf[i] = theta
        }
        
        // Compute ROC if DMU was efficient initially but not at final time
        if (theta_tk[i] != . && theta_tf[i] != .) {
            // Check if efficient at tk (theta = 1 for both orientations)
            was_efficient_tk = (abs(theta_tk[i] - 1) < 1e-6)
            
            // Check if inefficient at tf (depends on orientation)
            // Input-oriented: theta < 1 means inefficient (can reduce inputs)
            // Output-oriented: theta > 1 means inefficient (need to expand outputs)
            if (orientation == "input") {
                is_inefficient_tf = (theta_tf[i] < 1 - 1e-6)
            }
            else {
                // Output-oriented: theta > 1 means inefficient
                is_inefficient_tf = (theta_tf[i] > 1 + 1e-6)
            }
            
            if (was_efficient_tk && is_inefficient_tf) {
                // DMU was SOA initially, now inefficient
                // Compute t_eff: weighted average of reference times (Eq. 9)
                lambda_sum_roc = sum(lambda_roc)
                if (lambda_sum_roc > 1e-10) {
                    t_eff_roc = sum(lambda_roc :* time_t') / lambda_sum_roc
                }
                else {
                    t_eff_roc = tf
                }
                
                // ROC formula using t_eff (Eq. 10): γ = φ^(1/(t_eff - tk))
                elapsed = t_eff_roc - tk_val
                if (elapsed > 0 && theta_tf[i] > 0) {
                    roc_vec[i] = (theta_tf[i])^(1/elapsed)
                    tk_eff[i] = tk_val
                    t_eff_roc_vec[i] = t_eff_roc  // Store t_eff for TDW weighting
                }
            }
        }
        
        if (verbose && mod(i, 100) == 0) {
            printf("{txt}  Processed %g/%g observations\n", i, n)
        }
    }
    
    // Store ROC results
    st_store(., "_tfdea_theta_tk", theta_tk)
    st_store(., "_tfdea_theta_tf", theta_tf)
    st_store(., "_tfdea_tk_eff", tk_eff)
    st_store(., "_tfdea_roc", roc_vec)
    st_store(., "_tfdea_t_eff_roc", t_eff_roc_vec)
    
    // Compute mean ROC (Standard and TDW)
    real scalar mean_roc, mean_roc_tdw, n_roc
    real colvector valid_roc, valid_t_eff
    real colvector valid_idx
    
    valid_idx = selectindex((roc_vec :!= .) :& (roc_vec :> 0) :& (roc_vec :< 2))
    n_roc = rows(valid_idx)
    
    if (n_roc > 0) {
        valid_roc = roc_vec[valid_idx]
        valid_t_eff = t_eff_roc_vec[valid_idx]
        
        // Standard: Simple mean
        mean_roc = mean(valid_roc)
        
        // TDW: Time-decay weighted mean (if decay_delta > 0)
        if (decay_delta > 0) {
            real colvector tau, weights
            real scalar sum_weights
            
            // tau = time distance from frontier year
            // Newer observations (higher t_eff, closer to tf) get higher weight
            tau = tf :- valid_t_eff
            
            // Compute decay weights based on decay_type
            if (decay_type == "exponential") {
                weights = exp(-decay_delta :* tau)
            }
            else if (decay_type == "hyperbolic") {
                weights = 1 :/ (1 :+ decay_delta :* tau)
            }
            else if (decay_type == "linear") {
                weights = 1 :- decay_delta :* tau
                weights = weights :* (weights :> 0.01) + 0.01 :* (weights :<= 0.01)  // element-wise max with 0.01
            }
            else if (decay_type == "gaussian") {
                weights = exp(-decay_delta :* (tau:^2))
            }
            else if (decay_type == "logistic") {
                weights = 2 :/ (1 :+ exp(decay_delta :* tau))
            }
            else {
                weights = J(n_roc, 1, 1)
            }
            
            sum_weights = sum(weights)
            if (sum_weights > 1e-10) {
                mean_roc_tdw = sum(weights :* valid_roc) / sum_weights
            }
            else {
                mean_roc_tdw = mean_roc
            }
            
            printf("{txt}  Mean ROC (Standard): %9.6f (from %g DMUs)\n", mean_roc, n_roc)
            printf("{txt}  Mean ROC (TDW δ=%g): %9.6f\n", decay_delta, mean_roc_tdw)
        }
        else {
            mean_roc_tdw = mean_roc
            printf("{txt}  Mean ROC: %9.6f (from %g DMUs)\n", mean_roc, n_roc)
        }
        
        // Store mean ROC values in all observations for reference
        st_store(., "_tfdea_mean_roc_std", J(n, 1, mean_roc))
        st_store(., "_tfdea_mean_roc_tdw", J(n, 1, mean_roc_tdw))
    }
    else {
        printf("{err}  Warning: No valid ROC values computed\n")
        return
    }
    
    // Step 2: Forecasting for DMUs with time > tf
    real colvector se_vec, tf_eff_vec, tfexp_vec
    real scalar n_forecast, j, theta_se
    real matrix X_ref, Y_ref
    real colvector idx_forecast
    real colvector time_ref
    real scalar weighted_time, lambda_sum
    real matrix X_all, Y_all
    string scalar se_model
    
    idx_forecast = selectindex((time_var :> tf) :& touse)
    n_forecast = rows(idx_forecast)
    
    if (n_forecast > 0) {
        printf("{txt}\nStep 2: Computing Technology Forecasts...\n")
        printf("{txt}  Observations to forecast: %g\n", n_forecast)
        printf("{txt}  Orientation: %s, RTS: %s\n", orientation, rts)
        
        // Determine super-efficiency model
        if (model == "supersbm" || model == "sbm") {
            se_model = "supersbm"
            printf("{txt}  Super-efficiency model: Super-SBM (Tone 2002)\n")
        }
        else {
            se_model = "radial"
            printf("{txt}  Super-efficiency model: Radial\n")
        }
        
        se_vec = J(n, 1, .)
        tf_eff_vec = J(n, 1, .)
        tfexp_vec = J(n, 1, .)
        
        // Reference set: all DMUs at or before tf
        X_ref = select(X, (time_var :<= tf) :& touse)
        Y_ref = select(Y, (time_var :<= tf) :& touse)
        time_ref = select(time_var, (time_var :<= tf) :& touse)
        
        // Declare loop variables
        real rowvector lambda_se
        real scalar t_eff
        
        for (j = 1; j <= n_forecast; j++) {
            i = idx_forecast[j]
            
            if (se_model == "supersbm") {
                // Use Super-SBM via dea_lp (call Stata command)
                // Create temporary dataset with reference + new DMU
                theta_se = _sdea_tfdea_supersbm(X_ref, Y_ref, X[i,.], Y[i,.], 
                                                 inputs_str, outputs_str, rts)
                t_eff = tf  // For Super-SBM, use tf as approximation
            }
            else {
                // Use radial super-efficiency with lambda extraction
                // Create combined dataset: reference + new DMU
                X_all = X_ref \ X[i, .]
                Y_all = Y_ref \ Y[i, .]
                
                // Compute super-efficiency and get lambda vector
                theta_se = _sdea_solve_supereff_tfdea(X_all, Y_all, X[i,.], Y[i,.],
                                                       rows(X_all), orientation, rts, lambda_se)
                
                // Calculate t_eff: weighted average of reference times
                // Formula (9) from Inman et al. (2006): t_eff = Σ(λ_j * t_j) / Σλ_j
                lambda_sum = sum(lambda_se)
                if (lambda_sum > 1e-10) {
                    t_eff = sum(lambda_se :* time_ref') / lambda_sum
                }
                else {
                    t_eff = tf  // Fallback if no reference DMUs used
                }
            }
            
            if (theta_se != . && theta_se > 0) {
                se_vec[i] = theta_se
                tf_eff_vec[i] = t_eff  // Store t_eff, not tf!
                
                // Compute forecasted release date
                // Following Inman et al. (2006) Equation (13):
                // t_expected = t_eff + ln(1/φ) / ln(γ̄)
                //
                // Key insight for radial super-efficiency:
                // - Output-oriented VRS: φ < 1 means advanced (exceeds frontier)
                // - Input-oriented: θ < 1 means advanced
                
                real scalar log_roc, log_distance, roc_for_forecast
                
                // Use TDW ROC if decay is enabled, otherwise standard ROC
                if (decay_delta > 0) {
                    roc_for_forecast = mean_roc_tdw
                }
                else {
                    roc_for_forecast = mean_roc
                }
                log_roc = ln(roc_for_forecast)
                
                if (abs(log_roc) > 1e-6) {  // ROC != 1
                    if (se_model == "supersbm") {
                        // Tone (2002) Super-SBM: delta > 1 means more advanced
                        if (theta_se > 1) {
                            log_distance = ln(theta_se)
                            tfexp_vec[i] = t_eff + log_distance / abs(log_roc)
                        }
                        else {
                            tfexp_vec[i] = t_eff
                        }
                    }
                    else if (orientation == "output") {
                        // Output-oriented radial: φ < 1 means more advanced
                        // Formula: t_exp = t_eff + ln(1/φ) / ln(γ̄)
                        if (theta_se < 1) {
                            log_distance = ln(1/theta_se)
                            tfexp_vec[i] = t_eff + log_distance / log_roc
                        }
                        else {
                            tfexp_vec[i] = t_eff
                        }
                    }
                    else {
                        // Input-oriented radial: θ > 1 means more advanced
                        if (theta_se > 1) {
                            log_distance = ln(theta_se)
                            tfexp_vec[i] = t_eff + log_distance / abs(log_roc)
                        }
                        else {
                            tfexp_vec[i] = t_eff
                        }
                    }
                }
                
                printf("{txt}    DMU %g: SE = %9.4f, t_eff = %7.2f, tf_exp = %9.2f\n", 
                       i, se_vec[i], t_eff, tfexp_vec[i])
            }
            else {
                // Infeasible super-efficiency: technology dominates frontier
                // This can happen when the DMU is so far beyond the frontier
                // that no reference combination can match it
                printf("{txt}    DMU %g: Super-efficiency infeasible (extreme frontier exceedance)\n", i)
                
                se_vec[i] = .
                tf_eff_vec[i] = tf
                tfexp_vec[i] = .
            }
        }
        
        st_store(., "_tfdea_se", se_vec)
        st_store(., "_tfdea_tf_eff", tf_eff_vec)
        st_store(., "_tfdea_tfexp", tfexp_vec)
    }
    else {
        printf("{txt}\nNote: No observations with time > tf, skipping forecasting\n")
    }
}

// Super-SBM computation for TFDEA — Tone (2002) formulation
// Internal implementation via Charnes-Cooper linearization
real scalar _sdea_tfdea_supersbm(real matrix X_ref, real matrix Y_ref,
                                  real rowvector x0, real rowvector y0,
                                  string scalar inputs_str, string scalar outputs_str,
                                  string scalar rts)
{
    // Call internal Tone (2002) Super-SBM solver with k = n+1 (new DMU is last)
    real matrix X_all, Y_all
    real scalar n_ref, se_result
    
    n_ref = rows(X_ref)
    X_all = X_ref \ x0
    Y_all = Y_ref \ y0
    
    // Tone (2002) Super-SBM: delta >= 1 for super-efficient
    se_result = _sdea_solve_supersbm_internal(X_all, Y_all, x0, y0, n_ref + 1, rts)
    
    return(se_result)
}

//==============================================================================
// SUPER-SBM Internal Solver — Tone (2002, EJOR 143:32-41)
// Super-efficiency SBM model (excludes DMU k from reference set)
//
// Charnes-Cooper linearization of the ratio form:
//   delta* = min [(1/m) SUM_i (x_bar_i / x_i0)] / [(1/s) SUM_r (y_bar_r / y_r0)]
//   s.t. x_bar >= X*lambda,  y_bar <= Y*lambda
//        x_bar >= x0,  y_bar <= y0
//        SUM lambda = 1 (VRS), lambda >= 0
//
// Let t = 1 / [(1/s) SUM_r (Y_bar_r / y_r0)] > 0
// Substitutions: Lambda = t*lambda, X_bar = t*x_bar, Y_bar = t*y_bar
//
// LP: min SUM_i [1/(m*x_i0)] * X_bar_i
//     s.t. SUM_r [1/(s*y_r0)] * Y_bar_r = 1       (normalization)
//          X*Lambda - X_bar <= 0                     (X_bar >= X*Lambda)
//          -Y*Lambda + Y_bar <= 0                    (Y_bar <= Y*Lambda)
//          -X_bar + t*x0 <= 0                        (X_bar >= t*x0)
//          Y_bar - t*y0 <= 0                         (Y_bar <= t*y0)
//          SUM Lambda - t = 0                         (VRS: SUM lambda = 1)
//          Lambda >= 0, t > 0, X_bar >= 0, Y_bar >= 0
//
// Variables: [Lambda(n_ref), X_bar(m), Y_bar(s), t(1)]
// Result delta* >= 1 for super-efficient, = 1 for frontier, < 1 NOT feasible
// (should only be applied to SBM-efficient DMUs)
//
// For infeasible cases (very extreme DMUs), returns missing value.
//==============================================================================
real scalar _sdea_solve_supersbm_internal(real matrix X_all, real matrix Y_all,
                                          real rowvector x0, real rowvector y0,
                                          real scalar k, string scalar rts)
{
    real scalar n, m, s, i, r, delta
    real matrix X_ref, Y_ref
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub
    real scalar n_ref, total_vars
    class LinearProgram scalar lp
    
    n = rows(X_all)
    m = cols(X_all)
    s = cols(Y_all)
    
    if (n <= 1) return(.)
    
    // Exclude DMU k from reference set
    if (k == 1) {
        X_ref = X_all[2::n, .]
        Y_ref = Y_all[2::n, .]
    }
    else if (k == n) {
        X_ref = X_all[1::(n-1), .]
        Y_ref = Y_all[1::(n-1), .]
    }
    else {
        X_ref = X_all[1::(k-1), .] \ X_all[(k+1)::n, .]
        Y_ref = Y_all[1::(k-1), .] \ Y_all[(k+1)::n, .]
    }
    
    n_ref = rows(X_ref)
    if (n_ref == 0) return(.)
    
    // Check for zero values (avoid division by zero)
    for (i = 1; i <= m; i++) {
        if (x0[i] <= 1e-10) return(.)
    }
    for (r = 1; r <= s; r++) {
        if (y0[r] <= 1e-10) return(.)
    }
    
    // Variables: [Lambda(n_ref), X_bar(m), Y_bar(s), t(1)]
    total_vars = n_ref + m + s + 1
    
    // Objective: min SUM_i [1/(m*x_i0)] * X_bar_i
    c = J(1, total_vars, 0)
    for (i = 1; i <= m; i++) {
        c[n_ref + i] = 1 / (m * x0[i])
    }
    
    // --- Inequality constraints ---
    // (1) X_ref' * Lambda - I_m * X_bar <= 0   (m rows)
    A_ineq = (X_ref', -I(m), J(m, s, 0), J(m, 1, 0))
    b_ineq = J(m, 1, 0)
    
    // (2) -Y_ref' * Lambda + I_s * Y_bar <= 0  (s rows)
    A_ineq = A_ineq \ (-Y_ref', J(s, m, 0), I(s), J(s, 1, 0))
    b_ineq = b_ineq \ J(s, 1, 0)
    
    // (3) -X_bar + t*x0 <= 0  i.e. X_bar >= t*x0  (m rows)
    real matrix A3
    A3 = J(m, n_ref, 0)
    A3 = (A3, -I(m), J(m, s, 0), x0')
    A_ineq = A_ineq \ A3
    b_ineq = b_ineq \ J(m, 1, 0)
    
    // (4) Y_bar - t*y0 <= 0  i.e. Y_bar <= t*y0  (s rows)
    real matrix A4
    A4 = J(s, n_ref, 0)
    A4 = (A4, J(s, m, 0), I(s), -y0')
    A_ineq = A_ineq \ A4
    b_ineq = b_ineq \ J(s, 1, 0)
    
    // --- Equality constraints ---
    // (1) Normalization: SUM_r [1/(s*y_r0)] * Y_bar_r = 1
    real rowvector eq_norm
    eq_norm = J(1, total_vars, 0)
    for (r = 1; r <= s; r++) {
        eq_norm[n_ref + m + r] = 1 / (s * y0[r])
    }
    
    if (rts == "vrs") {
        // (2) SUM Lambda = t  ->  SUM Lambda - t = 0
        real rowvector eq_vrs
        eq_vrs = (J(1, n_ref, 1), J(1, m, 0), J(1, s, 0), -1)
        A_eq = eq_norm \ eq_vrs
        b_eq = (1 \ 0)
    }
    else {
        // CRS: no convexity constraint on lambda (but t still present)
        // For CRS, lambda is free-scaling: no SUM constraint
        A_eq = eq_norm
        b_eq = 1
    }
    
    // --- Bounds ---
    // Lambda >= 0, X_bar >= 0, Y_bar >= 0, t >= epsilon
    lb = (J(1, n_ref, 0), J(1, m, 0), J(1, s, 0), 1e-8)
    ub = J(1, total_vars, .)
    
    // --- Solve ---
    lp = LinearProgram()
    lp.setCoefficients(c)
    lp.setMaxOrMin("min")
    lp.setInequality(A_ineq, b_ineq)
    lp.setEquality(A_eq, b_eq)
    lp.setBounds(lb, ub)
    lp.setTrace("off")
    
    (void) lp.optimize()
    
    if (lp.converged() == 1 & lp.errorcode() == 0) {
        delta = lp.value()
        // delta >= 1 for super-efficient DMUs (Tone 2002)
    }
    else {
        // Infeasible: typically means DMU is extremely super-efficient
        // or the problem is degenerate
        delta = .
    }
    
    return(delta)
}

// Basic DEA solver for TFDEA - WITH DUAL FALLBACK
// Tries primal (envelopment) first, falls back to dual (multiplier) if primal fails
// Strong duality ensures same efficiency value from both formulations
real scalar _sdea_solve_dea_basic(real matrix X_ref, real matrix Y_ref,
                                   real rowvector x0, real rowvector y0,
                                   string scalar orientation, string scalar rts)
{
    real scalar eff
    
    // Try primal first
    eff = _sdea_solve_dea_primal(X_ref, Y_ref, x0, y0, orientation, rts)
    
    // If primal failed, try dual
    if (eff == .) {
        eff = _sdea_solve_dea_dual(X_ref, Y_ref, x0, y0, orientation, rts)
    }
    
    return(eff)
}

// ============================================================================
// PRIMAL (ENVELOPMENT) MODEL
// Input-oriented: min θ s.t. Xλ ≤ θx₀, Yλ ≥ y₀, [Σλ=1], λ ≥ 0
// Output-oriented: max φ s.t. Xλ ≤ x₀, Yλ ≥ φy₀, [Σλ=1], λ ≥ 0
// ============================================================================
real scalar _sdea_solve_dea_primal(real matrix X_ref, real matrix Y_ref,
                                    real rowvector x0, real rowvector y0,
                                    string scalar orientation, string scalar rts)
{
    real scalar n, m, s, eff
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub
    class LinearProgram scalar lp
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (n == 0) return(.)
    
    // Special case: n=1 with VRS - compare directly
    if (n == 1 & rts == "vrs") {
        if (orientation == "input") {
            if (all(x0 :>= X_ref) & all(Y_ref :>= y0)) {
                eff = min(X_ref :/ x0)
                return(eff)
            }
            else if (all(X_ref :<= x0) & all(Y_ref :>= y0)) {
                return(1)
            }
            else {
                eff = min(X_ref :/ x0)
                if (eff <= 1) return(eff)
                return(.)
            }
        }
        else {
            if (all(X_ref :<= x0) & all(y0 :<= Y_ref)) {
                eff = max(y0 :/ Y_ref)
                if (eff > 0) eff = 1 / eff
                return(eff)
            }
            else if (all(X_ref :<= x0)) {
                eff = max(y0 :/ Y_ref)
                if (eff > 0) return(1/eff)
                return(.)
            }
            else {
                return(.)
            }
        }
    }
    
    lp = LinearProgram()
    lp.setTrace("off")
    
    if (orientation == "input") {
        // min θ: θx0 >= X'λ, Y'λ >= y0
        c = (1, J(1, n, 0))
        
        // -θx0 + X'λ <= 0
        A_ineq = (-x0', X_ref')
        b_ineq = J(m, 1, 0)
        
        // -Y'λ <= -y0
        A_ineq = A_ineq \ (J(s, 1, 0), -Y_ref')
        b_ineq = b_ineq \ (-y0')
        
        lb = (0, J(1, n, 0))
        ub = (., J(1, n, .))
        
        if (rts == "vrs") {
            A_eq = (0, J(1, n, 1))
            b_eq = 1
        }
        else {
            A_eq = J(0, n+1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lp.setCoefficients(c)
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        lp.setMaxOrMin("min")
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = lp.value()
        }
        else {
            eff = .
        }
    }
    else {
        // Output orientation: max φ s.t. X'λ <= x0, Y'λ >= φy0, λ >= 0
        c = (J(1, n, 0), 1)
        
        A_ineq = (X_ref', J(m, 1, 0))
        b_ineq = x0'
        
        A_ineq = A_ineq \ (-Y_ref', y0')
        b_ineq = b_ineq \ J(s, 1, 0)
        
        lb = (J(1, n, 0), 1)
        ub = (J(1, n, .), .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
        }
        else {
            A_eq = J(0, n+1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lp.setCoefficients(c)
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        lp.setMaxOrMin("max")
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = lp.value()
        }
        else {
            eff = .
        }
    }
    
    return(eff)
}

// DEA solver with lambda vector extraction (for ROC t_eff calculation)
// Returns efficiency score and populates lambda vector
real scalar _sdea_solve_dea_with_lambda(real matrix X_ref, real matrix Y_ref,
                                         real rowvector x0, real rowvector y0,
                                         string scalar orientation, string scalar rts,
                                         real rowvector lambda)
{
    real scalar n, m, s, eff
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub, params
    class LinearProgram scalar lp
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    lambda = J(1, n, 0)  // Initialize lambda vector
    
    if (n == 0) return(.)
    
    lp = LinearProgram()
    lp.setTrace("off")
    
    if (orientation == "input") {
        // Input-oriented: min θ s.t. θx₀ ≥ Xλ, Yλ ≥ y₀, λ ≥ 0
        // Variables: [θ, λ₁, ..., λₙ]
        c = (1, J(1, n, 0))
        
        // -θx₀ + Xλ ≤ 0 (θx₀ ≥ Xλ)
        A_ineq = (-x0', X_ref')
        b_ineq = J(m, 1, 0)
        
        // -Yλ ≤ -y₀ (Yλ ≥ y₀)
        A_ineq = A_ineq \ (J(s, 1, 0), -Y_ref')
        b_ineq = b_ineq \ (-y0')
        
        lb = (1e-10, J(1, n, 0))
        ub = (1, J(1, n, .))  // θ ≤ 1
        
        if (rts == "vrs") {
            A_eq = (0, J(1, n, 1))
            b_eq = 1
        }
        else {
            A_eq = J(0, n+1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lp.setCoefficients(c)
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        lp.setMaxOrMin("min")
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = lp.value()
            params = lp.parameters()
            lambda = params[2..cols(params)]
        }
        else {
            eff = .
        }
    }
    else {
        // Output-oriented: max φ s.t. Xλ ≤ x₀, φy₀ ≤ Yλ, λ ≥ 0
        // Variables: [λ₁, ..., λₙ, φ]
        c = (J(1, n, 0), -1)  // min -φ
        
        // Xλ ≤ x₀
        A_ineq = (X_ref', J(m, 1, 0))
        b_ineq = x0'
        
        // -Yλ + φy₀ ≤ 0 (φy₀ ≤ Yλ)
        A_ineq = A_ineq \ (-Y_ref', y0')
        b_ineq = b_ineq \ J(s, 1, 0)
        
        lb = (J(1, n, 0), 1)
        ub = (J(1, n, .), .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
        }
        else {
            A_eq = J(0, n+1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lp.setCoefficients(c)
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        lp.setMaxOrMin("min")
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = -lp.value()  // Return φ (not -φ)
            params = lp.parameters()
            lambda = params[1..n]
        }
        else {
            eff = .
        }
    }
    
    return(eff)
}

// ============================================================================
// DUAL (MULTIPLIER) MODEL - Used as fallback when primal fails
// Input-oriented CRS: max u'y₀ s.t. v'x₀ = 1, u'Y - v'X ≤ 0, u,v ≥ 0
// Input-oriented VRS: max u'y₀ + u₀ s.t. v'x₀ = 1, u'Y - v'X + u₀ ≤ 0
// Output-oriented CRS: min v'x₀ s.t. u'y₀ = 1, u'Y - v'X ≤ 0 → eff = 1/value
// Output-oriented VRS: min v'x₀ - u₀ s.t. u'y₀ = 1, u'Y - v'X + u₀ ≤ 0
// ============================================================================
real scalar _sdea_solve_dea_dual(real matrix X_ref, real matrix Y_ref,
                                  real rowvector x0, real rowvector y0,
                                  string scalar orientation, string scalar rts)
{
    real scalar n, m, s, eff, total_vars
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub
    class LinearProgram scalar lp
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (n == 0) return(.)
    
    lp = LinearProgram()
    lp.setTrace("off")
    
    if (orientation == "input") {
        // Input-oriented Multiplier (Dual)
        // max u'y₀ [+ u₀] s.t. v'x₀ = 1, u'Y - v'X [+ u₀] ≤ 0
        // Variables: [u(s), v(m)] or [u(s), v(m), u₀]
        
        if (rts == "crs") {
            total_vars = s + m
            c = (-y0, J(1, m, 0))  // min -u'y₀ = max u'y₀
            A_eq = (J(1, s, 0), x0)  // v'x₀ = 1
            A_ineq = (Y_ref, -X_ref)  // u'Y - v'X ≤ 0
            lb = J(1, total_vars, 0)
        }
        else {
            // VRS
            total_vars = s + m + 1
            c = (-y0, J(1, m, 0), -1)  // min -(u'y₀ + u₀)
            A_eq = (J(1, s, 0), x0, 0)  // v'x₀ = 1
            A_ineq = (Y_ref, -X_ref, J(n, 1, 1))  // u'Y - v'X + u₀ ≤ 0
            lb = (J(1, s+m, 0), .)  // u,v ≥ 0, u₀ free
        }
        b_eq = 1
        b_ineq = J(n, 1, 0)
        ub = J(1, total_vars, .)
        
        lp.setCoefficients(c)
        lp.setMaxOrMin("min")
        lp.setEquality(A_eq, b_eq)
        lp.setInequality(A_ineq, b_ineq)
        lp.setBounds(lb, ub)
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = -lp.value()  // θ = u'y₀ [+ u₀]
        }
        else {
            eff = .
        }
    }
    else {
        // Output-oriented Multiplier (Dual)
        // min v'x₀ [- u₀] s.t. u'y₀ = 1, u'Y - v'X [+ u₀] ≤ 0
        
        if (rts == "crs") {
            total_vars = s + m
            c = (J(1, s, 0), x0)  // min v'x₀
            A_eq = (y0, J(1, m, 0))  // u'y₀ = 1
            A_ineq = (Y_ref, -X_ref)  // u'Y - v'X ≤ 0
            lb = J(1, total_vars, 0)
        }
        else {
            // VRS
            total_vars = s + m + 1
            c = (J(1, s, 0), x0, -1)  // min v'x₀ - u₀
            A_eq = (y0, J(1, m, 0), 0)  // u'y₀ = 1
            A_ineq = (Y_ref, -X_ref, J(n, 1, 1))  // u'Y - v'X + u₀ ≤ 0
            lb = (J(1, s+m, 0), .)  // u,v ≥ 0, u₀ free
        }
        b_eq = 1
        b_ineq = J(n, 1, 0)
        ub = J(1, total_vars, .)
        
        lp.setCoefficients(c)
        lp.setMaxOrMin("min")
        lp.setEquality(A_eq, b_eq)
        lp.setInequality(A_ineq, b_ineq)
        lp.setBounds(lb, ub)
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = lp.value()  // 1/φ
            // Convert to φ to match primal output
            if (eff > 0) {
                eff = 1 / eff
            }
            else {
                eff = .
            }
        }
        else {
            eff = .
        }
    }
    
    return(eff)
}

// Super-efficiency solver for TFDEA forecasting
// Returns φ directly (not 1/φ), where φ < 1 means technology exceeds frontier
// FIX: Separate solver for SE (φ<1 allowed) - from sdea_tf.ado
real scalar _sdea_solve_supereff_basic(real matrix X_all, real matrix Y_all,
                                        real rowvector x0, real rowvector y0,
                                        real scalar idx,
                                        string scalar orientation, string scalar rts)
{
    real scalar phi
    real rowvector lambda
    
    // Call extended solver, discard lambda
    phi = _sdea_solve_supereff_tfdea(X_all, Y_all, x0, y0, idx, orientation, rts, lambda)
    
    return(phi)
}

// Extended Super-efficiency solver for TFDEA - returns both φ and λ
// Following Inman et al. (2006) methodology
real scalar _sdea_solve_supereff_tfdea(real matrix X_all, real matrix Y_all,
                                        real rowvector x0, real rowvector y0,
                                        real scalar idx,
                                        string scalar orientation, string scalar rts,
                                        real rowvector lambda)
{
    real scalar n, n_ref, m, s, eff
    real matrix X_ref, Y_ref, A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub, params
    class LinearProgram scalar lp
    
    n = rows(X_all)
    
    // Exclude the evaluated DMU from reference set
    if (idx == 1) {
        X_ref = X_all[2::n, .]
        Y_ref = Y_all[2::n, .]
    }
    else if (idx == n) {
        X_ref = X_all[1::(n-1), .]
        Y_ref = Y_all[1::(n-1), .]
    }
    else {
        X_ref = X_all[1::(idx-1), .] \ X_all[(idx+1)::n, .]
        Y_ref = Y_all[1::(idx-1), .] \ Y_all[(idx+1)::n, .]
    }
    
    n_ref = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    lambda = J(1, n_ref, 0)
    
    if (n_ref == 0) return(.)
    
    lp = LinearProgram()
    lp.setTrace("off")
    
    if (orientation == "input") {
        // Input-oriented super-efficiency: min θ
        // s.t. θx0 >= X'λ, Y'λ >= y0
        c = (1, J(1, n_ref, 0))
        
        // -θx0 + X'λ <= 0
        A_ineq = (-x0', X_ref')
        b_ineq = J(m, 1, 0)
        
        // -Y'λ <= -y0
        A_ineq = A_ineq \ (J(s, 1, 0), -Y_ref')
        b_ineq = b_ineq \ (-y0')
        
        // θ can be < 1 for super-efficiency (technology exceeds frontier)
        lb = (1e-6, J(1, n_ref, 0))
        ub = (., J(1, n_ref, .))
        
        if (rts == "vrs") {
            A_eq = (0, J(1, n_ref, 1))
            b_eq = 1
        }
        else {
            A_eq = J(0, n_ref+1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lp.setCoefficients(c)
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        lp.setMaxOrMin("min")
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = lp.value()  // Return θ directly
            params = lp.parameters()
            lambda = params[2..cols(params)]
        }
        else {
            eff = .
        }
    }
    else {
        // Output-oriented super-efficiency: max φ
        // s.t. X'λ <= x0, φy0 <= Y'λ, [Σλ=1 for VRS]
        // For super-efficiency, φ can be < 1 (technology exceeds frontier)
        
        // Variables: [λ(n_ref), φ]
        c = (J(1, n_ref, 0), -1)  // min -φ = max φ
        
        // X'λ <= x0
        A_ineq = (X_ref', J(m, 1, 0))
        b_ineq = x0'
        
        // -Y'λ + φy0 <= 0  (equivalent to Y'λ >= φy0)
        A_ineq = A_ineq \ (-Y_ref', y0')
        b_ineq = b_ineq \ J(s, 1, 0)
        
        // φ can be < 1 for super-efficiency
        lb = (J(1, n_ref, 0), 1e-6)
        ub = (J(1, n_ref, .), .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n_ref, 1), 0)
            b_eq = 1
        }
        else {
            A_eq = J(0, n_ref+1, 0)
            b_eq = J(0, 1, 0)
        }
        
        lp.setCoefficients(c)
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        lp.setMaxOrMin("min")
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            // Return φ directly
            eff = -lp.value()
            params = lp.parameters()
            lambda = params[1..n_ref]
        }
        else {
            eff = .
        }
    }
    
    return(eff)
}

//==============================================================================
// FDH (Free Disposal Hull) - Non-convex frontier
//==============================================================================
real scalar _sdea_solve_fdh(real matrix X_ref, real matrix Y_ref,
                             real rowvector x0, real rowvector y0,
                             string scalar orientation)
{
    real scalar n, m, s, j, i, r
    real scalar theta_min, theta_j, phi_max, phi_j
    real scalar is_feasible
    real rowvector xj, yj
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (n == 0) return(.)
    
    if (orientation == "input") {
        // Input-oriented: find min θ such that exists j with x_j ≤ θ*x_0 and y_j ≥ y_0
        theta_min = .
        
        for (j = 1; j <= n; j++) {
            xj = X_ref[j, .]
            yj = Y_ref[j, .]
            
            // Check if j dominates in outputs: y_j ≥ y_0
            is_feasible = 1
            for (r = 1; r <= s; r++) {
                if (yj[r] < y0[r] - 1e-10) {
                    is_feasible = 0
                    break
                }
            }
            
            if (is_feasible) {
                // Find min contraction factor: max_i (x_ji / x_i0)
                theta_j = 0
                for (i = 1; i <= m; i++) {
                    if (x0[i] > 1e-10) {
                        theta_j = max((theta_j, xj[i] / x0[i]))
                    }
                }
                if (theta_min == . | theta_j < theta_min) {
                    theta_min = theta_j
                }
            }
        }
        
        return(theta_min)
    }
    else {
        // Output-oriented: find max φ such that exists j with x_j ≤ x_0 and y_j ≥ φ*y_0
        phi_max = 0
        
        for (j = 1; j <= n; j++) {
            xj = X_ref[j, .]
            yj = Y_ref[j, .]
            
            // Check if j dominates in inputs: x_j ≤ x_0
            is_feasible = 1
            for (i = 1; i <= m; i++) {
                if (xj[i] > x0[i] + 1e-10) {
                    is_feasible = 0
                    break
                }
            }
            
            if (is_feasible) {
                // Find max expansion factor: min_r (y_jr / y_r0)
                phi_j = .
                for (r = 1; r <= s; r++) {
                    if (y0[r] > 1e-10) {
                        if (phi_j == .) {
                            phi_j = yj[r] / y0[r]
                        }
                        else {
                            phi_j = min((phi_j, yj[r] / y0[r]))
                        }
                    }
                }
                if (phi_j != . & phi_j > phi_max) {
                    phi_max = phi_j
                }
            }
        }
        
        // Efficiency = 1/φ
        if (phi_max > 1e-10) {
            return(1 / phi_max)
        }
        return(.)
    }
}

//==============================================================================
// Imprecise DEA (Bounded Data)
// Returns optimistic efficiency (best-case scenario)
//==============================================================================
real scalar _sdea_solve_imprecise(real matrix X_lo, real matrix X_hi,
                                   real matrix Y_lo, real matrix Y_hi,
                                   real scalar k,
                                   string scalar orientation, string scalar rts)
{
    real scalar n, m, s, total_vars
    real rowvector x0, y0, c, lb, ub
    real matrix X_opt, Y_opt, A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar eff
    class LinearProgram scalar lp
    
    n = rows(X_lo)
    m = cols(X_lo)
    s = cols(Y_lo)
    total_vars = n + 1
    
    if (orientation == "input") {
        // Optimistic for k: k uses X_lo, Y_hi; others use X_hi, Y_lo
        x0 = X_lo[k, .]'
        y0 = Y_hi[k, .]'
        
        X_opt = X_hi
        Y_opt = Y_lo
        X_opt[k, .] = X_lo[k, .]
        Y_opt[k, .] = Y_hi[k, .]
        
        // min θ: X'λ ≤ θx0, Y'λ ≥ y0
        c = (J(1, n, 0), 1)
        A_ineq = (X_opt', -x0) \ (-Y_opt', J(s, 1, 0))
        b_ineq = J(m, 1, 0) \ (-y0)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
        }
        else {
            A_eq = J(0, total_vars, 0)
            b_eq = J(0, 1, 0)
        }
        
        lb = J(1, total_vars, 0)
        ub = J(1, total_vars, .)
        
        lp = LinearProgram()
        lp.setCoefficients(c)
        lp.setMaxOrMin("min")
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = lp.value()
        }
        else {
            eff = .
        }
    }
    else {
        // Output-oriented optimistic
        x0 = X_hi[k, .]'
        y0 = Y_lo[k, .]'
        
        X_opt = X_lo
        Y_opt = Y_hi
        X_opt[k, .] = X_hi[k, .]
        Y_opt[k, .] = Y_lo[k, .]
        
        // max φ: X'λ ≤ x0, Y'λ ≥ φy0
        c = (J(1, n, 0), -1)
        A_ineq = (X_opt', J(m, 1, 0)) \ (-Y_opt', y0)
        b_ineq = x0 \ J(s, 1, 0)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), 0)
            b_eq = 1
        }
        else {
            A_eq = J(0, total_vars, 0)
            b_eq = J(0, 1, 0)
        }
        
        lb = J(1, total_vars, 0)
        ub = J(1, total_vars, .)
        
        lp = LinearProgram()
        lp.setCoefficients(c)
        lp.setMaxOrMin("min")
        lp.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
        lp.setBounds(lb, ub)
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = -lp.value()
            if (eff > 0) eff = 1 / eff
        }
        else {
            eff = .
        }
    }
    
    return(eff)
}

//==============================================================================
// COST EFFICIENCY
// min w'x* s.t. Xλ ≤ x*, Yλ ≥ y₀, (Σλ = 1 for VRS)
// CE = w'x* / w'x₀
//==============================================================================
real scalar _sdea_solve_cost(real matrix X_ref, real matrix Y_ref,
                             real rowvector x0, real rowvector y0,
                             real rowvector w0, string scalar rts)
{
    real scalar n, m, s, eff
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub, sol
    real colvector xstar
    real scalar cost0, coststar
    class LinearProgram scalar lp
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (n == 0) return(.)
    
    // Variables: [λ(n), x*(m)]
    c = (J(1, n, 0), w0)
    
    // Inequalities:
    // X'λ - x* ≤ 0 (m rows)
    // -Y'λ ≤ -y₀ (s rows)
    A_ineq = (X_ref', -I(m)) \ (-Y_ref', J(s, m, 0))
    b_ineq = J(m, 1, 0) \ (-y0')
    
    if (rts == "vrs") {
        A_eq = (J(1, n, 1), J(1, m, 0))
        b_eq = 1
    }
    else {
        A_eq = J(0, n+m, 0)
        b_eq = J(0, 1, 0)
    }
    
    lb = J(1, n+m, 0)
    ub = J(1, n+m, .)
    
    lp = LinearProgram()
    lp.setCoefficients(c)
    lp.setMaxOrMin("min")
    lp.setInequality(A_ineq, b_ineq)
    if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
    lp.setBounds(lb, ub)
    lp.setTrace("off")
    
    (void) lp.optimize()
    
    if (lp.converged() == 1 & lp.errorcode() == 0) {
        sol = lp.parameters()
        xstar = sol[(n+1)::(n+m)]'
        coststar = w0 * xstar
        cost0 = w0 * x0'
        if (cost0 > 0) eff = coststar / cost0
        else eff = .
    }
    else eff = .
    
    return(eff)
}

//==============================================================================
// REVENUE EFFICIENCY
// max p'y* s.t. Xλ ≤ x₀, Yλ ≥ y*, (Σλ = 1 for VRS)
// RE = p'y₀ / p'y*
//==============================================================================
real scalar _sdea_solve_revenue(real matrix X_ref, real matrix Y_ref,
                                real rowvector x0, real rowvector y0,
                                real rowvector p0, string scalar rts)
{
    real scalar n, m, s, eff
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub, sol
    real colvector ystar
    real scalar rev0, revstar
    class LinearProgram scalar lp
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (n == 0) return(.)
    
    // Variables: [λ(n), y*(s)]
    c = (J(1, n, 0), p0)
    
    // Inequalities:
    // X'λ ≤ x₀ (m rows)
    // -Y'λ + y* ≤ 0 (s rows)  → y* ≤ Y'λ
    A_ineq = (X_ref', J(m, s, 0)) \ (-Y_ref', I(s))
    b_ineq = x0' \ J(s, 1, 0)
    
    if (rts == "vrs") {
        A_eq = (J(1, n, 1), J(1, s, 0))
        b_eq = 1
    }
    else {
        A_eq = J(0, n+s, 0)
        b_eq = J(0, 1, 0)
    }
    
    lb = J(1, n+s, 0)
    ub = J(1, n+s, .)
    
    lp = LinearProgram()
    lp.setCoefficients(c)
    lp.setMaxOrMin("max")
    lp.setInequality(A_ineq, b_ineq)
    if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
    lp.setBounds(lb, ub)
    lp.setTrace("off")
    
    (void) lp.optimize()
    
    if (lp.converged() == 1 & lp.errorcode() == 0) {
        sol = lp.parameters()
        ystar = sol[(n+1)::(n+s)]'
        revstar = p0 * ystar
        rev0 = p0 * y0'
        if (revstar > 0) eff = rev0 / revstar
        else eff = .
    }
    else eff = .
    
    return(eff)
}

//==============================================================================
// DIRECTIONAL DISTANCE FUNCTION (DDF)
// max β s.t. Xλ + β*gx ≤ x₀, Yλ - β*gy ≥ y₀, (Σλ = 1 for VRS)
//==============================================================================
real scalar _sdea_solve_ddf(real matrix X_ref, real matrix Y_ref,
                            real rowvector x0, real rowvector y0,
                            real rowvector gx0, real rowvector gy0,
                            string scalar rts)
{
    real scalar n, m, s, beta
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub
    class LinearProgram scalar lp
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (n == 0) return(.)
    
    // Variables: [λ(n), β]
    c = (J(1, n, 0), 1)
    
    // Inequalities:
    // X'λ + gx*β ≤ x₀ (m rows)
    // -Y'λ - gy*β ≤ -y₀ (s rows)
    A_ineq = (X_ref', gx0') \ (-Y_ref', -gy0')
    b_ineq = x0' \ (-y0')
    
    if (rts == "vrs") {
        A_eq = (J(1, n, 1), 0)
        b_eq = 1
    }
    else {
        A_eq = J(0, n+1, 0)
        b_eq = J(0, 1, 0)
    }
    
    lb = (J(1, n, 0), 0)
    ub = (J(1, n, .), .)
    
    lp = LinearProgram()
    lp.setCoefficients(c)
    lp.setMaxOrMin("max")
    lp.setInequality(A_ineq, b_ineq)
    if (rows(A_eq) > 0) lp.setEquality(A_eq, b_eq)
    lp.setBounds(lb, ub)
    lp.setTrace("off")
    
    (void) lp.optimize()
    
    if (lp.converged() == 1 & lp.errorcode() == 0) {
        beta = lp.value()
    }
    else beta = .
    
    return(beta)
}

//==============================================================================
// SLACKS-BASED MEASURE (SBM)
// Charnes-Cooper transformation of non-radial SBM model
// min t - (1/m) Σ(s⁻/x₀) s.t. X'λ̃ + s⁻ = t*x₀, Y'λ̃ - s⁺ = t*y₀,
//     t + (1/s) Σ(s⁺/y₀) = 1, (Σλ̃ = t for VRS)
//==============================================================================
real scalar _sdea_solve_sbm(real matrix X_ref, real matrix Y_ref,
                            real rowvector x0, real rowvector y0,
                            string scalar rts)
{
    real scalar n, m, s, i, r, rho
    real matrix A_eq
    real colvector b_eq
    real rowvector c, lb, ub
    real matrix Eq1, Eq2, Eq3, Eq4
    real colvector b1, b2, b3, b4
    class LinearProgram scalar lp
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (n == 0) return(.)
    
    // Variables: [λ̃(n), s⁻(m), s⁺(s), t]
    // Objective: min t - (1/m) Σ(s⁻ᵢ/x₀ᵢ)
    c = J(1, n+m+s+1, 0)
    c[1, n+m+s+1] = 1
    for (i = 1; i <= m; i++) {
        if (x0[i] > 0) c[1, n+i] = -(1/m) * (1/x0[i])
    }
    
    // Equalities:
    // X'λ̃ + s⁻ - t*x₀ = 0 (m)
    Eq1 = (X_ref', I(m), J(m, s, 0), -x0')
    b1 = J(m, 1, 0)
    
    // Y'λ̃ - s⁺ - t*y₀ = 0 (s)
    Eq2 = (Y_ref', J(s, m, 0), -I(s), -y0')
    b2 = J(s, 1, 0)
    
    // t + (1/s) Σ(s⁺ᵣ/y₀ᵣ) = 1
    Eq3 = (J(1, n, 0), J(1, m, 0), J(1, s, 0), 1)
    for (r = 1; r <= s; r++) {
        if (y0[r] > 0) Eq3[1, n+m+r] = (1/s) * (1/y0[r])
    }
    b3 = 1
    
    if (rts == "vrs") {
        // Σλ̃ - t = 0
        Eq4 = (J(1, n, 1), J(1, m, 0), J(1, s, 0), -1)
        b4 = 0
        A_eq = Eq1 \ Eq2 \ Eq3 \ Eq4
        b_eq = b1 \ b2 \ b3 \ b4
    }
    else {
        A_eq = Eq1 \ Eq2 \ Eq3
        b_eq = b1 \ b2 \ b3
    }
    
    lb = J(1, n+m+s+1, 0)
    ub = J(1, n+m+s+1, .)
    
    lp = LinearProgram()
    lp.setCoefficients(c)
    lp.setMaxOrMin("min")
    lp.setEquality(A_eq, b_eq)
    lp.setBounds(lb, ub)
    lp.setTrace("off")
    
    (void) lp.optimize()
    
    if (lp.converged() == 1 & lp.errorcode() == 0) {
        rho = lp.value()
    }
    else rho = .
    
    return(rho)
}

//==============================================================================
// ADDITIVE MODEL
// max Σs⁻ + Σs⁺ s.t. X'λ + s⁻ = x₀, Y'λ - s⁺ = y₀, (Σλ = 1 for VRS)
//==============================================================================
real scalar _sdea_solve_additive(real matrix X_ref, real matrix Y_ref,
                                 real rowvector x0, real rowvector y0,
                                 string scalar rts)
{
    real scalar n, m, s, total_slack, n_eq
    real matrix A_eq
    real colvector b_eq
    real rowvector c, lb, ub
    class LinearProgram scalar lp
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (n == 0) return(.)
    
    // Variables: λ(n) + s⁻(m) + s⁺(s)
    // Objective: max Σs⁻ + Σs⁺ → min -(Σs⁻ + Σs⁺)
    c = (J(1, n, 0), J(1, m, -1), J(1, s, -1))
    
    // Constraints
    n_eq = m + s
    if (rts == "vrs") n_eq = n_eq + 1
    
    A_eq = J(n_eq, n+m+s, 0)
    b_eq = J(n_eq, 1, 0)
    
    // Input: X'λ + s⁻ = x₀
    A_eq[1::m, 1::n] = X_ref'
    A_eq[1::m, (n+1)::(n+m)] = I(m)
    b_eq[1::m] = x0'
    
    // Output: Y'λ - s⁺ = y₀
    A_eq[(m+1)::(m+s), 1::n] = Y_ref'
    A_eq[(m+1)::(m+s), (n+m+1)::(n+m+s)] = -I(s)
    b_eq[(m+1)::(m+s)] = y0'
    
    // VRS: Σλ = 1
    if (rts == "vrs") {
        A_eq[n_eq, 1::n] = J(1, n, 1)
        b_eq[n_eq] = 1
    }
    
    lb = J(1, n+m+s, 0)
    ub = J(1, n+m+s, .)
    
    lp = LinearProgram()
    lp.setCoefficients(c)
    lp.setMaxOrMin("min")
    lp.setEquality(A_eq, b_eq)
    lp.setBounds(lb, ub)
    lp.setTrace("off")
    
    (void) lp.optimize()
    
    if (lp.converged() == 1 & lp.errorcode() == 0) {
        total_slack = -lp.value()
    }
    else total_slack = .
    
    return(total_slack)
}

//==============================================================================
// MULTIPLIER MODEL
// Input-oriented: min x₀'v s.t. y₀'u = 1, Y'u - X'v ≤ 0, (u₀ free for VRS)
// Output-oriented: max y₀'u s.t. x₀'v = 1, Y'u - X'v ≤ 0
//==============================================================================
real scalar _sdea_solve_multiplier(real matrix X_ref, real matrix Y_ref,
                                   real rowvector x0, real rowvector y0,
                                   string scalar orientation, string scalar rts)
{
    real scalar n, m, s, eff, total_vars
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector c, lb, ub
    class LinearProgram scalar lp
    
    n = rows(X_ref)
    m = cols(X_ref)
    s = cols(Y_ref)
    
    if (n == 0) return(.)
    
    // Variables: u(s) + v(m) [+ u₀ for VRS]
    if (rts == "crs") total_vars = s + m
    else total_vars = s + m + 1
    
    if (orientation == "input") {
        // min x₀'v s.t. y₀'u = 1, Yu - Xv ≤ 0
        if (rts == "crs") {
            c = (J(1, s, 0), x0)
            A_eq = (y0, J(1, m, 0))
            A_ineq = (Y_ref, -X_ref)
            lb = J(1, total_vars, 0)
        }
        else {
            c = (J(1, s, 0), x0, 0)
            A_eq = (y0, J(1, m, 0), 1)
            A_ineq = (Y_ref, -X_ref, J(n, 1, 1))
            lb = (J(1, s+m, 0), .)
        }
        b_eq = 1
        b_ineq = J(n, 1, 0)
        ub = J(1, total_vars, .)
        
        lp = LinearProgram()
        lp.setCoefficients(c)
        lp.setMaxOrMin("min")
        lp.setEquality(A_eq, b_eq)
        lp.setInequality(A_ineq, b_ineq)
        lp.setBounds(lb, ub)
        lp.setTrace("off")
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = lp.value()
            if (eff > 0) eff = 1 / eff
            else eff = .
        }
        else eff = .
    }
    else {
        // Output orientation Multiplier (Dual of Output Envelopment)
        // Envelopment: max φ s.t. Xλ ≤ x₀, Yλ ≥ φy₀, [Σλ=1 for VRS]
        // Multiplier:  min v'x₀ - u₀ s.t. u'y₀ = 1, v'Xⱼ - u₀ ≥ u'Yⱼ
        //              Rewrite: min v'x₀ - u₀ s.t. u'y₀ = 1, u'Yⱼ - v'Xⱼ + u₀ ≤ 0
        
        if (rts == "crs") {
            // CRS: min v'x₀ s.t. u'y₀ = 1, u'Yⱼ - v'Xⱼ ≤ 0
            // Variables: [u(s), v(m)]
            c = (J(1, s, 0), x0)              // min v'x₀
            A_eq = (y0, J(1, m, 0))           // u'y₀ = 1
            A_ineq = (Y_ref, -X_ref)          // u'Y - v'X ≤ 0
            lb = J(1, total_vars, 0)
            b_eq = 1
            b_ineq = J(n, 1, 0)
        }
        else {
            // VRS: min v'x₀ - u₀ s.t. u'y₀ = 1, u'Yⱼ - v'Xⱼ + u₀ ≤ 0
            // Variables: [u(s), v(m), u₀]
            c = (J(1, s, 0), x0, -1)          // min v'x₀ - u₀
            A_eq = (y0, J(1, m, 0), 0)        // u'y₀ = 1 (u₀ not in equality)
            A_ineq = (Y_ref, -X_ref, J(n, 1, 1))  // u'Y - v'X + u₀ ≤ 0
            lb = (J(1, s+m, 0), .)            // u,v ≥ 0, u₀ free
            b_eq = 1
            b_ineq = J(n, 1, 0)
        }
        ub = J(1, total_vars, .)
        
        lp = LinearProgram()
        lp.setCoefficients(c)
        lp.setMaxOrMin("min")                 // Minimization
        lp.setEquality(A_eq, b_eq)
        lp.setInequality(A_ineq, b_ineq)
        lp.setBounds(lb, ub)
        lp.setTrace("off")
        
        (void) lp.optimize()
        
        if (lp.converged() == 1 & lp.errorcode() == 0) {
            eff = lp.value()
            // Envelopment output returns 1/φ (≤1), so we do the same
            if (eff > 0) eff = 1 / eff
            else eff = .
        }
        else eff = .
    }
    
    return(eff)
}

//==============================================================================
// ALLOCATIVE EFFICIENCY (Cost-based)
// AE = CE / TE
// Where CE = cost efficiency, TE = technical efficiency
//==============================================================================
real scalar _sdea_solve_allocative_cost(real matrix X_ref, real matrix Y_ref,
                                        real rowvector x0, real rowvector y0,
                                        real rowvector w0, string scalar rts)
{
    real scalar te, ce, ae
    
    // Step 1: Technical Efficiency (input-oriented)
    te = _sdea_solve_dea_basic(X_ref, Y_ref, x0, y0, "input", rts)
    if (te == . | te <= 0) return(.)
    
    // Step 2: Cost Efficiency
    ce = _sdea_solve_cost(X_ref, Y_ref, x0, y0, w0, rts)
    if (ce == .) return(.)
    
    // Step 3: Allocative Efficiency = CE / TE
    ae = ce / te
    
    return(ae)
}

//==============================================================================
// ALLOCATIVE EFFICIENCY (Revenue-based)
// AE = RE / TE
// Where RE = revenue efficiency, TE = technical efficiency (output-oriented)
//==============================================================================
real scalar _sdea_solve_allocative_revenue(real matrix X_ref, real matrix Y_ref,
                                           real rowvector x0, real rowvector y0,
                                           real rowvector p0, string scalar rts)
{
    real scalar te, re, ae
    
    // Step 1: Technical Efficiency (output-oriented)
    te = _sdea_solve_dea_basic(X_ref, Y_ref, x0, y0, "output", rts)
    if (te == . | te <= 0) return(.)
    
    // Step 2: Revenue Efficiency
    re = _sdea_solve_revenue(X_ref, Y_ref, x0, y0, p0, rts)
    if (re == .) return(.)
    
    // Step 3: Allocative Efficiency = RE / TE
    // For output orientation, TE >= 1, RE <= 1
    ae = re * te  // RE = AE * (1/TE) for output
    
    return(ae)
}

//==============================================================================
// BATCH SOLVER - Solve DEA model for all DMUs
// Dispatches to appropriate solver based on model type
//==============================================================================
void _sdea_solve_model_batch(string scalar inputs_str, string scalar outputs_str,
                             string scalar eff_var, string scalar touse_str,
                             string scalar model, string scalar orientation,
                             string scalar rts, string scalar w_str,
                             string scalar p_str, string scalar gx_str,
                             string scalar gy_str)
{
    real matrix X, Y, W, P, GX, GY
    real colvector obs_idx
    real scalar n, k, eff, eff_idx
    real rowvector x0, y0, w0, p0, gx0, gy0
    
    // Load data
    X = st_data(., inputs_str, touse_str)
    Y = st_data(., outputs_str, touse_str)
    obs_idx = selectindex(st_data(., touse_str))
    eff_idx = st_varindex(eff_var)
    
    n = rows(X)
    if (n == 0) return
    
    // Load optional price/direction vectors
    if (w_str != "") W = st_data(., w_str, touse_str)
    if (p_str != "") P = st_data(., p_str, touse_str)
    if (gx_str != "") GX = st_data(., gx_str, touse_str)
    else GX = X  // Default: use own inputs as direction
    if (gy_str != "") GY = st_data(., gy_str, touse_str)
    else GY = Y  // Default: use own outputs as direction
    
    // Solve for each DMU
    for (k = 1; k <= n; k++) {
        x0 = X[k, .]
        y0 = Y[k, .]
        
        if (model == "te") {
            eff = _sdea_solve_dea_basic(X, Y, x0, y0, orientation, rts)
        }
        else if (model == "cost") {
            w0 = W[k, .]
            eff = _sdea_solve_cost(X, Y, x0, y0, w0, rts)
        }
        else if (model == "revenue") {
            p0 = P[k, .]
            eff = _sdea_solve_revenue(X, Y, x0, y0, p0, rts)
        }
        else if (model == "ddf") {
            gx0 = GX[k, .]
            gy0 = GY[k, .]
            eff = _sdea_solve_ddf(X, Y, x0, y0, gx0, gy0, rts)
        }
        else if (model == "sbm") {
            eff = _sdea_solve_sbm(X, Y, x0, y0, rts)
        }
        else if (model == "additive") {
            eff = _sdea_solve_additive(X, Y, x0, y0, rts)
        }
        else if (model == "multiplier") {
            eff = _sdea_solve_multiplier(X, Y, x0, y0, orientation, rts)
        }
        else if (model == "allocative") {
            if (w_str != "") {
                w0 = W[k, .]
                eff = _sdea_solve_allocative_cost(X, Y, x0, y0, w0, rts)
            }
            else if (p_str != "") {
                p0 = P[k, .]
                eff = _sdea_solve_allocative_revenue(X, Y, x0, y0, p0, rts)
            }
            else eff = .
        }
        else if (model == "supereff") {
            eff = _sdea_solve_supereff_basic(X, Y, x0, y0, k, orientation, rts)
        }
        else if (model == "supersbm") {
            eff = _sdea_solve_supersbm_internal(X, Y, x0, y0, k, rts)
        }
        else if (model == "fdh") {
            eff = _sdea_solve_fdh(X, Y, x0, y0, orientation)
        }
        else {
            // Default to TE
            eff = _sdea_solve_dea_basic(X, Y, x0, y0, orientation, rts)
        }
        
        st_store(obs_idx[k], eff_idx, eff)
    }
    
    st_numscalar("r(n_converged)", n - colmissing(st_data(., eff_var, touse_str)))
    st_numscalar("r(n_failed)", colmissing(st_data(., eff_var, touse_str)))
}


void _sdea_static_pareto(string scalar xs, string scalar ys, string scalar dest)
{
    real matrix Z
    real colvector perm, witnesses, keep
    real scalar i, j, k, n
    Z = (st_data(.,tokens(xs)), -st_data(.,tokens(ys)))
    perm = order((Z,(1::rows(Z))),1..(cols(Z)+1))
    witnesses = keep = J(rows(Z),1,0)
    n = 0
    for (i=1; i<=rows(Z); i++) {
        k = perm[i]
        for (j=1; j<=n; j++) {
            if (all(Z[witnesses[j],.] :<= Z[k,.])) break
        }
        if (j>n) {
            n++
            witnesses[n] = k
            keep[k] = 1
        }
    }
    st_store(.,dest,keep)
}

// Compare in the original units, relative to the target (including exact zeros).
real scalar _sdea_static_feasible(real colvector w, real matrix X, real matrix Y,
    real rowvector xmax, real rowvector ymin, real scalar tol)
{
    real rowvector xfit, yfit
    xfit = w' * X
    yfit = w' * Y
    if (hasmissing((xfit,yfit,xmax,ymin))) return(0)
    return(all(xfit :- xmax :<= tol:*xmax) & all(ymin :- yfit :<= tol:*ymin))
}

real scalar _sdea_static_score(real matrix X, real matrix Y, real matrix rawX,
    real matrix rawY, real rowvector x, real rowvector y, real rowvector rawx0,
    real rowvector rawy0, string scalar ort, string scalar rts)
{
    struct RSM_Prob scalar prob
    struct DEA_Info scalar info
    struct RSM_Result scalar sol
    real scalar radial, score
    real colvector w
    real rowvector xmax, ymin
    _dea_lp_stage1(prob,info,X,Y,x',y',ort,rts)
    sol = _rsm_two_phase(prob,1e-7,"")
    if (sol.status!=0) return(.)
    radial = _dea_extract_theta(sol,info)
    if (radial<=0 | radial>=.) return(.)
    w = sol.x[info.idx_lam]
    if (hasmissing(w) | min(w) < -1e-10) return(.)
    w = w :* (w:>0)
    if (rts=="vrs" & abs(sum(w)-1)>1e-7) return(.)
    xmax = rawx0
    ymin = rawy0
    if (ort=="input") {
        score = radial
        xmax = radial * xmax
    }
    else {
        score = 1/radial
        ymin = radial * ymin
    }
    if (score<=0 | score>1+1e-7) return(.)
    if (!_sdea_static_feasible(w,rawX,rawY,xmax,ymin,1e-7)) return(.)
    return(min((1,score)))
}

real colvector _sdea_score_rows(real matrix X, real matrix Y, real matrix rawX, real matrix rawY,
    real colvector fr, string scalar ort, string scalar rts, real scalar first,
    real scalar last, real scalar retries)
{
    real matrix FX,FY,origFX,origFY
    real colvector scores
    real scalar k
    FX = X[fr,.]
    FY = Y[fr,.]
    origFX = rawX[fr,.]
    origFY = rawY[fr,.]
    scores = J(last-first+1,1,.)
    retries = 0
    for (k=first; k<=last; k++) {
        scores[k-first+1] = _sdea_static_score(FX,FY,origFX,origFY,X[k,.],Y[k,.],
            rawX[k,.],rawY[k,.],ort,rts)
        if (missing(scores[k-first+1])) {
            retries++
            scores[k-first+1] = _sdea_static_score(X,Y,rawX,rawY,X[k,.],Y[k,.],
                rawX[k,.],rawY[k,.],ort,rts)
            if (missing(scores[k-first+1])) {
                errprintf("full-reference radial score failed for estimation-sample row %g\n",k)
                _error(430)
            }
        }
    }
    return(scores)
}

// All workers read the same normalized and original matrices and final frame.
void _sdea_worker_run(real scalar i, real scalar workers)
{
    real scalar f, first, last, retries
    real matrix X,Y,rawX,rawY
    real colvector fr,scores
    string rowvector paths,model
    f = fopen("../sample.bin","r")
    paths = fgetmatrix(f)
    model = fgetmatrix(f)
    X = fgetmatrix(f)
    Y = fgetmatrix(f)
    rawX = fgetmatrix(f)
    rawY = fgetmatrix(f)
    fr = fgetmatrix(f)
    fclose(f)
    first = floor((i-1)*rows(X)/workers)+1
    last = floor(i*rows(X)/workers)
    scores = _sdea_score_rows(X,Y,rawX,rawY,fr,model[1],model[2],first,last,retries)
    f = fopen("result.bin","w")
    fputmatrix(f,(first,last,retries))
    fputmatrix(f,scores)
    fclose(f)
}

string scalar _sdea_shell_quote(string scalar s)
{
    return("'" + subinstr(s,"'", "'"+char(92)+"''") + "'")
}

void _sdea_worker_cleanup(string scalar job, real scalar workers)
{
    real scalar i,j
    string scalar folder
    string colvector files
    for (i=1; i<=workers; i++) {
        folder = job+"/"+strofreal(i)
        if (!direxists(folder)) continue
        files = dir(folder,"files","*")
        for (j=1; j<=rows(files); j++) unlink(folder+"/"+files[j])
        rmdir(folder)
    }
    unlink(job+"/sample.bin")
    unlink(job+"/worker.do")
    unlink(job+"/launch.sh")
    rmdir(job)
}

real colvector _sdea_parallel_scores(real matrix X, real matrix Y, real matrix rawX, real matrix rawY,
    real colvector fr, string scalar ort, string scalar rts, real scalar workers, real scalar retries)
{
    real scalar f, script, i, first, last, rc
    real colvector scores, part
    real rowvector info
    string scalar job, folder, exe, driver
    job = st_local("job")
    exe = _sdea_shell_quote(st_local("executable"))
    driver = _sdea_shell_quote(pathjoin(pathgetparent(st_local("source")),"_sdea_worker.do"))
    f = fopen(job+"/sample.bin","w")
    fputmatrix(f,(st_local("source"),st_local("core")))
    fputmatrix(f,(ort,rts))
    fputmatrix(f,X)
    fputmatrix(f,Y)
    fputmatrix(f,rawX)
    fputmatrix(f,rawY)
    fputmatrix(f,fr)
    fclose(f)
    script = fopen(job+"/launch.sh","w")
    fput(script,"cd "+_sdea_shell_quote(job)+" || exit 1")
    fput(script,"cp "+driver+" worker.do || exit 1")
    for (i=1; i<=workers; i++) {
        folder = job+"/"+strofreal(i)
        mkdir(folder)
        f = fopen(folder+"/profile.do","w")
        fput(f,"sysdir set PERSONAL .")
        fput(f,"sysdir set PLUS .")
        fclose(f)
        fput(script,"(cd "+strofreal(i)+" && exec "+exe+" -b do ../worker.do "+
            strofreal(i)+" "+strofreal(workers)+") &")
    }
    fput(script,"wait")
    fclose(script)
    st_local("launcher",_sdea_shell_quote(job+"/launch.sh"))
    stata("shell /bin/sh "+char(96)+"macval(launcher)"+char(39))

    scores = J(rows(X),1,.)
    retries = 0
    for (i=1; i<=workers; i++) {
        folder = job+"/"+strofreal(i)
        if (!fileexists(folder+"/status")) {
            errprintf("Stata worker %g did not finish. Check the Stata/MP batch executable and ado installation.\n",i)
            _error(603)
        }
        f = fopen(folder+"/status","r")
        rc = strtoreal(fget(f))
        fclose(f)
        if (rc!=0) {
            errprintf("Stata worker %g failed with return code %g\n",i,rc)
            _error(rc)
        }
        f = fopen(folder+"/result.bin","r")
        info = fgetmatrix(f)
        part = fgetmatrix(f)
        fclose(f)
        first = floor((i-1)*rows(X)/workers)+1
        last = floor(i*rows(X)/workers)
        if (cols(info)!=3 | rows(info)!=1) _error(459)
        if (info[1]!=first | info[2]!=last | rows(part)!=last-first+1 | cols(part)!=1) _error(459)
        if (hasmissing(part) | min(part)<=0 | max(part)>1) _error(459)
        scores[|first\last|] = part
        retries = retries + info[3]
    }
    return(scores)
}

real rowvector _sdea_static_evaluate(string scalar xs, string scalar ys,
    string scalar nx, string scalar ny, string scalar candidate, string scalar member,
    string scalar dest, string scalar ort, string scalar rts, | real scalar workers)
{
    real matrix X, Y, rawX, rawY, FX, FY, origFX, origFY
    real colvector keep, fr, omitted, scores
    real rowvector w
    real scalar i,k,alpha,retained,retries
    if (args()==9) workers=1
    rawX = st_data(.,tokens(xs))
    rawY = st_data(.,tokens(ys))
    X = st_data(.,tokens(nx))
    Y = st_data(.,tokens(ny))
    keep = st_data(.,member)
    fr = vec(selectindex(keep:==1))
    omitted = vec(selectindex(st_data(.,candidate) :& !keep))
    // Certify each omitted Pareto candidate against the FINAL frame. This
    // also covers discarded starters and simultaneous prune decisions.
    retained = 0
    FX = X[fr,.]
    FY = Y[fr,.]
    origFX = rawX[fr,.]
    origFY = rawY[fr,.]
    for (i=1; i<=rows(omitted); i++) {
        k = omitted[i]
        if (rows(fr)) {
            w = J(1,0,.)
            alpha = _sdea_buildhull_eval(FX,FY,X[k,.]',Y[k,.]',"simplex",w)
            if (alpha<1e-9 & !hasmissing(w) & min(w)>=-1e-10) {
                w = w :* (w:>0)
                if (abs(sum(w)-1)<1e-9) {
                    w = w / sum(w)
                    if (_sdea_static_feasible(w',origFX,origFY,rawX[k,.],rawY[k,.],1e-10)) continue
                }
            }
        }
        keep[k] = 1
        retained++
    }
    fr = vec(selectindex(keep))
    if (workers==1) {
        scores = _sdea_score_rows(X,Y,rawX,rawY,fr,ort,rts,1,rows(X),retries)
    }
    else {
        scores = _sdea_parallel_scores(X,Y,rawX,rawY,fr,ort,rts,workers,retries)
    }
    st_store(.,member,keep)
    st_store(.,dest,scores)
    return((retained,retries))
}
end
