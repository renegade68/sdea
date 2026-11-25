{smcl}
{* *! version 5.8.0  26nov2026}{...}
{viewerjumpto "Syntax" "sdea##syntax"}{...}
{viewerjumpto "Description" "sdea##description"}{...}
{viewerjumpto "Options" "sdea##options"}{...}
{viewerjumpto "DEA Models" "sdea##models"}{...}
{viewerjumpto "Streaming Methods" "sdea##methods"}{...}
{viewerjumpto "Drift Detection" "sdea##drift"}{...}
{viewerjumpto "Stored Results" "sdea##results"}{...}
{viewerjumpto "Examples" "sdea##examples"}{...}
{viewerjumpto "Technical Notes" "sdea##technical"}{...}
{viewerjumpto "References" "sdea##references"}{...}
{viewerjumpto "Author" "sdea##author"}{...}
{title:Title}

{p2colset 5 13 15 2}{...}
{p2col:{bf:sdea} {hline 2}}Streaming Data Envelopment Analysis{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:sdea} {varlist} {ifin}{cmd:,} {opt inp:uts(varlist)} [{it:options}]

{pstd}
where {it:varlist} contains the output variables.


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt inp:uts(varlist)}}input variables (inputs to production process){p_end}

{syntab:Time/Panel}
{synopt:{opt time(varname)}}time period identifier (required for streaming methods){p_end}
{synopt:{opt dmu(varname)}}DMU identifier variable{p_end}

{syntab:Method}
{synopt:{opt meth:od(string)}}streaming method; default is {bf:buildhull}{p_end}
{synopt:{opt mod:el(string)}}DEA model; default is {bf:te}{p_end}
{synopt:{opt refs:et(full|buildhull)}}reference set for {cmd:method(standard) model(te)}; default {cmd:full}{p_end}

{syntab:DEA Specification}
{synopt:{opt ori:entation(string)}}orientation: {bf:input} or {bf:output} (default){p_end}
{synopt:{opt rts(string)}}returns to scale: {bf:vrs} (default) or {bf:crs}{p_end}

{syntab:Model-Specific Options}
{synopt:{opt w(varlist)}}input prices for cost/allocative efficiency{p_end}
{synopt:{opt p(varlist)}}output prices for revenue/allocative efficiency{p_end}
{synopt:{opt gx(varlist)}}input direction vector for DDF{p_end}
{synopt:{opt gy(varlist)}}output direction vector for DDF{p_end}
{synopt:{opt bounded(varlist)}}lower bound variables for imprecise DEA{p_end}
{synopt:{opt upper(varlist)}}upper bound variables for imprecise DEA{p_end}

{syntab:Method-Specific Options}
{synopt:{opt star:ters(#)}}initial DMUs for modified BuildHull; see details below{p_end}
{synopt:{opt bhvariant(string)}}BuildHull variant: {bf:modified} (default) or {bf:dula}{p_end}
{synopt:{opt solver(string)}}hull LP solver: {bf:ipm}, {bf:simplex}, or {bf:auto}; defaults depend on method{p_end}
{synopt:{opt win:dow(#)}}window size for window/moving methods; default is 3{p_end}
{synopt:{opt tf(#)}}frontier time cutoff for TFDEA{p_end}
{synopt:{opt decay(#)}}time-decay parameter for TDW-TFDEA; default is 0 (off){p_end}
{synopt:{opt decaytype(string)}}TDW decay function; default is {bf:exponential}{p_end}
{synopt:{opt semodel(string)}}super-efficiency model for TFDEA: {bf:supereff} or {bf:supersbm}{p_end}
{synopt:{opt global}}use global frontier for Malmquist{p_end}
{synopt:{opt prune(#)}}pruning interval for BuildHull (0=no pruning){p_end}
{synopt:{opt gate(newvar)}}create frontier update gate status for BuildHull{p_end}
{synopt:{opt gateeps(#)}}tolerance for enveloped/no-update gate; default is 1e-8{p_end}
{synopt:{opt mlscore(varname)}}user-supplied rule-based or ML score for accepting candidate updates{p_end}
{synopt:{opt mlthreshold(#)}}minimum {cmd:mlscore()} required to update frontier; default is 0.5{p_end}

{syntab:Drift Detection}
{synopt:{opt drifttest}}perform drift diagnostics{p_end}
{synopt:{opt driftalgo(string)}}drift algorithm: {bf:multi} (default), {bf:frontier}, {bf:cusum}, {bf:adwin}, {bf:ph}, {bf:ddm}, or {bf:eddm}{p_end}
{synopt:{opt driftmode(string)}}detector input mode: {bf:dea} (default) or {bf:original}{p_end}
{synopt:{opt drifterror(varname)}}0/1 error indicator for {bf:driftmode(original)} with DDM or EDDM{p_end}
{synopt:{opt by(varname)}}numeric group variable for heterogeneity diagnostic with {bf:driftalgo(frontier)}{p_end}
{synopt:{opt driftwin:dow(#)}}rolling window size for drift detection; default is 3{p_end}
{synopt:{opt thrvol(#)}}volatility threshold; default is 0.1{p_end}
{synopt:{opt threxp(#)}}expansion threshold; default is 0.2{p_end}
{synopt:{opt thrturn(#)}}turnover threshold; default is 0.3{p_end}
{synopt:{opt threshold(#)}}CUSUM threshold; default is 5{p_end}

{syntab:Output}
{synopt:{opt gen:erate(namelist)}}names for efficiency, frame, and drift variables{p_end}
{synopt:{opt replace}}replace existing variables{p_end}
{synopt:{opt verbose}}display detailed progress{p_end}
{synopt:{opt gr:aph}}display visualization graphs{p_end}
{synopt:{opt saving(string)}}save combined graph to file{p_end}
{synopt:{opt saveref(string)}}save reference set to file (note: may contain missing values if dual fallback was used){p_end}
{synopt:{opt savehistory(string)}}save frame history to file{p_end}
{synopt:{opt summarysave(string)}}save period-level monitoring summary dataset{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:sdea} performs Streaming Data Envelopment Analysis, combining traditional DEA
efficiency measurement with time-series and streaming data techniques. It supports
12 DEA models and 9 methods. The streaming {cmd:method(buildhull)} uses Stata's
{cmd:LinearProgram()} interior-point solver by default and can optionally use the
package simplex engine for comparison and smaller problems.

{pstd}
The command evaluates the relative efficiency of Decision Making Units (DMUs) using
various DEA formulations and can track how efficiency evolves over time using
streaming algorithms.

{pstd}
{cmd:sdea} 5.8.0 requires Stata 16 or later and the complete {cmd:ado} directory of
{browse "https://github.com/renegade68/DEA-Stata/archive/ad77287466b22c3edd1ae3cc074981225c95d424.zip":DEA-Stata 2.0.73 (tested ZIP)}.
Extract the ZIP's top-level folder beside {cmd:load.do}. In a fresh Stata session,
change to the folder containing {cmd:load.do} and run {cmd:do load.do}.
The README includes the Stata download command and instructions for an existing core installation.
{cmd:sdea.ado} adds the streaming-specific components:
BuildHull frame maintenance, update gates, drift diagnostics, and time-indexed
frontier-management quantities.

{pstd}
{bf:Key Features:}

{phang2}• 12 DEA models (technical, cost, revenue, SBM, DDF, etc.){p_end}
{phang2}• 9 time-aware or streaming methods (BuildHull, FullStream, Malmquist, TFDEA, etc.){p_end}
{phang2}• Drift diagnostics for input-output relationships and DEA frontier signals{p_end}
{phang2}• Malmquist Productivity Index with TFPCH/EFFCH/TECHCH/PECH/SECH decomposition{p_end}
{phang2}• Technology Forecasting DEA (TFDEA){p_end}
{phang2}• Interior-point optimization via LinearProgram() for robust convergence, with simplex comparison available for BuildHull{p_end}
{phang2}• Dual fallback: automatic primal-to-dual switching for numerical stability{p_end}
{phang2}• Prune option for computational efficiency on large datasets{p_end}


{title:Static radial scores with a reduced reference set}

{pstd}
{cmd:method(standard) model(te) refset(buildhull)} evaluates every selected firm against
a common final reference set for a single cross-section. Inputs and outputs must be
nonnegative with at least one positive input and output per observation. Both CRS and
VRS and both orientations are supported. Use {cmd:if} or {cmd:in} to choose a period
for a static run. With {cmd:refset(full)}, the ordinary standard estimator also supports
{cmd:time()} and estimates each period separately.

{pstd}
The first {cmd:generate()} variable contains input efficiency theta or inverse output
expansion 1/phi. The second equals one for retained reference firms, which can include
redundant generators. The third is zero (static drift). All three are missing outside
the sample. The score is radial stage 1. Slacks, peer weights and a scale decomposition
require separate estimation. The default streaming {cmd:method(buildhull)} returns its
usual arrival-time hull distance alpha.

{pstd}
The static path screens for coordinatewise dominance, builds a VRS frame and checks
omitted candidates for coverage. Its default is {cmd:bhvariant(modified) solver(simplex) prune(0)}.
{cmd:solver(ipm|simplex|auto)} selects the construction policy. Final coverage checks and
radial scoring use the DEA-Stata 2.0.73 revised simplex. {cmd:bhvariant(dula)} uses
hyperplane selection. Positive {cmd:starters()} and {cmd:prune()} apply to the modified
variant. A failed coverage check retains the candidate. A failed score or primal
residual check triggers a full-reference retry for that firm. If that retry also fails,
the command exits with r(430) before publishing the requested variables.

{pstd}
For example, using two output variables and five inputs:

{phang2}{cmd:. sdea y1 y2, inputs(x1 x2 x3 x4 x5) method(standard) refset(buildhull) rts(vrs) orientation(output) generate(te isref drift)}{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt inputs(varlist)} specifies the input variables. These represent resources
consumed in the production process.

{dlgtab:Time/Panel}

{phang}
{opt time(varname)} specifies the time period variable. Required for all streaming
methods except {bf:standard}.

{phang}
{opt dmu(varname)} specifies the DMU identifier variable. Required for {bf:malmquist}
method; optional for others.

{dlgtab:Method}

{phang}
{opt method(string)} specifies the streaming method. Options are:

{p2colset 9 25 27 2}{...}
{p2col:{bf:buildhull}}BuildHull VRS frame-maintenance algorithm (default){p_end}
{p2col:{bf:fullstream}}FullStream with self-inclusive evaluation{p_end}
{p2col:{bf:cumulative}}Cumulative expanding frontier{p_end}
{p2col:{bf:moving}}Moving window frontier{p_end}
{p2col:{bf:standard}}Standard cross-sectional DEA{p_end}
{p2col:{bf:time}}Time-based sequential evaluation{p_end}
{p2col:{bf:window}}Window DEA analysis{p_end}
{p2col:{bf:malmquist}}Malmquist Productivity Index{p_end}
{p2col:{bf:tfdea}}Technology Forecasting DEA{p_end}

{phang}
{opt model(string)} specifies the DEA model. See {help sdea##models:DEA Models} below.

{dlgtab:DEA Specification}

{phang}
{opt orientation(string)} specifies whether the DEA scoring layer is input-oriented
(minimize inputs) or output-oriented (maximize outputs). Default is {bf:output}.
For {bf:method(buildhull)}, the Dulá-style frame-maintenance step is non-oriented;
orientation applies only to methods that compute conventional DEA scores.

{phang}
{opt rts(string)} specifies returns to scale assumption. {bf:vrs} assumes variable
returns to scale (default); {bf:crs} assumes constant returns to scale.
The current {bf:method(buildhull)} implementation follows the VRS hull/frame
construction and does not change its hull LP when {bf:rts(crs)} is specified.

{dlgtab:Model-Specific Options}

{phang}
{opt w(varlist)} specifies input price variables. Required for {bf:cost} and
{bf:allocative} (cost-based) models. Must have same number of variables as inputs.

{phang}
{opt p(varlist)} specifies output price variables. Required for {bf:revenue} and
{bf:allocative} (revenue-based) models. Must have same number of variables as outputs.

{phang}
{opt gx(varlist)} specifies input direction vector for DDF model. If not specified,
own inputs are used as direction.

{phang}
{opt gy(varlist)} specifies output direction vector for DDF model. If not specified,
own outputs are used as direction.

{phang}
{opt bounded(varlist)} specifies lower bound variables for imprecise DEA.

{phang}
{opt upper(varlist)} specifies upper bound variables for imprecise DEA.

{dlgtab:Method-Specific Options}

{phang}
{opt starters(#)} specifies the number of initial DMUs used to build the initial
frame in the modified BuildHull method. The automatic size is the largest of 30, three times the number of inputs and outputs combined,
and the integer part of the square root of sample size. The streaming method caps its
automatic size at 80% of the sample. Static BuildHull applies the rule to the
Pareto candidates and caps the starter count at the number of candidates.

{phang}
{opt bhvariant(string)} specifies the BuildHull variant for {bf:method(buildhull)} or {bf:refset(buildhull)}.
{bf:modified} is the default streaming-oriented implementation: exterior arrivals
are appended to the maintained frame and optional pruning removes redundant frame
members. {bf:dula} uses a closer implementation of Dulá's original BuildHull phase-1
frame-identification procedure: it solves the hull-membership LP, obtains the
separating hyperplane from the dual, and selects the next frame element by
hyperplane translation over the remaining candidate set.

{phang}
{opt solver(string)} specifies the LP solver used by the BuildHull
hull-membership and Dulá-style dual-plane LPs. The streaming method defaults to {bf:ipm}, which
uses Stata's {cmd:LinearProgram()} interior-point solver and is recommended for
large or degenerate streaming problems.  {bf:simplex} uses the package's
two-phase revised simplex engine.  {bf:auto} tries the IPM solver first and
falls back to simplex if the BuildHull LP does not converge.  The BuildHull
classification is solver-independent in principle when the LP is solved to
optimality, but near-frontier observations can be sensitive to numerical
tolerances.

{phang}
{opt window(#)} specifies the window size for {bf:window} and {bf:moving} methods.
Default is 3.

{phang}
{opt tf(#)} specifies the frontier time cutoff for TFDEA. DMUs with time <= tf
form the frontier; DMUs with time > tf are evaluated against this frontier.
The cutoff may be a real value, which is useful for decimal time variables.

{phang}
{opt decay(#)} specifies the time-decay parameter for time-decay weighted
TFDEA aggregation. The default is 0, which disables time-decay weighting and
uses the standard unweighted aggregate ROC.

{phang}
{opt decaytype(string)} specifies the TDW decay function. Options are
{bf:exponential}, {bf:hyperbolic}, {bf:linear}, {bf:gaussian}, and
{bf:logistic}. The default is {bf:exponential}.

{phang}
{opt semodel(string)} specifies the super-efficiency model for TFDEA. Options are
{bf:supereff} (radial super-efficiency, default) or {bf:supersbm} (Super-SBM,
Tone 2002, EJOR 143:32-41). Super-SBM uses the Charnes-Cooper linearization
of the original ratio formulation.

{phang}
{opt global} uses a global frontier (pooled across all periods) for Malmquist
index calculation instead of adjacent period frontiers.

{phang}
{opt prune(#)} specifies the pruning interval for modified BuildHull. A positive value triggers pruning after
initialization, after every {it:#} accepted exterior updates, and at the end.
The default is zero.

{phang}
{opt gate(newvar)} creates a frontier update gate variable for {bf:method(buildhull)}.
For each incoming observation, {cmd:sdea} first checks whether the observation is
already enveloped by the currently maintained frontier frame. If it is enveloped,
the frontier is not updated. If it is not enveloped, the observation is treated as
a candidate frontier-changing observation. The generated gate variable takes four
values: 0 = enveloped/no update, 1 = accepted frontier update, 2 = candidate
rejected, and 3 = starter observation.

{phang}
{opt gateeps(#)} specifies the numerical tolerance used by the frontier update
gate. A BuildHull alpha less than or equal to {it:#} is treated as enveloped by
the maintained frontier. The default is 1e-8.

{phang}
{opt mlscore(varname)} supplies an optional externally computed rule-based or
machine-learning score for candidate frontier-changing observations. When
specified, a candidate updates the maintained frontier only if {it:varname} is
nonmissing and at least {cmd:mlthreshold()}. Candidates below the threshold are
rejected and do not enter the frame. This option does not estimate a machine
learning model internally; it provides an interface for externally computed
admissibility scores. If {cmd:mlscore()} is specified without {cmd:gate()}, the
admissibility rule still operates and gate counts are stored in {cmd:r()}, but no
observation-level gate status variable is saved. The gate interface is supported
with {cmd:bhvariant(modified)}, the default BuildHull implementation.

{phang}
{opt mlthreshold(#)} specifies the minimum score required for a candidate
frontier-changing observation to update the maintained frontier when
{cmd:mlscore()} is supplied. The default is 0.5.


{marker drift}{...}
{title:Drift Detection Options}

{phang}
{opt drifttest} performs drift diagnostics for input-output relationship shifts
and DEA-derived frontier signals over time.

{phang}
{opt driftalgo(string)} specifies the drift detection algorithm:

{p2colset 9 20 22 2}{...}
{p2col:{bf:multi}}Multi-indicator framework using volatility, expansion, and turnover (default){p_end}
{p2col:{bf:frontier}}Conditional frontier relationship drift; distinguishes support expansion from input-output relationship shifts{p_end}
{p2col:{bf:cusum}}CUSUM-based mean-shift detection on the DEA signal{p_end}
{p2col:{bf:adwin}}ADWIN-style subwindow mean comparison for the DEA signal{p_end}
{p2col:{bf:ph}}Page-Hinkley detection{p_end}
{p2col:{bf:ddm}}DDM-style deviation monitoring for the DEA signal{p_end}
{p2col:{bf:eddm}}EDDM-style distance-between-deviations monitoring{p_end}

{phang}
{opt driftwindow(#)} specifies the rolling window size for computing drift indicators.
Default is 3 periods.

{phang}
{bf:driftalgo(frontier)} estimates a rolling log-linear input-output relationship
from the preceding window and tests whether current-period output residuals shift
relative to that relationship. This option is intended to separate ordinary
frontier expansion caused by new support points from concept drift in the
conditional production relationship. It uses {opt thrvol()} as a z-score threshold;
values below 1.5 are treated as 1.5 for this algorithm.

{phang}
{opt by(varname)} is a heterogeneity diagnostic for {bf:driftalgo(frontier)}. It
fits the rolling conditional relationship separately within the numeric groups
defined by {it:varname} and reports whether drift remains after conditioning on
observable group composition. This option does not define a group-specific DEA
technology or alter the maintained BuildHull reference set; it is a sensitivity
check for interpreting drift signals under observable DMU heterogeneity.

{phang}
The ADWIN, DDM, and EDDM options are adaptations of standard data-stream
detectors to univariate DEA-derived signals. They are not literal supervised
classification-error implementations. Use them for robustness checks against
standard change-detection logic, and interpret them as signal detectors rather
than as complete concept-drift adaptation methods.

{phang}
All drift algorithms can be used for score-level monitoring when the DEA output
series is explicitly treated as the data stream generated by a fixed DEA
workflow. In that interpretation, drift means a change in the distribution of
DEA-derived benchmarking outcomes, such as efficiency scores or BuildHull alpha
hull-distance values. It should not be interpreted as production-technology
drift unless supported by input-output relationship or frontier-structure
evidence.

{phang}
{opt driftmode(string)} specifies the detector input mode. The default,
{bf:driftmode(dea)}, applies the selected detector to DEA-derived signals. With
{bf:driftalgo(ddm)}, the monitored input is |1 - alpha|. With
{bf:driftalgo(eddm)}, an error event is a DEA signal deviation from alpha = 1.
The option {bf:driftmode(original)} implements the original supervised
DDM/EDDM input convention by using a user-supplied binary error stream.

{phang}
{opt drifterror(varname)} supplies the 0/1 error indicator required by
{bf:driftmode(original)}. It is currently supported only with
{bf:driftalgo(ddm)} and {bf:driftalgo(eddm)}. This option is intended for
reproducibility checks against the original DDM/EDDM detectors or for settings
where the analyst has externally defined whether each arriving observation is an
error relative to a maintained prediction rule. It does not alter the DEA
frontier, reference set, or efficiency calculation.

{phang}
{opt thrvol(#)} specifies the threshold for efficiency volatility indicator (V_t).
Default is 0.1. Lower values increase sensitivity.

{phang}
{opt threxp(#)} specifies the threshold for frontier expansion rate indicator (D_t).
Default is 0.2. For {bf:method(buildhull)}, this indicator uses the observed
share of arrivals entering the maintained frame in each period.

{phang}
{opt thrturn(#)} specifies the threshold for reference set turnover indicator (T_t).
Default is 0.3. For {bf:method(buildhull)}, this indicator uses frame entries
and pruned frame members to approximate period-to-period frame composition change.

{phang}
{opt threshold(#)} specifies the CUSUM threshold when using {bf:driftalgo(cusum)}.
Default is 5.

{pstd}
{bf:Multi-indicator drift diagnostic:} A drift signal is raised when at least 2 of
the 3 indicators (V_t, D_t, T_t) exceed their respective thresholds simultaneously.
This reduces false positives compared to single-indicator methods.

{dlgtab:Output}

{phang}
{opt generate(namelist)} specifies names for generated variables. Provide up to
three names: main score or distance, frame indicator, and drift magnitude. Default names
are {bf:_sdea_alpha}, {bf:_sdea_frame}, and {bf:_sdea_drift}.

{phang}
{opt replace} allows replacing existing variables with the same names.

{phang}
{opt verbose} displays detailed progress information during computation.

{phang}
{opt graph} displays visualization graphs including evolution, drift detection,
and frame size panels.

{phang}
{opt saving(string)} saves the combined graph to the specified filename.

{phang}
{opt saveref(string)} saves the final reference set (frame) to the specified file.

{phang}
{opt savehistory(string)} saves the frame membership history to the specified file.

{phang}
{opt summarysave(string)} saves a period-level monitoring dataset. The saved file
contains one row per value of {cmd:time()} and includes the number of observations,
mean and maximum score or hull distance, frame members, frame entries, pruned frame
members, drift counts, and gate counts when {cmd:gate()} or {cmd:mlscore()} is used.
This option is intended for Stata monitoring workflows in which users want a compact
frontier-update audit trail that can be graphed, merged, exported, or appended across
runs.


{marker models}{...}
{title:DEA Models}

{pstd}
{cmd:sdea} supports 12 DEA models, all implemented internally:

{p2colset 5 18 20 2}{...}
{p2col:{bf:Model}}Description{p_end}
{p2line}
{p2col:{bf:te}}Technical Efficiency (radial){p_end}
{p2col:{bf:cost}}Cost Efficiency (requires {opt w()}){p_end}
{p2col:{bf:revenue}}Revenue Efficiency (requires {opt p()}){p_end}
{p2col:{bf:ddf}}Directional Distance Function{p_end}
{p2col:{bf:sbm}}Slacks-Based Measure (non-radial){p_end}
{p2col:{bf:additive}}Additive Model{p_end}
{p2col:{bf:supereff}}Super-efficiency (radial, orientation-dependent){p_end}
{p2col:{bf:supersbm}}Super-SBM (Tone 2002){p_end}
{p2col:{bf:allocative}}Allocative Efficiency (requires {opt w()} or {opt p()}){p_end}
{p2col:{bf:multiplier}}Multiplier (dual) formulation{p_end}
{p2col:{bf:fdh}}Free Disposal Hull{p_end}
{p2col:{bf:imprecise}}Imprecise DEA with bounded data{p_end}
{p2line}

{pstd}
{bf:Efficiency Score Interpretation:}

{p2colset 5 22 24 2}{...}
{p2col:{bf:Model}}Efficient Value{space 4}Range{p_end}
{p2line}
{p2col:{bf:te, cost, revenue}}1{space 14}(0, 1]{p_end}
{p2col:{bf:sbm}}1{space 14}(0, 1]{p_end}
{p2col:{bf:ddf, additive}}0{space 14}[0, inf){p_end}
{p2col:{bf:supereff (input)}}>=1{space 12}(0, inf){p_end}
{p2col:{bf:supereff (output)}}<=1{space 12}(0, inf){p_end}
{p2col:{bf:supersbm}}>=1{space 12}[1, inf){p_end}
{p2col:{bf:fdh}}1{space 14}(0, 1] or [1, inf){p_end}
{p2line}


{marker methods}{...}
{title:Streaming Methods}

{pstd}
{bf:BuildHull} - Maintains a VRS frame/reference set. The default
{bf:bhvariant(modified)} is a streaming-oriented modification in which each new
DMU is evaluated against the current frame, exterior arrivals are added, and
the {opt prune()} option can remove redundant frame members. The option
{bf:bhvariant(dula)} implements the original Dulá-style phase-1 frame search more
closely by using the dual separating hyperplane and hyperplane translation to
select the next frame element from the remaining candidate set. Cite Dulá (2011)
for the two-phase frame algorithm and Dulá and López (2013) for DEA with
streaming arrivals. See {help sdea##references:References}.

{pstd}
{bf:FullStream} - Self-inclusive evaluation. New DMU is added to frontier first,
then evaluated. Removed if inefficient.

{pstd}
{bf:Cumulative} - Expanding frontier. All historical DMUs form the reference set.
Frontier only grows, never shrinks.

{pstd}
{bf:Moving} - Moving window frontier. Reference set includes only DMUs within the
specified window of periods.

{pstd}
{bf:Standard} - Traditional cross-sectional DEA. Each period analyzed separately.

{pstd}
{bf:Time} - Time-based sequential. DMUs evaluated against all prior periods'
cumulative frontier.

{pstd}
{bf:Window} - Window DEA. Each DMU appears in multiple overlapping windows;
average efficiency reported.

{pstd}
{bf:Malmquist} - Malmquist Productivity Index with full decomposition:
{p_end}
{phang2}• TFPCH: Total Factor Productivity Change{p_end}
{phang2}• EFFCH: Efficiency Change (catch-up){p_end}
{phang2}• TECHCH: Technical Change (frontier shift){p_end}
{phang2}• PECH: Pure Efficiency Change{p_end}
{phang2}• SECH: Scale Efficiency Change{p_end}

{pstd}
{bf:TFDEA} - Technology Forecasting DEA. Evaluates future DMUs against historical
frontier using super-efficiency.


{marker results}{...}
{title:Stored Results}

{pstd}
{cmd:sdea} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}
{synopt:{cmd:r(N_efficient)}}number of efficient DMUs{p_end}
{synopt:{cmd:r(N_frame)}}number of frontier members{p_end}
{synopt:{cmd:r(refset)}}full or buildhull for standard radial TE{p_end}
{synopt:{cmd:r(N_candidates)}}static BuildHull candidates after dominance screening{p_end}
{synopt:{cmd:r(N_hull_ipm)}}IPM calls during frame construction (before coverage checks){p_end}
{synopt:{cmd:r(N_hull_simplex)}}simplex calls during construction, including retries with {cmd:solver(auto)}{p_end}
{synopt:{cmd:r(N_hull_failed)}}construction LPs unresolved after the requested solver policy{p_end}
{synopt:{cmd:r(N_uncertain)}}omitted candidates restored by the final coverage check{p_end}
{synopt:{cmd:r(N_full_retry)}}firms rescored against the full reference set{p_end}
{synopt:{cmd:r(hull_solver)}}requested static hull construction solver policy{p_end}
{synopt:{cmd:r(score_solver)}}radial score solver used by the static BuildHull path{p_end}
{synopt:{cmd:r(N_drifts)}}number of drift event periods detected{p_end}
{synopt:{cmd:r(N_gate_enveloped)}}number of BuildHull observations enveloped without update{p_end}
{synopt:{cmd:r(N_gate_updates)}}number of accepted event-triggered frontier updates{p_end}
{synopt:{cmd:r(N_gate_rejected)}}number of candidate updates rejected by gate rule{p_end}
{synopt:{cmd:r(N_gate_starters)}}number of starter observations used to initialize frame{p_end}
{synopt:{cmd:r(mean_eff)}}mean efficiency score{p_end}
{synopt:{cmd:r(window)}}drift detection window size{p_end}
{synopt:{cmd:r(thr_volatility)}}volatility threshold used{p_end}
{synopt:{cmd:r(thr_expansion)}}expansion threshold used{p_end}
{synopt:{cmd:r(thr_turnover)}}turnover threshold used{p_end}
{synopt:{cmd:r(n_vol_triggered)}}periods with volatility triggered{p_end}
{synopt:{cmd:r(n_exp_triggered)}}periods with expansion triggered{p_end}
{synopt:{cmd:r(n_turn_triggered)}}periods with turnover triggered{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(method)}}streaming method used{p_end}
{synopt:{cmd:r(model)}}DEA model used{p_end}
{synopt:{cmd:r(orientation)}}orientation (input/output){p_end}
{synopt:{cmd:r(rts)}}returns to scale (vrs/crs){p_end}
{synopt:{cmd:r(algo_name)}}drift detection algorithm name{p_end}

{pstd}
{bf:Generated Variables:}

{p2colset 5 25 27 2}{...}
{p2col:{bf:Variable}}Description{p_end}
{p2line}
{p2col:{bf:_sdea_alpha}}Efficiency score for standard DEA methods; BuildHull alpha/hull distance for {bf:method(buildhull)}{p_end}
{p2col:{bf:_sdea_frame}}Frontier membership indicator (0/1){p_end}
{p2col:{bf:_sdea_drift}}Drift diagnostic value or magnitude (if {opt drifttest} specified){p_end}
{p2col:{bf:_sdea_frame_size}}BuildHull frame size after each arrival{p_end}
{p2col:{bf:_sdea_frame_entry}}Indicator that an arriving observation entered the BuildHull frame{p_end}
{p2col:{bf:_sdea_frame_pruned}}Number of BuildHull frame members pruned at that arrival{p_end}
{p2col:{bf:_sdea_drift_V}}Multi-indicator volatility diagnostic{p_end}
{p2col:{bf:_sdea_drift_D}}Multi-indicator frontier-expansion diagnostic{p_end}
{p2col:{bf:_sdea_drift_T}}Multi-indicator frame-turnover diagnostic{p_end}
{p2col:{bf:_sdea_drift_frontier}}Conditional frontier-shift diagnostic for {bf:driftalgo(frontier)}{p_end}
{p2col:{bf:gate(newvar)}}Frontier update status: 0=enveloped, 1=update, 2=rejected, 3=starter{p_end}
{p2col:{bf:_tfpch}}Total Factor Productivity Change (Malmquist){p_end}
{p2col:{bf:_effch}}Efficiency Change (Malmquist){p_end}
{p2col:{bf:_techch}}Technical Change (Malmquist){p_end}
{p2col:{bf:_pech}}Pure Efficiency Change (Malmquist){p_end}
{p2col:{bf:_sech}}Scale Efficiency Change (Malmquist){p_end}
{p2line}


{marker examples}{...}
{title:Examples}

{pstd}Setup: Load Swedish pharmacy data{p_end}
{phang2}{cmd:. use pharmacy_yearly_1980_1989.dta, clear}{p_end}
{phang2}{cmd:. xtset dmu year}{p_end}

{pstd}Basic BuildHull streaming DEA frame construction{p_end}
{phang2}{cmd:. sdea y1 y2 y3 y4, inputs(x1 x2 x3 x4) time(year) method(buildhull) replace}{p_end}

{pstd}BuildHull with Multi-indicator Drift Diagnostics{p_end}
{phang2}{cmd:. sdea y1 y2 y3 y4, inputs(x1 x2 x3 x4) time(year) dmu(dmu) ///}{p_end}
{phang2}{cmd:      method(buildhull) drifttest driftalgo(multi) ///}{p_end}
{phang2}{cmd:      thrvol(0.1) threxp(0.2) thrturn(0.3) graph replace}{p_end}

{pstd}BuildHull with Pruning for Large Datasets{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) method(buildhull) prune(50) replace}{p_end}

{pstd}BuildHull with simplex solver for comparison{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) method(buildhull) solver(simplex) replace}{p_end}

{pstd}BuildHull with frontier update gate{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) method(buildhull) ///}{p_end}
{phang2}{cmd:      gate(update_gate) gateeps(1e-8) replace}{p_end}

{pstd}BuildHull with a saved period-level monitoring table{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) method(buildhull) ///}{p_end}
{phang2}{cmd:      gate(update_gate) summarysave(monitoring_summary.dta) replace}{p_end}

{pstd}BuildHull with an externally supplied admissibility score{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) method(buildhull) ///}{p_end}
{phang2}{cmd:      gate(update_gate) mlscore(frontier_score) mlthreshold(0.5) replace}{p_end}

{pstd}Malmquist Productivity Index{p_end}
{phang2}{cmd:. sdea y1 y2 y3 y4, inputs(x1 x2 x3 x4) time(year) dmu(dmu) ///}{p_end}
{phang2}{cmd:      method(malmquist) rts(crs) orientation(output) verbose replace}{p_end}

{pstd}Compare Multi-indicator vs CUSUM drift diagnostics{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) method(buildhull) drifttest driftalgo(multi) replace}{p_end}
{phang2}{cmd:. local multi_drifts = r(N_drifts)}{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) method(buildhull) drifttest driftalgo(cusum) replace}{p_end}
{phang2}{cmd:. local cusum_drifts = r(N_drifts)}{p_end}
{phang2}{cmd:. di "Multi-indicator: `multi_drifts' vs CUSUM: `cusum_drifts'"}{p_end}

{pstd}Cost Efficiency with Input Prices{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) model(cost) w(w1 w2) replace}{p_end}

{pstd}Window DEA with 3-period window{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) method(window) window(3) replace}{p_end}

{pstd}TFDEA with frontier cutoff at period 7{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) method(tfdea) tf(7) replace}{p_end}

{pstd}Custom variable names{p_end}
{phang2}{cmd:. sdea y1 y2, inputs(x1 x2) time(year) generate(eff frontier drift) replace}{p_end}


{marker technical}{...}
{title:Technical Notes}

{pstd}
{bf:Dual Fallback Strategy}

{pstd}
Version 5.5.0 includes automatic primal-to-dual switching (dual fallback) to improve
numerical stability. When the primal (envelopment) LP fails to converge due to
cycling, degeneracy, or other numerical issues, the solver automatically attempts
the dual (multiplier) formulation.

{pstd}
By the Strong Duality Theorem, both formulations return identical efficiency values.
However, reference set weights (lambda values) are only available from the primal
formulation. If dual fallback is triggered:

{phang2}• Efficiency scores remain accurate{p_end}
{phang2}• Reference set (lambda) values will be missing (.){p_end}

{pstd}
Version 5.7.0 adds {cmd:summarysave()} for period-level monitoring summaries.
Version 5.6.0 adds an event-triggered frontier update gate for BuildHull and an
optional {cmd:mlscore()} interface for externally supplied rule-based or
machine-learning admissibility scores. Version 5.5.0 also updates
super-efficiency handling, implements the Tone (2002)
Super-SBM formulation, returns missing values for infeasible super-efficiency
cases, permits real-valued {opt tf()} cutoffs, and adds TDW-TFDEA aggregation via
{opt decay()} and {opt decaytype()}.

{pstd}
This is most relevant when using {opt saveref()} with large datasets or unusual
data structures.

{pstd}
If {cmd:sdea} is maintained as a published package, future code organization should
factor shared DEA model and solver routines more explicitly into the common
{cmd:dea} layer, while retaining in {cmd:sdea} only the streaming frame,
gate, drift, and time-indexed frontier-management code.


{marker references}{...}
{title:References}

{phang}
Banker, R. D., A. Charnes, and W. W. Cooper. 1984. Some models for estimating
technical and scale inefficiencies in data envelopment analysis.
{it:Management Science} 30: 1078-1092.

{phang}
Charnes, A., W. W. Cooper, and E. Rhodes. 1978. Measuring the efficiency of
decision making units. {it:European Journal of Operational Research} 2: 429-444.

{phang}
Dulá, J. H. 2011. An Algorithm for Data Envelopment Analysis.
{it:INFORMS Journal on Computing} 23(2), 284-296.
{browse "https://doi.org/10.1287/ijoc.1100.0400":doi:10.1287/ijoc.1100.0400}.

{phang}
Dulá, J. H., and F. J. López. 2013. DEA with streaming data.
{it:Omega} 41(1), 41-47.
{browse "https://doi.org/10.1016/j.omega.2011.07.010":doi:10.1016/j.omega.2011.07.010}.

{phang}
Fare, R., S. Grosskopf, B. Lindgren, and P. Roos. 1992. Productivity changes
in Swedish pharmacies 1980-1989: A non-parametric Malmquist approach.
{it:Journal of Productivity Analysis} 3: 85-101.

{phang}
Fare, R., S. Grosskopf, M. Norris, and Z. Zhang. 1994. Productivity growth,
technical progress, and efficiency change in industrialized countries.
{it:American Economic Review} 84: 66-83.

{phang}
Ji, Y., and C. Lee. 2010. Data envelopment analysis. {it:Stata Journal} 10: 267-280.

{phang}
Lee, C. 2022. {it:A Handbook of Data Envelopment Analysis using Stata}. KDP Press.

{phang}
Lee, C. 2025. {it:A Handbook of Data Envelopment Analysis using Stata}. 2nd ed. KDP Press.

{phang}
Lee, C., and N.-S. Cho. 2019. The matching problem in Stata. Presented at the 2019 Stata
Conference.

{phang}
Lee, C., and B.-K. Jung. 2019. Technology forecasting using data envelopment analysis in Stata.
Presented at the Oceania Stata Conference 2019.

{phang}
Lee, C., K. R. Lee, and B. I. Leem. 2010. An efficient data envelopment analysis with large data
set in Stata. Presented at the 2010 Stata Conference, Boston.

{phang}
Lee, C., K. R. Lee, and B. I. Leem. 2011. Malmquist productivity index using DEA frontier in Stata.
Presented at the 2011 Stata Conference (CHI11).

{phang}
Lee, C., K. R. Lee, and B. I. Leem. 2012. Allocative efficiency analysis using DEA in Stata.
Presented at the 2012 Stata Conference (SAN12).

{phang}
Lee, C., K. R. Lee, and B. I. Leem. 2013. Mathematical optimization in Stata: LP and MILP.
Presented at the 2013 Stata Conference.

{phang}
Lim, D. J., T. R. Anderson, and T. Shott. 2014. Technological forecasting based
on segmented rate of change. {it:Technological Forecasting and Social Change}
85: 182-193.

{phang}
Tone, K. 2001. A slacks-based measure of efficiency in data envelopment analysis.
{it:European Journal of Operational Research} 130: 498-509.

{phang}
Tone, K. 2002. A slacks-based measure of super-efficiency in data envelopment analysis.
{it:European Journal of Operational Research} 143: 32-41.


{marker author}{...}
{title:Author}

{pstd}
Choonjoo Lee{break}
National Security College{break}
Korea National Defense University{break}
Nonsan-si, Republic of Korea{break}
bloom.rampike@gmail.com

{pstd}
For bug reports and feature requests, please contact the author.

{pstd}
Copyright (c) 2026 Choonjoo Lee.{break}
Licensed under the {browse "https://spdx.org/licenses/AFL-3.0.html":Academic Free License version 3.0} (AFL-3.0).


{title:Also See}

{psee}
Online: {helpb dea_lp} (if installed)
{p_end}
