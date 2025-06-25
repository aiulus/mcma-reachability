function [completed, results, R_id, R_val] = flexGrayBoxConform(varargin)
% flexGrayBoxConform  –  Reach‑set conformance synthesis for **gray‑box** models.
%
%     • "graySim" – simultaneous optimisation of system parameters *and*
%                   uncertainty sets (outer‑loop fmincon).
%     • "graySeq" – sequential optimisation (first parameters, then sets).
%     • "grayLS"  – least‑squares cost instead of max/∞‑norm.
%
%   Name–Value pairs (all optional):
%     'dynamics'      – char/string, system identifier for loadDynamics()  (default "Square")
%     'testSuites'    – 1×3 cell {id , train , val}                        (default: auto‑generate)
%     'plot_settings' – struct as in getConfig()                           (default: cfg.plot_settings)
%     'sysparams'     – struct forwarded to custom_loadDynamics            (default: empty)
%     'grayAlg'       – "graySim" | "graySeq" | "grayLS"                   (default "graySeq")
%
%   Outputs:
%     • completed – logical flag (true if everything ran)
%     • results   – cell {trueModel , grayModel}
%     • R_id      – reach‑sets for identification data
%     • R_val     – reach‑sets for validation data
%
%   --------------------------------------------------------------------
%   Example (sequential identification):
%       flexGrayBoxConform('dynamics','platoon','grayAlg','graySeq');
%
%   --------------------------------------------------------------------
%   © 2025 – Adapted from:
%     • contDynamics/@contDynamics/conform.m
%     • contDynamics/@contDynamics/private/priv_conform_gray.m
% ---------------------------------------------------------------------

%% 0)  Parse options ------------------------------------------------------

p = inputParser;
addParameter(p,'plot_settings',[]);
addParameter(p,'dynamics',"Square");
addParameter(p,'testSuites',[]);
addParameter(p,'sysparams',struct());
addParameter(p,'grayAlg',"graySeq");
parse(p,varargin{:});
plot_settings = p.Results.plot_settings;
DYN           = p.Results.dynamics;
TS_in         = p.Results.testSuites;
sysparams     = p.Results.sysparams;
modes = {'graySim','graySeq','grayLS'};
grayAlg = validatestring(p.Results.grayAlg, modes); 

% Obtain global config (single source‑of‑truth)
if isfield(sysparams,'cfg')
    cfg = sysparams.cfg;
else
    cfg = getConfig();
end
settings   = cfg.settings;
optsReach  = cfg.options_reach;

if isempty(plot_settings)
    plot_settings = cfg.plot_settings;
end

% Tunable knobs (override here or expose via new NVP):
cost_norm   = "interval";   % "interval" | "frob" – only relevant inside conform_white sub‑call
constraints = "half";       % "half" | "gen"

%% 1)  Load system --------------------------------------------------------

[sys, params_true.R0, params_true.U, params_true.p_true] = custom_loadDynamics(DYN,"rand",sysparams);
params_true.tFinal = sys.dt * settings.n_k - sys.dt;

%% 2)  Build / import test suites ----------------------------------------

if isempty(TS_in)
    params_true.testSuite       = createTestSuite(sys, params_true, settings.n_k      , settings.n_m      , settings.n_s      , cfg.options_testS);
    params_true.testSuite_train = createTestSuite(sys, params_true, settings.n_k_train, settings.n_m_train, settings.n_s_train, cfg.options_testS);
    params_true.testSuite_val   = createTestSuite(sys, params_true, settings.n_k_val  , settings.n_m_val  , settings.n_s_val  , cfg.options_testS);
else
    params_true.testSuite       = TS_in{1};
    params_true.testSuite_train = TS_in{2};
    params_true.testSuite_val   = TS_in{3};
end

%% 3)  Assemble solver options ------------------------------------------

options               = getConformanceOptions(optsReach, cost_norm, constraints, sys);
options.cs.verbose     = false;
options.cs.cost        = cost_norm;
options.cs.constraints = constraints;

% Gray‑box specific: initial parameter vector p0 and setter
p0_centres = [center(params_true.R0); center(params_true.U)];
options.cs.p0    = p0_centres;
options.cs.cp_lim = inf;
options.cs.p_min = p0_centres - options.cs.cp_lim;
options.cs.p_max = p0_centres + options.cs.cp_lim;
options.cs.set_p = @(p,params) set_p_centers(p, params, sys);

% Timeout safeguard (optional)
options.cs.timeout = 600;  % seconds – can be overridden by caller

% Weighting vector (for LS / max error)
options.cs.w      = ones(1, settings.n_k);

%% 4)  Initial disturbance‑set guesses -----------------------------------

c_R0 = zeros(size(center(params_true.R0)));
c_U  = zeros(size(center(params_true.U)));
params_id_init       = params_true;
params_id_init.R0    = zonotope([c_R0 eye(sys.nrOfStates)]);
params_id_init.U     = zonotope([c_U  eye(sys.nrOfInputs)]);

%% 5)  GRAY‑BOX identification -------------------------------------------

fprintf("\n[flexGrayBoxConform]  Identification using %s …\n", grayAlg);
tmr = tic;
[params_gray, resGray] = conform(sys, params_id_init, options, grayAlg); 
T_ident = toc(tmr);
fprintf("      finished in %.2f s\n", T_ident);

% Bundle outputs for downstream analysis
results          = cell(2,1);
results{1}       = struct('sys',sys, 'params',params_true , 'options',optsReach, 'name','true');
results{2}       = struct('sys',resGray.sys, 'params',params_gray, 'options',optsReach, 'name',grayAlg);

%% 6)  Reachability on identification data -------------------------------

R_id = cell(size(params_true.testSuite));
for m = 1:numel(params_true.testSuite)
    [R_id{m}, ~] = validateReach(params_true.testSuite{m}, results, 1); 
end

%% 7)  Reachability on validation data -----------------------------------

TS_val = params_true.testSuite_val;
R_val  = cell(size(TS_val));
for m = 1:numel(TS_val)
    [R_val{m}, ~] = validateReach(TS_val{m}, results, 1, plot_settings);
end

% Optional containment stats
validateReachableSets(TS_val, results, settings.n_k_val, {"true",grayAlg}, ...
    'plot_settings', plot_settings, 'label', 'VALIDATION DATA');

completed = true;

end  

%% -----------------------------------------------------------------------
% Helper – default parameter setter (centre‑vector update only)
% ------------------------------------------------------------------------
function [sys_out, params_out] = set_p_centers(p, params_in, sys_in)
    c_R0 = p(1:dim(params_in.R0));
    c_U  = p(dim(params_in.R0)+1:end);
    params_out      = params_in;
    params_out.R0   = zonotope(c_R0, generators(params_in.R0));
    params_out.U    = zonotope(c_U , generators(params_in.U));
    sys_out         = sys_in;  % no structural change to the model
end

function validateReachableSets(testSuite, configs, n_k_val, methods, varargin)
%  Simple containment statistics + optional plotting
%
%  INPUTS (same as call site)
%    testSuite   – cell array of testCase objects (validation data)
%    configs     – cell array with fields {sys,params,options,name}
%    n_k_val     – # of timesteps in each validation test case
%    methods     – cellstr or string array naming each config (order must
%                  match `configs`)
%
%  OPTIONAL name–value pairs:
%      'plot_settings'        – struct (same layout as in validateReach)
%      'label'                – char   (header shown in console)
%      'require_full_containment' – logical (default = false)
% -------------------------------------------------------------------------

    p = inputParser;
    addOptional(p,'plot_settings',[]);
    addOptional(p,'label','VALIDATION DATA');
    addOptional(p,'require_full_containment',false);
    parse(p,varargin{:});
    plot_settings          = p.Results.plot_settings;
    label                  = p.Results.label;
    require_full_contain   = p.Results.require_full_containment;

    num_out = zeros(numel(methods),1);
    num_in  = zeros(numel(methods),1);

    for m = 1:numel(testSuite)
        if isempty(plot_settings)
            [~, eval] = validateReach(testSuite{m}, configs, true);
        else
            [~, eval] = validateReach(testSuite{m}, configs, true, plot_settings);
        end
        num_out = num_out + eval.num_out;
        num_in  = num_in  + eval.num_in;
    end

    samplesTotal = numel(testSuite) * n_k_val * size(testSuite{1}.y,3);
    fprintf("\n%s:\n", label);
    for i = 1:numel(methods)
        pctContained = 100 * num_in(i) / (num_in(i)+num_out(i));
        notValid     = 100 * (samplesTotal - (num_out(i)+num_in(i))) / samplesTotal;

        if require_full_contain
            suffix = " (must be 100 %)";
        else
            suffix = "";
        end

        fprintf("  %-10s: %.2f %% samples contained%s\n", methods{i}, pctContained, suffix);
        fprintf("              %.2f %% samples invalid\n", notValid);
    end
end
