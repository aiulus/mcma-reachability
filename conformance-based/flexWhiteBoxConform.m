function [completed, results, R_id, R_val] = flexWhiteBoxConform(varargin)
% flexWhiteBoxConform  –  Reach‑set conformance identification for **white‑box** models.
%
%   Name–Value pairs:
%     'dynamics'      – String   : system identifier for loadDynamics()   (default "Square")
%     'testSuites'    – 1×3 cell : {train+id , train , val} pre‑computed  (default = auto‑generate)
%     'plot_settings' – struct   : plotting flags (see getConfig())        (default = cfg.plot_settings)
%     'sysparams'     – struct   : extra args forwarded to custom_loadDynamics (default = struct())
%
%   Outputs:
%     completed – logical flag (true if run finished)
%     results   – cell array with entries {trueModel , identifiedModel}
%     R_id      – reach‑sets for identification test data
%     R_val     – reach‑sets for validation data
%
%   --------------------------------------------------------------------
%   Example (minimal):
%
%       [~, res] = flexWhiteBoxConform('dynamics','platoon');
%
%   --------------------------------------------------------------------
%   © 2025 – Adapted from:
%     • contDynamics/@contDynamics/conform.m
%     • contDynamics/@contDynamics/private/priv_conform_white.m
% ----------------------------------------------------------------------

p = inputParser;
addParameter(p,'plot_settings',[]);
addParameter(p,'dynamics',"Square");
addParameter(p,'testSuites',[]);
addParameter(p,'sysparams',struct());
parse(p,varargin{:});
plot_settings = p.Results.plot_settings;
DYN = p.Results.dynamics;
TS_in = p.Results.testSuites;
sysparams = p.Results.sysparams;

% Global configuration (single source of truth)
if isfield(sysparams,'cfg')
    cfg = sysparams.cfg;
else
    cfg = getConfig();
end
settings = cfg.settings;
opts_reach = cfg.options_reach;

% Fall‑back to default plot settings
if isempty(plot_settings)
    plot_settings = cfg.plot_settings;
end

% Hard‑coded conformance knobs (override by editing or extending signature)
cost_norm   = "interval";   % "interval" | "frob"
constraints = "half";       % "half"     | "gen"

% ---------------------------------------------------------------------
% 1)  Load / build system
% ---------------------------------------------------------------------

if isfield(sysparams,'cfg'); end  %# (placeholder – keeps linter happy)

[sys, params_true.R0, params_true.U, params_true.p_true] = custom_loadDynamics(DYN,"rand",sysparams);
params_true.tFinal = sys.dt * settings.n_k - sys.dt;

% ---------------------------------------------------------------------
% 2)  Prepare test suites (generate or use provided)
% ---------------------------------------------------------------------

if isempty(TS_in)
    % auto‑generate via CORA createTestSuite
    params_true.testSuite = createTestSuite(sys, params_true, settings.n_k, settings.n_m, settings.n_s, cfg.options_testS);
    params_true.testSuite_train = createTestSuite(sys, params_true, settings.n_k_train, settings.n_m_train, settings.n_s_train, cfg.options_testS);
    params_true.testSuite_val  = createTestSuite(sys, params_true, settings.n_k_val, settings.n_m_val, settings.n_s_val, cfg.options_testS);
else
    params_true.testSuite = TS_in{1};
    params_true.testSuite_train = TS_in{2};
    params_true.testSuite_val = TS_in{3};
end

% ---------------------------------------------------------------------
% 3)  Build option struct for conformance solver
% ---------------------------------------------------------------------

options = getConformanceOptions(opts_reach, cost_norm, constraints, sys);
options.cs.verbose = false;            % suppress solver chatter
options.cs.cost = cost_norm;        % propagate cost norm
options.cs.constraints = constraints;      % propagate constraint type

% ---------------------------------------------------------------------
% 4)  Initial guesses for uncertainty sets (centred unit zonotopes)
% ---------------------------------------------------------------------

c_R0 = zeros(size(center(params_true.R0)));
c_U = zeros(size(center(params_true.U)));
params_id_init = params_true;
params_id_init.R0 = zonotope([c_R0 eye(sys.nrOfStates)]);
params_id_init.U = zonotope([c_U  eye(sys.nrOfInputs)]);

% ---------------------------------------------------------------------
% 5)  **WHITE‑BOX** conformance identification
% ---------------------------------------------------------------------

fprintf("\n[flexWhiteBoxConform]  Identification (white‑box) …\n");
tmr = tic;
[params_white, ~] = conform(sys, params_id_init, options, "white");
T_ident = toc(tmr);
fprintf("finished in %.2f s\n", T_ident);

% Collect configs for downstream reach/validation
results = cell(2,1);
results{1}.sys = sys;
results{1}.params = params_true;
results{1}.options = opts_reach;
results{1}.name = "true";

results{2}.sys = sys;
results{2}.params = params_white;
results{2}.options = opts_reach;
results{2}.name = "white";

% ---------------------------------------------------------------------
% 6)  Reachability on identification data
% ---------------------------------------------------------------------

R_id = cell(size(params_true.testSuite));
for m = 1:numel(params_true.testSuite)
    [R_id{m}, ~] = validateReach(params_true.testSuite{m}, results, 1); 
end

% ---------------------------------------------------------------------
% 7)  Reachability on *validation* data (optional)
% ---------------------------------------------------------------------

TS_val = params_true.testSuite_val;
R_val  = cell(size(TS_val));
for m = 1:numel(TS_val)
    [R_val{m}, ~] = validateReach(TS_val{m}, results, 1, plot_settings);
end

% ---------------------------------------------------------------------
% 8)  Optional containment statistics (console only)
% ---------------------------------------------------------------------

validateReachableSets(TS_val, results, settings.n_k_val, {"true","white"}, ...
    'plot_settings', plot_settings, 'label', 'VALIDATION DATA');

% ---------------------------------------------------------------------
completed = true;

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
