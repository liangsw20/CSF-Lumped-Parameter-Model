%% run_all_scenarios.m -- reproduce every scenario reported in Tables 2 and 3.
%
% Run from anywhere; the paths below are resolved from this file's location, so
% the folder containing this file only has to sit next to SimulinkModel/.
%
% ORDER MATTERS: the baseline B is simulated first, because every other runner
% reads its reference values from Results/H_series/B.mat through base_ref.m.
% The whole batch takes roughly 7 h on a current desktop; each runner prints its
% own wall clock and can be run on its own.
%
% ASCII only.
here = fileparts(mfilename('fullpath'));
root = here;
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(fullfile(root,'SimulinkModel','Batch'));
addpath(fullfile(root,'SimulinkModel'));

fprintf('=== run_all_scenarios ===\n');
fprintf('root: %s\n', root);

% quick pre-flight: every scenario code must produce its expected parameters
check_codes;

groups = {'run_base_final', ...   % A0  baseline            ~6 min
          'run_h_series',   ...   % A1  hypertension        ~25 min
          'run_c_series',   ...   % A2  outflow impairment  ~100 min
          'run_cref_series',...   % A2-ref single-outflow   ~10 min
          'run_a3_series',  ...   % A3  combined pathology  ~50 min
          'run_s_series'};        % Table 3 sensitivity     ~110 min

tAll = tic;
for k = 1:numel(groups)
    fprintf('\n##############################################################\n');
    fprintf('###  %s   ( %d / %d )\n', groups{k}, k, numel(groups));
    fprintf('##############################################################\n');
    tG = tic;
    feval(groups{k});
    fprintf('\n###  %s done in %.1f min\n', groups{k}, toc(tG)/60);
end

% master table: one row per scenario (Results/collate_all.csv)
collate_all;

fprintf('\n=== run_all_scenarios: all groups done in %.2f h ===\n', toc(tAll)/3600);
