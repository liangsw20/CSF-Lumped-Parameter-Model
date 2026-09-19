%% base_ref.m -- B-scenario reference values, read from the saved B run.
%
% Four runners (run_venr_series, run_a3_series, run_cref_series, run_s_series) need a
% set of B reference values to express their results as percentage changes.
% Those were hard-coded literals, so every time the model changed they silently
% went stale: after one revision the per-cycle swing baselines were off by 20-50%,
% which would have mis-stated the pulsatility columns of every sensitivity table
% without any error being raised.
%
% This function removes the failure mode by reading the values from the B run
% itself.  Call it AFTER running run_base_final, which writes
% Results/H_series/B.mat.
%
%   ref = base_ref('ref')
%       struct with the scalar references the O / A4 / Cref_series runners use, plus
%       the three per-cycle swings run_cref_series prints for context.
%
%   [K, S] = base_ref('KS', KEY, SW)
%       the reference vectors in the CALLER's own signal orders, for run_s_series
%       (KEY = steady-state means, SW = per-cycle swings).  A signal that is not
%       logged in the current model yields NaN, which the caller's delta table
%       then shows as an empty cell rather than a wrong number.
%
% Peak-to-peak swings are averaged over the last N cycles and steady-state means
% over the last W cycles, matching every runner and collate_all.
%
% ASCII only.

function varargout = base_ref(what, varargin)

CYCLE = 0.8;   % s, cardiac period (T in bp_compute)
W     = 400;   % cycles averaged for the steady-state means
N     = 100;   % cycles averaged for the peak-to-peak swings

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
f = fullfile(root,'SimulinkModel','Results','H_series','B.mat');
if ~isfile(f)
    error(['base_ref: %s not found.\n' ...
           'Run  matlab -batch "run_base_final"  first: the reference values are ' ...
           'read from that run so they cannot go stale.'], f);
end
S = load(f,'ls');  ls = S.ls;
has = @(s) any(strcmp(ls.getElementNames(), s));
M = @(s) ssm(ls, s, CYCLE, W);
P = @(s) swp(ls, s, CYCLE, N);

switch lower(what)

    case 'ref'
        r = struct();
        r.P_ECS    = M('P_ECS');
        r.P_VEN    = M('P_VEN');
        r.P_SAS    = M('P_SAS');
        r.grad     = r.P_VEN - r.P_SAS;
        r.V_VEN    = M('V_VEN');
        r.V_SPI    = M('V_SPI');
        r.V_ECS    = M('V_ECS');
        r.Q_out = M('Q_out');
        r.Q_VEN    = M('Q_VEN');
        r.Q_SPI    = M('Q_SPI');
        r.Q_SLS   = M('Q_SLS');
        r.Q_GS   = M('Q_GS');
        r.Q_BBB    = M('Q_BBB');
        r.Q_A      = M('Q_A');
        r.alloc    = 100*[M('Q_SSS') M('Q_OLS') M('Q_SLS')] / M('Q_out');
        % per-cycle swings, for the run_cref_series context line
        r.sw_V_SPI = P('V_SPI');
        r.sw_P_ECS = P('P_ECS');
        r.sw_Q_SPI = P('Q_SPI');
        r.sw_V_Vessel = P('V_Vessel');
        r.V_T      = M('V_T');      % logged since the 2026-09-16 revision
        varargout{1} = r;

    case 'ks'
        KEY = varargin{1};  SW = varargin{2};
        K = nan(1, numel(KEY));
        for q = 1:numel(KEY)
            if has(KEY{q}), K(q) = M(KEY{q}); end
        end
        Sv = nan(1, numel(SW));
        for q = 1:numel(SW)
            if has(SW{q}), Sv(q) = P(SW{q}); end
        end
        varargout{1} = K;  varargout{2} = Sv;

    otherwise
        error('base_ref: unknown request "%s" (use "ref" or "KS")', what);
end
end

% ---------------------------------------------------------------- helpers
function v = ssm(ls, name, CYCLE, W)
[y,t] = getYT(ls,name);
n = floor(t(end)/CYCLE);
v = mean(y(t > (n-W)*CYCLE));
end
function d = swp(ls, name, CYCLE, N)
[y,t] = getYT(ls,name);
n = floor(t(end)/CYCLE);  R = nan(N,1);
for c = (n-N+1):n
    m = (t > (c-1)*CYCLE) & (t <= c*CYCLE);
    if any(m), R(c-(n-N)) = max(y(m)) - min(y(m)); end
end
d = mean(R,'omitnan');
end
function [y,t] = getYT(ls, name)
e = ls.getElement(name);
if isa(e,'Simulink.SimulationData.Dataset'), ts = e.getElement(1).Values; else, ts = e.Values; end
y = double(ts.Data(:)); t = double(ts.Time(:));
end
