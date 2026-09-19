%% compare_runs.m -- what changed between two runs of the same scenario?
%
%   compare_runs('Results/H_series/B.mat', 'Batch/_B_pre_IEGH.mat')
%
% Compares the two logsout datasets signal by signal: sample count, time grid,
% and the mean over the last W cycles plus the per-cycle swing.  Used to quantify
% the effect of a parameter change (e.g. R_IEGH 0.64 -> 1.14) without re-running
% the whole batch, and to prove that a pure rename changed nothing.
%
% ASCII only.  matlab -batch "compare_runs('a.mat','b.mat')"
function compare_runs(fNew, fOld)
if nargin < 2 || isempty(fOld)
    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    fNew = fullfile(root,'SimulinkModel','Results','H_series','B.mat');
    fOld = fullfile(here,'_B_pre_IEGH.mat');
end
N = load(fNew); O = load(fOld);
nA = N.ls.getElementNames(); nB = O.ls.getElementNames();
CYCLE = 0.8; W = 400;

fprintf('=== compare_runs ===\n  new: %s\n  old: %s\n', fNew, fOld);
fprintf('%d vs %d signals\n\n', numel(nA), numel(nB));
missing = setdiff(nB, nA);
if ~isempty(missing), fprintf('*** only in the old run: %s\n', strjoin(missing, ', ')); end

fprintf('%-14s %14s %14s %10s   %14s %14s %10s\n', ...
        'signal','new mean','old mean','delta %%','new swing','old swing','delta %%');
for k = 1:numel(nA)
    s = nA{k};
    if ~any(strcmp(nB, s)), continue; end
    a = local_ts(N.ls, s);  b = local_ts(O.ls, s);
    ta = double(a.Time(:));  tb = double(b.Time(:));
    va = double(a.Data(:));  vb = double(b.Data(:));
    na = floor(ta(end)/CYCLE); nb = floor(tb(end)/CYCLE);
    ma = ta > (na-W)*CYCLE;    mb = tb > (nb-W)*CYCLE;
    mA = mean(va(ma));         mB = mean(vb(mb));
    sA = local_swing(ta, va, na); sB = local_swing(tb, vb, nb);
    d1 = 100*(mA-mB)/max(abs(mB),eps);
    d2 = 100*(sA-sB)/max(abs(sB),eps);
    flag = '';
    if abs(d1) > 1 || abs(d2) > 1, flag = '  <<<'; end
    fprintf('%-14s %14.6g %14.6g %+9.2f%%   %14.6g %14.6g %+9.2f%%%s\n', ...
            s, mA, mB, d1, sA, sB, d2, flag);
end
fprintf('\n=== DONE ===\n');
end
% =====================================================================
function ts = local_ts(ls, name)
e = ls.getElement(name);
if isa(e,'Simulink.SimulationData.Dataset'), ts = e.getElement(1).Values;
else, ts = e.Values; end
end
function R = local_swing(t, v, n)
R = nan(1,100);
for c = (n-99):n
    m = (t > (c-1)*0.8) & (t <= c*0.8);
    if any(m), R(c-(n-100)) = max(v(m)) - min(v(m)); end
end
R = mean(R,'omitnan');
end
