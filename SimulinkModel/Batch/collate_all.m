%% collate_all.m -- one master table over every scenario .mat in Results/
%
% Every batch runner saves the whole logsout Dataset (-v7.3), so each file
% carries all 46 logged signals regardless of what that runner printed.  This
% script walks Results/, computes the same steady-state statistics the runners
% use (cycle mean over the last W cardiac cycles), and writes ONE table with a
% row per scenario.  It replaces nothing: the per-series compare_*.m scripts
% still make their own figures, but every number in the write-ups can come from
% here, and a stale file is obvious from its timestamp column.
%
% The 2026-09-15 model revision renamed/dropped three signals:
%   P_CC   -> P_C
%   V_PNS -> no longer logged, but still a compartment: recovered here from
%             the Monro-Kellie closure V_PNS = 1500 - (sum of the seven
%             logged cranial compartments), which was verified to hold exactly.
%   V_All  -> no longer logged; the craniospinal total is conserved, so its
%             only dynamics are dV_SPI and it carries no extra information.
%
% ASCII only.  matlab -batch "collate_all"

fprintf('=== collate_all (MATLAB %s) ===\n', version);
tAll = tic;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
resDir = fullfile(root,'SimulinkModel','Results');

CYCLE   = 0.8;   % s, cardiac period (T in bp_compute)
W       = 400;   % cycles averaged, same as the runners
V_T0  = 1500;  % mL, invariant cranial compartment volume

MEAN = {'P_ECS','P_SAS','P_SPI','P_VEN','P_PAS','P_PNS','P_A','P_C','P_V','E_A', ...
        'V_SPI','V_VEN','V_ECS','V_SAS','V_PAS','V_brain','V_Vessel','V_ICF','V_ICS', ...
        'Q_SPI','Q_SAS','Q_SLS','Q_SSS','Q_OLS','Q_out','Q_GS','Q_ECS', ...
        'Q_BBB','Q_A','Q_VEN','Q_PAS','Q_SASE','Q_PNS','R_IEG','dV_Vessel'};
SWING = {'P_ECS','P_SAS','P_VEN','P_A','V_SPI','V_Vessel','V_VEN','V_ECS','V_SAS', ...
         'Q_SPI','Q_A'};
VOL7  = {'V_A','V_C','V_V','V_VEN','V_SAS','V_PAS','V_brain'};

files = dir(fullfile(resDir,'**','*.mat'));
files = files(~startsWith({files.name},'~'));
fprintf('found %d .mat files under Results/\n\n', numel(files));

rows = {};   % each row: a struct
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    series = relseries(resDir, files(k).folder);
    code   = erase(files(k).name,'.mat');
    fprintf('[%2d/%2d] %-28s %s', k, numel(files), code, series);
    try
        S  = load(f,'ls');
        ls = S.ls;
        nms = ls.getElementNames();
        has = @(s) any(strcmp(nms,s));
        R = struct();
        R.series = series; R.code = code; R.file = f;
        R.mat_time = char(datetime(files(k).datenum,'ConvertFrom','datenum', ...
                        'Format','yyyy-MM-dd HH:mm:ss'));
        tEnd = NaN;
        for q = 1:numel(MEAN)
            s = MEAN{q};
            if ~has(s), continue; end
            [v,t] = rd(ls,s);
            n = floor(t(end)/CYCLE);
            m = t > (n-W)*CYCLE;
            if ~any(m), m = true(size(t)); end
            R.(s) = mean(v(m));
            tEnd = max(tEnd, t(end));
        end
        for q = 1:numel(SWING)
            s = SWING{q};
            if ~has(s), continue; end
            [v,t] = rd(ls,s);
            n = floor(t(end)/CYCLE); n100 = min(100,n-1);
            cy = nan(max(n100,1),1);
            for c = (n-n100+1):n
                m2 = (t > (c-1)*CYCLE) & (t <= c*CYCLE);
                if any(m2), cy(c-(n-n100)) = max(v(m2)) - min(v(m2)); end
            end
            R.(['sw_' s]) = mean(cy,'omitnan');
        end
        R.StopTime = tEnd;
        R.ncycles  = floor(tEnd/CYCLE);
        % derived: cranial closure and outflow split
        if all(cellfun(@(s) isfield(R,s), VOL7))
            cr = 0; for q = 1:numel(VOL7), cr = cr + R.(VOL7{q}); end
            R.V_vPVS_derived = V_T0 - cr;
        end
        if isfield(R,'Q_out') && R.Q_out ~= 0
            for p = {'Q_OLS','Q_SLS','Q_SSS'}
                if isfield(R,p{1}), R.(['pct_' p{1}]) = 100*R.(p{1})/R.Q_out; end
            end
        end
        if isfield(R,'P_VEN') && isfield(R,'P_SAS'), R.grad_VEN_SAS = R.P_VEN - R.P_SAS; end
        if isfield(R,'P_ECS') && isfield(R,'P_SPI'),  R.grad_ECS_SPI = R.P_ECS - R.P_SPI; end
        rows{end+1} = R; %#ok<SAGROW>
        fprintf('  ok  P_ECS=%.4f  V_SPI=%.3f  Q_out=%.5f\n', ...
            getf(R,'P_ECS'), getf(R,'V_SPI'), getf(R,'Q_out'));
    catch ME
        if contains(ME.message,'ls')
            % CentralWaveforms/*.mat are not logsout files (they hold harmonic
            % fits), so they have no 'ls' variable and are not scenarios.
            fprintf('  skipped (not a logsout file)\n');
        else
            fprintf('  *** FAILED: %s\n', ME.message);
        end
    end
end

%% ---------------- write the master CSV --------------------------------
if isempty(rows)
    fprintf('\nno data -- nothing written\n'); return;
end
fields = {};
for k = 1:numel(rows)
    fields = union(fields, fieldnames(rows{k}));
end
% fieldnames/union both return COLUMN cell arrays; normalise to rows so the
% concatenation below cannot mismatch orientation.
fields = fields(:)';
% stable, readable column order: identity, time, then signals alphabetically
ident  = {'series','code','mat_time','StopTime','ncycles'};
ident  = ident(:)';
fields = [ident, setdiff(fields, ident)];
% drop the absolute file path from the CSV (keep it only in the printout)
fields = setdiff(fields, {'file'}, 'stable');
fields = fields(:)';

outCsv = fullfile(resDir,'collate_all.csv');
fid = fopen(outCsv,'w');
fprintf(fid,'%s', strjoin(fields,','));
fprintf(fid,'\n');
for k = 1:numel(rows)
    for q = 1:numel(fields)
        v = getf(rows{k}, fields{q});
        if ischar(v)
            fprintf(fid,'%s', v);
        elseif isnan(v)
            fprintf(fid,'');
        else
            fprintf(fid,'%.6g', v);
        end
        if q < numel(fields), fprintf(fid,','); end
    end
    fprintf(fid,'\n');
end
fclose(fid);
fprintf('\nwrote %s  (%d scenarios x %d columns)\n', outCsv, numel(rows), numel(fields));

%% ---------------- per-series overview ---------------------------------
% cellfun rather than {rows{1:end}.series}: chaining further indexing onto a
% comma-separated list is illegal ("Intermediate '{}' indexing produced a
% comma-separated list ... but must produce a single value when followed by
% another indexing operation"), which is why the first version wrote the CSV and
% then died on this line.
allseries = unique(cellfun(@(r) r.series, rows, 'UniformOutput', false));
for q = 1:numel(allseries)
    sel = rows(cellfun(@(r) strcmp(r.series, allseries{q}), rows));
    fprintf('\n===== %s (%d) =====\n', allseries{q}, numel(sel));
    fprintf('%-22s %10s %10s %10s %12s %10s %8s\n', ...
        'code','P_ECS','P_VEN','V_SPI','Q_out','V_Vessel','sw_V_SPI');
    for k = 1:numel(sel)
        fprintf('%-22s %10.4f %10.4f %10.4f %12.5f %10.4f %8.4f\n', ...
            sel{k}.code, getf(sel{k},'P_ECS'), getf(sel{k},'P_VEN'), ...
            getf(sel{k},'V_SPI'), getf(sel{k},'Q_out'), ...
            getf(sel{k},'V_Vessel'), getf(sel{k},'sw_V_SPI'));
    end
    if isfield(sel{1},'pct_Q_OLS')
        fprintf('  outflow split (lymph / spinal / AG):');
        for k = 1:numel(sel)
            fprintf(' %s %.2f/%.2f/%.2f', sel{k}.code, ...
                getf(sel{k},'pct_Q_OLS'), getf(sel{k},'pct_Q_SLS'), getf(sel{k},'pct_Q_SSS'));
        end
        fprintf('\n');
    end
end

fprintf('\n=== total wall clock %.1f s ===\n=== DONE ===\n', toc(tAll));

% ---------------------------------------------------------------- helpers
function s = relseries(resDir, folder)
s = strrep(folder, [resDir filesep], '');
if isempty(s), s = '(root)'; end
end

function v = getf(R, name)
if isfield(R,name), v = R.(name); else, v = NaN; end
end

function [v,t] = rd(ls,name)
e = ls.getElement(name);
if isa(e,'Simulink.SimulationData.Dataset'), ts = e.getElement(1).Values; else, ts = e.Values; end
v = double(ts.Data(:));
t = double(ts.Time(:));
end
