%% run_s_series.m -- B-series one-at-a-time sensitivity runs
%
% Every scenario is "B + exactly one parameter change" (see check_s_codes.m
% for the parameter-level verification).  Combine with '+' to change the base,
% e.g. 'S_bEA_0+H_G2A'.
%
% StopTime: R_eff = Rout_eff + (R_VEN - R_VEN_base), P_eq from the outflow
% balance, C = k_m/P_eq, tau = 0.78*R_eff*C, StopTime = max(3000,ceil(6*tau*60/500)*500).
%
% NON-DESTRUCTIVE: no save_system, no logging reconfiguration.
% ASCII only.  matlab -batch "run_s_series"
%   optional selection:  S_CODES={'S_PVI_30'}; run_s_series
%
% NOTE: variable names i / S / fn are avoided -- Brain_Parameters clears them.

fprintf('=== run_s_series (MATLAB %s) ===\n', version);
tAll = tic;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
mdl  = 'CSF_Model';
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(here); addpath(fullfile(root,'SimulinkModel'));
outDir = fullfile(root,'SimulinkModel','Results','S_series');
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~bdIsLoaded(mdl), load_system(fullfile(root,'SimulinkModel',[mdl '.slx'])); end
CYCLE = 0.8;

codes = { ...
  'S_WAVE_RADIAL','S_WAVE_CENTRAL','S_WAVE_AMP120','S_WAVE_AMP080', ...
  'S_PVI_20','S_PVI_30', ...
  'S_HR_60','S_HR_90', ...
  'S_RFM_0.002','S_RFM_0.05', ...
  'S_PSPI0_7.5','S_PSPI0_10.5', ...
  'S_RBBB_0.25','S_RBBB_4', ...
  'S_bEA_0','S_bEA_0.02','S_bEA_0+H_G2A', ...
  'S_SWAP_AGMAX','S_VALVE_OFF','S_VALVE_HIGH', ...
  'S_RECS_0.5','S_RECS_2'};
if exist('S_CODES','var') && ~isempty(S_CODES), codes = S_CODES; end

% ---- B reference values, refreshed 2026-09-15 for the current model ----
% The revision changed E_VEN (13.3 -> 3) and E_SAS (9.5 -> 17.1).  The
% steady-state MEANS moved <0.4%, but the per-cycle SWINGS moved 20-50%, so the
% swing rows below are the ones that mattered.  Revision-pinned literals; the
% sensitivity deltas are all relative to these.
KEY = {'P_ECS','P_VEN','P_SAS','P_SPI','V_VEN','V_SPI','V_ECS','V_All', ...
       'V_A','V_C','V_V','V_Vessel','V_SAS','V_PAS','V_brain','Q_VEN','Q_SPI', ...
       'Q_SAS','Q_SLS','Q_SSS','Q_OLS','Q_out','Q_GS','Q_ECS', ...
       'Q_BBB','Q_A','P_A'};
SW = {'P_ECS','P_SAS','P_VEN','V_VEN','V_SPI','V_Vessel','V_A','V_C','V_V', ...
      'Q_VEN','Q_SPI','Q_A','P_A'};
% B references read from Results/H_series/B.mat (base_ref.m) in the
% orders declared by KEY and SW, so the sensitivity deltas below cannot go stale
% when the model changes.  A signal the current model does not log returns NaN
% and its delta cell simply comes out empty.
[base.K, base.S] = base_ref('KS', KEY, SW);

sumCode = {}; sumMat = []; sumSW = []; sumStop = []; sumTau = []; sumRatio = [];

for kk = 1:numel(codes)
    code = codes{kk};
    fprintf('\n==================================================\n');
    fprintf('###  %s   ( %d / %d )\n', code, kk, numel(codes));
    fprintf('==================================================\n');
    tSc = tic;
    try
        SCENARIO_CODE = code; %#ok<NASGU>
        Brain_Parameters;
        ap = evalin('base','sApply'); sv = evalin('base','sValue');
        pv = evalin('base','PVI');  km = evalin('base','k_m');
        be = evalin('base','b_EA'); rf = evalin('base','R_FM');
        ps = evalin('base','P_SPI0'); rb = evalin('base','R_BBB');
        re = evalin('base','R_ECS'); reI = evalin('base','r_IEG');
        rvi = evalin('base','R_IEGH');
        Rag = evalin('base','R_AG'); Rcr = evalin('base','R_CP');
        Rsp = evalin('base','R_SPI'); RVen = evalin('base','R_VEN');
        RVenb = evalin('base','R_VEN_base');
        Qcsf = evalin('base','Q_CSF'); Psss = evalin('base','P_SSS');
        Pspi = evalin('base','P_SPI0'); Pcr = evalin('base','P_OLS');
        cc = evalin('base','central_cycle'); tf = evalin('base','t_fit');
        G = 1/Rag + 1/Rcr + 1/Rsp;
        Rout = 1/G;
        Peq  = (Qcsf + Psss/Rag + Pspi/Rsp + Pcr/Rcr)/G;
        Ccap = km/Peq;
        Reff = Rout + (RVen - RVenb);
        tauMin = 0.78*Reff*Ccap;
        TSTOP = max(3000, ceil(6*tauMin*60/500)*500);
        HRimpl = 60/((floor(tf(end)/0.01)+1)*0.01);
        fprintf('  override        %s = %s\n', ap, sv);
        fprintf('  PVI=%g k_m=%.4f b_EA=%.5f R_FM=%.4f P_SPI0=%.2f R_BBB=%.2f R_ECS=%.3f\n', ...
            pv, km, be, rf, ps, rb, re);
        fprintf('  r_IEG=%g R_IEGH=%.2f   R_AG=%.2f R_CP=%.2f R_SPI=%.2f\n', ...
            reI, rvi, Rag, Rcr, Rsp);
        fprintf('  waveform        mean %.4f  pulse %.4f mmHg  n=%d  implied HR %.2f bpm\n', ...
            mean(cc), max(cc)-min(cc), numel(cc), HRimpl);
        fprintf('  Rout=%.4f P_eq=%.4f C=%.4f tau=%.2f min -> StopTime %d s\n', ...
            Rout, Peq, Ccap, tauMin, TSTOP);

        tS = tic;
        out = sim(mdl,'StopTime',num2str(TSTOP),'SaveOutput','off','SignalLogging','on', ...
                      'SignalLoggingName','logsout','ReturnWorkspaceOutputs','on');
        ls = out.logsout;
        tvec = getT(ls,'P_ECS');
        fprintf('  sim %.1f s  (samples %d, cycles %d)\n', ...
            toc(tS), numel(tvec), floor(tvec(end)/min(1.0,60/HRimpl)));

        f = fullfile(outDir,[strrep(code,'+','_') '.mat']);
        save(f,'ls','-v7.3');
        dd = dir(f);
        fprintf('  saved %s (%.1f MB)\n', f, dd.bytes/1e6);

        nms = ls.getElementNames();
        bad = {};
        for q = 1:numel(nms)
            y = gv(ls, nms{q});
            if any(~isfinite(y)), bad{end+1} = nms{q}; end %#ok<SAGROW>
        end
        if isempty(bad)
            fprintf('  all %d logged signals finite   OK\n', numel(nms));
        else
            fprintf('  *** NON-FINITE in: %s\n', strjoin(bad, ', '));
        end

        % ---- steady state, window in cycles ----------------------------
        cyc = min(1.0, 60/HRimpl);         % cycle length to fold the means over
        W = 400;
        vals = nan(1,numel(KEY));
        for q = 1:numel(KEY)
            if any(strcmp(nms, KEY{q}))
                vals(q) = ssm(ls, KEY{q}, cyc, W);
            end
        end
        fprintf('\n  --- steady state (last %d cycles) vs B ---\n', W);
        fprintf('    %-10s %14s %14s %9s\n','signal','value','B','delta %');
        for q = 1:numel(KEY)
            if isnan(vals(q)), continue; end
            b = base.K(q);
            d = NaN; if abs(b) > 1e-12, d = (vals(q)-b)/b*100; end
            fprintf('    %-10s %14.6f %14.6f %9.2f\n', KEY{q}, vals(q), b, d);
        end

        fprintf('\n    ---- outflow allocation (%% of Q_out) ----\n');
        iT = strcmp(KEY,'Q_out');
        tt = vals(iT);
        for s = {'Q_SSS','Q_OLS','Q_SLS'}
            q = strcmp(KEY, s{1});
            if ~isnan(vals(q)), fprintf('      %-8s %7.2f%%\n', s{1}, vals(q)/tt*100); end
        end
        fprintf('      (B  AG:crib:spine = 31.14 : 46.29 : 22.56)\n');

        % ---- swings ----------------------------------------------------
        swv = nan(1,numel(SW));
        for q = 1:numel(SW)
            if any(strcmp(nms, SW{q})), swv(q) = swing(ls, SW{q}, cyc, 100); end
        end
        fprintf('\n  --- per-cycle swing (last 100 cycles) vs B ---\n');
        fprintf('    %-10s %14s %14s %9s\n','signal','value','B','delta %');
        for q = 1:numel(SW)
            if isnan(swv(q)), continue; end
            b = base.S(q);
            d = NaN; if abs(b) > 1e-12, d = (swv(q)-b)/b*100; end
            fprintf('    %-10s %14.6f %14.6f %9.2f\n', SW{q}, swv(q), b, d);
        end

        % ---- convergence ----------------------------------------------
        fprintf('\n  --- convergence diagnostic (12 blocks) ---\n');
        rat = NaN(1,3);
        sl = {'P_ECS','V_VEN','Q_GS'};
        for q = 1:numel(sl)
            if ~any(strcmp(nms, sl{q})), continue; end
            v = gv(ls, sl{q}); t = getT(ls, sl{q});
            n = floor(t(end)/cyc); cm = nan(n,1);
            for c = 1:n
                m = (t > (c-1)*cyc) & (t <= c*cyc);
                if any(m), cm(c) = mean(v(m)); end
            end
            perB = floor(n/12); BM = nan(12,1);
            for b = 1:12
                r = (b-1)*perB + (1:perB);
                BM(b) = mean(cm(r));
            end
            d = diff(BM); rat(q) = d(end)/d(end-1);
            fprintf('    %-8s last diff %+.6f  ratio %+.3f\n', sl{q}, d(end), rat(q));
        end

        sumCode{end+1} = code; %#ok<SAGROW>
        sumMat(end+1,:) = vals; %#ok<SAGROW>
        sumSW(end+1,:) = swv; %#ok<SAGROW>
        sumStop(end+1) = TSTOP; sumTau(end+1) = tauMin; %#ok<SAGROW>
        sumRatio(end+1) = max(abs(rat)); %#ok<SAGROW>
        clear out ls
        fprintf('\n  wall clock %.1f s\n', toc(tSc));
    catch ME
        fprintf('\n  *** %s FAILED ***\n   id = %s\n', code, ME.identifier);
        m = ME.message; asc = repmat('.',1,numel(m));
        for q = 1:numel(m)
            c = double(m(q)); if c >= 32 && c < 127, asc(q) = char(c); end
        end
        fprintf('   msg = %s\n', asc);
        for s = 1:numel(ME.stack)
            fprintf('     at %s line %d\n', ME.stack(s).name, ME.stack(s).line);
        end
    end
end

% ---- summary ----------------------------------------------------------
fprintf('\n==================================================\nSUMMARY: S series (%d done)\n==================================================\n', ...
    numel(sumCode));
csv = fullfile(outDir,'S_series_summary.csv');
% Overwrite, never append.  This script used to open with 'a' and write the
% header only when the file did not exist, so every re-run left the previous
% run's rows in place on top of the new ones.  After the 2026-09-15 re-run the
% file held 32 rows for 21 scenarios, the stale ones recognisable by their
% non-empty V_All column -- V_All stopped being logged in that revision.  Every
% other runner already opens its summary with 'w'.
fid = fopen(csv,'w');
fprintf(fid,'code,StopTime,tau_min,conv_ratio'); fprintf(fid,',%s',KEY{:});
fprintf(fid,',SW_%s',SW{:}); fprintf(fid,'\n');
for kk = 1:numel(sumCode)
    fprintf(fid,'%s,%d,%.3f,%.3f', sumCode{kk}, sumStop(kk), sumTau(kk), sumRatio(kk));
    fprintf(fid,',%.8g', sumMat(kk,:));
    fprintf(fid,',%.8g', sumSW(kk,:));
    fprintf(fid,'\n');
end
fclose(fid);
fprintf('wrote %s\n', csv);

csv2 = fullfile(outDir,'S_series_delta_pct.csv');
fid = fopen(csv2,'w');   % overwrite for the same reason as the summary above
fprintf(fid,'code');
fprintf(fid,',%s',KEY{:}); fprintf(fid,',SW_%s',SW{:}); fprintf(fid,'\n');
for kk = 1:numel(sumCode)
    fprintf(fid,'%s', sumCode{kk});
    dv = (sumMat(kk,:) - base.K)./base.K*100;
    fprintf(fid,',%.4f', dv);
    ds = (sumSW(kk,:) - base.S)./base.S*100;
    fprintf(fid,',%.4f', ds);
    fprintf(fid,'\n');
end
fclose(fid);
fprintf('wrote %s\n', csv2);

fprintf('\n=== total wall clock %.1f s (%.2f hours) ===\n=== DONE ===\n', ...
    toc(tAll), toc(tAll)/3600);

% =====================================================================
function v = gv(ls, name)
    e = ls.getElement(name);
    if isa(e,'Simulink.SimulationData.Dataset'), ts = e.getElement(1).Values; else, ts = e.Values; end
    v = double(ts.Data(:));
end
function t = getT(ls, name)
    e = ls.getElement(name);
    if isa(e,'Simulink.SimulationData.Dataset'), ts = e.getElement(1).Values; else, ts = e.Values; end
    t = double(ts.Time(:));
end
function v = ssm(ls, name, CYCLE, W)
    % Average over EXACTLY W whole cycles.  The upper bound matters: Q_SPI has a
    % pulsatile amplitude ~4000x its mean, so including even a fraction of an
    % extra cycle biases its mean badly.  Every previously recorded scenario had
    % a StopTime that was an integer multiple of the cycle (3000/0.8 = 3750 etc.)
    % and so was never affected; S_HR_90 (cycle 0.6700 s, 3000 s) is.
    y = gv(ls,name); t = getT(ls,name);
    n = floor(t(end)/CYCLE);
    m = (t > (n-W)*CYCLE) & (t <= n*CYCLE);
    v = mean(y(m));
end
function d = swing(ls, name, CYCLE, N)
    y = gv(ls,name); t = getT(ls,name);
    n = floor(t(end)/CYCLE); R = nan(N,1);
    for c = (n-N+1):n
        m = (t > (c-1)*CYCLE) & (t <= c*CYCLE);
        R(c-(n-N)) = max(y(m)) - min(y(m));
    end
    d = mean(R);
end
