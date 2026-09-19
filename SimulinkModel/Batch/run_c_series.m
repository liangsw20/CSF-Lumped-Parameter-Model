%% run_c_series.m -- CSF outflow resistance impairment (communicating) series
%
%   7 pathway combinations x 2 alpha levels = 14 scenarios
%     L = lymphatic/cribriform (R_CP), V = arachnoid granulations (R_AG),
%     S = spinal nerve sheaths (R_SPI)
%     alpha = 2.5 (mean of the NPH literature) and 5 (model's own extreme)
%
% The StopTime is derived per scenario from the steady-state flow balance:
%   Rout_eff = 1/(1/R_AG + 1/R_CP + 1/R_SPI)
%   P_eq     = (Q_CSF + P_SSS/R_AG + P_SPI0/R_SPI) / (1/Rout_eff)
%   C        = k_m / P_eq                        (exponential P-V relation)
%   tau      = 0.78 * Rout_eff * C   [minutes]   (0.78 fitted on B/H_G2A)
%   StopTime = max(3000, ceil(6*tau*60 / 500)*500)
%
% NON-DESTRUCTIVE: no save_system, no logging reconfiguration.
% ASCII only.  matlab -batch "run_c_series"
%
% NOTE: variable names i / S / fn are avoided -- Brain_Parameters clears them.

fprintf('=== run_c_series (MATLAB %s) ===\n', version);
tAll = tic;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
mdl  = 'CSF_Model';
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(here);
addpath(fullfile(root,'SimulinkModel'));
outDir = fullfile(root,'SimulinkModel','Results','C_series');
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~bdIsLoaded(mdl), load_system(fullfile(root,'SimulinkModel',[mdl '.slx'])); end
CYCLE = 0.8;

codes = {'C_O2.5','C_V2.5','C_S2.5','C_OV2.5','C_OS2.5','C_VS2.5','C_OVS2.5', ...
         'C_O5'  ,'C_V5'  ,'C_S5'  ,'C_OV5'  ,'C_OS5'  ,'C_VS5'  ,'C_OVS5'};
% Optional selection, matched to run_venr_series / run_a3_series / run_s_series:
%   matlab -batch "C_CODES={'C_OVS5'}; C_STOP_OVERRIDE=15000; run_c_series"
%
% This hook did NOT exist before the 2026-09-15 fix-up pass, although the
% override comment further down referred to it.  Without it C_CODES was silently
% ignored and the full 14-scenario sweep ran: the first attempt at a single
% C_OVS5 re-run therefore started C_O2.5 and, combined with a forced StopTime of
% 15000 s, would have taken over five hours instead of 23 minutes.  If a new
% runner is added, give it the same <CODE>_CODES hook.
if exist('C_CODES','var') && ~isempty(C_CODES), codes = C_CODES; end

KEYP = {'P_ECS','P_SAS','P_SPI','P_VEN','V_SPI','V_All','V_VEN','V_ECS'};
KEYF = {'Q_SPI','Q_SAS','Q_SLS','Q_SSS','Q_OLS','Q_out','Q_GS','Q_BBB','Q_A','Q_VEN'};

sumCode = {}; sumRow = []; sumStop = [];

for kk = 1:numel(codes)
    code = codes{kk};
    fprintf('\n==================================================\n');
    fprintf('###  %s   ( %d / %d )\n', code, kk, numel(codes));
    fprintf('==================================================\n');
    tSc = tic;
    try
        SCENARIO_CODE = code; %#ok<NASGU>
        Brain_Parameters;
        Rag = evalin('base','R_AG'); Rcr = evalin('base','R_CP');
        Rsp = evalin('base','R_SPI'); Qcsf = evalin('base','Q_CSF');
        Psss = evalin('base','P_SSS'); Pspi = evalin('base','P_SPI0');
        kkm  = evalin('base','k_m');
        G = 1/Rag + 1/Rcr + 1/Rsp;
        Rout = 1/G;
        Peq  = (Qcsf + Psss/Rag + Pspi/Rsp) / G;
        Ccap = kkm / Peq;
        tauMin = 0.78 * Rout * Ccap;
        TSTOP = max(3000, ceil(6*tauMin*60/500)*500);
        % The analytic tau above uses the CRANIAL compliance k_m/P_eq, but in the
        % stiffest scenarios the slowest mode is the spinal/whole-system volume
        % equilibration, so 6*tau under-estimates the StopTime needed.  C_OVS5
        % (all three resistances x5) came out at 7000 s with a last-block P_ECS
        % drift of 0.139 mmHg per 250 s and a block ratio of 0.80, i.e. not
        % converged: its inflow-outflow residual was still 6.9% of Q_out,
        % while every other C scenario was below 2%.  C_STOP_OVERRIDE forces a
        % longer StopTime for such a re-run:
        %   matlab -batch "C_CODES={'C_OVS5'}; C_STOP_OVERRIDE=15000; run_c_series"
        if exist('C_STOP_OVERRIDE','var') && ~isempty(C_STOP_OVERRIDE)
            TSTOP = C_STOP_OVERRIDE;
            fprintf('  [override] StopTime forced to %d s\n', TSTOP);
        end
        fprintf('  R_CP=%.1f  R_AG=%.1f  R_SPI=%.1f  Rout_eff=%.2f\n', Rcr, Rag, Rsp, Rout);
        fprintf('  P_eq(est)=%.2f mmHg  tau=%.1f min (%.0f s)  -> StopTime %d s\n', ...
            Peq, tauMin, tauMin*60, TSTOP);

        tS = tic;
        out = sim(mdl,'StopTime',num2str(TSTOP),'SaveOutput','off','SignalLogging','on', ...
                      'SignalLoggingName','logsout','ReturnWorkspaceOutputs','on');
        ls = out.logsout;
        tvec = getT(ls,'P_ECS');
        fprintf('  sim %.1f s   (samples %d, cycles %d)\n', ...
            toc(tS), numel(tvec), floor(tvec(end)/CYCLE));

        f = fullfile(outDir,[code '.mat']);
        save(f,'ls','-v7.3');
        dd = dir(f);
        fprintf('  saved %s (%.1f MB)\n', f, dd.bytes/1e6);

        W = 400;
        fprintf('\n  --- steady state, mean over last %d cycles ---\n', W);
        vals = struct();
        for s = [KEYP, KEYF]
            if ~any(strcmp(ls.getElementNames(), s{1})), continue; end
            v = gv(ls, s{1}); t = getT(ls, s{1});
            n = floor(t(end)/CYCLE);
            m = t > (n-W)*CYCLE;
            vals.(s{1}) = mean(v(m));
            fprintf('    %-10s %16.6f\n', s{1}, vals.(s{1}));
        end
        if isfield(vals,'Q_out')
            tt = vals.Q_out;
            fprintf('    ---- outflow allocation (%% of Q_out) ----\n');
            for s = {'Q_OLS','Q_SLS','Q_SSS'}
                if isfield(vals,s{1}), fprintf('      %-8s %6.2f%%\n', s{1}, vals.(s{1})/tt*100); end
            end
        end

        fprintf('\n  --- per-cycle cardiac swing (last 100 cycles) ---\n');
        for s = {'V_SPI','P_ECS','V_Vessel','P_SAS','Q_SPI','Q_A'}
            if ~any(strcmp(ls.getElementNames(), s{1})), continue; end
            v = gv(ls, s{1}); t = getT(ls, s{1});
            n = floor(t(end)/CYCLE); R = nan(100,1);
            for c = (n-99):n
                m = (t > (c-1)*CYCLE) & (t <= c*CYCLE);
                R(c-(n-100)) = max(v(m)) - min(v(m));
            end
            fprintf('    %-10s %14.6f +/- %.6f\n', s{1}, mean(R), std(R));
        end

        fprintf('\n  --- convergence diagnostic (12 blocks) ---\n');
        for s = {'V_VEN','V_SPI','P_ECS'}
            if ~any(strcmp(ls.getElementNames(), s{1})), continue; end
            v = gv(ls, s{1}); t = getT(ls, s{1});
            n = floor(t(end)/CYCLE); cm = nan(n,1);
            for c = 1:n
                m = (t > (c-1)*CYCLE) & (t <= c*CYCLE);
                if any(m), cm(c) = mean(v(m)); end
            end
            perB = floor(n/12); BM = nan(12,1);
            for b = 1:12
                r = (b-1)*perB + (1:perB);
                BM(b) = mean(cm(r));
            end
            d = diff(BM);
            fprintf('    %-8s last diff %+.6f  ratio %+.3f\n', s{1}, d(end), d(end)/d(end-1));
        end

        sumCode{end+1} = code; %#ok<SAGROW>
        sumStop(end+1) = TSTOP; %#ok<SAGROW>
        sumRow(end+1,:) = [vals.P_ECS, vals.V_SPI, vals.Q_SPI, vals.Q_SLS, ...
            vals.Q_out, vals.Q_GS, vals.Q_BBB, vals.Q_A]; %#ok<SAGROW>
        clear out ls
        fprintf('\n  scenario wall clock %.1f s\n', toc(tSc));
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

fprintf('\n==================================================\nSUMMARY\n==================================================\n');
hdr = {'code','StopTime','P_ECS','V_SPI','Q_SPI','Q_SLS','Q_out','Q_GS','Q_BBB','Q_A'};
fprintf('%-10s %9s %11s %11s %11s %11s %11s %11s %11s %11s\n', hdr{:});
for kk = 1:numel(sumCode)
    fprintf('%-10s %9d %11.5f %11.5f %11.6f %11.6f %11.6f %11.6f %11.6f %11.3f\n', ...
        sumCode{kk}, sumStop(kk), sumRow(kk,:));
end
csv = fullfile(outDir,'C_series_summary.csv');
fid = fopen(csv,'w');
fprintf(fid,'%s\n', strjoin(hdr,','));
for kk = 1:numel(sumCode)
    fprintf(fid,'%s,%d,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.4f\n', ...
        sumCode{kk}, sumStop(kk), sumRow(kk,:));
end
fclose(fid);
fprintf('\nwrote %s\n', csv);
fprintf('\n=== total wall clock %.1f s (%.1f hours) ===\n=== DONE ===\n', ...
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
