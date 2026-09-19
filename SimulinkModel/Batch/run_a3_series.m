%% run_a3_series.m -- A3 combined pathology series
%
%   H+C   : Grade-2 hypertension + all three outflow pathways at 500 %
%   HA+C  : acute Grade-2 hypertension (baseline vasculature, PVI 25) + the
%           same three-pathway impairment
%   (off-table extras kept for the record: H_G2+Ven_R7.5,
%    H_G2+C_OVS2.5+Ven_E0.25)
%
% A combined code is split on '+': the H_* part supplies the vasculature
% (pulseType / PVI / r_R_Hyper / r_E / r_BBB) and the C_* part supplies the CSF
% side (outflow-resistance multiples).  The Table 2 shorthands H+C and HA+C are
% expanded inside bp_compute.  See check_codes.m for the parameter-level check.
%
% StopTime: the compliance is k_m/P_eq with k_m = 0.4343*PVI, so the estimate
% below is built from the CRANIAL compliance only and UNDER-estimates the real
% settling time by roughly a factor of two (the spinal compartment carries 69 %
% of the craniospinal compliance).  Measured: H+C needs 9000 s where the formula
% says 4000 s.  Use A3_STOP_OVERRIDE when in doubt:
%   R_eff    = Rout_eff + (R_VEN - R_VEN_base)
%   P_eq     = (Q_CSF + P_SSS/R_AG + P_SPI0/R_SPI + P_OLS/R_CP) / (1/Rout_eff)
%   C        = k_m / P_eq
%   tau      = 0.78 * R_eff * C   [min]
%   StopTime = max(3000, ceil(6*tau*60/500)*500)
%
% NON-DESTRUCTIVE: no save_system, no logging reconfiguration.
% ASCII only.  matlab -batch "run_a3_series"
%
% NOTE: variable names i / S / fn are avoided -- Brain_Parameters clears them.

fprintf('=== run_a3_series (MATLAB %s) ===\n', version);
tAll = tic;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
mdl  = 'CSF_Model';
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(here);
addpath(fullfile(root,'SimulinkModel'));
outDir = fullfile(root,'SimulinkModel','Results','A3_series');
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~bdIsLoaded(mdl), load_system(fullfile(root,'SimulinkModel',[mdl '.slx'])); end
CYCLE = 0.8;

codes = {'H+C','HA+C'};
if exist('A3_CODES','var') && ~isempty(A3_CODES), codes = A3_CODES; end

% B reference values, READ FROM THE BASE RUN via base_ref.m (see the note in
% run_venr_series.m): reading them from Results/H_series/B.mat means they
% cannot go stale when the model changes.
ref = base_ref('ref');

KEYP = {'P_ECS','P_VEN','P_SAS','P_SPI','V_VEN','V_SPI','V_ECS','V_All'};
KEYF = {'Q_VEN','Q_SPI','Q_SAS','Q_SLS','Q_SSS','Q_OLS','Q_out', ...
        'Q_GS','Q_BBB','Q_A'};

sumCode = {}; sumRow = []; sumStop = []; sumAlloc = [];

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
        kkm  = evalin('base','k_m');   PVIv = evalin('base','PVI');
        re   = evalin('base','r_E');   rbb = evalin('base','r_BBB');
        rrh  = evalin('base','r_R_Hyper');
        RVen = evalin('base','R_VEN'); RVenb = evalin('base','R_VEN_base');
        G = 1/Rag + 1/Rcr + 1/Rsp;
        Rout = 1/G;
        Peq  = (Qcsf + Psss/Rag + Pspi/Rsp) / G;
        Ccap = kkm / Peq;
        Reff = Rout + (RVen - RVenb);
        tauMin = 0.78 * Reff * Ccap;
        TSTOP = max(3000, ceil(6*tauMin*60/500)*500);
        % Same under-estimate as in run_c_series: tau is built from the CRANIAL
        % compliance, but when all three outflow resistances are raised the
        % slowest mode is the spinal/whole-system volume equilibration.  In the
        % first pass H+C was given 4000 s and came out with a last-block
        % ratio of 0.786/0.787 and a P_ECS drift of 0.12 mmHg per 250 s, i.e. not
        % converged (its partner H_G2+Ven_R7.5 converged at 0.42-0.46 with 6500 s).
        % A3_STOP_OVERRIDE forces a longer StopTime for the re-run:
        %   matlab -batch "A3_CODES={'H+C'}; A3_STOP_OVERRIDE=12000; run_a3_series"
        if exist('A3_STOP_OVERRIDE','var') && ~isempty(A3_STOP_OVERRIDE)
            TSTOP = A3_STOP_OVERRIDE;
            fprintf('  [override] StopTime forced to %d s\n', TSTOP);
        end
        fprintf('  vasculature: PVI=%g  r_R=%.2f  r_E=%.2f  r_BBB=%.2f  k_m=%.4f\n', ...
            PVIv, rrh, re, rbb, kkm);
        fprintf('  CSF side   : R_AG=%.2f R_CP=%.2f R_SPI=%.2f  Rout_eff=%.4f\n', ...
            Rag, Rcr, Rsp, Rout);
        fprintf('               R_VEN=%.4f (baseline %.4f, +%.4f)\n', ...
            RVen, RVenb, RVen-RVenb);
        fprintf('  R_eff=%.3f  P_eq=%.4f  C=%.4f  tau=%.2f min (%.0f s)  -> StopTime %d s\n', ...
            Reff, Peq, Ccap, tauMin, tauMin*60, TSTOP);

        tS = tic;
        out = sim(mdl,'StopTime',num2str(TSTOP),'SaveOutput','off','SignalLogging','on', ...
                      'SignalLoggingName','logsout','ReturnWorkspaceOutputs','on');
        ls = out.logsout;
        tvec = getT(ls,'P_ECS');
        fprintf('  sim %.1f s   (samples %d, cycles %d)\n', ...
            toc(tS), numel(tvec), floor(tvec(end)/CYCLE));

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

        W = 400;
        fprintf('\n  --- steady state, mean over last %d cycles ---\n', W);
        vals = struct();
        for s = [KEYP, KEYF]
            if ~any(strcmp(nms, s{1})), continue; end
            v = gv(ls, s{1}); t = getT(ls, s{1});
            n = floor(t(end)/CYCLE);
            m = t > (n-W)*CYCLE;
            vals.(s{1}) = mean(v(m));
            fprintf('    %-10s %16.6f\n', s{1}, vals.(s{1}));
        end

        grad = vals.P_VEN - vals.P_SAS;
        fprintf('\n  --- CSF pressure structure ---\n');
        fprintf('    P_VEN            %14.6f   (B %.6f)\n', vals.P_VEN, ref.P_VEN);
        fprintf('    P_SAS            %14.6f   (B %.6f)\n', vals.P_SAS, ref.P_SAS);
        fprintf('    P_VEN - P_SAS    %14.6f   (B %.6f, delta %+.6f)\n', ...
            grad, ref.grad, grad - ref.grad);
        fprintf('    V_VEN            %14.6f   (B %.6f, delta %+.6f mL)\n', ...
            vals.V_VEN, ref.V_VEN, vals.V_VEN - ref.V_VEN);
        fprintf('    V_ECS            %14.6f   (B %.6f, delta %+.6f mL)\n', ...
            vals.V_ECS, ref.V_ECS, vals.V_ECS - ref.V_ECS);
        fprintf('    Monro-Kellie check: dV_VEN + dV_ECS = %+.6f mL\n', ...
            (vals.V_VEN - ref.V_VEN) + (vals.V_ECS - ref.V_ECS));

        if isfield(vals,'Q_out')
            tt = vals.Q_out;
            fprintf('\n    ---- outflow allocation (%% of Q_out) ----\n');
            al = zeros(1,3); q = 0;
            for s = {'Q_SSS','Q_OLS','Q_SLS'}
                q = q + 1;
                if isfield(vals,s{1})
                    al(q) = vals.(s{1})/tt*100;
                    fprintf('      %-8s %7.2f%%\n', s{1}, al(q));
                end
            end
            sumAlloc(end+1,:) = al; %#ok<SAGROW>
            fprintf('      (B  AG:crib:spine = %.2f : %.2f : %.2f)\n', ref.alloc);
        end

        fprintf('\n  --- per-cycle cardiac swing (last 100 cycles) ---\n');
        for s = {'P_VEN','P_SAS','P_ECS','V_VEN','V_SPI','Q_VEN','Q_SPI','Q_A'}
            if ~any(strcmp(nms, s{1})), continue; end
            v = gv(ls, s{1}); t = getT(ls, s{1});
            n = floor(t(end)/CYCLE); R = nan(100,1);
            for c = (n-99):n
                m = (t > (c-1)*CYCLE) & (t <= c*CYCLE);
                R(c-(n-100)) = max(v(m)) - min(v(m));
            end
            fprintf('    %-10s %14.6f +/- %.6f\n', s{1}, mean(R), std(R));
        end

        fprintf('\n  --- convergence diagnostic (12 blocks) ---\n');
        for s = {'P_VEN','P_SAS','P_ECS','V_VEN'}
            if ~any(strcmp(nms, s{1})), continue; end
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
            fprintf('    %-8s blocks [%s]\n', s{1}, num2str(BM(:)', '%.4f '));
            fprintf('    %-8s last diff %+.6f  ratio %+.3f  (first diff %+.6f)\n', ...
                s{1}, d(end), d(end)/d(end-1), d(1));
        end

        sumCode{end+1} = code; %#ok<SAGROW>
        sumStop(end+1) = TSTOP; %#ok<SAGROW>
        sumRow(end+1,:) = [vals.P_VEN, vals.P_SAS, grad, vals.V_VEN, vals.P_ECS, ...
            vals.Q_VEN, vals.Q_out, vals.Q_GS, vals.Q_BBB]; %#ok<SAGROW>
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

fprintf('\n==================================================\nSUMMARY: A3 combined\n==================================================\n');
hdr = {'code','StopTime','P_VEN','P_SAS','grad','V_VEN','P_ECS','Q_VEN','Q_out','Q_GS','Q_BBB'};
fprintf('%-14s %9s %11s %11s %11s %11s %11s %11s %11s %11s %11s\n', hdr{:});
fprintf('%-14s %9d %11.5f %11.5f %11.5f %11.5f %11.5f %11.6f %11.6f %11.6f %11.6f\n', ...
    'B', 0, ref.P_VEN, ref.P_SAS, ref.grad, ref.V_VEN, ref.P_ECS, ref.Q_VEN, ...
    ref.Q_out, ref.Q_GS, ref.Q_BBB);
for kk = 1:numel(sumCode)
    fprintf('%-14s %9d %11.5f %11.5f %11.5f %11.5f %11.5f %11.6f %11.6f %11.6f %11.6f\n', ...
        sumCode{kk}, sumStop(kk), sumRow(kk,:));
end

csv = fullfile(outDir,'A3_series_summary.csv');
fid = fopen(csv,'w');
fprintf(fid,'%s\n', strjoin(hdr,','));
fprintf(fid,'B,0,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f\n', ...
    ref.P_VEN, ref.P_SAS, ref.grad, ref.V_VEN, ref.P_ECS, ref.Q_VEN, ...
    ref.Q_out, ref.Q_GS, ref.Q_BBB);
for kk = 1:numel(sumCode)
    fprintf(fid,'%s,%d,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f\n', ...
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
