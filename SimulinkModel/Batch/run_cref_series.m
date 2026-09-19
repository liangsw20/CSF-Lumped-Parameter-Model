%% run_cref_series.m -- single-outflow reference model (configuration A)
%
% PURPOSE (Reviewer 1: multi-pathway vs traditional single-outflow model)
%
%   Reference model: two of the three CSF outflow pathways are BLOCKED and the
%   surviving pathway carries all drainage.  Its resistance is set to the
%   PARALLEL EQUIVALENT of the three baseline pathways, so the total equivalent
%   outflow resistance is identical to the three-pathway model and the
%   comparison isolates the effect of pathway STRUCTURE from the effect of the
%   overall resistance LEVEL:
%
%       1/Rout_base = 1/(25*r_R) + 1/(65*r_R) + 1/(20*r_R)
%       Rout_base   = 9.489051 mmHg*min/mL
%       C_Ref : R_AG   = Rout_base , R_CP = R_BLOCK , R_SPI = R_BLOCK
%       C_Leq : R_CP = Rout_base , R_AG   = R_BLOCK , R_SPI = R_BLOCK
%       C_Seq : R_SPI  = Rout_base , R_AG   = R_BLOCK , R_CP = R_BLOCK
%
%   R_BLOCK = 1e9 (finite, NOT Inf -- Inf produces Inf/Inf -> NaN in the
%   nonlinear resistance subsystems and aborts with AlgStateNotFinite).
%   R_FM keeps its physiological value: craniospinal communication stays open
%   and the spinal SAS acts as a compliant extension of the SAS, which is what
%   a lumped single-chamber model assumes.
%
% StopTime is derived per scenario from the steady-state balance:
%   Rout_eff = 1/(1/R_AG + 1/R_CP + 1/R_SPI)
%   P_eq     = (Q_CSF + P_SSS/R_AG + P_SPI0/R_SPI) / (1/Rout_eff)
%   C        = k_m / P_eq
%   tau      = 0.78 * Rout_eff * C   [min]
%   StopTime = max(3000, ceil(6*tau*60/500)*500)
%
% NON-DESTRUCTIVE: no save_system, no logging reconfiguration.
% ASCII only.  matlab -batch "run_cref_series"
%
% NOTE: variable names i / S / fn are avoided -- Brain_Parameters clears them.

fprintf('=== run_cref_series (MATLAB %s) ===\n', version);
tAll = tic;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
mdl  = 'CSF_Model';
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(here);
addpath(fullfile(root,'SimulinkModel'));
outDir = fullfile(root,'SimulinkModel','Results','Cref_series');
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~bdIsLoaded(mdl), load_system(fullfile(root,'SimulinkModel',[mdl '.slx'])); end
CYCLE = 0.8;

% The reference model is AG-only BY DEFINITION: the classical CSF absorption
% model treats the arachnoid granulations as the sole absorption site, while the
% cribriform-lymphatic and spinal nerve-sheath pathways are exactly what the
% three-pathway model adds.  C_Leq / C_Seq correspond to no traditional model
% and are therefore NOT part of the study plan; the implementation stays in
% bp_compute only for future mechanistic sensitivity work.
codes = {'C_Ref'};

KEYP = {'P_ECS','P_SAS','P_SPI','P_VEN','V_SPI','V_All','V_VEN','V_ECS'};
KEYF = {'Q_SPI','Q_SAS','Q_SLS','Q_SSS','Q_OLS','Q_out','Q_GS','Q_BBB', ...
        'Q_A','Q_VEN'};

% B reference values, READ FROM THE BASE RUN via base_ref.m (see the note in
% run_venr_series.m): they cannot go stale when the model changes.
ref = base_ref('ref');
% % of Q_out, ordered AG(V) : cribriform(L) : spinal(S)
% NOTE the labelling: cribriform/lymphatic is the LARGEST baseline route
% because it drains against atmospheric pressure while AG drains against P_SSS,
% not because its resistance is lowest.

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
        kkm  = evalin('base','k_m');  Rbase = evalin('base','Rout_base');
        Rblk = evalin('base','R_BLOCK'); Rfm = evalin('base','R_FM');
        G = 1/Rag + 1/Rcr + 1/Rsp;
        Rout = 1/G;
        Peq  = (Qcsf + Psss/Rag + Pspi/Rsp) / G;
        Ccap = kkm / Peq;
        tauMin = 0.78 * Rout * Ccap;
        TSTOP = max(3000, ceil(6*tauMin*60/500)*500);
        fprintf('  CONFIG A: R_AG=%.6f  R_CP=%.6g  R_SPI=%.6g  (R_BLOCK=%.0g)\n', ...
            Rag, Rcr, Rsp, Rblk);
        fprintf('  R_FM=%.4g (physiological, craniospinal coupling kept open)\n', Rfm);
        fprintf('  parallel check: 1/(1/Rag+1/Rcr+1/Rsp) = %.9f   vs Rout_base = %.9f', ...
            Rout, Rbase);
        if abs(Rout-Rbase) < 1e-6, fprintf('   MATCH\n'); else, fprintf('   *** MISMATCH ***\n'); end
        fprintf('  P_eq(est)=%.4f mmHg   C=%.4f mL/mmHg   tau=%.2f min (%.0f s)\n', ...
            Peq, Ccap, tauMin, tauMin*60);
        fprintf('  StopTime %d s  (~%.1f min wall)\n', TSTOP, TSTOP/10000*10.1);

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
        if isfield(vals,'Q_out')
            tt = vals.Q_out;
            fprintf('    ---- outflow allocation (%% of Q_out) ----\n');
            al = zeros(1,3); q = 0;
            for s = {'Q_OLS','Q_SLS','Q_SSS'}
                q = q + 1;
                if isfield(vals,s{1})
                    al(q) = vals.(s{1})/tt*100;
                    fprintf('      %-8s %7.2f%%\n', s{1}, al(q));
                end
            end
            sumAlloc(end+1,:) = al; %#ok<SAGROW>
            fprintf('      (B allocation  AG:crib:spine = %.2f : %.2f : %.2f)\n', ...
                ref.alloc(1), ref.alloc(2), ref.alloc(3));
        end

        fprintf('\n  --- per-cycle cardiac swing (last 100 cycles) ---\n');
        for s = {'V_SPI','P_ECS','V_Vessel','P_SAS','Q_SPI','Q_A','V_VEN'}
            if ~any(strcmp(nms, s{1})), continue; end
            v = gv(ls, s{1}); t = getT(ls, s{1});
            n = floor(t(end)/CYCLE); R = nan(100,1);
            for c = (n-99):n
                m = (t > (c-1)*CYCLE) & (t <= c*CYCLE);
                R(c-(n-100)) = max(v(m)) - min(v(m));
            end
            fprintf('    %-10s %14.6f +/- %.6f\n', s{1}, mean(R), std(R));
        end
        if isfield(vals,'V_SPI')
            fprintf('    (B per-cycle V_SPI %.6f, P_ECS %.6f, Q_SPI %.6f)\n', ...
                ref.sw_V_SPI, ref.sw_P_ECS, ref.sw_Q_SPI);
        end

        fprintf('\n  --- convergence diagnostic (12 blocks) ---\n');
        for s = {'V_SPI','P_ECS','V_VEN','Q_SPI'}
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

fprintf('\n==================================================\nSUMMARY: config-A single-outflow vs B\n==================================================\n');
hdr = {'code','StopTime','P_ECS','V_SPI','Q_SPI','Q_SLS','Q_out','Q_GS','Q_BBB','Q_A'};
fprintf('%-8s %9s %11s %11s %11s %11s %11s %11s %11s %11s\n', hdr{:});
fprintf('%-8s %9s %11.5f %11.5f %11.6f %11.6f %11.6f %11.6f %11.6f %11.3f\n', ...
    'B', 0, ref.P_ECS, ref.V_SPI, ref.Q_SPI, ref.Q_SLS, ref.Q_out, ...
    ref.Q_GS, ref.Q_BBB, ref.Q_A);
for kk = 1:numel(sumCode)
    fprintf('%-8s %9d %11.5f %11.5f %11.6f %11.6f %11.6f %11.6f %11.6f %11.3f\n', ...
        sumCode{kk}, sumStop(kk), sumRow(kk,:));
end
fprintf('\n%-8s %11s %11s %11s   (B: %.5f)\n', 'ICP err', 'dP_ECS', 'dP_ECS %', '', ref.P_ECS);
for kk = 1:numel(sumCode)
    dp = sumRow(kk,1) - ref.P_ECS;
    fprintf('%-8s %11.5f %11.3f\n', sumCode{kk}, dp, dp/ref.P_ECS*100);
end

csv = fullfile(outDir,'SingleEq_summary.csv');
fid = fopen(csv,'w');
fprintf(fid,'%s\n', strjoin(hdr,','));
fprintf(fid,'B,0,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.4f\n', ...
    ref.P_ECS, ref.V_SPI, ref.Q_SPI, ref.Q_SLS, ref.Q_out, ...
    ref.Q_GS, ref.Q_BBB, ref.Q_A);
for kk = 1:numel(sumCode)
    fprintf(fid,'%s,%d,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.4f\n', ...
        sumCode{kk}, sumStop(kk), sumRow(kk,:));
end
fclose(fid);
fprintf('\nwrote %s\n', csv);

if ~isempty(sumAlloc)
    csv2 = fullfile(outDir,'SingleEq_allocation.csv');
    fid = fopen(csv2,'w');
    fprintf(fid,'code,AG_V,crib_L,spinal_S\n');
    fprintf(fid,'B,%12.6f,%12.6f,%12.6f\n', ref.alloc);
    for kk = 1:size(sumAlloc,1)
        fprintf(fid,'%s,%12.6f,%12.6f,%12.6f\n', sumCode{kk}, sumAlloc(kk,:));
    end
    fclose(fid);
    fprintf('wrote %s\n', csv2);
end

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
