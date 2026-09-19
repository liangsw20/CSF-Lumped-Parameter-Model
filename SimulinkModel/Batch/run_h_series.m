%% run_h_series.m -- run the H (hypertension) series and store the raw results
%
%   H_HM  (132/84, high-normal)   H_G1 (144/88, grade 1)
%   H_G2  (164/92, grade 2)       H_G2A (164/92, acute -- pressure only)
%
% Stores out.logsout for each scenario in Results/H_series/<code>.mat and
% prints steady-state values plus a convergence diagnostic, so a scenario that
% needs a longer run can be identified without redoing the others.
%
% NON-DESTRUCTIVE: no save_system, no logging reconfiguration.
% ASCII only.  matlab -batch "run_h_series"
%
% NOTE: variable names i / S / fn are avoided on purpose -- Brain_Parameters
% ends with `clear S fn i`.

fprintf('=== run_h_series (MATLAB %s) ===\n', version);
tAll = tic;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
mdl  = 'CSF_Model';
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(here);
addpath(fullfile(root,'SimulinkModel'));
outDir = fullfile(root,'SimulinkModel','Results','H_series');
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~bdIsLoaded(mdl), load_system(fullfile(root,'SimulinkModel',[mdl '.slx'])); end

codes = {'H_HM','H_G1','H_G2','H_G2A'};
TSTOP = 3000;
CYCLE = 0.8;

% summary accumulators
sumCode = {}; sumRow = [];
KEYP = {'P_ECS','P_SAS','P_SPI','P_VEN','V_SPI','V_All','V_VEN','V_ECS'};
KEYF = {'Q_SPI','Q_SAS','Q_SLS','Q_SSS','Q_OLS','Q_out','Q_GS','Q_BBB','Q_A','Q_VEN'};

for kk = 1:numel(codes)
    code = codes{kk};
    fprintf('\n==================================================\n');
    fprintf('###  %s   ( %d / %d )   StopTime %g s\n', code, kk, numel(codes), TSTOP);
    fprintf('==================================================\n');
    tSc = tic;
    try
        SCENARIO_CODE = code; %#ok<NASGU>
        Brain_Parameters;
        ccv = evalin('base','central_cycle');
        tfv = evalin('base','t_fit');
        fprintf('  grid: central_cycle %s  t_fit %s  span [%.4f, %.4f]  step %.4f\n', ...
            mat2str(size(ccv)), mat2str(size(tfv)), tfv(1), tfv(end), mean(diff(tfv)));
        fprintf('  PVI=%.0f r_R=%.2f r_E=%.2f r_BBB=%.2f  R_AG=%.2f R_CP=%.2f R_SPI=%.2f\n', ...
            evalin('base','PVI'), evalin('base','r_R_Hyper'), evalin('base','r_E'), ...
            evalin('base','r_BBB'), evalin('base','R_AG'), evalin('base','R_CP'), ...
            evalin('base','R_SPI'));

        tS = tic;
        out = sim(mdl,'StopTime',num2str(TSTOP),'SaveOutput','off','SignalLogging','on', ...
                      'SignalLoggingName','logsout','ReturnWorkspaceOutputs','on');
        ls = out.logsout;
        fprintf('  sim %.1f s\n', toc(tS));

        f = fullfile(outDir,[code '.mat']);
        save(f,'ls','-v7.3');
        dd = dir(f);
        fprintf('  saved %s (%.1f MB)\n', f, dd.bytes/1e6);

        %% steady state (last 400 cycles)
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

        %% per-cycle cardiac swings (last 100 cycles)
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

        %% convergence diagnostic: block means of the slowest states
        fprintf('\n  --- convergence diagnostic (12 blocks of 250 s) ---\n');
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
            fprintf('    %-8s block diffs:', s{1});
            for b = 1:numel(d), fprintf(' %+.6f', d(b)); end
            fprintf('\n             last diff %+.6f  ratio(last/prev) %+.3f\n', ...
                d(end), d(end)/d(end-1));
        end

        %% record a summary row
        sumCode{end+1} = code; %#ok<SAGROW>
        sumRow(end+1,:) = [vals.P_ECS, vals.V_SPI, vals.V_VEN, ...
            vals.Q_SPI, vals.Q_SLS, vals.Q_out, vals.Q_A]; %#ok<SAGROW>
        clear out ls
        fprintf('\n  scenario wall clock %.1f s\n', toc(tSc));
    catch ME
        fprintf('\n  *** %s FAILED ***\n   id  = %s\n', code, ME.identifier);
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

%% summary table + CSV
fprintf('\n==================================================\n');
fprintf('SUMMARY\n');
fprintf('==================================================\n');
hdr = {'code','P_ECS','V_SPI','V_VEN','Q_SPI','Q_SLS','Q_out','Q_A'};
fprintf('%-8s %12s %12s %12s %12s %12s %12s %12s\n', hdr{:});
for kk = 1:numel(sumCode)
    fprintf('%-8s %12.5f %12.5f %12.5f %12.6f %12.6f %12.6f %12.3f\n', ...
        sumCode{kk}, sumRow(kk,:));
end
csv = fullfile(outDir,'H_series_summary.csv');
fid = fopen(csv,'w');
fprintf(fid,'%s\n', strjoin(hdr,','));
for kk = 1:numel(sumCode)
    fprintf(fid,'%s,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.6f,%12.4f\n', ...
        sumCode{kk}, sumRow(kk,:));
end
fclose(fid);
fprintf('\nwrote %s\n', csv);
fprintf('\n=== total wall clock %.1f s (%.1f min) ===\n=== DONE ===\n', ...
    toc(tAll), toc(tAll)/60);

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
