%% run_base_final.m -- B scenario: simulate once, then draw the figures
%
% NON-DESTRUCTIVE: no save_system, no logging reconfiguration.
% Follows the conventions of baselineResultDraw.m (2x2 tiledlayout, x axis in
% percent of the cardiac cycle, same colour scheme, exportgraphics 300 dpi).
%
% Signal names are mapped to the CURRENT model:
%   P_BR -> P_ECS,  P_aPVS -> P_PAS,  Q_aPVS -> Q_PAS,  meanQ_GS -> computed
%
% ASCII only.  matlab -batch "run_base_final"

fprintf('=== run_base_final (MATLAB %s) ===\n', version);
tAll = tic;
here   = fileparts(mfilename('fullpath'));
root   = fileparts(fileparts(here));
mdl    = 'CSF_Model';
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(here);
addpath(fullfile(root,'SimulinkModel'));
dataDir = fullfile(root,'SimulinkModel','Results','H_series');
figDir  = fullfile(root,'SimulinkModel','Results','B');
if ~exist(dataDir,'dir'), mkdir(dataDir); end
if ~exist(figDir,'dir'),  mkdir(figDir);  end

CYCLE = 0.8;          % cardiac period (s) -- now exact
TSTOP = 3000;         % s

%% ================= STEP 1: parameters + simulation =================
fprintf('\n########## STEP 1: simulate B, StopTime %g s ##########\n', TSTOP);
SCENARIO_CODE = 'B'; %#ok<NASGU>
Brain_Parameters;
cc_tmp = evalin('base','central_cycle');
tf_tmp = evalin('base','t_fit');
fprintf('  central_cycle %s   t_fit %s   span [%.4f, %.4f]  step %.4f\n', ...
    mat2str(size(cc_tmp)), mat2str(size(tf_tmp)), tf_tmp(1), tf_tmp(end), mean(diff(tf_tmp)));
if ~bdIsLoaded(mdl), load_system(fullfile(root,'SimulinkModel',[mdl '.slx'])); end
tS = tic;
out = sim(mdl,'StopTime',num2str(TSTOP),'SaveOutput','off','SignalLogging','on', ...
              'SignalLoggingName','logsout','ReturnWorkspaceOutputs','on');
fprintf('  sim %.1f s\n', toc(tS));
ls = out.logsout;
f = fullfile(dataDir,'B.mat');
save(f,'ls','-v7.3');
fprintf('  saved %s\n', f);

has = @(s) any(strcmp(ls.getElementNames(), s));
gv  = @(s) getVec(ls, s);

%% ================= STEP 2: cardiac-cycle time series figure =========
fprintf('\n########## STEP 2: cardiac-cycle figure ##########\n');
nCycles  = 2;
prePhase = 0.10;

tV  = gv('P_A');  tvec = getT(ls,'P_A');  dt = mean(diff(tvec)); fs = round(1/dt);
nPtsCycle  = round(CYCLE * fs);
nPreOffset = round(prePhase * nPtsCycle);
nShow      = nCycles * nPtsCycle;
fprintf('  fs %d Hz, points/cycle %d, offset %d, show %d\n', fs, nPtsCycle, nPreOffset, nShow);

sig_PA   = gv('P_A');
sig_PVEN = gv('P_VEN');
sig_PSAS = gv('P_SAS');
sig_PECS = gv('P_ECS');
sig_PPAS = gv('P_PAS');
sig_QPAS = gv('Q_PAS');
sig_QSPI = gv('Q_SPI');
sig_QSAS = gv('Q_SAS');
sig_GSin = gv('Q_GS');
% one-cycle moving average of Q_GS (the model no longer logs meanQ_GS)
nTotal0  = numel(sig_GSin);
kern     = ones(nPtsCycle,1)/nPtsCycle;
meanQ_GS = filter(kern, 1, sig_GSin);
meanQ_GS(1:nPtsCycle-1) = meanQ_GS(nPtsCycle);

nTotal = numel(sig_PA);
searchLen = min((nCycles+4)*nPtsCycle, nTotal);
searchSeg = sig_PA(end-searchLen+1:end);
searchOff = nTotal - searchLen;
[~, troughLoc] = findpeaks(-searchSeg, 'MinPeakDistance', round(nPtsCycle*0.60), ...
                                         'MinPeakProminence', max(searchSeg)*0.05);
troughIdx = searchOff + troughLoc;
fprintf('  troughs found %d\n', numel(troughIdx));
if numel(troughIdx) >= nCycles+1
    anchorTrough = troughIdx(end-nCycles);
else
    anchorTrough = troughIdx(1);
end
startIdx = max(anchorTrough - nPreOffset, 1);
endIdx   = min(startIdx + nShow - 1, nTotal);
fprintf('  window %d..%d (%.3f s)\n', startIdx, endIdx, (endIdx-startIdx+1)/fs);

X  = @(s) s(startIdx:endIdx);
tPct = linspace(0, nCycles*100, endIdx-startIdx+1);
troughXall = prePhase*100 + (0:nCycles-1)*100;

clr_PA = [0.84 0.18 0.15];  clr_PVEN = [0.13 0.47 0.71];
clr_PSAS = [0.17 0.63 0.17]; clr_PECS = [0.20 0.20 0.20];
clr_PPAS = [0.80 0.40 0.00]; clr_QSPI = [0.13 0.47 0.71];
clr_QSAS = [0.17 0.63 0.17]; clr_GSin = [0.58 0.15 0.68];
clr_mGSin = [0.20 0.20 0.20]; lw = 1.5;
xLimAll = [0, nCycles*100];
xtv = 0:25:nCycles*100;
xtlbl = arrayfun(@(v) sprintf('%d%%',v), xtv, 'UniformOutput', false);

figTS = figure('Name','Fig_B_Timeseries','NumberTitle','off', ...
    'Units','centimeters','Position',[5 2 30 24],'Color','white','PaperPositionMode','auto');
tl = tiledlayout(figTS, 2, 2, 'TileSpacing','compact','Padding','compact');

% (a) P_A
ax = nexttile(tl); hold(ax,'on');
plot(ax, tPct, X(sig_PA), '-', 'Color', clr_PA, 'LineWidth', lw);
padY(ax, X(sig_PA));
cycleLines(ax, nCycles, troughXall);
set(ax,'XLim',xLimAll,'XTick',xtv,'XTickLabel',xtlbl,'YGrid','on','Box','off','GridAlpha',0.25);
ylabel(ax,'Pressure (mmHg)'); xlabel(ax,'Cardiac Cycle (%)');
title(ax,'(a)  Arterial Pressure  P_A','FontWeight','bold','FontSize',9);
legend(ax,{'P_A'},'Location','northeast','Box','off','FontSize',8);

% (b) CSF pressures
ax = nexttile(tl); hold(ax,'on');
plot(ax, tPct, X(sig_PVEN), '-', 'Color', clr_PVEN, 'LineWidth', lw);
plot(ax, tPct, X(sig_PSAS), '-', 'Color', clr_PSAS, 'LineWidth', lw);
plot(ax, tPct, X(sig_PECS), '-', 'Color', clr_PECS, 'LineWidth', lw);
plot(ax, tPct, X(sig_PPAS), '-', 'Color', clr_PPAS, 'LineWidth', lw);
padY(ax, [X(sig_PVEN); X(sig_PSAS); X(sig_PECS); X(sig_PPAS)]);
cycleLines(ax, nCycles, troughXall);
set(ax,'XLim',xLimAll,'XTick',xtv,'XTickLabel',xtlbl,'YGrid','on','Box','off','GridAlpha',0.25);
ylabel(ax,'Pressure (mmHg)'); xlabel(ax,'Cardiac Cycle (%)');
title(ax,'(b)  CSF/ISF Compartment Pressures','FontWeight','bold','FontSize',9);
legend(ax,{'P_{VEN}','P_{SAS}','P_{ECS} (ICP)','P_{PAS}'}, ...
    'Location','northeast','Box','off','FontSize',8,'NumColumns',2);

% (c) CSF flows
ax = nexttile(tl); hold(ax,'on');
plot(ax, tPct, X(sig_QPAS), '-', 'Color', clr_PPAS, 'LineWidth', lw);
plot(ax, tPct, X(sig_QSPI), '-', 'Color', clr_QSPI, 'LineWidth', lw);
plot(ax, tPct, X(sig_QSAS), '-', 'Color', clr_QSAS, 'LineWidth', lw);
yline(ax, 0, '-', 'Color',[0.75 0.75 0.75],'LineWidth',0.7,'HandleVisibility','off');
padY(ax, [X(sig_QPAS); X(sig_QSPI); X(sig_QSAS)]);
cycleLines(ax, nCycles, troughXall);
set(ax,'XLim',xLimAll,'XTick',xtv,'XTickLabel',xtlbl,'YGrid','on','Box','off','GridAlpha',0.25);
ylabel(ax,'Flow (mL/min)'); xlabel(ax,'Cardiac Cycle (%)');
title(ax,'(c)  CSF Flow Rates','FontWeight','bold','FontSize',9);
legend(ax,{'Q_{PAS}','Q_{SPI}','Q_{SAS}'},'Location','northeast','Box','off','FontSize',8,'NumColumns',3);

% (d) glymphatic inflow
ax = nexttile(tl); hold(ax,'on');
plot(ax, tPct, X(sig_GSin),    '-',  'Color', clr_GSin,  'LineWidth', lw);
plot(ax, tPct, X(meanQ_GS),    '--', 'Color', clr_mGSin, 'LineWidth', 1.2);
yline(ax, 0, '-', 'Color',[0.75 0.75 0.75],'LineWidth',0.7,'HandleVisibility','off');
segM = X(meanQ_GS);
padY(ax, [X(sig_GSin); segM]);
steadyGS = mean(segM(end-round(nPtsCycle*0.5):end));
text(ax, xLimAll(2)*0.97, steadyGS + (max([X(sig_GSin);segM])-min([X(sig_GSin);segM]))*0.10, ...
    sprintf('mean = %.4f mL/min', steadyGS), 'HorizontalAlignment','right', ...
    'FontSize',7.5,'Color',clr_mGSin,'FontWeight','bold');
cycleLines(ax, nCycles, troughXall);
set(ax,'XLim',xLimAll,'XTick',xtv,'XTickLabel',xtlbl,'YGrid','on','Box','off','GridAlpha',0.25);
ylabel(ax,'Flow (mL/min)'); xlabel(ax,'Cardiac Cycle (%)');
title(ax,'(d)  Glymphatic Inflow','FontWeight','bold','FontSize',9);
legend(ax,{'Q_{GSin}','cycle mean'},'Location','northeast','Box','off','FontSize',8,'NumColumns',2);

exportgraphics(figTS, fullfile(figDir,'Fig_B_Timeseries.pdf'),'ContentType','vector','Resolution',300);
exportgraphics(figTS, fullfile(figDir,'Fig_B_Timeseries.png'),'Resolution',300,'BackgroundColor','white');
fprintf('  wrote Fig_B_Timeseries.pdf/.png\n');

%% ================= STEP 3: convergence figure ======================
fprintf('\n########## STEP 3: convergence figure ##########\n');
SIGC = {'P_ECS','V_SPI','Q_SPI','Q_SAS','Q_SLS','Q_SSS','Q_OLS','Q_out'};
CM = [];
for i = 1:numel(SIGC)
    v = gv(SIGC{i});  t = getT(ls, SIGC{i});
    n = floor(t(end)/CYCLE);
    col = nan(n,1);
    for c = 1:n
        m = (t > (c-1)*CYCLE) & (t <= c*CYCLE);
        if any(m), col(c) = mean(v(m)); end
    end
    CM = [CM, col]; %#ok<AGROW>
end
tc = ((1:size(CM,1))' - 0.5) * CYCLE;      % cycle-centre time (s)
fprintf('  cycles %d\n', size(CM,1));

figCV = figure('Name','Fig_B_Convergence','NumberTitle','off', ...
    'Units','centimeters','Position',[5 2 30 20],'Color','white','PaperPositionMode','auto');
tl2 = tiledlayout(figCV, 2, 2, 'TileSpacing','compact','Padding','compact');

ax = nexttile(tl2);
plot(ax, tc, CM(:,1), '-', 'Color', clr_PECS, 'LineWidth', 1.2);
padY(ax, CM(:,1));
grid(ax,'on'); set(ax,'Box','off','GridAlpha',0.25);
xlabel(ax,'Time (s)'); ylabel(ax,'P_{ECS} cycle mean (mmHg)');
title(ax,'(a)  ICP convergence','FontWeight','bold','FontSize',9);

ax = nexttile(tl2);
plot(ax, tc, CM(:,2), '-', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.2);
padY(ax, CM(:,2));
grid(ax,'on'); set(ax,'Box','off','GridAlpha',0.25);
xlabel(ax,'Time (s)'); ylabel(ax,'V_{SPI} cycle mean (mL)');
title(ax,'(b)  Spinal CSF volume convergence','FontWeight','bold','FontSize',9);

ax = nexttile(tl2); hold(ax,'on');
plot(ax, tc, CM(:,3), '-', 'Color', clr_QSPI, 'LineWidth', 1.2);
plot(ax, tc, CM(:,4), '-', 'Color', clr_QSAS, 'LineWidth', 1.2);
yline(ax, 0, '-', 'Color',[0.75 0.75 0.75],'LineWidth',0.7,'HandleVisibility','off');
padY(ax, [CM(:,3); CM(:,4)]);
grid(ax,'on'); set(ax,'Box','off','GridAlpha',0.25);
xlabel(ax,'Time (s)'); ylabel(ax,'Flow cycle mean (mL/min)');
title(ax,'(c)  Craniospinal exchange','FontWeight','bold','FontSize',9);
legend(ax,{'Q_{SPI}','Q_{SAS}'},'Location','east','Box','off','FontSize',8);

ax = nexttile(tl2); hold(ax,'on');
plot(ax, tc, CM(:,5), '-', 'Color', [0.17 0.63 0.17], 'LineWidth', 1.2);
plot(ax, tc, CM(:,6), '-', 'Color', [0.13 0.47 0.71], 'LineWidth', 1.2);
plot(ax, tc, CM(:,7), '-', 'Color', [0.84 0.18 0.15], 'LineWidth', 1.2);
padY(ax, [CM(:,5); CM(:,6); CM(:,7)]);
grid(ax,'on'); set(ax,'Box','off','GridAlpha',0.25);
xlabel(ax,'Time (s)'); ylabel(ax,'Flow cycle mean (mL/min)');
title(ax,'(d)  Outflow pathways','FontWeight','bold','FontSize',9);
legend(ax,{'Q_{SLS} (spinal)','Q_{SSS} (AG)','Q_{OLS} (olfactory lymph)'}, ...
    'Location','east','Box','off','FontSize',8);

exportgraphics(figCV, fullfile(figDir,'Fig_B_Convergence.pdf'),'ContentType','vector','Resolution',300);
exportgraphics(figCV, fullfile(figDir,'Fig_B_Convergence.png'),'Resolution',300,'BackgroundColor','white');
fprintf('  wrote Fig_B_Convergence.pdf/.png\n');

%% ================= STEP 4: steady-state summary ====================
fprintf('\n########## STEP 4: steady state (mean over last 400 cycles) ##########\n');
W = 400;
fprintf('%-10s %16s\n','signal','mean');
for s = {'P_ECS','P_SAS','P_SPI','P_VEN','P_PAS','V_SPI','V_All','V_VEN','V_ECS', ...
         'Q_SPI','Q_SAS','Q_SLS','Q_SSS','Q_OLS','Q_out','Q_GS','Q_BBB','Q_A','Q_VEN'}
    if ~has(s{1}), continue; end
    v = gv(s{1}); t = getT(ls,s{1});
    n = floor(t(end)/CYCLE);
    m = t > (n-W)*CYCLE;
    fprintf('%-10s %16.6f\n', s{1}, mean(v(m)));
end

fprintf('\n--- per-cycle cardiac swings (last 100 cycles) ---\n');
for s = {'V_SPI','P_ECS','V_Vessel','P_SAS','Q_A','Q_SPI'}
    if ~has(s{1}), continue; end
    v = gv(s{1}); t = getT(ls,s{1});
    n = floor(t(end)/CYCLE); R = nan(100,1);
    for c = (n-99):n
        m = (t > (c-1)*CYCLE) & (t <= c*CYCLE);
        R(c-(n-100)) = max(v(m)) - min(v(m));
    end
    fprintf('  %-10s %14.6f +/- %.6f\n', s{1}, mean(R), std(R));
end

fprintf('\n  Q_GS cycle mean (from fig d) = %.6f mL/min\n', steadyGS);
fprintf('\n=== total wall clock %.1f s ===\n=== DONE ===\n', toc(tAll));

% =====================================================================
function v = getVec(ls, name)
    e = ls.getElement(name);
    if isa(e,'Simulink.SimulationData.Dataset'), ts = e.getElement(1).Values; else, ts = e.Values; end
    v = double(ts.Data(:));
end
function t = getT(ls, name)
    e = ls.getElement(name);
    if isa(e,'Simulink.SimulationData.Dataset'), ts = e.getElement(1).Values; else, ts = e.Values; end
    t = double(ts.Time(:));
end
function padY(ax, v)
    yl = [min(v(:)), max(v(:))];
    s  = yl(2) - yl(1);
    if s < 1e-9, s = max(abs(yl(2)),1); end
    ylim(ax, [yl(1)-s*0.15, yl(2)+s*0.20]);
end
function cycleLines(ax, nCyc, ~)
    yl = ylim(ax);
    for ci = 1:nCyc-1
        xv = ci*100;
        plot(ax, [xv xv], yl, '-', 'Color',[0.55 0.55 0.55], ...
            'LineWidth',1.0,'HandleVisibility','off');
    end
end
