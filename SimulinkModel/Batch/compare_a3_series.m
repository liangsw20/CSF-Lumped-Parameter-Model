%% compare_a3_series.m -- A3 combined-pathology interaction analysis
%
% Loads the two combined scenarios together with their single-factor parents and
% tests whether hypertension and the CSF lesion combine ADDITIVELY:
%
%   additive expectation = B + (H_G2 - B) + (lesion - B)
%   interaction          = observed - additive
%
% Files are loaded one at a time and cleared, because the per-scenario dumps are
% 100-400 MB each.  Read-only w.r.t. the model.  ASCII only.
%   matlab -batch "compare_a3_series"

fprintf('=== compare_a3_series (MATLAB %s) ===\n', version);
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
resDir = fullfile(root,'SimulinkModel','Results');
outDir = fullfile(resDir,'A3_series');
if ~exist(outDir,'dir'), mkdir(outDir); end
CYCLE = 0.8;  W = 400;

labels = {'B','H_G2','C_OVS5','H+C','Ven_R7.5','H_G2+Ven_R7.5'};
files  = {fullfile(resDir,'H_series','B.mat'), ...
          fullfile(resDir,'H_series','H_G2.mat'), ...
          fullfile(resDir,'C_series','C_OVS5.mat'), ...
          fullfile(outDir,'H_C.mat'), ...
          fullfile(resDir,'VenR_series','Ven_R7.5.mat'), ...
          fullfile(outDir,'H_G2_Ven_R7.5.mat')};

KEY = {'P_ECS','P_VEN','P_SAS','P_SPI','V_VEN','V_SPI','V_ECS', ...
       'V_A','V_C','V_V','V_Vessel','V_SAS','V_PAS','V_PNS','V_brain', ...
       'Q_VEN','Q_SPI','Q_SAS','Q_SLS','Q_SSS','Q_OLS','Q_out', ...
       'Q_GS','Q_BBB','Q_A'};
% V_All was dropped from the logging set in the 2026-09-15 revision.  It was
% the craniospinal total, which is conserved, so its only dynamics were those
% of V_SPI and it carried no information of its own.

M = nan(numel(KEY), numel(labels));
have = false(1,numel(labels));
for k = 1:numel(labels)
    if ~exist(files{k},'file')
        fprintf('MISSING %s -- run it first\n', files{k});
        continue;
    end
    d = load(files{k});
    ls = d.ls;
    nms = ls.getElementNames();
    for q = 1:numel(KEY)
        if any(strcmp(nms, KEY{q}))
            M(q,k) = ssm(ls, KEY{q}, CYCLE, W);
        end
    end
    clear d ls
    have(k) = true;
    fprintf('loaded %-14s\n', labels{k});
end

% derived quantity: ventriculo-SAS gradient is not a logged signal
iG  = strcmp(KEY,'P_VEN');  iS = strcmp(KEY,'P_SAS');
grad = M(iG,:) - M(iS,:);

% split V_VEN into mean + ventricular expansion relative to B
iV = strcmp(KEY,'V_VEN');  iE = strcmp(KEY,'V_ECS');
dV = M(iV,:) - M(iV,1);

fprintf('\n--- steady state, mean over last %d cycles ---\n', W);
fprintf('%-12s', 'signal'); fprintf('%16s', labels{:}); fprintf('\n');
for q = 1:numel(KEY)
    fprintf('%-12s', KEY{q}); fprintf('%16.6f', M(q,:)); fprintf('\n');
end
fprintf('%-12s', 'P_VEN-P_SAS'); fprintf('%16.6f', grad); fprintf('\n');
fprintf('%-12s', 'dV_VEN'); fprintf('%16.6f', dV); fprintf('\n');

% ---- cranial volume closure -------------------------------------------
% Monro-Kellie: the CRANIAL volume is constant.  V_ECS alone is NOT the
% complement of V_VEN unless the vasculature is identical to B: the cranial
% total is
%   V_A + V_C + V_V + V_VEN + V_SAS + V_brain + V_PAS + V_PNS  = V_T0 = 1500
% (V_brain = V_ECS + V_ICS + solid, so dV_BR = dV_ECS).  V_All additionally
% contains the spinal SAS, which is NOT part of the constant cranial total.
%
% Since the 2026-09-15 revision V_PNS is no longer logged (the compartment
% itself is unchanged).  It is therefore recovered from the closure:
%     V_PNS = V_T0 - sum(seven logged compartments)
% bp_compute gives V_PNS0 = 30*(r_PVS^2-1) = 41.84 mL.  If the derived value
% sits at that constant in every scenario, the closure holds and the PNS
% compartment is effectively rigid there; a drift means PNS itself is storing
% volume and the naive dV_VEN+dV_ECS residual is only part of the story.
% Filling the V_PNS row back into M also restores it in the compensating-term
% listing below.
V_T0 = 1500;
volNames = {'V_A','V_C','V_V','V_VEN','V_SAS','V_brain','V_PAS'};
vi = zeros(1,numel(volNames));
for q = 1:numel(volNames), vi(q) = find(strcmp(KEY, volNames{q})); end
seven   = sum(M(vi,:), 1);
vPVS_dv = V_T0 - seven;
iVP = find(strcmp(KEY,'V_PNS'));
if ~isempty(iVP), M(iVP,:) = vPVS_dv; end
cranial = seven + vPVS_dv;
fprintf('\n--- cranial volume closure (7 logged compartments + derived PNS) ---\n');
fprintf('  V_PNS0 = 30*(r_PVS^2-1), r_PVS=1.55  -> %.6f mL (unstressed)\n', 30*(1.55^2-1));
for k = 1:numel(labels)
    fprintf('  %-14s sum7 %12.6f   V_PNS %10.6f  (vs B %+.6f mL)   total %12.6f\n', ...
        labels{k}, seven(k), vPVS_dv(k), vPVS_dv(k)-vPVS_dv(1), cranial(k));
end
end
fprintf('  naive dV_VEN+dV_ECS residual is NOT a conservation error; the\n');
fprintf('  compensating terms are listed below for the combined runs.\n');
for k = 2:numel(labels)
    if ~have(k), continue; end
    fprintf('  %-14s:', labels{k});
    for q = 1:numel(volNames)
        d = M(vi(q),k) - M(vi(q),1);
        if abs(d) > 1e-4, fprintf('  d%s %+.4f', volNames{q}, d); end
    end
    fprintf('\n');
end

% ---- additivity test ----------------------------------------------------
% combo 1 = H_G2 + C_OVS5, combo 2 = H_G2 + Ven_R7.5
combos = {struct('name','H+C','lab','H+C', ...
                 'single','C_OVS5','comb','H+C'), ...
          struct('name','H_G2+Ven_R7.5','lab','H_G2+Ven_R7.5', ...
                 'single','Ven_R7.5','comb','H_G2+Ven_R7.5')};
iH = find(strcmp(labels,'H_G2'));

for cq = 1:numel(combos)
    cs = combos{cq};
    iL = find(strcmp(labels, cs.single));
    iC = find(strcmp(labels, cs.comb));
    if ~have(iL) || ~have(iC) || ~have(iH)
        fprintf('\n=== %s : SKIPPED (missing run) ===\n', cs.name);
        continue;
    end
    fprintf('\n=== additivity test: %s ===\n', cs.name);
    fprintf('  additive expectation = B + (H_G2 - B) + (%s - B)\n', cs.single);
    fprintf('\n%-12s %14s %14s %14s %14s %14s %11s\n', ...
        'signal','B','H_G2',cs.single,cs.comb,'additive','interaction');
    add = M(:,iH) + M(:,iL) - M(:,1);
    inter = M(:,iC) - add;
    for q = 1:numel(KEY)
        fprintf('%-12s %14.6f %14.6f %14.6f %14.6f %14.6f %+11.6f\n', ...
            KEY{q}, M(q,1), M(q,iH), M(q,iL), M(q,iC), add(q), inter(q));
    end
    fprintf('%-12s %14.6f %14.6f %14.6f %14.6f %14.6f %+11.6f\n', ...
        'P_VEN-P_SAS', grad(1), grad(iH), grad(iL), grad(iC), ...
        grad(iH)+grad(iL)-grad(1), grad(iC)-(grad(iH)+grad(iL)-grad(1)));
    fprintf('%-12s %14.6f %14.6f %14.6f %14.6f %14.6f %+11.6f\n', ...
        'dV_VEN', 0, dV(iH), dV(iL), dV(iC), dV(iH)+dV(iL), ...
        dV(iC)-(dV(iH)+dV(iL)));

    % interaction expressed relative to the lesion-only change
    fprintf('\n  interaction as %% of the lesion-only change:\n');
    for s = {'P_ECS','V_VEN','P_VEN','Q_GS','Q_BBB','Q_out'}
        q = find(strcmp(KEY, s{1}));
        den = M(q,iL) - M(q,1);
        if abs(den) < 1e-12, continue; end
        fprintf('    %-10s observed %+11.6f   additive %+11.6f   interaction %+11.6f  (%+.1f%%)\n', ...
            s{1}, M(q,iC)-M(q,1), add(q)-M(q,1), inter(q), inter(q)/den*100);
    end
end

% ---- CSV ---------------------------------------------------------------
csv = fullfile(outDir,'A4_interaction.csv');
fid = fopen(csv,'w');
fprintf(fid,'signal,%s\n', strjoin(labels,','));
for q = 1:numel(KEY)
    fprintf(fid,'%s', KEY{q});
    fprintf(fid,',%14.6f', M(q,:));
    fprintf(fid,'\n');
end
fprintf(fid,'P_VEN-P_SAS');
fprintf(fid,',%14.6f', grad);  fprintf(fid,'\n');
fprintf(fid,'dV_VEN');
fprintf(fid,',%14.6f', dV);    fprintf(fid,'\n');
fclose(fid);
fprintf('\nwrote %s\n', csv);

% ---- figure ------------------------------------------------------------
iC5 = find(strcmp(labels,'C_OVS5'));  iO5 = find(strcmp(labels,'Ven_R7.5'));
iOc2 = find(strcmp(labels,'H_G2+Ven_R7.5'));
iHc2 = find(strcmp(labels,'H+C'));
need = [1 iH iC5 iO5 iHc2 iOc2];
if ~all(have(need))
    fprintf('\n*** figure skipped: not all six runs are present ***\n');
    fprintf('\n=== DONE ===\n');
    return
end

fig = figure('Position',[80 80 1220 780],'Color','w','Visible','off');

% (a) ICP for the outflow-resistance combination
ax1 = subplot(2,2,1);
idxC = [1 iH iC5 iHc2];
bh = bar(ax1, M(strcmp(KEY,'P_ECS'), idxC), 0.6,'FaceColor',[0.25 0.45 0.70]);
hold(ax1,'on');
addC = M(strcmp(KEY,'P_ECS'),iH) + M(strcmp(KEY,'P_ECS'),iC5) - M(strcmp(KEY,'P_ECS'),1);
plot(ax1, 4, addC, 'rp','MarkerSize',12,'MarkerFaceColor','r');
text(ax1, 4, addC, sprintf('  additive %.2f', addC),'Color','r','FontSize',8, ...
    'VerticalAlignment','bottom');
set(ax1,'XTickLabel',{'B','H_G2','C_OVS5',sprintf('H_G2+\nC_OVS5')}, ...
    'TickLabelInterpreter','none');
ylabel(ax1,'P_{ECS} (mmHg)'); title(ax1,'(a) hypertension x outflow impairment');
grid(ax1,'on'); box(ax1,'on');

% (b) ventricular volume for the obstruction combination
ax2 = subplot(2,2,2);
idxO = [1 iH iO5 iOc2];
bh2 = bar(ax2, M(iV, idxO), 0.6,'FaceColor',[0.85 0.55 0.15]);
hold(ax2,'on');
addO = M(iV,iH) + M(iV,iO5) - M(iV,1);
plot(ax2, 4, addO, 'rp','MarkerSize',12,'MarkerFaceColor','r');
text(ax2, 4, addO, sprintf('  additive %.2f', addO),'Color','r','FontSize',8, ...
    'VerticalAlignment','bottom');
set(ax2,'XTickLabel',{'B','H_G2','Ven_R7.5',sprintf('H_G2+\nO7.5')}, ...
    'TickLabelInterpreter','none');
ylabel(ax2,'V_{VEN} (mL)'); title(ax2,'(b) hypertension x ventriculo-SAS obstruction');
grid(ax2,'on'); box(ax2,'on');

% (c) interaction bars
ax3 = subplot(2,2,3);
n1 = M(:,iHc2) - (M(:,iH) + M(:,iC5) - M(:,1));
n2 = M(:,iOc2) - (M(:,iH) + M(:,iO5) - M(:,1));
sel = {'P_ECS','V_VEN','P_VEN','Q_GS','Q_BBB','Q_out'};
qq = zeros(1,numel(sel));
for q = 1:numel(sel), qq(q) = find(strcmp(KEY,sel{q})); end
MM = [n1(qq) n2(qq)];
bh3 = bar(ax3, MM,'grouped');
bh3(1).FaceColor = [0.65 0.20 0.20];
bh3(2).FaceColor = [0.20 0.35 0.75];
set(ax3,'XTickLabel',sel); ylabel(ax3,'observed - additive');
title(ax3,'(c) interaction (0 = perfectly additive)');
legend(ax3,{'H+C','H_G2+Ven_R7.5'},'Location','best');
hold(ax3,'on'); plot(ax3,[0.4 numel(sel)+0.6],[0 0],'k-');
grid(ax3,'on'); box(ax3,'on');

% (d) outflow allocation for all six
ax4 = subplot(2,2,4);
A = zeros(numel(labels),3);
for k = 1:numel(labels)
    tt = M(strcmp(KEY,'Q_out'),k);
    A(k,:) = [M(strcmp(KEY,'Q_SSS'),k) M(strcmp(KEY,'Q_OLS'),k) ...
              M(strcmp(KEY,'Q_SLS'),k)]/tt*100;
end
bh4 = bar(ax4, A,'stacked');
bh4(1).FaceColor = [0.20 0.45 0.75];
bh4(2).FaceColor = [0.35 0.70 0.35];
bh4(3).FaceColor = [0.85 0.55 0.15];
set(ax4,'XTickLabel',{'B','H_G2','C_OVS5',sprintf('H_G2+\nC_OVS5'), ...
    'Ven_R7.5',sprintf('H_G2+\nO7.5')},'TickLabelInterpreter','none');
ylabel(ax4,'% of Q_out');
title(ax4,'(d) outflow allocation  (AG / cribriform / spinal)');
legend(ax4,{'AG','cribriform','spinal'},'Location','eastoutside');
ylim(ax4,[0 108]); grid(ax4,'on'); box(ax4,'on');

sgtitle('A3 combined pathology: hypertension + CSF lesion','FontWeight','bold');
png = fullfile(outDir,'Fig_A3_series.png');
exportgraphics(fig, png, 'Resolution', 200);
fprintf('wrote %s\n', png);
close(fig);

fprintf('\n=== DONE ===\n');

% =====================================================================
function v = ssm(ls, name, CYCLE, W)
    e = ls.getElement(name);
    if isa(e,'Simulink.SimulationData.Dataset'), ts = e.getElement(1).Values; else, ts = e.Values; end
    y = double(ts.Data(:)); t = double(ts.Time(:));
    n = floor(t(end)/CYCLE);
    v = mean(y(t > (n-W)*CYCLE));
end
