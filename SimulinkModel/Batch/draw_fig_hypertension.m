%% ========================================================
%  FIGURE 4: Hypertension impact on intracranial fluid dynamics
%  Five states: Baseline / High-normal / Grade 1 / Grade 2 / Acute G2
%
%  Style follows MATLAB_script/DrawResult0520/ResultDraw.m (FIGURE 1):
%  same C palette, same bar conventions, same %-annotation and reference line.
%  Font sizes carry the same 1.5x scale as Fig3, with the figure widened to
%  34 x 20 cm so that the printed text size matches the other figures
%  (both land at ~6.5-7 pt when scaled to \textwidth).
%
%  Panel order follows the causal chain:
%    (a) vascular volume pulsatility  -- the driving force
%    (b) glymphatic inflow Q_GS             -- the impaired transport
%    (c) PAS reverse flow Q_PAS       -- the mechanism (baseline vs chronic vs acute)
%    (d) BBB permeation Q_BBB         -- the compensatory influx
%    (e) intracranial pressure P_ECS  -- the resulting pressure
%    (f) total efflux distribution    -- the downstream redistribution
%
%  Output: PDF only, into <root>/新插图/, together with a copy of this script.
%
%  Usage: matlab -batch "draw_fig_hypertension"
% ========================================================
function draw_fig_hypertension(matDir, outPdf)
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(here); addpath(fullfile(root,'SimulinkModel'));

if nargin < 1 || isempty(matDir)
    matDir = fullfile(root,'SimulinkModel','Results','H_series');
end
if nargin < 2 || isempty(outPdf)
    outPdf = fullfile(root,'新插图','Fig4_Hypertension.pdf');
end
if ~exist(fileparts(outPdf),'dir'), mkdir(fileparts(outPdf)); end

STATES = {'B','H_HM','H_G1','H_G2','H_G2A'};
LABELS = STATES;           % 横轴直接用论文 Table 2 的实验代号
LABELS_TEX = strrep(LABELS,'_','\_');   % tex 解释器会把 _ 当上下标，显示时转义
WAVE   = [1 4 5];          % states overlaid in panel (c): B, Grade 2, Acute G2
CYCLE  = 0.8;              % cardiac period (s)
TAVG   = 80;               % s averaged for the steady-state means
FONT_SCALE = 1.5;
fs_tick = 12; fs_axlab = 15; fs_title = 9*FONT_SCALE;
fs_leg  = 12; fs_pct = 7.5*FONT_SCALE; fs_layer = 7*FONT_SCALE;

fprintf('=== draw_fig_hypertension ===\n  数据目录: %s\n  输出: %s\n', matDir, outPdf);

%% ── 颜色（与 DrawResult0520/ResultDraw.m 完全一致）─────
C.white = [1 1 1];
C.base  = [0.35 0.35 0.35];
C.baseBar = [0.45 0.45 0.45];
barLightness = 0.2;
lighten = @(rgb,f) rgb + (C.white - rgb) * f;

C.htn.highNormal = lighten([0.92 0.78 0.42], barLightness);   % 暖金黄
C.htn.light      = lighten([0.88 0.55 0.45], barLightness);   % 浅红/珊瑚
C.htn.main       = lighten([0.75 0.28 0.24], barLightness);   % 深红 (Grade 2)
C.htn.acute      = lighten([0.85 0.52 0.25], barLightness);   % 深橙 (Acute)
C.series.htn = [C.base; C.htn.highNormal; C.htn.light; C.htn.main; C.htn.acute];

pathSaturation = 0.4;
C.path.lymph    = hsv2rgb([0.33, pathSaturation, 0.75]);
C.path.ag       = hsv2rgb([0.58, pathSaturation, 0.75]);
C.path.spinal   = hsv2rgb([0.78, pathSaturation, 0.75]);

C.bar.edge  = [0.35 0.35 0.35];  C.bar.edgeWidth = 0.6;
C.label.percent = [0.1 0.1 0.1];
C.line.ref  = [0.40 0.40 0.40];
C.line.zero = [0.78 0.78 0.78];

%% ── 读取五个状态并计算指标 ────────────────────────────
nS = numel(STATES);
M = struct();
for k = 1:nS
    f = fullfile(matDir, [STATES{k} '.mat']);
    if ~exist(f,'file'), error('draw_fig_hypertension:noFile','缺少 %s', f); end
    S = load(f, 'ls');   ls = S.ls;
    t = local_sig(ls,'P_A');  T = local_time(ls,'P_A');
    iAvg = T >= T(end) - TAVG;

    M(k).code    = STATES{k};
    M(k).qGSin   = local_mean(ls,'Q_GS', iAvg);
    M(k).qBBB    = local_mean(ls,'Q_BBB',  iAvg);
    M(k).pECS    = local_mean(ls,'P_ECS',  iAvg);
    M(k).qOLS    = local_mean(ls,'Q_OLS',  iAvg);
    M(k).qSSS    = local_mean(ls,'Q_SSS',  iAvg);
    M(k).qSLS    = local_mean(ls,'Q_SLS',  iAvg);
    M(k).qOut    = M(k).qOLS + M(k).qSSS + M(k).qSLS;

    % 搏动幅度（末 40 s）
    iSw = T >= T(end) - 40;
    M(k).swVVes  = local_swing(ls,'V_Vessel', iSw);
    M(k).swVT    = local_swing(ls,'V_T',      iSw);
    M(k).swPPAS  = local_swing(ls,'P_PAS',    iSw);

    % 单个心动周期的 Q_PAS 波形（与 Fig3 相同的谷点对齐）
    nPtsCycle = round(CYCLE * round(1/mean(diff(T))));
    tr = local_troughs(t, nPtsCycle);
    a  = tr(max(1,end-1));
    i0 = max(a - round(0.10*nPtsCycle), 1);
    i1 = min(i0 + nPtsCycle - 1, numel(t));
    M(k).qPASc   = local_sig(ls,'Q_PAS');
    M(k).qPASc   = M(k).qPASc(i0:i1);
    M(k).qPASmin = min(M(k).qPASc);
    fprintf('  %-7s Q_GS %.5f  Q_BBB %.5f  P_ECS %8.4f  dV_Vessel %.4f  dV_T %.4f  minQ_PAS %+.4f  swP_PAS %.3f\n', ...
        STATES{k}, M(k).qGSin, M(k).qBBB, M(k).pECS, M(k).swVVes, M(k).swVT, M(k).qPASmin, M(k).swPPAS);
end
base = M(1);
xpos = 1:nS;
xCyc = linspace(0, 100, numel(M(1).qPASc));

%% ── 图形 ──────────────────────────────────────────────
set(groot,'DefaultAxesFontName','Arial'); set(groot,'DefaultTextFontName','Arial');
fig = figure('Name','Fig4_Hypertension','NumberTitle','off', ...
    'Units','centimeters','Position',[5 2 34 20], ...
    'Color',C.white,'PaperPositionMode','auto');
set(fig,'DefaultAxesFontSize',fs_tick,'DefaultTextFontSize',fs_tick);
tl = tiledlayout(fig, 2, 3, 'TileSpacing','compact','Padding','compact');

%% ── (a) 血管容积搏动 ─────────────────────────────────
% 审稿意见：条形图从 0 起（本量 0 = 无搏动），避免截断夸张
ax = nexttile;  vals = [M.swVVes]';
bar(ax, xpos, vals, 'FaceColor','flat','EdgeColor',C.bar.edge, ...
    'LineWidth',C.bar.edgeWidth,'BarWidth',0.65,'CData',C.series.htn);
local_barAxis(ax, vals, xpos, LABELS_TEX, C, fs_tick, fs_pct, 0);
ylabel(ax,'\Delta V_{Vessel} (mL)');
title(ax,'(a)  Vascular Volume Pulsatility','FontWeight','bold','FontSize',fs_title);

%% ── (b) 类淋巴流入 ───────────────────────────────────
% 审稿意见：面板 a 类的流量条形允许偏移轴（沿用原图 a 的处理）
ax = nexttile;  vals = [M.qGSin]';
bar(ax, xpos, vals, 'FaceColor','flat','EdgeColor',C.bar.edge, ...
    'LineWidth',C.bar.edgeWidth,'BarWidth',0.65,'CData',C.series.htn);
local_barAxis(ax, vals, xpos, LABELS_TEX, C, fs_tick, fs_pct, NaN);
ylabel(ax,'Flow Rate (mL/min)');
title(ax,'(b)  Glymphatic Inflow (Q_{GS})','FontWeight','bold','FontSize',fs_title);

%% ── (c) Q_PAS 倒流（单周期波形叠加）──────────────────
ax = nexttile;  hold(ax,'on');
yline(ax, 0, '-', 'Color',C.line.zero, 'LineWidth',0.8, 'HandleVisibility','off');
h = gobjects(1,numel(WAVE));
for j = 1:numel(WAVE)
    k = WAVE(j);
    h(j) = plot(ax, xCyc, M(k).qPASc, '-', 'Color', C.series.htn(k,:), 'LineWidth', 1.5);
end
allV = [M(WAVE).qPASc];
local_ylim(ax, allV(:));
set(ax,'XTick',0:20:100,'XTickLabel',arrayfun(@(v) sprintf('%d%%',v),0:20:100,'UniformOutput',false), ...
    'YGrid','off','XGrid','off','Box','off');
xtickangle(ax,0);
xlabel(ax,'Cardiac Cycle (%)');
ylabel(ax,'Q_{PAS} (mL/min)');
% 图例括号内给出该状态的正向、逆向平均流量
legTxt = cell(1,numel(WAVE));
for j = 1:numel(WAVE)
    v = M(WAVE(j)).qPASc;
    legTxt{j} = sprintf('%s  (+%.2f / %.2f)', LABELS_TEX{WAVE(j)}, mean(v(v>0)), mean(v(v<0)));
end
legend(ax, h, legTxt, 'Location','best','Box','off','FontSize',fs_leg);   % 与 Fig5(b) 图例同字号
title(ax,'(c)  PAS Reverse Flow (Q_{PAS})','FontWeight','bold','FontSize',fs_title);

%% ── (d) BBB 通透 ─────────────────────────────────────
% 审稿意见明确要求：这一格必须在 0 处相交
ax = nexttile;  vals = [M.qBBB]';
bar(ax, xpos, vals, 'FaceColor','flat','EdgeColor',C.bar.edge, ...
    'LineWidth',C.bar.edgeWidth,'BarWidth',0.65,'CData',C.series.htn);
local_barAxis(ax, vals, xpos, LABELS_TEX, C, fs_tick, fs_pct, 0);
ylabel(ax,'Flow Rate (mL/min)');
title(ax,'(d)  BBB Permeation (Q_{BBB})','FontWeight','bold','FontSize',fs_title);

%% ── (e) 颅内压 ───────────────────────────────────────
% 审稿意见：非零交点必须"有助于突出各列差异"。ICP 全部落在 10.3–10.8 mmHg，
% 从 0 起三根柱子几乎等高；这里取整 10 mmHg 作为交点并在图注中说明。
ax = nexttile;  vals = [M.pECS]';
bar(ax, xpos, vals, 'FaceColor','flat','EdgeColor',C.bar.edge, ...
    'LineWidth',C.bar.edgeWidth,'BarWidth',0.65,'CData',C.series.htn);
local_barAxis(ax, vals, xpos, LABELS_TEX, C, fs_tick, fs_pct, 10);
ylabel(ax,'ICP (mmHg)');
title(ax,'(e)  Intracranial Pressure (P_{ECS})','FontWeight','bold','FontSize',fs_title);

%% ── (f) 流出路径堆叠图 ───────────────────────────────
ax = nexttile;  hold(ax,'on');
stackMat = [[M.qOLS]', [M.qSSS]', [M.qSLS]'];
layerLabels = {'Q_{OLS}','Q_{SSS}','Q_{SLS}'};
layerColors = [C.path.lymph; C.path.ag; C.path.spinal];
bh = bar(ax, xpos, stackMat, 'stacked', 'EdgeColor',C.bar.edge, ...
         'LineWidth',C.bar.edgeWidth, 'BarWidth',0.65);
for li = 1:numel(bh), bh(li).FaceColor = layerColors(li,:); end
totH = sum(stackMat,2)';
yline(ax, totH(1), '--', 'Color',C.line.ref, 'LineWidth',1.2);
ys = max(totH) - min(totH);  if ys < 1e-9, ys = max(totH)*0.3 + 0.01; end
ylim(ax, [0, max(totH) + ys*2.6]);            % 上方多留一点给居中图例
dyChar = local_charData(ax, fs_pct);          % 半个字符的抬升量（data 单位）
for k = 2:nS
    text(ax, xpos(k), totH(k) + 0.2*ys + 0.5*dyChar, ...
        sprintf('%+.0f%%', 100*(totH(k)/totH(1)-1)), ...
        'HorizontalAlignment','center','FontSize',fs_pct,'Color',C.label.percent,'FontWeight','bold');
end
cumBot = zeros(1,nS);
for li = 1:size(stackMat,2)
    layH = stackMat(:,li)';
    for k = 1:nS
        if layH(k) > totH(k)*0.10
            text(ax, xpos(k), cumBot(k)+layH(k)/2, sprintf('%.0f%%', 100*layH(k)/totH(k)), ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontSize',fs_layer,'Color',C.white,'FontWeight','bold');
        end
    end
    cumBot = cumBot + layH;
end
set(ax,'XTick',xpos,'XTickLabel',LABELS_TEX,'XTickLabelRotation',0, ...
    'YGrid','off','Box','off');
ylabel(ax,'Total Outflow Rate (mL/min)');
title(ax,'(f)  Total Efflux Pathway Distribution','FontWeight','bold','FontSize',fs_title);
legend(ax, bh, layerLabels, 'Location','north','FontSize',fs_leg,'Box','off', ...
       'Orientation','horizontal','NumColumns',3);

%% ── 字号自检 ─────────────────────────────────────────
fprintf('  字号(pt): 刻度 %.1f | 轴标题 %.1f | 面板标题 %.1f | %%标注 %.1f\n', ...
    fs_tick, fs_axlab, fs_title, fs_pct);

%% ── 导出 + 复制脚本 ──────────────────────────────────
try, java.lang.System.setProperty('user.home', getenv('TEMP')); catch, end
exportgraphics(fig, outPdf, 'ContentType','vector','Resolution',300);
fprintf('✓ 图已生成: %s\n', outPdf);
try
    thisFile = [mfilename('fullpath') '.m'];
    if ~strcmpi(fileparts(thisFile), fileparts(outPdf))
        copyfile(thisFile, fullfile(fileparts(outPdf),'draw_fig_hypertension.m'));
        fprintf('  脚本副本: %s\n', fullfile(fileparts(outPdf),'draw_fig_hypertension.m'));
    end
catch err
    fprintf('  脚本复制失败: %s\n', err.message);
end

%% ── 控制台汇总（供写图注）────────────────────────────
fprintf('\n--- 相对基线(%%) ---\n');
fprintf('  %-10s %8s %8s %8s %8s %8s %8s\n','state','Q_GS','Q_BBB','P_ECS','dV_Ves','dV_T','Q_out');
for k = 2:nS
    fprintf('  %-10s %+8.1f %+8.1f %+8.1f %+8.1f %+8.1f %+8.1f\n', STATES{k}, ...
        100*(M(k).qGSin/base.qGSin-1), 100*(M(k).qBBB/base.qBBB-1), ...
        100*(M(k).pECS/base.pECS-1),   100*(M(k).swVVes/base.swVVes-1), ...
        100*(M(k).swVT/base.swVT-1),   100*(M(k).qOut/base.qOut-1));
end
fprintf('\n--- 流出占比 (%%) ---\n');
for k = 1:nS
    fprintf('  %-10s OLS %5.2f | SSS %5.2f | SLS %5.2f  (总 %.5f mL/min)\n', ...
        STATES{k}, 100*M(k).qOLS/M(k).qOut, 100*M(k).qSSS/M(k).qOut, ...
        100*M(k).qSLS/M(k).qOut, M(k).qOut);
end
fprintf('\n--- (c) 面板统计 ---\n');
for k = WAVE
    fprintf('  %-10s min Q_PAS %+.4f mL/min | P_PAS 振荡 %.3f mmHg | 倒流时长占比 %.1f%%\n', ...
        STATES{k}, M(k).qPASmin, M(k).swPPAS, 100*mean(M(k).qPASc<0));
end
fprintf('=== DONE ===\n');
end

%% ════════════════════════════════════════════════════════
%  辅助函数
%% ════════════════════════════════════════════════════════
function v = local_sig(ls, name)
ts = ls.getElement(name).Values;
v  = double(ts.Data(:));
end

function t = local_time(ls, name)
ts = ls.getElement(name).Values;
t  = double(ts.Time(:));
end

function m = local_mean(ls, name, idx)
v = local_sig(ls, name);
m = mean(v(idx));
end

function s = local_swing(ls, name, idx)
v = local_sig(ls, name); v = v(idx);
s = max(v) - min(v);
end

function tr = local_troughs(PA, nPtsCycle)
% 与 Fig3 相同的谷点检测
searchLen = min(6*nPtsCycle, numel(PA));
seg = PA(end-searchLen+1:end);
[~, loc] = findpeaks(-seg, 'MinPeakDistance', round(nPtsCycle*0.60), ...
                     'MinPeakProminence', max(seg)*0.05);
tr = (numel(PA)-searchLen) + loc;
end

function local_ylim(ax, vals)
vals = vals(~isnan(vals));
sp = max(vals) - min(vals);  if sp < 1e-9, sp = max(abs(vals))*0.1+0.01; end
ylim(ax, [min(vals)-sp*0.12, max(vals)+sp*0.18]);
end

function dy = local_charData(ax, fs)
% 一个字符的高度折算成该坐标轴的 data 单位
fig = ancestor(ax,'figure');  figCm = fig.Position(3:4);
axCm = ax.Position(4) * figCm(2);
dy = (fs*2.54/72) / axCm * diff(ylim(ax));
end

function local_barAxis(ax, vals, xpos, labels, C, fs_tick, fs_pct, yStart)
% 与 DrawResult0520 fig1 一致的条形样式。
%   yStart = 0   条形从 0 起（面积可比）
%   yStart = 数值 指定交点（需在图注中说明）
%   yStart = NaN  自动（数据最小值下方留 10% 跨度）
% 参考虚线先于文字绘制，保证它落在文字下层。
if nargin < 8, yStart = NaN; end
hold(ax,'on');
base_v = vals(1);
sp = max(vals) - min(vals);  if sp < 1e-9, sp = abs(base_v)*0.2+0.01; end
yline(ax, base_v, '--', 'Color',C.line.ref, 'LineWidth',1);
if isnan(yStart)
    y0 = min(vals) - sp*0.10;
elseif yStart == 0
    y0 = 0;
else
    y0 = yStart;
end
y1 = max(vals) + max(sp*0.40, (max(vals)-y0)*0.16);
ylim(ax, [y0, y1]);
off = 0.5 * local_charData(ax, fs_pct);       % 标注离柱顶半个字符
for k = 2:numel(vals)
    text(ax, xpos(k), vals(k) + off, sprintf('%+.0f%%', 100*(vals(k)/base_v-1)), ...
        'HorizontalAlignment','center','FontSize',fs_pct, ...
        'Color',C.label.percent,'FontWeight','bold');
end
set(ax,'XTick',xpos,'XTickLabel',labels,'XTickLabelRotation',0, ...
    'YGrid','off','XGrid','off','Box','off');
xtickangle(ax, 0);
end
