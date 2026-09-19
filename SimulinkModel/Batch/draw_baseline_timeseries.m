%% ========================================================
%  Baseline Time-series Visualization
%  Last N cardiac cycles, x-axis = % of cardiac cycle
%  Cycle start = 10% before an arterial pressure trough
%
%  Style follows MATLAB_script/DrawResult0520/baselineResultDraw.m
%  (40 x 24 cm, tiledlayout 2x3 'loose', hsv-based palette, % x axis).
%
%  Panel content (paper figure, baseline state):
%    (a) radial (measured) and central (model inlet P_A) arterial pressure
%    (b) the five CSF/ISF compartment pressures
%    (c) cranial SAS flow Q_SAS and spinal flow Q_SPI
%    (d) change of the compartment volumes about their cycle mean (dV_SPI, dV_SAS, dV_T)
%    (e) glymphatic flow Q_PAS, Q_GS, Q_ECS, Q_BBB（图中显示名；模型信号名为 Q_GS / Q_ECS）
%    (f) steady-state distribution of the three outflow pathways
%
%  Output: PDF only.
%
%  Usage
%    matlab -batch "draw_baseline_timeseries"
%    matlab -batch "draw_baseline_timeseries('in.mat','out.pdf')"
% ========================================================
function draw_baseline_timeseries(matFile, outPdf)
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(here); addpath(fullfile(root,'SimulinkModel'));

%% ── 参数设置 ──────────────────────────────────────────
nCycles  = 2;      % 展示的周期数
cycleDur = 0.8;    % 周期时长 (s)，75 BPM
prePhase = 0.10;   % 谷点前偏移比例

% Panels drawn on a split left/right y axis; empty = every panel a single axis.
% e.g. SPLIT_ON = {'b','d','e'};  the largest swing stays left, anything more
% than 8x smaller moves right.
SPLIT_ON = {};
SHOW_CYCLE_LINES = false;

% 全局字号倍率（1.5 = 放大 1.5 倍；标题/图例/饼图标签与坐标轴字号一起缩放）
FONT_SCALE = 1.5;
% 面板标题 / 图例 / 轴标题按需求再放大 1.2 倍（只改这三类，刻度与标注字号不动）
FONT_ZOOM  = 1.2;
fs_axis = 10 * FONT_SCALE;            % 坐标轴与刻度（MATLAB 默认 10）
fs_t    = 9  * FONT_SCALE * FONT_ZOOM;   % 面板标题
fs_l    = 8  * FONT_SCALE * FONT_ZOOM;   % 图例
fs_p    = 7.5* FONT_SCALE * FONT_ZOOM;   % 饼图标签（按需求一并放大）
% 刻度标签与轴标题字号（原始大小 10 pt 为基准）：
%   横/纵轴刻度统一用 fs_xtick（10*1.2*1.2），横轴刻度按 50% 间隔只保留 5 个；
%   横轴标题 Cardiac Cycle (%) 10*1.5*1.2，纵轴标题与横轴标题同尺寸
fs_xtick = 12 * FONT_ZOOM;            % 横/纵轴刻度标签（两轴统一）
fs_axlab = 15 * FONT_ZOOM;            % 横/纵轴标题
PAD_TOP_A = 0.25;            % (a) 曲线以上额外留白，给居中的图例让位
PAD_TOP_B = 0.45;            % (b) 曲线以上额外留白，给图例让位（其余面板 0.20）
PAD_TOP_C = 0.20;            % (c) 同上
% (a)/(b) 图例水平微调：单位 = 字符宽度（正 = 右移）。
% 居中之后按版面观感微调，数值越小越接近正中。
LEG_SHIFT_A =  0.5;          % (a) 右移半个字符
LEG_SHIFT_B = -1.0;          % (b) 左移一个字符
PIE_INSIDE = 0.62;           % (f) 标注沿半径向圆心缩放的系数：1 = 饼图外(默认)，0.6 左右 = 饼图内
% (d) 图例平移：LEG_SHIFT_D 为 figure 归一化单位的右移量；LEG_DROP_D 为"字符高度"倍数（负 = 下移）
LEG_SHIFT_D = 0.012;
LEG_DROP_D  = -0.5;
% (f) 饼图内文字微调，单位 = 字符宽度/字符高度倍数
PIE_NUDGE_OLS_X = 2;         % Q_{OLS} 右移
PIE_NUDGE_SSS_Y = 2;         % Q_{SSS} 上移

if nargin < 1 || isempty(matFile), matFile = localFindBaseline(root); end
if nargin < 2 || isempty(outPdf)
    outPdf = fullfile(root,'新插图','Fig3_Baseline_Timeseries.pdf');
end
if ~exist(fileparts(outPdf),'dir'), mkdir(fileparts(outPdf)); end
fprintf('=== draw_baseline_timeseries ===\n  数据: %s\n  输出: %s\n', matFile, outPdf);

%% ── 读取时间轴和采样率 ────────────────────────────────
S = load(matFile);
if ~isfield(S,'ls'), error('draw_baseline_timeseries:noLs','%s 中没有 logsout 数据集 ls', matFile); end
ds = S.ls;

tVec = getTime(ds, 'P_A');
dt   = mean(diff(tVec));
fs   = round(1 / dt);
fprintf('检测到采样率: %d Hz，采样间隔: %.4f s\n', fs, dt);

nPtsCycle  = round(cycleDur * fs);
nPreOffset = round(prePhase * nPtsCycle);
nShow      = nCycles * nPtsCycle;

%% ── 读取各信号 ────────────────────────────────────────
sig_PA   = getSignal(ds, 'P_A');
sig_PVEN = getSignal(ds, 'P_VEN');
sig_PSAS = getSignal(ds, 'P_SAS');
sig_PECS = getSignal(ds, 'P_ECS');
sig_PPAS = getSignal(ds, 'P_PAS');
sig_PSPI = getSignal(ds, 'P_SPI');
sig_QSAS = getSignal(ds, 'Q_SAS');
sig_QSPI = getSignal(ds, 'Q_SPI');
sig_QPAS = getSignal(ds, 'Q_PAS');
sig_QGin = getSignal(ds, 'Q_GS');
sig_QGot = getSignal(ds, 'Q_ECS');
sig_QBBB = getSignal(ds, 'Q_BBB');
sig_VT   = getSignal(ds, 'V_T');
sig_VSAS = getSignal(ds, 'V_SAS');
sig_VSPI = getSignal(ds, 'V_SPI');
sig_QOLS = getSignal(ds, 'Q_OLS');
sig_QSSS = getSignal(ds, 'Q_SSS');
sig_QSLS = getSignal(ds, 'Q_SLS');

nTotal = numel(sig_PA);
fprintf('信号总长度: %d 点 (%.2f s)\n', nTotal, nTotal/fs);

%% ── 谷点检测：在末尾足够长片段内找 PA 极小值 ─────────
searchLen = min((nCycles + 4) * nPtsCycle, nTotal);
searchSeg = sig_PA(end - searchLen + 1 : end);
searchOff = nTotal - searchLen;

[~, troughLoc] = findpeaks(-searchSeg, ...
    'MinPeakDistance',   round(nPtsCycle * 0.60), ...
    'MinPeakProminence', max(searchSeg) * 0.05);

troughIdx = searchOff + troughLoc;
fprintf('检测到 %d 个谷点\n', numel(troughIdx));

%% ── 以谷点前 10% 为起始，截取 nCycles 个完整周期 ──────
if numel(troughIdx) >= nCycles + 1
    anchorTrough = troughIdx(end - nCycles);
else
    warning('谷点数量不足，使用最早检测到的谷点');
    anchorTrough = troughIdx(1);
end

startIdx = max(anchorTrough - nPreOffset, 1);
endIdx   = min(startIdx + nShow - 1, nTotal);

fprintf('锚点谷点位置: %d，起始(谷前 %.0f%%): %d，结束: %d\n', ...
        anchorTrough, prePhase*100, startIdx, endIdx);
fprintf('截取共 %d 点 (%.3f s)\n', endIdx-startIdx+1, (endIdx-startIdx+1)/fs);

%% ── 截取各信号片段 ────────────────────────────────────
exSeg = @(s) s(startIdx:endIdx);

seg_PVEN = exSeg(sig_PVEN);   seg_PSAS = exSeg(sig_PSAS);
seg_PECS = exSeg(sig_PECS);   seg_PPAS = exSeg(sig_PPAS);
seg_PSPI = exSeg(sig_PSPI);
seg_QSAS = exSeg(sig_QSAS);   seg_QSPI = exSeg(sig_QSPI);
seg_QPAS = exSeg(sig_QPAS);   seg_QGin = exSeg(sig_QGin);
seg_QGot = exSeg(sig_QGot);   seg_QBBB = exSeg(sig_QBBB);
seg_VT   = exSeg(sig_VT);
seg_VSAS = exSeg(sig_VSAS);
seg_VSPI = exSeg(sig_VSPI);
seg_QOLS = exSeg(sig_QOLS);   seg_QSSS = exSeg(sig_QSSS);
seg_QSLS = exSeg(sig_QSLS);

% 体积变化量：相对所显示窗口的周期均值
seg_VSPI = seg_VSPI - mean(seg_VSPI);
seg_VSAS = seg_VSAS - mean(seg_VSAS);
seg_VT   = seg_VT   - mean(seg_VT);

nPts = numel(seg_PVEN);

%% ── x 轴：0 ~ nCycles*100% ───────────────────────────
tPct = linspace(0, nCycles * 100, nPts);
troughXpct = prePhase * 100;
troughXall = troughXpct + (0:nCycles-1)*100;

%% ── 颜色：按"生理实体"统一，同一符号在六个子图里同色 ────
% 实体 -> 颜色：
%   红  0.02  动脉侧    P_A(中心) / 桡动脉(虚线) / P_PAS / Q_PAS
%   蓝  0.58  静脉侧    P_VEN / Q_ECS (ECS->PNS) / (f) 的 SSS 静脉窦
%   绿  0.33  颅侧 SAS  P_SAS / Q_SAS / V_SAS / (f) 的 OLS 淋巴通路
%   紫  0.78  脊髓侧    P_SPI / Q_SPI / V_SPI / (f) 的 SLS 脊髓通路
%   亮橙 0.08 类淋巴流入 Q_GS   （比动脉红明显更亮更黄，同框也不混）
%   深青 0.45 BBB       Q_BBB
%   黑  0.10  ECS/ICP   P_ECS 与 V_T —— 二者共用同一条指数 P-V 关系，故同色
pathSaturation = 0.4;      % 饼图
hue_lymph   = 0.33;        % 绿
hue_ag      = 0.58;        % 蓝
hue_spinal  = 0.78;        % 紫
hue_red     = 0.02;        % 红
hue_gs      = 0.08;        % 亮橙 —— 类淋巴流入
hue_bbb     = 0.45;        % 青 —— BBB

clr_Lymph  = hsv2rgb([hue_lymph,  pathSaturation, 0.75]);
clr_SSS    = hsv2rgb([hue_ag,     pathSaturation, 0.75]);
clr_Spinal = hsv2rgb([hue_spinal, pathSaturation, 0.75]);

timeSeriesSat = 0.7;       % 时序图：提高饱和度、降低明度
timeSeriesVal = 0.7;
ts = @(h) hsv2rgb([h, timeSeriesSat, timeSeriesVal]);

clr_PA    = [0.85 0.20 0.15];   % 动脉压（中心实线 / 桡动脉虚线，同色同粗细）
clr_PArad = clr_PA;             % 桡动脉：与中心动脉压完全同色，靠虚线区分
clr_mean  = [0.25 0.25 0.25];   % 灰色（参考线）

clr_PVEN  = ts(hue_ag);         % 静脉侧 深蓝
clr_PSAS  = ts(hue_lymph);      % SAS   深绿
clr_ECS   = [0.10 0.10 0.10];   % ECS/ICP 黑（中性色，避免与橙/橄榄/红混淆）
clr_PECS  = clr_ECS;            % ICP = P_ECS
clr_PPAS  = ts(hue_red);        % PAS   深红（与动脉侧一致）
clr_PSPI  = ts(hue_spinal);     % 脊髓   深紫（与 (f) 的 Q_SLS 一致）

clr_QSAS  = ts(hue_lymph);      % Q_SAS 绿（与 P_SAS、V_SAS 一致）
clr_QSPI  = ts(hue_spinal);     % Q_SPI 紫（与 P_SPI、V_SPI、(f) Q_SLS 一致）
clr_QPAS  = ts(hue_red);        % Q_PAS 红（与 P_PAS 一致）
clr_QGin  = hsv2rgb([hue_gs, 0.90, 0.95]);    % Q_GS 亮橙（比动脉红更亮更黄）
clr_QGot  = ts(hue_ag);         % Q_ECS 蓝（流入静脉侧 PNS，取静脉色）
clr_QBBB  = hsv2rgb([hue_bbb, 0.85, 0.60]);   % Q_BBB 深青（与蓝拉开明度）

clr_VSPI  = ts(hue_spinal);     % 深紫（同 Q_SPI、P_SPI）
clr_VSAS  = ts(hue_lymph);      % 深绿（同 Q_SAS、P_SAS）
clr_VT    = clr_ECS;            % V_T：与 ICP 共用同一条指数 P-V 公式，故同色

lw = 1.5;

%% ── 创建图形：2×3 ────────────────────────────────────
figTS = figure('Name','Fig_Baseline_Waveforms','NumberTitle','off', ...
    'Units','centimeters','Position',[5 2 40 24], ...
    'Color','white','PaperPositionMode','auto');
set(figTS, 'DefaultAxesFontSize', fs_axis, 'DefaultTextFontSize', fs_axis);
tl = tiledlayout(figTS, 2, 3, 'TileSpacing','loose','Padding','compact');

xLimAll = [0, nCycles * 100];
xtv    = 0:50:nCycles*100;    % 刻度间隔 50%（原 25% 太密，放大字号后取消）
xtlbl  = arrayfun(@(v) sprintf('%d%%',v), xtv, 'UniformOutput', false);

%% ── (a) 桡动脉与中心动脉压 ───────────────────────────
ax_a = nexttile(tl);   hold(ax_a,'on');
abpFile = fullfile(root,'SimulinkModel','ABP_Data','120.txt');
if exist(abpFile,'file')
    [radial, tFit] = ABP_Get(abpFile, false);
    central = ABP_to_Central_GTF(radial, tFit, false);
    xr = linspace(0, 100, numel(radial) + 1);   xr(end) = [];
    % 桡动脉与中心动脉压：同色同粗细，桡动脉用虚线区分
    plot(ax_a, [xr, xr+100], [radial(:); radial(:)],   '--', 'Color', clr_PA, 'LineWidth', lw);
    plot(ax_a, [xr, xr+100], [central(:); central(:)], '-',  'Color', clr_PA, 'LineWidth', lw);
    vals_a = [radial(:); central(:)];
else
    plot(ax_a, tPct, exSeg(sig_PA), '-', 'Color', clr_PA, 'LineWidth', lw);
    vals_a = exSeg(sig_PA);
end
localYlim(ax_a, vals_a, PAD_TOP_A);
if SHOW_CYCLE_LINES, addCycleLines(ax_a, nCycles, troughXall); end
set(ax_a, 'XLim',xLimAll, 'XTick',xtv, 'XTickLabel',xtlbl, ...
    'YGrid','off', 'Box','off', 'GridAlpha',0.25);
ylabel(ax_a,'Pressure (mmHg)');
xlabel(ax_a,'Cardiac Cycle (%)');
title(ax_a,'(a)  Arterial Pressure','FontWeight','bold','FontSize',fs_t);
legend(ax_a,{'Radial (measured)','Central (input)'}, ...
    'Location','north','Box','off','FontSize',fs_l, ...
    'Orientation','horizontal','NumColumns',2);
localShiftLegend(ax_a.Legend, LEG_SHIFT_A * localCharW(figTS, fs_l), 0);

%% ── (b) 五个腔室压力 ─────────────────────────────────
ax_b = nexttile(tl);   hold(ax_b,'on');
Yb = [seg_PVEN, seg_PSAS, seg_PECS, seg_PPAS, seg_PSPI];
labB = {'P_{VEN}','P_{SAS}','P_{ECS}','P_{PAS}','P_{SPI}'};
split_b = localSplitFor('b', Yb, SPLIT_ON);
[hb1,hb2] = drawSeries(ax_b, tPct, Yb, {clr_PVEN,clr_PSAS,clr_PECS,clr_PPAS,clr_PSPI}, split_b, lw);
localYlimSplit(ax_b, Yb, split_b, PAD_TOP_B);
if SHOW_CYCLE_LINES, addCycleLines(ax_b, nCycles, troughXall); end
set(ax_b, 'XLim',xLimAll, 'XTick',xtv, 'XTickLabel',xtlbl, ...
    'YGrid','off', 'Box','off', 'GridAlpha',0.25);
localYlabel(ax_b, split_b, 'Pressure (mmHg)', 'Pressure (mmHg)');
xlabel(ax_b,'Cardiac Cycle (%)');
title(ax_b,'(b)  CSF/ISF Compartment Pressures','FontWeight','bold','FontSize',fs_t);
hb = [hb1(:); hb2(:)]';
legend(hb, labB(1:numel(hb)), ...
    'Location','northeast','Box','off','FontSize',fs_l,'NumColumns',3);
localShiftLegend(ax_b.Legend, LEG_SHIFT_B * localCharW(figTS, fs_l), 0);
localReport('b', labB, Yb, split_b);

%% ── (c) Q_SAS 与 Q_SPI ───────────────────────────────
ax_c = nexttile(tl);   hold(ax_c,'on');
Yc = [seg_QSAS, seg_QSPI];
labC = {'Q_{SAS}','Q_{SPI}'};
split_c = localSplitFor('c', Yc, SPLIT_ON);
[hc1,hc2] = drawSeries(ax_c, tPct, Yc, {clr_QSAS, clr_QSPI}, split_c, lw);
localYlimSplit(ax_c, Yc, split_c, PAD_TOP_C);
yline(ax_c, 0, '-', 'Color',[0.4 0.4 0.4],'LineWidth',0.7,'HandleVisibility','off');
if SHOW_CYCLE_LINES, addCycleLines(ax_c, nCycles, troughXall); end
set(ax_c, 'XLim',xLimAll, 'XTick',xtv, 'XTickLabel',xtlbl, ...
    'YGrid','off', 'Box','off', 'GridAlpha',0.25);
localYlabel(ax_c, split_c, 'Flow Rate (mL/min)', 'Flow Rate (mL/min)');
xlabel(ax_c,'Cardiac Cycle (%)');
title(ax_c,'(c)  Cranial SAS and Spinal Flow','FontWeight','bold','FontSize',fs_t);
hc = [hc1(:); hc2(:)]';
legend(hc, labC(1:numel(hc)), ...
    'Location','north','Box','off','FontSize',fs_l,'NumColumns',2);
localReport('c', {'Q_{SAS}','Q_{SPI}'}, Yc, split_c);

%% ── (d) 体积变化量：ΔV_SPI, ΔV_SAS, ΔV_T ─────────────
ax_d = nexttile(tl);   hold(ax_d,'on');
Yd = [seg_VSPI, seg_VSAS, seg_VT];
labD = {'\Delta V_{SPI}','\Delta V_{SAS}','\Delta V_{T}'};
split_d = localSplitFor('d', Yd, SPLIT_ON);
[hd1,hd2] = drawSeries(ax_d, tPct, Yd, {clr_VSPI, clr_VSAS, clr_VT}, split_d, lw);
localYlimSplit(ax_d, Yd, split_d);
yline(ax_d, 0, '-', 'Color',[0.4 0.4 0.4],'LineWidth',0.7,'HandleVisibility','off');
if SHOW_CYCLE_LINES, addCycleLines(ax_d, nCycles, troughXall); end
set(ax_d, 'XLim',xLimAll, 'XTick',xtv, 'XTickLabel',xtlbl, ...
    'YGrid','off', 'Box','off', 'GridAlpha',0.25);
localYlabel(ax_d, split_d, '\Delta V (mL)', '\Delta V (mL)');
xlabel(ax_d,'Cardiac Cycle (%)');
title(ax_d,'(d)  Compartment Volume Change','FontWeight','bold','FontSize',fs_t);
hd = [hd1(:); hd2(:)]';
legend(hd, labD(1:numel(hd)), ...
    'Location','northeast','Box','off','FontSize',fs_l,'NumColumns',3);
localReport('d', labD, Yd, split_d);

%% ── (e) 血管旁 / 类淋巴流量 ──────────────────────────
ax_e = nexttile(tl);   hold(ax_e,'on');
Ye = [seg_QPAS, seg_QGin, seg_QGot, seg_QBBB];
% 显示名与论文统一：Q_GS -> Q_{GS}，Q_ECS -> Q_{ECS}
% （模型信号名仍是 Q_GS / Q_ECS，见上面的 getSignal 调用，不要改）
labE = {'Q_{PAS}','Q_{GS}','Q_{ECS}','Q_{BBB}'};
split_e = localSplitFor('e', Ye, SPLIT_ON);
[he1,he2] = drawSeries(ax_e, tPct, Ye, {clr_QPAS, clr_QGin, clr_QGot, clr_QBBB}, split_e, lw);
localYlimSplit(ax_e, Ye, split_e);
yline(ax_e, 0, '-', 'Color',[0.4 0.4 0.4],'LineWidth',0.7,'HandleVisibility','off');
if SHOW_CYCLE_LINES, addCycleLines(ax_e, nCycles, troughXall); end
set(ax_e, 'XLim',xLimAll, 'XTick',xtv, 'XTickLabel',xtlbl, ...
    'YGrid','off', 'Box','off', 'GridAlpha',0.25);
localYlabel(ax_e, split_e, 'Flow Rate (mL/min)', 'Flow Rate (mL/min)');
xlabel(ax_e,'Cardiac Cycle (%)');
title(ax_e,'(e)  Glymphatic flow','FontWeight','bold','FontSize',fs_t);
he = [he1(:); he2(:)]';
legend(he, labE(1:numel(he)), ...
    'Location','northeast','Box','off','FontSize',fs_l, ...
    'Orientation','horizontal','NumColumns',4);
localShiftLegend(ax_e.Legend, LEG_SHIFT_D, LEG_DROP_D * localCharH(figTS, fs_l));
localReport('e', labE, Ye, split_e);

%% ── (f) 流出路径饼图 ─────────────────────────────────
ax_f = nexttile(tl);

% 窗口内均值（两个完整周期，即周期均值）
avgLymph = mean(seg_QOLS);
avgSSS   = mean(seg_QSSS);
avgSPI   = mean(seg_QSLS);

pieData   = abs([avgLymph, avgSSS, avgSPI]);
pieLabels = {'Q_{OLS}','Q_{SSS}','Q_{SLS}'};
pieColors = [clr_Lymph; clr_SSS; clr_Spinal];

if any(pieData <= 0) || any(isnan(pieData)) || any(isinf(pieData))
    warning('饼图数据存在无效值');
    pieData = [0.1, 0.2, 0.3];
end

p = pie(ax_f, pieData*10);
for i = 1:2:length(p)
    idx = (i+1)/2;
    set(p(i), 'FaceColor', pieColors(idx,:), 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.8);
end
total    = sum(pieData);
percents = pieData / total * 100;
for i = 2:2:length(p)
    idx = i/2;
    p(i).String = sprintf('%s\n%.1f%%\n(%.3f mL/min)', ...
        pieLabels{idx}, percents(idx), pieData(idx));
    p(i).FontSize   = fs_p;
    p(i).FontWeight = 'bold';
    p(i).Color      = [1 1 1];        % 扇区内文字用白色（压在彩色扇区上更清晰）
    % 标注移到饼图内部：沿半径向圆心缩放（饼图圆心在原点）
    pos = p(i).Position;
    p(i).Position = [pos(1)*PIE_INSIDE, pos(2)*PIE_INSIDE, pos(3)];
    p(i).HorizontalAlignment = 'center';
    p(i).VerticalAlignment   = 'middle';
end
title(ax_f,'(f)  Total Efflux Pathways Distribution','FontWeight','bold','FontSize',fs_t);

% 饼图内文字微调：Q_{OLS} 右移、Q_{SSS} 上移（单位=字符尺寸，换算成饼图坐标单位）
if numel(p) >= 4
    axPos = ax_f.Position;   figCm = figTS.Position(3:4);
    cmX = axPos(3)*figCm(1) / diff(xlim(ax_f));   % 1 个数据单位 = 多少 cm
    cmY = axPos(4)*figCm(2) / diff(ylim(ax_f));
    charHcm = fs_p * 2.54/72;   charWcm = 0.5 * charHcm;
    pos = p(2).Position;
    p(2).Position = [pos(1) + PIE_NUDGE_OLS_X*charWcm/cmX, pos(2), pos(3)];
    pos = p(4).Position;
    p(4).Position = [pos(1), pos(2) + PIE_NUDGE_SSS_Y*charHcm/cmY, pos(3)];
end

%% ── 刻度标签与轴标题字号（横/纵统一），横轴不旋转 ────
for axk = [ax_a ax_b ax_c ax_e ax_d ax_f]
    axk.XAxis.FontSize  = fs_xtick;
    axk.XLabel.FontSize = fs_axlab;
    xtickangle(axk, 0);
    for kk = 1:numel(axk.YAxis)
        axk.YAxis(kk).FontSize       = fs_xtick;
        axk.YAxis(kk).Label.FontSize = fs_axlab;
    end
end

%% ── 字号自检（读回实际生效的字号）────────────────────
if isempty(ax_b.Legend), fs_leg = NaN; else, fs_leg = ax_b.Legend.FontSize; end
fprintf(['  字号(pt): 刻度 横%.1f/纵%.1f | 轴标题 横%.1f/纵%.1f | ' ...
         '面板标题 %.1f | 图例 %.1f | 饼图标签 %.1f\n'], ...
    ax_b.XAxis.FontSize, ax_b.YAxis(1).FontSize, ...
    ax_b.XLabel.FontSize, ax_b.YAxis(1).Label.FontSize, ...
    ax_b.Title.FontSize, fs_leg, p(2).FontSize);

%% ── 导出（仅 PDF）────────────────────────────────────
% MATLAB 把 FOP 字体缓存写到 %USERPROFILE%\.fop；本会话该目录不可写，故临时把
% user.home 指向平台临时目录。
try
    java.lang.System.setProperty('user.home', getenv('TEMP'));
catch
end
exportgraphics(figTS, outPdf, 'ContentType','vector','Resolution',300);
fprintf('\n✓ 图已生成: %s\n', outPdf);

%% ── 把本脚本复制一份到插图目录，使该目录自带可复现的绘图代码 ──
try
    thisFile = [mfilename('fullpath') '.m'];
    if ~strcmpi(fileparts(thisFile), fileparts(outPdf))
        copyfile(thisFile, fullfile(fileparts(outPdf), 'draw_baseline_timeseries.m'));
        fprintf('  脚本副本: %s\n', fullfile(fileparts(outPdf),'draw_baseline_timeseries.m'));
    end
catch err
    fprintf('  脚本复制失败: %s\n', err.message);
end

%% ── 控制台汇总 ────────────────────────────────────────
fprintf('  谷点前偏移: %.0f%% 周期 (%d 点)\n', prePhase*100, nPreOffset);
fprintf('  流出路径分布（窗口均值）:\n');
fprintf('    Lymphatic (Q_OLS): %.4f mL/min (%.2f%%)\n', avgLymph, percents(1));
fprintf('    Sinus     (Q_SSS): %.4f mL/min (%.2f%%)\n', avgSSS,   percents(2));
fprintf('    Spinal    (Q_SLS): %.4f mL/min (%.2f%%)\n', avgSPI,   percents(3));
fprintf('  体积搏动幅度 (mL): SPI %.4f | SAS %.4f | T %.4f\n', ...
    max(seg_VSPI)-min(seg_VSPI), max(seg_VSAS)-min(seg_VSAS), ...
    max(seg_VT)-min(seg_VT));
fprintf('=== DONE ===\n');
end

%% ════════════════════════════════════════════════════════
%  辅助函数
%% ════════════════════════════════════════════════════════

function vec = getSignal(ds, name)
    ts  = getElem(ds, name).Values;
    vec = double(ts.Data(:));
end

function tvec = getTime(ds, name)
    ts   = getElem(ds, name).Values;
    tvec = double(ts.Time(:));
end

function e = getElem(ds, name)
    % logsout 是 Simulink.SimulationData.Dataset；参考代码里的 baselineDataset
    % 用的是 .get，这里两种都兼容
    try
        e = ds.get(name);
    catch
        e = ds.getElement(name);
    end
end

function localShiftLegend(lg, dx, dy)
    % 平移图例（figure 归一化单位）：dx 右为正，dy 上为正
    if isempty(lg), return; end
    if nargin < 3 || isempty(dy), dy = 0; end
    lg.Units = 'normalized';
    pos = lg.Position;
    lg.Position = [pos(1) + dx, pos(2) + dy, pos(3), pos(4)];
end

function h = localCharH(fig, fs)
    % 一个字符的高度换算成 figure 高度的比例
    h = (fs * 2.54/72) / fig.Position(4);
end

function w = localCharW(fig, fs)
    % 一个字符的宽度换算成 figure 宽度的比例（Arial 平均字宽 ≈ 0.5 字高）
    w = 0.5 * (fs * 2.54/72) / fig.Position(3);
end

function addCycleLines(ax, nCyc, troughXall)
    % 周期分界线（默认关闭，SHOW_CYCLE_LINES = true 可开启）
    yl = ylim(ax);
    for k = 1:numel(troughXall)
        xv = troughXall(k);
        plot(ax, [xv xv], yl, '--', 'Color',[0.6 0.6 0.6], ...
            'LineWidth',0.8, 'HandleVisibility','off');
    end
end

function localYlim(ax, vals, padTop)
    if nargin < 3 || isempty(padTop), padTop = 0.20; end
    vals = vals(~isnan(vals));
    yspan = max(vals) - min(vals);
    if yspan < 1e-9, yspan = max(abs(max(vals)), 1) * 0.1; end
    ylim(ax, [min(vals) - yspan*0.15, max(vals) + yspan*padTop]);
end

function localYlimSplit(ax, Y, split, padTop)
    if nargin < 4, padTop = []; end
    if isempty(ax.YAxis) || numel(ax.YAxis) < 2 || ~strcmp(ax.YAxis(2).Visible,'on')
        localYlim(ax, Y(:), padTop);
    else
        yyaxis(ax,'left');  localYlim(ax, Y(:,split), padTop);
        yyaxis(ax,'right'); localYlim(ax, Y(:,~split), padTop);
        yyaxis(ax,'left');
    end
end

function [hL, hR] = drawSeries(ax, x, Y, colours, split, lw)
% 画曲线并返回句柄（hL = 左轴/单轴系列，hR = 右轴系列），供 legend 显式引用
hL = gobjects(1,0);  hR = gobjects(1,0);
    if all(split)
        for k = 1:size(Y,2)
            hL(end+1) = plot(ax, x, Y(:,k), '-', 'Color', colours{k}, 'LineWidth', lw); %#ok<AGROW>
        end
    else
        yyaxis(ax,'left');
        for k = find(split),  hL(end+1) = plot(ax, x, Y(:,k), '-', 'Color', colours{k}, 'LineWidth', lw); end %#ok<AGROW>
        yyaxis(ax,'right');
        for k = find(~split), hR(end+1) = plot(ax, x, Y(:,k), '-', 'Color', colours{k}, 'LineWidth', lw); end %#ok<AGROW>
        kL = find(split,1);  kR = find(~split,1);
        ax.YAxis(1).Color = colours{kL};   ax.YAxis(2).Color = colours{kR};
        ax.YAxis(1).Label.Color = colours{kL};  ax.YAxis(2).Label.Color = colours{kR};
        yyaxis(ax,'left');
    end
end

function localYlabel(ax, split, txtL, txtR)
    if all(split)
        ylabel(ax, txtL);
    else
        ylabel(ax.YAxis(1), txtL);
        ylabel(ax.YAxis(2), txtR);
    end
end

function split = localSplitFor(tag, Y, splitOn)
    % 单 y 轴，除非该面板字母列在 SPLIT_ON 中
    if any(strcmp(splitOn, tag))
        rng = max(Y,[],1,'omitnan') - min(Y,[],1,'omitnan');
        rng(isnan(rng)) = 0;
        split = rng * 8 >= max(rng);
        if all(split) || ~any(split), split = true(1,size(Y,2)); end
    else
        split = true(1, size(Y,2));
    end
end

function localReport(tag, names, Y, split)
    % 每条曲线占坐标轴高度的百分比——占比越小则越接近直线
    fprintf('  面板 (%s) 各曲线占轴高:\n', tag);
    for side = 1:2
        if side == 1, m = split; else, m = ~split; end
        idx = find(m);  if isempty(idx), continue; end
        if all(split)
            v = Y(:,idx);  v = v(~any(isnan(v),2), :);
            lo = min(v(:)); hi = max(v(:));
            span = hi - lo;  if span < 1e-9, span = 1; end
            h = span*1.35;
        else
            v = Y(:,idx);  v = v(~any(isnan(v),2), :);
            lo = min(v(:)); hi = max(v(:));
            span = hi - lo;  if span < 1e-9, span = 1; end
            h = span*1.35;
        end
        for k = idx
            if all(isnan(Y(:,k))), continue; end
            sw = max(Y(:,k)) - min(Y(:,k));
            fprintf('    %-10s %+9.4f .. %+9.4f  swing %8.4f  (%4.1f%%)\n', ...
                names{k}, min(Y(:,k)), max(Y(:,k)), sw, 100*sw/h);
        end
    end
end

function f = localFindBaseline(root)
    cand = { fullfile(root,'SimulinkModel','Results','H_series','B.mat'), ...
             fullfile(root,'SimulinkModel','Results','H_series','B.mat'), ...
             fullfile(root,'SimulinkModel','Results','H_series','B.mat'), ...
             fullfile(root,'SimulinkModel','Results','B_series','B.mat') };
    for k = 1:numel(cand)
        if exist(cand{k},'file'), f = cand{k}; return; end
    end
    d = dir(fullfile(root,'SimulinkModel','Results','**','*.mat'));
    for k = 1:numel(d)
        if ~isempty(regexpi(d(k).name,'^(b|base)[_\-\.]','once'))
            f = fullfile(d(k).folder, d(k).name); return
        end
    end
    error('draw_baseline_timeseries:noData','Results/ 下找不到基线 .mat');
end
