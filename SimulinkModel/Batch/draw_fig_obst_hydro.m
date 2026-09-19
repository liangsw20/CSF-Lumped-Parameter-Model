%% ========================================================
%  FIGURE 6: Obstructive hydrocephalus (ventricular mechanics) together with
%  the two combined scenarios that contain a ventricular setting.
%
%  Six scenarios on one common x-axis, grouped into two colour families
%  (the colour coding alone marks the grouping -- no separator rule is drawn):
%    resistance family (warm)   R2        = Ven_R2
%                               R7.5      = Ven_R7.5
%                               H+R7.5    = H_G2 + Ven_R7.5
%    elastance family (cool)    E0.5      = Ven_E0.5
%                               E0.25     = Ven_E0.25
%                               H+C+E     = H_G2 + C_OVS2.5 + Ven_E0.25
%  The baseline B is not drawn as a column: it is the dashed grey reference line
%  of the pressure panel and the denominator of every percentage.
%
%  Two panels:
%    (a) ventricular pressure P_VEN                    -- the driving pressure
%    (b) volume change of the ventricle and the ECS    -- where the volume goes
%
%  Style follows MATLAB_script/DrawResult0520/ResultDraw.m: same %-annotation and
%  reference line, same 1.5x font scale as Fig3/Fig4/Fig5.
%
%  Data: SimulinkModel/Results/collate_all.csv (scenario means, from
%        Batch/collate_all.m).  Output: PDF only into <root>/新插图/, plus a
%        copy of this script.
%
%  Usage: matlab -batch "draw_fig_obst_hydro"
% ========================================================
function draw_fig_obst_hydro(csvFile, outPdf)
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(here);

if nargin < 1 || isempty(csvFile)
    csvFile = fullfile(root,'SimulinkModel','Results','collate_all.csv');
end
if nargin < 2 || isempty(outPdf)
    outPdf = fullfile(root,'新插图','Fig6_ObstHydro.pdf');
end
if ~exist(fileparts(outPdf),'dir'), mkdir(fileparts(outPdf)); end

%% ── 场景定义（顺序：R 族三根，再 E 族三根）────────────
CODES  = {'Ven_R2','Ven_R7.5','H_G2_Ven_R7.5', ...        % 阻力族（暖色系）
          'Ven_E0.5','Ven_E0.25','H_G2_C_OVS2.5_Ven_E0.25'}; % 弹性族（冷色系）
SHORT  = {'R2','R7.5','H+R','E0.5','E0.25','H+C+E'};   % x 轴显示名（完整代号见图注）
BASE   = 'B';
NS     = numel(CODES);
% 两族只靠配色区分，不画竖线

% 配色：族内同一色系，族间冷暖分明；同族内随严重度/叠加加重
PAL = [ 0.96 0.80 0.50;    % R2        浅金
        0.88 0.60 0.22;    % R7.5      金
        0.72 0.33 0.20;    % H+R7.5    砖红（同暖色系）
        0.70 0.84 0.94;    % E0.5      浅蓝
        0.42 0.66 0.87;    % E0.25     中蓝
        0.20 0.40 0.68];   % H+C+E     深蓝（同冷色系）

C.white = [1 1 1];
C.ecsBar = [0.45 0.45 0.45];           % ECS/总量族：ΔV_ECS 统一灰
C.bar.edge = [0.35 0.35 0.35];  C.bar.edgeWidth = 0.6;
C.line.ref = [0.40 0.40 0.40];
C.label.percent = [0.10 0.10 0.10];

FONT_SCALE = 1.5;
% 横/纵轴刻度字号按需求放大到约 1.2 倍（原为 12 / 11）
fs_tick = 14.4; fs_xtick = 13.2; fs_ax = 15; fs_title = 9*FONT_SCALE;
fs_pct = 7.5*FONT_SCALE;

fprintf('=== draw_fig_obst_hydro ===\n  数据: %s\n  输出: %s\n', csvFile, outPdf);

%% ── 读表 ──────────────────────────────────────────────
T = readtable(csvFile, 'VariableNamingRule','preserve');
codesAll = string(T.code);
iB = find(codesAll == BASE, 1);
if isempty(iB), error('draw_fig_obst_hydro:noBase','表里找不到基线行 %s', BASE); end
idx = zeros(1,NS+1);  idx(1) = iB;              % 第 1 行始终是基线，仅作参考不画柱
for k = 1:NS
    i = find(codesAll == string(CODES{k}), 1);
    if isempty(i), error('draw_fig_obst_hydro:noRow','表里缺少场景 %s', CODES{k}); end
    idx(k+1) = i;
end
fprintf('  基线 %s + %d 个场景全部找到\n', BASE, NS);

g = @(col) local_col(T, col, idx);              % 1 x (NS+1)，含基线
qGSin = g('Q_GS');  qBBB = g('Q_BBB');  pECS = g('P_ECS');
pVEN  = g('P_VEN');   vVEN = g('V_VEN');  vECS = g('V_ECS');
qOLS  = g('Q_OLS');   qSSS = g('Q_SSS');  qSLS = g('Q_SLS');
tot   = qOLS + qSSS + qSLS;
dVVEN = vVEN - vVEN(1);   dVECS = vECS - vECS(1);
kS    = 2:NS+1;                                 % 要画柱的行（不含基线）
xpos  = 1:NS;
fprintf('  B: P_VEN %.4f | V_VEN %.4f | V_ECS %.3f | 总出口 %.5f\n', ...
    pVEN(1), vVEN(1), vECS(1), tot(1));

%% ── 图形：1 x 2 ───────────────────────────────────────
set(groot,'DefaultAxesFontName','Arial'); set(groot,'DefaultTextFontName','Arial');
fig = figure('Name','Fig6_ObstHydro','NumberTitle','off', ...
    'Units','centimeters','Position',[3 2 26 12.5], ...
    'Color',C.white,'PaperPositionMode','auto');
set(fig,'DefaultAxesFontSize',fs_tick,'DefaultTextFontSize',fs_tick);
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact','Padding','compact');

%% ── (a) 心室压 P_VEN ──────────────────────────────────
% 压强远离 0：取整 10 mmHg 作交点（审稿意见允许，须在图注写明）
ax = nexttile;
local_barPanel(ax, xpos, pVEN(kS), pVEN(1), PAL, SHORT, C, fs_xtick, fs_pct, 10);
ylabel(ax,'P_{VEN} (mmHg)');
title(ax,'(a)  Ventricular Pressure (P_{VEN})','FontWeight','bold','FontSize',fs_title);
ax_a = ax;

%% ── (b) 心室与 ECS 的容积变化（ΔV_VEN / ΔV_ECS）───────
% 两个量同单位（mL）、同量级（约 ±2.7 mL），共用一个从 0 起的对称轴；
% 每组左柱 = ΔV_VEN（沿用 (a) 的场景配色），右柱 = ΔV_ECS（灰，ECS/总量族）。
% 四项单独场景里两者严格互为镜像（1:1 体积互换），叠加高血压后不再守恒。
ax = nexttile;  hold(ax,'on');
axMag = ceil(2*max(abs([dVVEN(kS), dVECS(kS)]))*1.12)/2;   % 对称圆整幅度（此前为 3.5）
axLo  = -axMag;
axHi  = axMag + 0.5;                            % 上界再抬 0.5，给轴内顶部图例留位
barTop = max([dVVEN(kS), dVECS(kS)]);           % 竖线只画到柱顶，不穿图例
bvol = bar(ax, xpos, [dVVEN(kS)', dVECS(kS)'], 'grouped', ...
    'EdgeColor',C.bar.edge, 'LineWidth',C.bar.edgeWidth, 'BarWidth',0.8);
bvol(1).FaceColor = 'flat';  bvol(1).CData = PAL;   % 左柱按场景配色，与 (a) 一致
bvol(2).FaceColor = C.ecsBar;                       % 右柱统一灰
ylim(ax, [axLo, axHi]);
xlim(ax, [0.5 NS+0.5]);            % 与 (a) 相同取景；否则 0 参考线会把轴撑出额外留白
plot(ax, [0.5 NS+0.5], [0 0], '-', 'Color',[0.55 0.55 0.55], 'LineWidth',0.8, 'HandleVisibility','off');
set(ax,'XTick',xpos,'XTickLabel',SHORT,'XTickLabelRotation',0, ...
    'TickLabelInterpreter','none','YGrid','off','XGrid','off','Box','off');
ax.XAxis.FontSize = fs_xtick;
ylabel(ax,'\Delta V (mL)');   ax.YLabel.FontSize = fs_ax;
title(ax,'(b)  Volume Change of the Ventricle and the ECS','FontWeight','bold','FontSize',fs_title);
legend(ax, bvol, {'\Delta V_{VEN}  (left)','\Delta V_{ECS}  (right)'}, 'Location','north', ...
       'FontSize',fs_xtick,'Box','off','Orientation','horizontal','NumColumns',2);
ylim(ax, [axLo, axHi]);
ax_b = ax;

%% ── 版面自检 ─────────────────────────────────────────
drawnow;
% 先确认两格等宽：tiledlayout 会按各格装饰微调，若不等宽则刻度无法对齐
% （tiledlayout 的子轴禁止手写 Position/PositionConstraint，只能报出来人工处理）
axs = [ax_a ax_b];
pw  = reshape([axs.Position], 4, []).';
figW = fig.Position(3);  figH = fig.Position(4);
if max(pw(:,3)) - min(pw(:,3)) > 1e-6
    fprintf('  ⚠ 两格轴宽不一致 %s，需手工调整\n', mat2str(round(pw(:,3)',4)));
else
    fprintf('  两格轴宽一致 %.4f 图窗宽 (= %.2f cm)\n', pw(1,3), pw(1,3)*figW);
end
colw = pw(1,3)*figW/NS;
wl = zeros(1,NS);
for r = 1:NS
    h = text(0,0,SHORT{r},'FontName','Arial','FontSize',fs_xtick, ...
             'Units','centimeters','Visible','off','Interpreter','none');
    wl(r) = h.Extent(3);  delete(h);
end
gap  = colw - (wl(1:end-1)+wl(2:end))/2;      % 相邻代号
gapE = 0.5*colw - wl/2;                       % 首/末代号到面板边缘
[gmin, imin] = min(gap);
if min([gmin, gapE]) >= 0.10, st = 'OK'; else, st = '间距不足，请降低 fs_xtick!'; end
fprintf('  字号(pt): 刻度 %.1f | 横轴 %.1f | 轴标题 %.1f | 面板标题 %.1f（印到 17.4 cm 宽时 ~%.1f pt）\n', ...
    fs_tick, fs_xtick, fs_ax, fs_title, fs_tick*17.4/figW);
fprintf('  横轴代号: 最宽 %s (%.2f cm) | 列宽 %.2f cm | 最紧相邻 %s|%s 余量 %.2f cm | 边缘余量 %.2f/%.2f cm %s\n', ...
    SHORT{find(wl==max(wl),1)}, max(wl), colw, SHORT{imin}, SHORT{imin+1}, gmin, gapE(1), gapE(end), st);
fprintf('  (a) ylim [%.1f %.2f] xlim %s | (b) ylim [%.1f %.1f] xlim %s，柱顶 %+.3f（不画分组竖线）\n', ...
    ylim(ax_a), mat2str(xlim(ax_a)), axLo, axHi, mat2str(xlim(ax_b)), barTop);
fprintf('  (b) ΔV_VEN+ΔV_ECS: %s（单独场景≈0，叠加高血压不再守恒）\n', ...
    mat2str(round(dVVEN(kS)+dVECS(kS),3)));

%% ── 导出 + 复制脚本 ──────────────────────────────────
try, java.lang.System.setProperty('user.home', getenv('TEMP')); catch, end
exportgraphics(fig, outPdf, 'ContentType','vector','Resolution',300);
fprintf('✓ 图已生成: %s\n', outPdf);
try
    thisFile = [mfilename('fullpath') '.m'];
    if ~strcmpi(fileparts(thisFile), fileparts(outPdf))
        copyfile(thisFile, fullfile(fileparts(outPdf),'draw_fig_obst_hydro.m'));
        fprintf('  脚本副本: %s\n', fullfile(fileparts(outPdf),'draw_fig_obst_hydro.m'));
    end
catch err
    fprintf('  脚本复制失败: %s\n', err.message);
end

%% ── 控制台汇总（供写图注）────────────────────────────
fprintf('\n--- 相对基线 (%%) ---\n');
fprintf('  %-26s %8s %8s %8s %8s %8s %8s %8s\n', ...
    'code','P_VEN','P_ECS','Q_GS','Q_BBB','dV_VEN','dV_ECS','Q_out');
for k = kS
    fprintf('  %-26s %+8.1f %+8.1f %+8.1f %+8.1f %+8.1f %+8.1f %+8.1f\n', CODES{k-1}, ...
        100*(pVEN(k)/pVEN(1)-1), 100*(pECS(k)/pECS(1)-1), ...
        100*(qGSin(k)/qGSin(1)-1), 100*(qBBB(k)/qBBB(1)-1), ...
        100*(vVEN(k)/vVEN(1)-1), 100*(vECS(k)/vECS(1)-1), 100*(tot(k)/tot(1)-1));
end
fprintf('\n--- 体积 (mL) ---\n');
fprintf('  %-26s %10s %10s %10s %10s %12s\n','code','V_VEN','dV_VEN','V_ECS','dV_ECS','sum');
for k = [1 kS]
    nm = 'B';  if k > 1, nm = CODES{k-1}; end
    fprintf('  %-26s %10.4f %10.4f %10.3f %10.4f %12.4f\n', ...
        nm, vVEN(k), dVVEN(k), vECS(k), dVECS(k), dVVEN(k)+dVECS(k));
end
fprintf('\n--- 流出占比 (%%) ---\n');
for k = [1 kS]
    nm = 'B';  if k > 1, nm = CODES{k-1}; end
    fprintf('  %-26s OLS %5.2f | SSS %5.2f | SLS %5.2f  (总 %.5f mL/min)\n', ...
        nm, 100*qOLS(k)/tot(k), 100*qSSS(k)/tot(k), 100*qSLS(k)/tot(k), tot(k));
end
fprintf('=== DONE ===\n');
end

%% ════════════════════════════════════════════════════════
%  辅助函数
%% ════════════════════════════════════════════════════════
function v = local_col(T, name, idx)
if ~ismember(name, T.Properties.VariableNames)
    error('draw_fig_obst_hydro:noCol','表里缺少列 %s', name);
end
col = T.(name);  v = col(idx);
if iscell(v), v = cell2mat(v); end
v = double(v(:))';
end

function local_barPanel(ax, xpos, vals, baseVal, pal, labels, C, fs_xtick, fs_pct, yStart)
% 条形面板：基线参考虚线 + 柱顶 % 变化。基线不画柱，只以虚线 + 百分比分母的形式出现。
%   yStart = 0    从 0 起（有物理零点的量，审稿意见）
%   yStart = 10   压强类量：取整 10 mmHg 作交点（审稿意见允许，须在图注写明）
if nargin < 10 || isempty(yStart), yStart = 0; end
hold(ax,'on');
bar(ax, xpos, vals, 'FaceColor','flat', 'EdgeColor',C.bar.edge, ...
    'LineWidth',C.bar.edgeWidth, 'BarWidth',0.65, 'CData',pal);
sp = max(vals) - min(vals);  if sp < 1e-9, sp = max(abs(vals))*0.2 + 0.01; end
yTop = max(vals) + max(sp*0.40, (max(vals)-yStart)*0.16);
ylim(ax, [yStart, yTop]);
xlim(ax, [0.5 numel(vals)+0.5]);      % 与 (b) 相同的取景，避免自动留白不一致
yline(ax, baseVal, '--', 'Color',C.line.ref, 'LineWidth',1.0, 'HandleVisibility','off');
set(ax,'XTick',xpos,'XTickLabel',labels,'XTickLabelRotation',0, ...
    'TickLabelInterpreter','none','YGrid','off','XGrid','off','Box','off');
ax.XAxis.FontSize = fs_xtick;
dy = local_charData(ax, fs_pct);
for k = 1:numel(vals)
    d = 100*(vals(k)/baseVal-1);
    if abs(d) < 0.05, txt = '0%'; else, txt = sprintf('%+.0f%%', d); end
    text(ax, xpos(k), vals(k) + 0.5*dy, txt, 'HorizontalAlignment','center', ...
        'FontSize',fs_pct,'Color',C.label.percent,'FontWeight','bold');
end
end

function dy = local_charData(ax, fs)
% 一个字符的高度折算成该坐标轴的 data 单位
fig = ancestor(ax,'figure');  figCm = fig.Position(3:4);
axCm = ax.Position(4) * figCm(2);
dy = (fs*2.54/72) / axCm * diff(ylim(ax));
end
