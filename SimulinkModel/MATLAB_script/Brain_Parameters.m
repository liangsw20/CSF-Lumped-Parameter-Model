% 单位基准mmHg, mL, min
% 交互式运行：默认场景 'B'；批量驱动可先设定 SCENARIO_CODE 再调用本脚本。
% 全部计算在 bp_compute 函数内完成（单一数据源），本脚本只负责把结果
% 推送到 base 工作区，供 Simulink 解析模块参数。
if ~exist('SCENARIO_CODE','var') || isempty(SCENARIO_CODE)
    SCENARIO_CODE = 'B';
end
close all;

S  = bp_compute(SCENARIO_CODE);
fn = fieldnames(S);
for i = 1:numel(fn)
    assignin('base', fn{i}, S.(fn{i}));
end
fprintf('Brain_Parameters: 场景 %s 已载入，%d 个变量写入 base 工作区。\n', SCENARIO_CODE, numel(fn));
clear S fn i
