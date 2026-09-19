function [mean_cycle,t_fit] = ABP_Get(fileName, ifPlot)
    fs = 125;  
    fid = fopen(fileName, 'r', 'n', 'UTF-8');
    rawText = fread(fid, inf, '*char')';
    fclose(fid);
    data = jsondecode(rawText);
    
    sbp_target = data.SBP;   
    dbp_target = data.DBP;   
    raw_signal = data.signal;
    raw_min = min(raw_signal);
    raw_max = max(raw_signal);
    abp_mmHg = ((raw_signal - raw_min) / (raw_max - raw_min)) * (sbp_target - dbp_target) + dbp_target;
    
    % 4. 识别波峰
    [pks, locs] = findpeaks(abp_mmHg, 'MinPeakHeight', dbp_target + (sbp_target-dbp_target)*0.6, ...
                             'MinPeakDistance', fs*0.5);
    
    % 5. 周期提取、重采样并调整起始点
    num_cycles = length(locs) - 1; 
    T_period   = 0.8;    % 【设定】一个心动周期 0.8 s（心率 75 bpm）
    dt_fit     = 0.001;  % 【设定】输出波形的时间步长 (s)
    % 波形表【不得包含】周期终点 T_period。
    % Repeating Sequence Interpolated 的有效周期 =（跨度内 SampleTime 网格点数）
    %   x SampleTime。模型输入块 SampleTime = 0.01 s：
    %     跨度含 0.80 端点 -> 81 个网格点 -> 周期 0.81000000 s（实测）
    %       与 0.8 s 折叠会产生 64.8 s 摩尔纹，是【数值假振荡，非物理】
    %     跨度 [0, 0.799]  -> 80 个网格点 -> 周期 0.80000000 s（实测）
    target_len = round(T_period/dt_fit);   % = 800 点，跨度 [0, 0.799]
    cycle_matrix = zeros(num_cycles, target_len);
    
    for i = 1:num_cycles
        % 提取从峰值到峰值的原始段
        segment = abp_mmHg(locs(i) : locs(i+1)-1);
        
        % 线性插值统一长度
        resampled_seg = interp1(1:length(segment), segment, ...
                                 linspace(1, length(segment), target_len), 'linear');
        
        % 【关键修改】：寻找该段内的最小值位置（即波谷）
        [~, min_idx] = min(resampled_seg);
        
        % 循环移位：将最小值点移动到数组的第一个位置
        % 这样波形顺序就变成了：舒张压 -> 收缩压 -> 舒张压
        cycle_matrix(i, :) = circshift(resampled_seg, -(min_idx-1));
    end
    
    % 6. 计算平均值
    mean_cycle = mean(cycle_matrix, 1);
    std_cycle = std(cycle_matrix, 0, 1);

    % --- 6b.【归一化回标签血压】平均步骤本身不作任何改动，只把平均后的波形
    % 重新映射回本记录标签的 SBP/DBP。
    %
    % 为什么必须补这一步：上面的 min/max 定标是在【整段记录】上取的全局极值，
    % 单搏确实能到 120/80，但各搏的幅度与波形形态并不相同。averaged cycle
    % 于是只剩 116.20/83.02（实测），脉压 33.18 而非标签的 40，比标签少 17%；
    % 再经 ABP_to_Central_GTF（桡动脉→中心主动脉，正确地把脉压乘 ~0.81）后，
    % 模型实际看到的动脉脉压只有 26.8 mmHg，而生理上 120/80 对应的中心动脉/
    % 颈动脉脉压应为 30-35 mmHg。
    %
    % 归一化的代价可以忽略：平均压仅从 95.64 变到 95.22 mmHg（-0.4%），
    % 稳态标定不受影响；收益是把中心动脉脉压从 26.8 抬到 32.2 mmHg，
    % 血管每搏容积摆幅从 0.67 升到约 0.81 mL（后者已由 S_WAVE_AMP080 /
    % AMP120 两个现成实验证明对输入脉压线性：x0.80 -> x0.800, x1.20 -> x1.199）。
    %
    % 注意：若某记录波形平坦（极差为 0），跳过以免除零。
    rng_cycle = max(mean_cycle) - min(mean_cycle);
    if rng_cycle > 0
        k_renorm   = (sbp_target - dbp_target) / rng_cycle;
        mean_cycle = (mean_cycle - min(mean_cycle)) * k_renorm + dbp_target;
        std_cycle  = std_cycle * k_renorm;   % 同步缩放，保持 ±std 图一致
    end
    
    % --- 9. 傅里叶解析与前 5 个谐波拟合 ---
    Y = fft(mean_cycle);
    L = length(mean_cycle);
    k = 4; 
    mask = zeros(1, L);
    mask(1:k+1) = 1;           
    mask(end-k+1:end) = 1;     
    Y_filtered = Y .* mask;
    fit_cycle = real(ifft(Y_filtered)); 
    
    % --- 10. 提取参数 ---
    amplitudes = abs(Y(1:k+1)) / (L/2);
    amplitudes(1) = amplitudes(1) / 2; 
    phases = angle(Y(1:k+1));
    
    % --- 11. 引入真实时间尺度 ---
    t_fit = (0:target_len-1) * dt_fit;   % 0 .. T_period-dt_fit，步长 = dt_fit
    if(ifPlot)
        % --- 12. 绘图展示 (确保窗口不重叠) ---
        close all;
        % 窗口 1: 原始数据与检测
        figure(1); 
        set(gcf, 'Color', 'w', 'Name', '原始信号分析', 'Position', [100, 500, 600, 400]);
        subplot(2,1,1);
        plot((0:length(abp_mmHg)-1)/fs, abp_mmHg, 'Color', [0.7 0.7 0.7]); hold on;
        plot(locs/fs, pks, 'rv', 'MarkerFaceColor', 'r');
        title('原始多周期信号及峰值检测');
        xlabel('时间 (s)'); ylabel('压力 (mmHg)');
        grid on;
        
        subplot(2,1,2);
        fill([t_fit, fliplr(t_fit)], [mean_cycle+std_cycle, fliplr(mean_cycle-std_cycle)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
        plot(t_fit, mean_cycle, 'r', 'LineWidth', 2);
        title(['重排后的平均周期 (时间轴还原, n=', num2str(num_cycles), ')']);
        xlabel('时间 (s)'); ylabel('压力 (mmHg)');
        grid on;
        
        % 窗口 2: 傅里叶拟合对比
        figure(2);
        set(gcf, 'Color', 'w', 'Name', '傅里叶拟合', 'Position', [750, 500, 600, 400]);
        plot(t_fit, mean_cycle, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', '原始平均数据');
        hold on;
        plot(t_fit, fit_cycle, 'b', 'LineWidth', 2, 'DisplayName', '前5个谐波拟合');
        title(['单个心动周期的傅里叶重构 (T = ', num2str(T_period), 's)']);
        xlabel('时间 (s)');
        ylabel('压力 (mmHg)');
        legend('Location', 'best');
        grid on;
    end        
    % --- 13. 打印结果 ---
    fprintf('\n处理完成！\n');
    fprintf('设定周期长度: %.2f s (心率约 %.1f bpm)\n', T_period, 60/T_period);
    fprintf('输出时间轴: %d 点, 步长 %.4f s\n', target_len, dt_fit);
    fprintf('傅里叶拟合参数已提取，可用于公式构建。\n');
end
