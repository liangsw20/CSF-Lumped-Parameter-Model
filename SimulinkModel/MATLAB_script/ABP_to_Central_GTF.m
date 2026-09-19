function central_cycle = ABP_to_Central_GTF(mean_cycle, t_fit, ifPlot)
% ABP_TO_CENTRAL_GTF  Converts the averaged radial pulse waveform (as
% returned by ABP_Get) into a central aortic pressure waveform using a
% generalized transfer function.
%
%   [central_cycle, t_fit] = ABP_to_Central_GTF(mean_cycle, t_fit, ifPlot)
%
%   Inputs:
%       mean_cycle  -  averaged radial pulse waveform over one cardiac
%                      cycle (mmHg), i.e. the first output of ABP_Get.
%                      Row (1xN) or column (Nx1) vector.
%       t_fit       -  time axis of mean_cycle (s), i.e. the second output
%                      of ABP_Get.  It must NOT contain the cycle endpoint
%                      T_period (ABP_Get returns e.g. (0:799)*0.001, i.e.
%                      span [0, 0.799]); see the note in section 3 below.
%       ifPlot      -  optional logical: plot the radial vs central
%                      waveforms (default false).
%
%   Outputs (same format as ABP_Get):
%       central_cycle  -  1 x N central aortic pressure waveform over one
%                         cardiac cycle (mmHg), aligned like mean_cycle
%       t_fit          -  1 x N time axis (same grid as the input t_fit)
%
%   Pipeline: averaged radial cycle -> harmonic truncation + resample
%             -> per_cen_TF (Karamanoglu radial -> central aortic GTF)
%
%   Dependencies: per_cen_TF.m (P.H. Charlton, MIT license,
%                 https://github.com/peterhcharlton/transfer_functions)
%   References:   Karamanoglu et al. Eur Heart J 1993;14:160-167
%                 Karamanoglu & Feneley. Am J Physiol 1996;271:H2399-H2404
%                 Chen et al. Circulation 1997;95:1827-1836

if nargin < 3, ifPlot = false; end
if nargin < 2
    error('ABP_to_Central_GTF requires mean_cycle and t_fit (outputs of ABP_Get).');
end

% ensure the dependency function is on the path
funcDir = fileparts(mfilename('fullpath'));
if isempty(funcDir), funcDir = pwd; end
addpath(funcDir);

v_rad    = mean_cycle(:);            % mmHg, column
t_old    = t_fit(:);                 % s, column
T_period = max(t_old) - min(t_old);  % cardiac cycle duration (e.g. 0.8 s)
n_out    = numel(t_fit);             % output points per cycle (e.g. 801)
fs_target = 200;                     % resampling frequency (Hz) for the GTF

% ---- 1. harmonic truncation + resample to fs_target ---------------------
% The input is the AVERAGED radial cycle from ABP_Get. Truncate to the
% first 12 harmonics (<= 15 Hz) to remove residual noise before
% re-sampling (the GTF Bode data are only defined up to ~9-12 Hz), then
% interpolate onto the target grid.
n_rs  = round(T_period * fs_target);        % 160 pts @ 200 Hz
t_new = (0 : n_rs-1)' / fs_target;          % s, within one cycle
Y = fft(v_rad);
L = numel(v_rad);
nKeep = min(12, floor(L/2));                % keep harmonics 1..nKeep
Y(nKeep+2 : end-nKeep) = 0;
v_filt = real(ifft(Y));
v_rs = interp1(t_old, v_filt, t_new, 'pchip');
v_rs = v_rs(:);

% ---- 2. averaged radial cycle -> central aortic via GTF ------------------
S.v  = v_rs;
S.fs = fs_target;
opts.do_plot           = false;
opts.retain_onset_time = true;    % keep the same time origin as the input
transformed = per_cen_TF(S, 'radiaABP', opts);
v_cent = transformed.centrABP.v(:);
if numel(v_cent) ~= n_rs            % per_cen_TF may drop one endpoint
    v_cent = interp1(linspace(0,1,numel(v_cent)), v_cent, ...
                     linspace(0,1,n_rs), 'pchip');
    v_cent = v_cent(:);             % interp1 returns a row for row xi
end

% ---- 3. resample back onto the caller's grid ----------------------------
% The output must sit on EXACTLY the grid that was passed in.  Do NOT append
% the cycle endpoint: the input table must not contain t = T_period, because
% the Repeating Sequence Interpolated block then counts one extra sample time
% per period (0.81 s instead of 0.80 s, measured in this model), which folds
% against the nominal 0.8 s cycle and produces a spurious 64.8 s oscillation.
t_out = t_old;                                % same grid as the input
t_ext = [t_new; t_new(1) + T_period];         % one period of extension
v_ext = [v_cent; v_cent(1)];
central_cycle = interp1(t_ext, v_ext, t_out, 'pchip')';   % 1 x n_out row vector
t_fit = t_old;                                % same grid (interface compatibility)

% ---- 4. optional plot ----------------------------------------------------
if ifPlot
    figure('Name','Radial vs Central (GTF)','NumberTitle','off','Color','w');
    plot(t_old, v_rad,  '-', 'LineWidth', 1.4); hold on;
    plot(t_fit, central_cycle, '-', 'LineWidth', 1.6);
    legend({'Radial (averaged input)', 'Central (GTF)'}, 'Box','off');
    xlabel('Time (s)'); ylabel('Pressure (mmHg)');
    title(sprintf('radial SBP %.1f -> central SBP %.1f mmHg', ...
        max(v_rad), max(central_cycle)));
    grid on;
end

end
