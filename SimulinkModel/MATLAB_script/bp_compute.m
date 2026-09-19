function S = bp_compute(SCENARIO_CODE)
% BP_COMPUTE  Compute every model workspace parameter for one scenario.
%
%   S = bp_compute('B')
%
% Returns a struct whose fields are exactly the variable names the Simulink
% model resolves from the base workspace (k_m, R_AG, P_ECS0, central_cycle,
% t_fit, ...).
%
% This is the single source of truth for the model parameters.  The script
% Brain_Parameters.m is a thin wrapper around it, so the interactive workflow
% is unchanged while batch drivers can call this function directly.
%
% NOTE: no `clear` is used anywhere here -- the function scope already
% guarantees that no stale variable leaks into the computation.

    if nargin < 1 || isempty(SCENARIO_CODE)
        SCENARIO_CODE = 'B';
    end

    %% Table 2 shorthands
    %   The A3 combined scenarios are written as H+C and HA+C in the manuscript,
    %   i.e. grade-2 hypertension -- chronic or acute -- together with all three
    %   outflow pathways at 500 %.  They are expanded here to the explicit form
    %   so that one canonical code reaches the rest of the parser and the
    %   runners; the file name keeps whatever the caller passed in.
    switch upper(SCENARIO_CODE)
        case 'H+C',  SCENARIO_CODE = 'H_G2+C_OVS5';
        case 'HA+C', SCENARIO_CODE = 'H_G2A+C_OVS5';
    end

    %% Scenario codes -- Table 2 of the manuscript
    %   A combined code is split on '+' and each part supplies one aspect of the
    %   parameter set:
    %       H_*          vasculature (pulseType, PVI, r_R_Hyper, r_E, r_BBB)
    %       C_<OVS><a>   CSF outflow resistance multiples; the pathway letters are
    %                    those of Table 2:
    %                        O = olfactory lymphatic    (R_CP)
    %                        V = arachnoid granulation  (R_AG)
    %                        S = spinal lymphatic       (R_SPI)
    %       C_Ref        single-outflow reference scenario (AG alone)
    %       Ven_R<n>     ventriculo-SAS resistance multiple
    %       Ven_E<n>     ventricular elastance multiple
    %       S_*          sensitivity override (supplementary series)
    %   A single code is simply a one-element split.
    parts = strsplit(SCENARIO_CODE, '+');
    vascCode = 'B'; cCode = ''; venRCode = ''; venECode = ''; sCode = '';
    for q = 1:numel(parts)
        p = strtrim(parts{q});
        if isempty(p), continue; end
        if ~isempty(regexp(p, '^C_', 'once'))
            cCode = p;
        elseif ~isempty(regexp(p, '^Ven_R', 'once'))
            venRCode = p;
        elseif ~isempty(regexp(p, '^Ven_E', 'once'))
            venECode = p;
        elseif ~isempty(regexp(p, '^S_', 'once'))
            sCode = p;             % sensitivity override
        else
            vascCode = p;          % H_* (or B)
        end
    end

    %% Public Parameters
    [pulseType, PVI, r_R_Hyper, r_E, r_BBB] = getSimulinkParamsFromSBP(vascCode);
    abpFile = fullfile(fileparts(mfilename('fullpath')), '..', 'ABP_Data', pulseType);
    [mean_cycle, t_fit] = ABP_Get(abpFile, false);
    central_cycle = ABP_to_Central_GTF(mean_cycle, t_fit, false);
    k_m = 0.4343 * PVI; % proportional constant(Rekate, 1988)
    r_R = 1;  % resistence ratio to control exit flow rate
    P_ECS0 = 10;
    T = 0.8;
    %% CSF outflow configuration codes (A2, A2-ref)
    %   C_<pathways><alpha> : e.g. C_OS2.5, C_OVS5 -- multiply the resistance of
    %                         the named pathways by alpha (2.5 or 5 in Table 2).
    %   C_Ref               : single-outflow reference scenario: the arachnoid
    %                         granulation route ALONE, with its resistance set to
    %                         the PARALLEL EQUIVALENT of the three baseline
    %                         pathways and the other two occluded,
    %                             R_AG = Rout_base
    %                                  = 1/(1/25 + 1/65 + 1/20) = 9.489051
    %                         The total equivalent outflow resistance therefore
    %                         equals the baseline three-pathway value, so the
    %                         comparison isolates the effect of the pathway
    %                         STRUCTURE from the effect of the resistance level.
    %                         This is the conventional AG-only model and the
    %                         analogue of Vinje et al. model 0.
    %                         R_FM is left at its physiological value, i.e. the
    %                         craniospinal communication stays open and the
    %                         spinal SAS acts as a compliant extension of the
    %                         SAS -- exactly what a lumped single-chamber model
    %                         assumes.
    kO = 1; kV = 1; kS = 1; kAlpha = 1;
    singleOnly = '';
    % "infinite" resistance for a blocked pathway.  A finite 1e9 is used
    % instead of Inf because Inf propagates as Inf/Inf -> NaN through the
    % nonlinear resistance subsystems (e.g. SASE/Divide1) and aborts the
    % simulation with Simulink:Engine:AlgStateNotFinite.  1e9 mmHg*min/mL
    % corresponds to a conductance of 1e-9, i.e. leakage < 1e-8 of the open
    % pathway, so it is a true no-flow condition numerically.
    R_BLOCK = 1e9;
    tok = regexp(cCode, '^C_([OVS]+)(\d+(?:\.\d+)?)$', 'tokens', 'once');
    if ~isempty(tok)
        kAlpha = str2double(tok{2});
        if any(tok{1} == 'O'), kO = kAlpha; end
        if any(tok{1} == 'V'), kV = kAlpha; end
        if any(tok{1} == 'S'), kS = kAlpha; end
    elseif strcmp(cCode, 'C_Ref')
        singleOnly = 'V';               % arachnoid granulation alone
    end
    blockedO = ~isempty(singleOnly) && ~strcmp(singleOnly,'O');
    blockedV = ~isempty(singleOnly) && ~strcmp(singleOnly,'V');
    blockedS = ~isempty(singleOnly) && ~strcmp(singleOnly,'S');
    % baseline equivalent outflow resistance of the three-pathway model
    Rout_base = 1/(1/(25*r_R) + 1/(65*r_R) + 1/(20*r_R));

    %% Ventricular configuration codes (A3-R obstruction, A3-E elastance)
    %   Ven_R<n> : multiply the ventriculo-SAS (aqueduct / ventricular outlet)
    %   resistance R_VEN by the nominal factor n.  The two study levels are
    %   pinned to independent clinical anchors, so their actual multipliers are
    %   2.0560 and 7.4638 rather than the round numbers in the code.
    %   The three absorption pathways stay at baseline: Bech-Azeddine 2007
    %   reports lumbar Rout within the normal range (11.0, 6.7-14.0) in
    %   aqueductal stenosis, i.e. SAS absorption is preserved.
    %   Ven_E<n> : multiply the ventricular elastance by n (0.5, 0.25).  Table 2
    %   requires the unstressed ventricular volume to be held at its baseline
    %   value, so V_VEN_free is always computed from E_VEN_base, never from the
    %   reduced E_VEN.
    %   Both codes were extracted by the '+' split above; kVEN / kEVEN are
    %   resolved next to R_VEN (they need Q_CSF, defined with the capillary
    %   parameters).

    %% Arteries parameters (Zagzoule and Marc-Vergnes, 1986, Segment number 1~21)
    %  The vessel resistances are proportionally adjusted to a level approaching physiological blood flow value)
    P_A0 = 95.5;                   % mmHg, Arteriolar pressure  (Doron et al, 2021)
    b_EA = 0.004;
    E_A_ref = r_E * 0.4 * 10000/133.32; % mmHg/mL, Vessel elastance
    E_A0 = E_A_ref / exp(109*b_EA);
    R_A0 = 0.091 * r_R_Hyper;                % mmHg*min/mL, Poiseuille flow resistance, R = 8 * miu * L / (pi * R^4)
    V_A0 = 45;                % mL, Volume of arteries
    V_Afree = V_A0 - (P_A0 - P_ECS0)/(E_A_ref / r_E);

    %% Capillaries (Zagzoule and Marc-Vergnes, 1986, Segment number 22~24)
    E_C = r_E  * 0.25 * 10000/133.32;
    R_C0 = 0.0182 * r_R_Hyper;
    V_C0 = 20;
    R_BBB = 128/r_BBB; % (Vinje et al, 2020)
    Q_CSF = 0.3;
    P_C0 = 30;
    pi_C0 = 25; % (Vinje et al, 2020)， we distinguish crystalloid osmotic pressure and colloid osmotic pressure
    V_Cfree = V_C0 - (P_C0 - P_ECS0)/(E_C / r_E);

    %% Veins (Zagzoule and Marc-Vergnes, 1986, Segment number 25~31)
    E_V = r_E * 0.2* 10000/133.32;
    R_V0 = 0.013 * r_R_Hyper;
    V_V0 = 75;
    P_V0 = 17;
    V_Vfree = V_V0 - (P_V0 - P_ECS0)/(E_V / r_E);

    %% Ventricles(Doron et al, 2021, elastance refer to Rekate et al, 1988)
    R_VEN_base = 250/60;            % mmHg*min/mL, baseline 4.1667
    switch venRCode
        case 'Ven_R2'
            % Bech-Azeddine 2007 (n=11): VRout - LRout = 4.4 mmHg*min/mL in
            % aqueductal stenosis.  A ventricular infusion must traverse the
            % obstruction while a lumbar infusion reaches the SAS directly, so
            % the obstruction ADDS 4.4 to the ventriculo-SAS resistance.
            kVEN = 1 + 4.4/R_VEN_base;                  % = 2.0560
        case 'Ven_R7.5'
            % Tabibkhooei 2025: ETV lowers ventricular pressure 19.33 -> 10.00
            % mmHg, i.e. the obstruction carries dP = 9.33 mmHg at Q_CSF.
            kVEN = (9.33/Q_CSF)/R_VEN_base;             % = 7.4638
        case ''
            kVEN = 1;
        otherwise
            tk = regexp(venRCode, '^Ven_R(\d+(?:\.\d+)?)$', 'tokens', 'once');
            if isempty(tk)
                error('无法解析心室梗阻代号 %s（应为 Ven_R<n>）', venRCode);
            end
            kVEN = str2double(tk{1});                   % free multiple
    end
    V_VEN0 = 25;
    R_VEN = R_VEN_base * kVEN;

    % ventricular elastance: Ven_E<n> multiplies E_VEN_base by n
    E_VEN_base = 3; % mmHg/mL
    kEVEN = 1;
    if ~isempty(venECode)
        tk = regexp(venECode, '^Ven_E(\d+(?:\.\d+)?)$', 'tokens', 'once');
        if isempty(tk)
            error('无法解析心室弹性代号 %s（应为 Ven_E<n>）', venECode);
        end
        kEVEN = str2double(tk{1});
    end
    E_VEN = E_VEN_base * kEVEN;
    % Table 2: the unstressed ventricular volume is held at its baseline value,
    % so it is derived from E_VEN_base and not from the reduced E_VEN.
    V_VEN_free = V_VEN0 - (1.7/E_VEN_base);

    %% Subarachnoid Space(Doron et al, 2021)
    R_SAS = 280/60;
    V_SAS0 = 35; % not include spinal SAS
    E_SAS = 10;

    %% CMC
    R_FM = 0.01;

    %% spinal subarachnoid space(Vinje, 2020)
    P_SPI0 = 10.4;
    P_SLS = 9;
    if blockedS
        R_SPI = R_BLOCK;           % blocked in the single-outflow reference model
    elseif strcmp(singleOnly,'S')
        R_SPI = Rout_base;         % sole pathway -> parallel equivalent of the three
    else
        R_SPI = 20 * r_R * kS;
    end
    V_SPI0 = 80;
    f_SPI = 0.69;

        %% Olfactory lymphatic pathway (Vinje, 2020) -- cribriform plate
    P_OLS = 0;
    if blockedO
        R_CP = R_BLOCK;          % blocked in the single-outflow reference model
    elseif strcmp(singleOnly,'O')
        R_CP = Rout_base;        % sole pathway -> parallel equivalent of the three
    else
        R_CP = 65 * r_R * kO;
    end

    %% SSS
    P_SSS = 7.5;  % (Doron et al, 2021)
    if blockedV
        R_AG = R_BLOCK;            % blocked in the single-outflow reference model
    elseif ~isempty(singleOnly)
        R_AG = Rout_base;          % sole pathway -> parallel equivalent of the three
    else
        R_AG = 25 * r_R * kV; % (Vinje et al, 2020)
    end

    %% Paravascular Space(Vinje, 2020 cite from Faghih, 2018)
    R_PAS = 1.14;
    R_PNS = 1.75/1000;
    R_PCS = 32.24;
    R_IEGL = 0.57; % (Vinje, 2020)
    R_IEGH = 0.57*2; % (Vinje, 2020)
    r_PVS = 1.55; % average radius of PVS compared to its correspond vessel
    V_PAS0 = 15 * (r_PVS^2 - 1);
    V_PNS0 = 30 * (r_PVS^2 - 1);
    E_IEG = 10000/133.32/V_PAS0; % (Gan et al, 2023)
    V_PAS_free = V_PAS0 - (P_ECS0/E_IEG);

    %% Extracellular Space & Intracellular space
    % The volumes of these two spaces are calculated by subtracting the other volumes from the constant total volume.
    % It is considered that 80% of the rest of volume is ICS while the other 20% is ECS.
    % Thus the elasticity and volume of them is included in other spaces.
    V_brain0 = 1300;
    V_brain0  = V_brain0 - V_PAS0 - V_PNS0;
    alpha_w = 0.8; % water proportion in brain tissue
    alpha_ECS = 0.2;                                % volume ratio  of ECS in Whole brain, (Yousefnezhad, 2016)
    R_ECS = 0.57 + 0.64; % (Vinje, 2020)
    V_ECS0 = 0.2 * alpha_w * V_brain0;
    V_ICS0 = 0.8 * alpha_w * V_brain0;
    V_BT_solid = V_brain0 * (1 - alpha_w); % solid part of brain
    pi_ECS0 = 10; % (Vinje et al, 2020)， we distinguish crystalloid osmotic pressure and colloid osmotic pressure
    f_ECS = 1- f_SPI;

    %% other public parameters
    V_T0 = V_VEN0 + V_SAS0 + (V_A0 + V_C0 + V_V0) + V_brain0 + V_PAS0 + V_PNS0; % Total cranial volume, constant

    %% ---- B-series sensitivity overrides: S_<param>_<value> --------------
    %   Applied AFTER every parameter is computed, so that an override can also
    %   refresh the quantities derived from it.  Combine with '+' to change the
    %   base configuration, e.g. 'S_bEA_0+H_G2A'.
    %
    %   WAVE  RADIAL | CENTRAL | AMP<pp>   input pressure waveform
    %   HR    <bpm>                        heart rate (cycle period)
    %   PVI   <mL>                         pressure-volume index
    %   bEA   <mmHg^-1>                    arterial pressure-sensitivity coeff.
    %   RBBB  <multiple>                   r_BBB (R_BBB = 128/r_BBB)
    %   RFM   <mmHg*min/mL>                craniospinal (foramen magnum) R
    %   PSPI0 <mmHg>                       spinal outflow sink pressure
    %   RECS  <multiple>                   R_ECS = 0.57 * multiple
    %   SWAP  AGMAX                        swap the AG and cribriform resistances
    %   VALVE OFF                          make R_IEG constant (no rectification)
    sApply = ''; sValue = '';
    if ~isempty(sCode)
        tk = regexp(sCode, '^S_([A-Za-z][A-Za-z0-9]*)_(.+)$', 'tokens', 'once');
        if isempty(tk)
            error('无法解析敏感性代号 %s（应为 S_<param>_<value>）', sCode);
        end
        sApply = tk{1}; sValue = tk{2};
    end
    HR_eff = NaN;    % reported for the HR scenarios

    switch sApply
        case ''
            % no override
        case 'PVI'
            PVI = str2double(sValue);
            k_m = 0.4343 * PVI;                 % Rekate 1988
            % Since the 2026-09-15 revision E_VEN (Ventricles) and E_SAS (SAS)
            % are explicit, literature-anchored parameters, so a PVI override must
            % NOT recompute them.  It only rescales k_m, which is what the cranial
            % (f_ECS) and spinal (f_SPI) pressure-volume relations use.
        case 'bEA'
            b_EA = str2double(sValue);
            % E_A0 is defined so that E_A = E_A_ref exactly at P_A = 116 mmHg
            % (the Grade-2 MAP of Table 2), which is why b_EA changes the
            % BASELINE arterial elastance but not the Grade-2 reference.
            E_A0 = E_A_ref / exp(116*b_EA);
        case 'RBBB'
            R_BBB = 128 / str2double(sValue);
        case 'RFM'
            R_FM = str2double(sValue);
        case 'PSPI0'
            P_SPI0 = str2double(sValue);
        case 'RECS'
            R_ECS = 0.57 * str2double(sValue);
        case 'fSPI'
            % f_SPI is ACTIVE in the spinal pressure subsystem (Int P_SPI): the
            % spinal P-V relation is built from k_m and f_SPI, i.e. the
            % craniospinal compliance constant k_m = 0.4343*PVI is PARTITIONED
            % between the intracranial (f_ECS) and spinal (f_SPI) compartments.
            % f_SPI = 0.69 therefore makes the spinal compartment the more
            % compliant of the two, which is what the MRI-based partition asks
            % for.  There is no spinal elastance parameter in the model: the
            % spinal pressure comes from the exponential P-V relation, exactly
            % as the cranial pressure does.
            f_SPI = str2double(sValue);
            f_ECS = 1 - f_SPI;
        case 'SWAP'
            if ~strcmpi(sValue,'AGMAX')
                error('S_SWAP 只支持 AGMAX');
            end
            % Put the LARGEST of the three resistances on the AG pathway by
            % swapping AG and cribriform.  The multiset {25,65,20} and hence the
            % total Rout are unchanged, so this isolates the ASSIGNMENT of the
            % resistances from their LEVEL -- which is exactly Reviewer 1 #4.
            tmpSwap = R_AG; R_AG = R_CP; R_CP = tmpSwap;
        case 'WAVE'
            u = upper(sValue);
            if strcmp(u,'RADIAL')
                central_cycle = mean_cycle;     % measured radial waveform
            elseif strcmp(u,'CENTRAL')
                % unchanged: central_cycle from ABP_to_Central_GTF
            else
                tk2 = regexp(u, '^AMP(\d+)$', 'tokens', 'once');
                if isempty(tk2), error('未知 S_WAVE 档位 %s', sValue); end
                amp = str2double(tk2{1})/100;   % AMP120 -> x1.20
                mc = mean(central_cycle);
                central_cycle = mc + (central_cycle - mc)*amp;   % scale pulse about the mean
            end
        case 'HR'
            hr = str2double(sValue);
            % The Repeating Sequence Interpolated block samples the table on its
            % own 0.01 s grid, so the EFFECTIVE cycle is
            %     (floor(span/0.01) + 1) * 0.01
            % (measured: span [0,0.799] -> 0.800 s; span [0,0.800] -> 0.810 s).
            % The span is therefore chosen so the effective period lands EXACTLY
            % on a 0.01 s multiple and never on a table endpoint -- this is the
            % rule that prevented the spurious 64 s moire oscillation.
            T_target = 60/hr;
            n100  = max(2, round(T_target/0.01));   % 0.01 s samples per cycle
            T_eff = n100*0.01;
            span  = (n100-1)*0.01 + 0.009;          % inside [., T_eff)
            dtf   = t_fit(2) - t_fit(1);
            T_old = t_fit(end) + dtf;               % period of the table we have
            nNew  = round(span/dtf) + 1;
            tNew  = (0:nNew-1)*dtf;
            src   = tNew * (T_old/T_eff);
            cc    = interp1([t_fit(:); T_old], [central_cycle(:); central_cycle(1)], ...
                            src(:), 'linear');
            t_fit = tNew(:);
            central_cycle = cc(:);
            HR_eff = 60/T_eff;
        case 'VALVE'
            % R_IEG is a two-valued Switch in the PAS subsystem: the open
            % (low) branch R_IEGL when P_PAS >= P_ECS and the closed (high)
            % branch R_IEGH otherwise.  These two settings remove the switching
            % altogether and hold the resistance at one of its two values for
            % the whole run:
            %   S_VALVE_OFF  : R_IEG = R_IEGL everywhere (valve held open)
            %   S_VALVE_HIGH : R_IEG = R_IEGH everywhere (valve held closed)
            if strcmpi(sValue,'OFF')
                R_IEGH = R_IEGL;      % both branches low
            elseif strcmpi(sValue,'HIGH')
                R_IEGL = R_IEGH;      % both branches high
            else
                error('S_VALVE 只支持 OFF 或 HIGH');
            end
            r_IEG  = 1;               % legacy ratio; the model no longer uses it
        otherwise
            error('未知的敏感性参数 %s（见 bp_compute 头部的表）', sApply);
    end

    %% ---- pack every local variable into a struct ----
    names = who;
    S = struct();
    for i = 1:numel(names)
        S.(names{i}) = eval(names{i});
    end
end

% =====================================================================
function [pulseType, PVI, r_R_Hyper, r_E, r_BBB] = getSimulinkParamsFromSBP(groupID)
    % C_* codes (CSF outflow resistance impairment) keep the baseline
    % vasculature; only the outflow resistances differ, and bp_compute scales
    % those separately.
    if ischar(groupID) && ~isempty(regexp(groupID, '^C_', 'once'))
        pulseType   = '120.txt';
        PVI         = 25;
        r_R_Hyper   = 1;
        r_E         = 1;
        r_BBB       = 1;
        return
    end

    % Ven_* codes (A3: ventriculo-SAS resistance or ventricular elastance
    % scaled) likewise keep the baseline vasculature; only R_VEN / E_VEN
    % differ, both scaled in bp_compute.
    if ischar(groupID) && ~isempty(regexp(groupID, '^Ven_', 'once'))
        pulseType   = '120.txt';
        PVI         = 25;
        r_R_Hyper   = 1;
        r_E         = 1;
        r_BBB       = 1;
        return
    end

    switch groupID
        case {1, 120, 'Baseline', 'B'}
            pulseType   = '120.txt';
            PVI        = 25;
            r_R_Hyper  = 1;
            r_E        = 1;
            r_BBB      = 1;

        case {2, 132, 'HighNormal', 'H_HM'}
            pulseType   = '132.txt';
            PVI        = 22;
            r_R_Hyper  = 1.06;
            r_E        = 1.24;
            r_BBB      = 1.05;

        case {3, 144, 'Grade 1', 'H_G1'}
            pulseType   = '144.txt';
            PVI        = 20;
            r_R_Hyper  = 1.11;
            r_E        = 1.51;
            r_BBB      = 1.10;

        case {4, 164, 'Grade 2', 'H_G2'}
            pulseType   = '164.txt';
            PVI        = 15;
            r_R_Hyper  = 1.22;
            r_E        = 1.96;
            r_BBB      = 1.20;

        case {5, 'Acute Grade 2', 'H_G2A'}
            pulseType   = '164.txt';
            PVI        = 25;
            r_R_Hyper  = 1.;
            r_E        = 1.;
            r_BBB      = 1.;

        otherwise
            error(['无效的 scenario 代号。允许: B, H_HM, H_G1, H_G2, H_G2A, ' ...
                   'C_<OVS...><alpha>, C_Ref, Ven_R<n>, Ven_E<n>, ' ...
                   '以及用 + 组合的 A4 联用码（如 H_G2+Ven_R7.5、H_G2+C_OVS5）']);
    end
end
