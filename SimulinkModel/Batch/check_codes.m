%% check_codes.m -- smoke test of every scenario code in Table 2.
%
% Calls bp_compute for each code and checks the parameters it must produce.  No
% simulation, so it takes seconds and catches a broken code parser before a
% multi-hour batch.  Supersedes check_a4_codes / check_c_codes / check_o_codes.
%
% The expectation table is padded to a constant width because a MATLAB cell
% literal cannot have rows of different lengths.
%
% ASCII only.  matlab -batch check_codes
function check_codes()
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(here);

C = {
 'B'              'R_CP' 65        'R_AG' 25        'R_SPI' 20        'R_VEN' 4.166667   'E_VEN' 3
 'H_HM'           'R_VEN' 4.166667 'E_VEN' 3        'PVI' 22          '' ''              '' ''
 'H_G1'           'R_VEN' 4.166667 'E_VEN' 3        'PVI' 20          '' ''              '' ''
 'H_G2'           'R_VEN' 4.166667 'E_VEN' 3        'PVI' 15          '' ''              '' ''
 'H_G2A'          'R_VEN' 4.166667 'E_VEN' 3        'PVI' 25          '' ''              '' ''
 'C_O2.5'         'R_CP' 162.5     'R_AG' 25        'R_SPI' 20        '' ''              '' ''
 'C_O5'           'R_CP' 325       'R_AG' 25        'R_SPI' 20        '' ''              '' ''
 'C_S5'           'R_CP' 65        'R_AG' 25        'R_SPI' 100       '' ''              '' ''
 'C_V2.5'         'R_CP' 65        'R_AG' 62.5      'R_SPI' 20        '' ''              '' ''
 'C_OS2.5'        'R_CP' 162.5     'R_AG' 25        'R_SPI' 50        '' ''              '' ''
 'C_OV5'          'R_CP' 325       'R_AG' 125       'R_SPI' 20        '' ''              '' ''
 'C_VS5'          'R_CP' 65        'R_AG' 125       'R_SPI' 100       '' ''              '' ''
 'C_OVS2.5'       'R_CP' 162.5     'R_AG' 62.5      'R_SPI' 50        '' ''              '' ''
 'C_OVS5'         'R_CP' 325       'R_AG' 125       'R_SPI' 100       '' ''              '' ''
 'C_Ref'          'R_CP' 1e9       'R_AG' 9.489051  'R_SPI' 1e9       '' ''              '' ''
 'Ven_R2'         'R_VEN' 8.566667 'E_VEN' 3        '' ''             '' ''              '' ''
 'Ven_R7.5'       'R_VEN' 31.100   'E_VEN' 3        '' ''             '' ''              '' ''
 'Ven_E0.5'       'R_VEN' 4.166667 'E_VEN' 1.5      'V_VEN_free' 24.433333 '' ''         '' ''
 'Ven_E0.25'      'R_VEN' 4.166667 'E_VEN' 0.75     'V_VEN_free' 24.433333 '' ''         '' ''
 'H_G2+Ven_R7.5'  'R_VEN' 31.100   'E_VEN' 3        'PVI' 15          '' ''              '' ''
 'H+C'    'R_VEN' 4.166667 'R_CP' 325       'R_AG' 125        'R_SPI' 100        'PVI' 15
 'HA+C'   'R_CP' 325       'R_AG' 125       'R_SPI' 100       'PVI' 25           '' ''
 'H_G2+C_OVS2.5+Ven_E0.25' 'R_VEN' 4.166667 'R_CP' 162.5 'R_AG' 62.5   'R_SPI' 50         'E_VEN' 0.75
};

fprintf('=== check_codes: %d scenario codes ===\n', size(C,1));
nBad = 0;
for k = 1:size(C,1)
    code = C{k,1};
    S = bp_compute(code);
    bad = {};
    for j = 2:2:size(C,2)
        f = C{k,j};
        if isempty(f), continue; end
        want = C{k,j+1};
        if ~isfield(S,f)
            bad{end+1} = sprintf('%s missing', f); %#ok<AGROW>
            continue
        end
        got = S.(f);
        if ~isscalar(got) || ~isnumeric(got)
            bad{end+1} = sprintf('%s not scalar', f); %#ok<AGROW>
            continue
        end
        if abs(got - want) > 1e-6*max(1,abs(want))
            bad{end+1} = sprintf('%s = %.6f, expected %.6f', f, got, want); %#ok<AGROW>
        end
    end
    if isempty(bad)
        fprintf(['  %-14s OK   R_CP %-9.4g R_AG %-9.5g R_SPI %-9.4g ' ...
                 'R_VEN %-9.4g E_VEN %-6.3g V_VEN_free %.4f\n'], ...
            code, S.R_CP, S.R_AG, S.R_SPI, S.R_VEN, S.E_VEN, S.V_VEN_free);
    else
        nBad = nBad + numel(bad);
        fprintf('  %-14s *** %s\n', code, strjoin(bad, ' ; '));
    end
end

fprintf('\n=== check_codes: %s (%d problem(s)) ===\n', ...
        ternary(nBad==0,'PASS','FAIL'), nBad);
if nBad > 0, error('check_codes:FAIL','%d code expectation(s) failed', nBad); end
end
% =====================================================================
function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
