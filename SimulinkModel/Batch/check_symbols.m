%% check_symbols.m -- pre-flight: does the model's symbol set match sym_rename_map.txt?
%
% Reports, for the model's four symbol stores (block parameters, Goto tags,
% signal line names, per-port signal-logging names):
%   * any OLD symbol still present          -> FAIL
%   * which NEW symbols are present         -> OK
%   * which NEW symbols are missing         -> NOTE (some are script-only)
% and then runs a 2 s simulation and checks the logsout labels themselves.
%
% Run this BEFORE a batch: a signal whose label silently reverts to the old name
% would only show up as a blank column after hours of compute.
%
% ASCII only.  matlab -batch check_symbols
function check_symbols()
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
mdl  = 'CSF_Model';
addpath(genpath(fullfile(root,'SimulinkModel','MATLAB_script')));
addpath(here); addpath(fullfile(root,'SimulinkModel'));
if ~bdIsLoaded(mdl), load_system(fullfile(root,'SimulinkModel',[mdl '.slx'])); end

[old, new, kind] = local_map(fullfile(here,'sym_rename_map.txt'));
ok = true;

%% ---- collect every symbol the model stores ------------------------------
WL = {'Value','Gain','GotoTag','Amplitude','InitialCondition','UpperLimit', ...
      'OutValues','TimeValues','BreakpointsForDimension1','DataSpecification'};
blks = find_system(mdl,'FindAll','on','LookUnderMasks','all','FollowLinks','on', ...
                      'Type','block');
sym = {};                       % {'store','symbol'}
for i = 1:numel(blks)
    b = blks(i);
    try, nm = get_param(b,'Name'); catch, nm = ''; end
    if ischar(nm) && ~isempty(nm), sym{end+1} = {'block:Name', nm}; end %#ok<AGROW>
    for k = 1:numel(WL)
        try, v = get_param(b, WL{k}); catch, continue; end
        if ~ischar(v) || isempty(v), continue; end
        t = regexp(v, '[A-Za-z_]\w*', 'match');
        for j = 1:numel(t), sym{end+1} = {['param:' WL{k}], t{j}}; end %#ok<AGROW>
    end
    try, ph = get_param(b,'PortHandles'); catch, continue; end
    for fld = {'Outport','Inport'}
        if ~isfield(ph, fld{1}), continue; end
        for k = 1:numel(ph.(fld{1}))
            for pf = {'DataLoggingName','Name'}
                try, v = get_param(ph.(fld{1})(k), pf{1}); catch, continue; end
                if ischar(v) && ~isempty(v)
                    sym{end+1} = {['port:' pf{1}], v}; %#ok<AGROW>
                end
            end
        end
    end
end
lh = find_system(mdl,'FindAll','on','LookUnderMasks','all','FollowLinks','on', ...
                    'Type','line');
for i = 1:numel(lh)
    try, nm = get_param(lh(i),'Name'); catch, continue; end
    if ischar(nm) && ~isempty(nm), sym{end+1} = {'line', nm}; end %#ok<AGROW>
end
fprintf('=== check_symbols: %d symbol occurrences collected ===\n', numel(sym));

%% ---- 1. any OLD symbol left? -------------------------------------------
fprintf('\n--- old symbols that should be gone ---\n');
nOld = 0;
seen = {};
% E_A0 is both an old key (row "E_A0 -> E_A_ref") and the target name of the
% row above it, so it must not be reported as a leftover.
for k = 1:numel(old)
    if any(strcmp(new, old{k})), skipOld(k) = true; else, skipOld(k) = false; end %#ok<AGROW>
end
for i = 1:numel(sym)
    s = sym{i}{2};
    for k = 1:numel(old)
        if skipOld(k), continue; end
        if strcmp(s, old{k}) || ~isempty(regexp(s, ['\<' old{k} '\>'], 'once'))
            if ~any(strcmp(seen, [old{k} '|' s])), seen{end+1} = [old{k} '|' s]; end %#ok<AGROW>
            nOld = nOld + 1;
        end
    end
end
if nOld == 0
    fprintf('  none  OK\n');
else
    ok = false;
    u = unique(seen);
    for i = 1:numel(u), fprintf('  *** %s\n', u{i}); end
    fprintf('  *** %d occurrence(s) of old symbols -- run rename_symbols\n', nOld);
end

%% ---- 2. which NEW symbols are present ----------------------------------
fprintf('\n--- new symbols present ---\n');
miss = {};
for k = 1:numel(new)
    found = false;
    for i = 1:numel(sym)
        if involves(sym{i}{2}, new{k}), found = true; break; end
    end
    if found
        fprintf('  %-14s (%s)  present\n', new{k}, kind{k});
    else
        miss{end+1} = new{k}; %#ok<AGROW>
    end
end
if ~isempty(miss)
    fprintf('  not referenced by the model (script-only or unused): %s\n', ...
            strjoin(miss, ', '));
end

%% ---- 3. the live labels ------------------------------------------------
fprintf('\n--- live logsout labels (2 s sim) ---\n');
SCENARIO_CODE = 'B'; %#ok<NASGU>
Brain_Parameters;
out = sim(mdl,'StopTime','2','SaveOutput','off','SignalLogging','on', ...
              'SignalLoggingName','logsout','ReturnWorkspaceOutputs','on');
nms = out.logsout.getElementNames(); nms = sort(nms(:));
bad = {};
for i = 1:numel(nms)
    lbl = nms{i};
    for k = 1:numel(old)
        if skipOld(k), continue; end
        if strcmp(lbl, old{k}) || ~isempty(regexp(lbl, ['\<' old{k} '\>'], 'once'))
            bad{end+1} = lbl; %#ok<AGROW>
            break
        end
    end
end
fprintf('  %d labels: %s\n', numel(nms), strjoin(nms', ', '));
if isempty(bad)
    fprintf('  no old symbol among the labels  OK\n');
else
    ok = false;
    fprintf('  *** old symbols still logged: %s\n', strjoin(unique(bad), ', '));
end

fprintf('\n=== check_symbols: %s ===\n', ternary(ok,'PASS','FAIL'));
if ~ok, error('check_symbols:FAIL','symbol set is not consistent with the map'); end
end

% =====================================================================
function tf = involves(store, key)
tf = strcmp(store, key) || ~isempty(regexp(store, ['\<' key '\>'], 'once'));
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end

function [old, new, kind] = local_map(f)
old = {}; new = {}; kind = {};
fid = fopen(f,'r');
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    ln = strtrim(ln);
    if isempty(ln) || ln(1) == '#', continue; end
    c = strsplit(ln, sprintf('\t'));
    if numel(c) < 3, continue; end
    old{end+1} = c{1}; new{end+1} = c{2}; kind{end+1} = c{3}; %#ok<AGROW>
end
fclose(fid);
end
