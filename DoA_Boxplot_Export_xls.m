clc; clear; close all;

%% ============================================================
%  PHASE 1 — EXPORT DoA DATA TO A FLAT EXCEL TABLE
%
%  Output table columns (one row per Date × Mouse):
%    Date | Group | ID |
%    medPD_INC_Push1 ... medPD_INC_Push4  medPD_INC_AvgPush |
%    medPD_HMI_Push1 ... medPD_HMI_Push4  medPD_HMI_AvgPush
%
%  "medPD ratio" = the chosen statChoice value of Peak Displacement
%  (fI=1) for ROI INC_rect (roiIdx=3) and ROI HMI (roiIdx=1).
%  Push3/4 fall back to Push2 when fewer pushes were acquired.
%% ============================================================

%% ---- USER INPUT --------------------------------------------
rootDir   = 'F:\Mouse2';
outXlsx   = 'F:\Mouse2\DoA_FlatTable2.xlsx';
statChoice = 'meanP2P';     % 'medP2P' 'meanP2P' 'iqrP2P' 'stdP2P'

%% ---- Group code map ----------------------------------------
GroupCodeMap = containers.Map('KeyType','double','ValueType','char');
GroupCodeMap(37) = 'Control';
GroupCodeMap(38) = 'Chemotherapy';
GroupCodeMap(39) = 'Immunotherapy';

%% ---- Mouse ID alias map ------------------------------------
IDMap = containers.Map('KeyType','double','ValueType','double');
IDMap(101) = 62;  IDMap(122) = 66;  IDMap(43)  = 53;
IDMap(121) = 63;  IDMap(102) = 61;  IDMap(120) = 125;

%% ---- Push / ROI layout -------------------------------------
%  We export each individual push (with fallback) and the avg push,
%  for both ROIs. fI=1 (Peak Displacement).
pushOptions  = [1 2 3 4 0];          % 0 = average row
pushOptNames = {'Push1','Push2','Push3','Push4','AvgPush'};
fallbackMap  = containers.Map([1 2 3 4 0],[1 2 2 2 0]);
roiList      = [3 1];                % INC_rect=3, HMI=1
roiNames     = {'INC','HMI'};        % matched to roiList order
fI           = 1;                    % Peak Displacement

%% ============================================================
%  DISK SCAN
%% ============================================================
fprintf('=== Scanning disk under %s ===\n', rootDir);

rows = {};   % cell rows, assembled then converted to table

dateDirs = dir(rootDir);
dateDirs = dateDirs([dateDirs.isdir] & ~ismember({dateDirs.name},{'.','..'}));

for di = 1:numel(dateDirs)
    dateStr   = dateDirs(di).name;
    mouseDirs = dir(fullfile(rootDir, dateStr));
    mouseDirs = mouseDirs([mouseDirs.isdir] & ~ismember({mouseDirs.name},{'.','..'}));

    for mi = 1:numel(mouseDirs)
        folderName = mouseDirs(mi).name;            % '37_102_Red_..._D19'
        parts = strsplit(folderName, '_');
        if numel(parts) < 2; continue; end

        groupCode = str2double(parts{1});
        rawID     = str2double(parts{2});
        if isnan(groupCode) || isnan(rawID);  continue; end
        if ~isKey(GroupCodeMap, groupCode);   continue; end

        % Resolve alias → canonical ID
        if isKey(IDMap, rawID); canonID = IDMap(rawID); else; canonID = rawID; end

        doaFile = fullfile(rootDir, dateStr, folderName, 'xCor', 'DoA_ImQ.mat');
        if ~exist(doaFile,'file'); continue; end

        try
            D = load(doaFile);
        catch ME
            warning('Cannot load %s: %s', doaFile, ME.message); continue;
        end
        if ~isfield(D, statChoice)
            warning('%s not in %s', statChoice, doaFile); continue;
        end
        statMat = D.(statChoice);   % (nPush+1) x 2 x 3

        % ---- Build one row ----------------------------------
        rowVals = nan(1, numel(roiList)*numel(pushOptions));
        col = 0;
        for ri = 1:numel(roiList)
            roiIdx = roiList(ri);
            for pi = 1:numel(pushOptions)
                col = col + 1;
                pReq = pushOptions(pi);
                pFal = fallbackMap(pReq);
                rowVals(col) = doaGetVal(statMat, pReq, pFal, fI, roiIdx);
            end
        end

        rows(end+1,:) = [{dateStr, groupCode, canonID}, num2cell(rowVals)]; %#ok<AGROW>
        fprintf('  %s | grp %d | ID %3d | %s\n', dateStr, groupCode, canonID, folderName);
    end
end

if isempty(rows)
    error('No data found under %s', rootDir);
end

%% ---- Build variable names ----------------------------------
varNames = {'Date','Group','ID'};
for ri = 1:numel(roiList)
    for pi = 1:numel(pushOptions)
        varNames{end+1} = sprintf('medPD_%s_%s', roiNames{ri}, pushOptNames{pi}); %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', varNames);

% Sort for readability: Group, ID, Date
T = sortrows(T, {'Group','ID','Date'});

writetable(T, outXlsx, 'Sheet', 'DoA');
fprintf('\nExported %d rows to:\n  %s\n', height(T), outXlsx);
fprintf('Statistic used: %s | fI=%d (Peak Displacement)\n', statChoice, fI);

%% ============================================================
%  LOCAL FUNCTION
%% ============================================================
function val = doaGetVal(medP2P, pushReq, fallback, fI, roiIdx)
    nPush = size(medP2P,1) - 1;
    if pushReq == 0
        pIdx = size(medP2P,1);
    elseif pushReq <= nPush
        pIdx = pushReq;
    else
        pIdx = fallback;
    end
    val = medP2P(pIdx, fI, roiIdx);
end