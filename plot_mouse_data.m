clc;
clear;

% ============================================================
% USER INPUT: ROOT DIRECTORY
% ============================================================
rootDir = 'F:\Mouse2';

% ============================================================
% FIND ALL filtDisp2D DIRECTORIES
% ============================================================
filtDirs = dir(fullfile(rootDir, '**', 'Process_Fun', 'filtDisp2D'));

fprintf('Found %d filtDisp2D directories\n', numel(filtDirs));

% ============================================================
% PREALLOCATE RESULT STRUCT
% ============================================================
Results = struct( ...
    'Date', {}, ...
    'ID', {}, ...
    'Loupass', {}, ...
    'xCor', {} );

row = 0;

% ============================================================
% LOOP OVER EACH filtDisp2D DIRECTORY
% ============================================================
for k = 1:numel(filtDirs)

    filtPath = fullfile(filtDirs(k).folder, filtDirs(k).name);
    fprintf('[%d/%d] Processing: %s\n', k, numel(filtDirs), filtPath);

    % --------------------------------------------------------
    % PARSE DATE AND ID FROM PATH
    % --------------------------------------------------------
    parts = strsplit(filtPath, filesep);

    dateIdx = find(~cellfun(@isempty, regexp(parts, '^\d{4}\.\d{2}\.\d{2}$')), 1);
    if isempty(dateIdx) || dateIdx + 1 > numel(parts)
        warning('Skipping (date/ID not found): %s', filtPath);
        continue;
    end

    DateStr = parts{dateIdx};
    IDStr   = parts{dateIdx + 1};

    % ========================================================
    % LOAD xCor (MANDATORY - ROBUST)
    % ========================================================
    xcorFile = fullfile(filtPath, 'Results_xCor', 'ImQ_P2PD.mat');

    if isfile(xcorFile)

        Sx = load(xcorFile);

    else
        % fallback → search for any xCor-related .mat
        tmp = dir(fullfile(filtPath, '*xCor*.mat'));

        if isempty(tmp)
            warning('Skipping (no xCor file): %s', filtPath);
            continue;
        end

        xcorFile = fullfile(tmp(1).folder, tmp(1).name);
        Sx = load(xcorFile);
    end

    % ---- Extract xCor features safely
    if isfield(Sx,'bMed') && isfield(Sx,'iMed')
        %xcorFeat = Sx.bMed.rect ./ Sx.iMed.rect;
        xcorFeat = Sx.bMed.rectAbs ./ Sx.iMed.rectAbs;
        xcorFeat = xcorFeat(:).';

        % enforce size consistency
        if numel(xcorFeat) ~= 13
            warning('Invalid xCor size (%d): %s', numel(xcorFeat), xcorFile);
            continue;
        end

    else
        warning('Invalid xCor structure: %s', xcorFile);
        continue;
    end

    % ========================================================
    % LOAD Loupass (OPTIONAL)
    % ========================================================
    louFile = fullfile(filtPath, 'Results_Lou', 'ImQ_P2PD.mat');

    if isfile(louFile)

        Sl = load(louFile);

        if isfield(Sl,'bMed') && isfield(Sl,'iMed')
            louFeat = Sl.bMed.rectAbs ./ Sl.iMed.rectAbs;
            louFeat = louFeat(:).';

            % enforce size
            if numel(louFeat) ~= 13
                warning('Invalid Loupass size (%d): %s', numel(louFeat), louFile);
                louFeat = nan(1,13);
            end
        else
            louFeat = nan(1,13);
        end

    else
        % missing Loupass → fill with NaN
        louFeat = nan(1,13);
    end

    % ========================================================
    % STORE RESULTS
    % ========================================================
    row = row + 1;

    Results(row).Date    = DateStr;
    Results(row).ID      = IDStr;
    Results(row).Loupass = louFeat;
    Results(row).xCor    = xcorFeat;

end

% ============================================================
% SAVE OUTPUT
% ============================================================
save('F:\Mouse2\Mouse_Lou_xCor_Features_Abs.mat', 'Results');

% ============================================================
% OPTIONAL: CONVERT TO TABLE (for ML / CSV)
% ============================================================
if ~isempty(Results)
    T = struct2table(Results);
    writetable(T, 'F:\Mouse2\Mouse_Features_Abs.csv');
end

fprintf('Extraction complete. %d experiments saved.\n', row);
%%
% ============================================================
% LOAD EXISTING FILE
% ============================================================
load('F:\Mouse2\Mouse_Lou_xCor_Features_Abs.mat','Results');

% ============================================================
% FUNCTION: Extract numeric ID
% ============================================================
getIDnum = @(str) str2double(regexp(str, '(?<=_)\d+(?=_)', 'match', 'once'));

% ============================================================
% USER INPUT (same as before)
% ============================================================
TimePoints = {
    {'2026.02.16','2026.02.17'}
    {'2026.02.23','2026.02.24'}
    {'2026.03.02','2026.03.03','2026.03.04'}
    {'2026.03.08','2026.03.09','2026.03.10'}
    {'2026.03.11','2026.03.12','2026.03.14'}
    {'2026.03.15','2026.03.16'}
    {'2026.03.17','2026.03.18'}
};

TimeLabels = {'TP1','TP2','TP3','TP4','TP5','TP6','TP7'};

Groups.Control        = [60,56,124,61,54];
Groups.Treatment      = [62,66,125,53,63];
Groups.Immunotherapy  = [57,59,58,67,65];

GroupNames = fieldnames(Groups);

% ---- ID mapping
IDMap = containers.Map('KeyType','double','ValueType','any');
IDMap(62)  = [62 101];
IDMap(66)  = [66 122];
IDMap(53)  = [53 43];
IDMap(63)  = [63 121];
IDMap(61)  = [61 102];
IDMap(125) = [125 120];

% ============================================================
% MODIFY STRUCT (ADD FIELDS)
% ============================================================

for r = 1:numel(Results)

    % ---- Extract numeric ID
    oldID = getIDnum(Results(r).ID);

    % ---- Assign NEW ID
    if isKey(IDMap, oldID)
        newID = IDMap(oldID);
        newID = newID(end);
    else
        newID = oldID;
    end

    % ---- Find TimePoint
    tpLabel = '';
    for t = 1:numel(TimePoints)
        if ismember(Results(r).Date, TimePoints{t})
            tpLabel = TimeLabels{t};
            break;
        end
    end

    % ---- Find Group
    grpName = '';
    for g = 1:numel(GroupNames)
        if ismember(oldID, Groups.(GroupNames{g}))
            grpName = GroupNames{g};
            break;
        end
    end

    % ---- Add new fields
    Results(r).NewID     = newID;
    Results(r).TimePoint = tpLabel;
    Results(r).Group     = grpName;

end

% ============================================================
% SAVE UPDATED FILE
% ============================================================
save('F:\Mouse2\Mouse_Lou_xCor_Features_Abs_WithMeta.mat','Results');

fprintf('✅ Updated MAT file saved with TimePoint, Group, and NewID.\n');
%%
clc; clear;close all


% ============================================================
% LOAD DATA
% ============================================================
load('F:\Mouse2\Mouse_Lou_xCor_Features_Abs_WithMeta.mat','Results');

% ============================================================
% FUNCTION: Extract numeric ID
% ============================================================
getIDnum = @(str) str2double(regexp(str, '(?<=_)\d+(?=_)', 'match', 'once'));

% ============================================================
% USER INPUT
% ============================================================

TimePoints = {
    {'2026.02.16','2026.02.17'},      % TP1
    {'2026.02.23','2026.02.24'},      % TP2
    {'2026.03.02','2026.03.03','2026.03.04'}, % TP3
    {'2026.03.08','2026.03.09','2026.03.10'}, % TP4
    {'2026.03.11','2026.03.12','2026.03.14'}, % TP5
    {'2026.03.15','2026.03.16'},      % TP6
    {'2026.03.17','2026.03.18'}       % TP7
};

TimeLabels = {'TP1','TP2','TP3','TP4','TP5','TP6','TP7'};

Groups.Control        = [60,56,124,61,54];
Groups.Treatment      = [62,66,125,53,63];
Groups.Immunotherapy  = [57,59,58,67,65];

GroupNames = fieldnames(Groups);

% ---- Alias Mapping
IDMap = containers.Map('KeyType','double','ValueType','any');
IDMap(62)  = [62 101];
IDMap(66)  = [66 122];
IDMap(53)  = [53 43];
IDMap(63)  = [63 121];
IDMap(61)  = [61 102];
IDMap(125) = [125 120];

% ============================================================
% CONSTANTS
% ============================================================
nFeat = 13;

FeatureNames = { ...
    '100','200','300','400','500',...
    '600','700','800','900','1000',...
    'InvSlope','PD','PV'};

outDir = 'E:\Mouse2\FinalPlots';
if ~exist(outDir,'dir'); mkdir(outDir); end

colors = [ ...
    0 0.447 0.741;      % Control
    0.85 0.325 0.098;   % Treatment
    0.466 0.674 0.188]; % Immunotherapy

offsets = [-0.25, 0, 0.25];

% ============================================================
% PLOTTING
% ============================================================
figure('Color','w','Position',[50 100 2200 900]);

for feat = 1:nFeat

    subplot(2,7,feat); hold on;

    data = [];
    tp   = {};
    grp  = {};

    % ---------------- GROUP LOOP ----------------
    for g = 1:numel(GroupNames)

        groupName = GroupNames{g};
        IDs = Groups.(groupName);

        for id = 1:numel(IDs)

            currID = IDs(id);

            % ---- Alias handling
            if isKey(IDMap, currID)
                currAliases = IDMap(currID);
            else
                currAliases = currID;
            end

            % ---------- BASELINE (TP1) ----------
            baseVals = [];

            for r = 1:numel(Results)
                if ismember(Results(r).Date, TimePoints{1}) && ...
                   ismember(getIDnum(Results(r).ID), currAliases)

                    baseVals(end+1) = Results(r).xCor(feat);
                end
            end

            if isempty(baseVals), continue; end
            baseVal = mean(baseVals);

            % ---------- TIMEPOINT LOOP ----------
            for t = 1:numel(TimePoints)

                vals = [];

                for r = 1:numel(Results)
                    if ismember(Results(r).Date, TimePoints{t}) && ...
                       ismember(getIDnum(Results(r).ID), currAliases)

                        vals(end+1) = Results(r).xCor(feat);
                    end
                end

                if isempty(vals), continue; end
                val = mean(vals);

                % ✅ CORRECT Δ (TP1 = 0)
                delta = (val);% - baseVal);%*100/baseVal;

                data(end+1) = delta;
                tp{end+1}   = TimeLabels{t};
                grp{end+1}  = groupName;

            end
        end
    end

    % ---------------- BOX PLOTS (SIDE-BY-SIDE) ----------------
    cats = categorical(tp, TimeLabels, 'Ordinal', true);
    x = double(cats);

    for g = 1:numel(GroupNames)

        groupName = GroupNames{g};
        idx = strcmp(grp, groupName);

        if any(idx)
            boxchart(x(idx) + offsets(g), data(idx), ...
                'BoxWidth', 0.2, ...
                'BoxFaceColor', colors(g,:), ...
                'WhiskerLineColor', colors(g,:), ...
                'MarkerStyle','.');
        end
    end

    title(FeatureNames{feat});
    grid on;
    xtickangle(45);

    if feat == 1
        ylabel('TPx');
        %ylabel('\Delta xCor (TP - TP1)');
    end
end

% ============================================================
% LEGEND
% ============================================================
hold on;
h(1) = plot(nan,nan,'s','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
h(2) = plot(nan,nan,'s','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
h(3) = plot(nan,nan,'s','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');

legend(h, {'Control','Treatment','Immunotherapy'}, ...
    'Position',[0.92 0.4 0.05 0.1]);

% ============================================================
% SAVE
% ============================================================
saveas(gcf, fullfile(outDir,'xCor_AllFeatures_DeltaTP1.png'));

fprintf('✅ Done: TP1-normalized grouped boxplots generated.\n');