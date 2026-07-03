clc; clear; close all;

%% ============================================================
%  USER INPUT — edit this section only
%% ============================================================

rootDir = 'F:\Mouse2';
outDir  = 'F:\Mouse2\Boxplots\DoA_PD_Boxplots';

% ---- Statistic to use from DoA_ImQ.mat ----------------------
%  Choose ONE of: 'medP2P'  'meanP2P'  'iqrP2P'  'stdP2P'
statChoice = 'medP2P';

%% ---- Time point definitions --------------------------------
TimePoints = {
    {'2026.02.16','2026.02.17'}                                              % TP1
    {'2026.02.23','2026.02.24'}                                              % TP2
    {'2026.03.02','2026.03.03','2026.03.04'}                                 % TP4
    {'2026.03.11','2026.03.12','2026.03.13','2026.03.14'}
    {'2026.03.15','2026.03.16'}
    {'2026.03.17','2026.03.18','2026.03.19','2026.03.22'} % TP6
};
TimeLabels = {'TP1','TP2','TP3','TP4','TP5','TP6'};
DayLabels  = [4 12 19 24 28 30];

%% ---- Group code map ----------------------------------------
GroupCodeMap = containers.Map('KeyType','double','ValueType','char');
GroupCodeMap(37) = 'Control';
GroupCodeMap(38) = 'Chemotherapy';
GroupCodeMap(39) = 'Immunotherapy';

GroupNames = {'Control','Chemotherapy','Immunotherapy'};
colors  = [0.000 0.447 0.741;
           0.850 0.325 0.098;
           0.466 0.674 0.188];
offsets = [-0.25, 0, 0.25];

%% ---- Mouse ID alias map ------------------------------------
%  IDMap(aliasID) = canonicalID
%  Both "37_aliasID_..." and "37_canonicalID_..." → same mouse
IDMap = containers.Map('KeyType','double','ValueType','double');
IDMap(101) = 62;
IDMap(122) = 66;
IDMap(43)  = 53;
IDMap(121) = 63;
IDMap(102) = 61;
IDMap(120) = 125;

%% ============================================================
%  FEATURE DEFINITIONS
%
%  medP2P: (nPush+1) x 2 x 3
%    Dim1: push index — 1..nPush = individual pushes, end = average
%    Dim2: 1=Peak(PD), 2=PeakV(PV)
%    Dim3: 1=roiHMI,   2=roiBKD_rect,  3=roiINC_rect
%
%  pushReq=0 → use end (average row)
%  Push3/4 fall back to Push2 when fewer pushes were acquired
%% ============================================================
%          Label                  pushReq  fallback  fI  roiIdx
FeatDef = {
    'Push1_PD_HMI',               1,       1,        1,  1;
    'Push2_PD_HMI',               2,       2,        1,  1;
    'Push3_PD_HMI',               3,       2,        1,  1;
    'Push4_PD_HMI',               4,       2,        1,  1;
    'AvgPush_PD_HMI',             0,       0,        1,  1;
    'Push1_PD_INC_rect',          1,       1,        1,  3;
    'Push2_PD_INC_rect',          2,       2,        1,  3;
    'Push3_PD_INC_rect',          3,       2,        1,  3;
    'Push4_PD_INC_rect',          4,       2,        1,  3;
    'AvgPush_PD_INC_rect',        0,       0,        1,  3;
};

nFeat    = size(FeatDef,1);
featName = FeatDef(:,1);
featPush = cell2mat(FeatDef(:,2));
featFall = cell2mat(FeatDef(:,3));
featFI   = cell2mat(FeatDef(:,4));
featROI  = cell2mat(FeatDef(:,5));

%% ============================================================
%  OUTPUT DIRECTORY
%% ============================================================
if ~exist(outDir,'dir'); mkdir(outDir); end

%% ============================================================
%  DISK SCAN
%  Path: rootDir\<date>\<mouseFolder>\xCor\DoA_ImQ.mat
%  Parse: parts{1}=groupCode, parts{2}=mouseID  (split on '_')
%% ============================================================
fprintf('=== Scanning disk under %s ===\n', rootDir);

diskRegistry = containers.Map('KeyType','double','ValueType','any');

dateDirs = dir(rootDir);
dateDirs = dateDirs([dateDirs.isdir] & ~ismember({dateDirs.name},{'.','..',}));

for di = 1:numel(dateDirs)
    dateStr   = dateDirs(di).name;
    mouseDirs = dir(fullfile(rootDir, dateStr));
    mouseDirs = mouseDirs([mouseDirs.isdir] & ~ismember({mouseDirs.name},{'.','..',}));

    for mi = 1:numel(mouseDirs)
        folderName = mouseDirs(mi).name;            % e.g. '37_102_Red_2p3mm_D19'

        parts = strsplit(folderName, '_');
        if numel(parts) < 2; continue; end

        groupCode = str2double(parts{1});           % 37 / 38 / 39
        rawID     = str2double(parts{2});           % e.g. 102

        if isnan(groupCode) || isnan(rawID); continue; end
        if ~isKey(GroupCodeMap, groupCode);  continue; end

        groupName = GroupCodeMap(groupCode);

        % Resolve alias → canonical ID
        if isKey(IDMap, rawID)
            canonID = IDMap(rawID);
        else
            canonID = rawID;
        end

        % Only register if DoA_ImQ.mat actually exists
        doaFile = fullfile(rootDir, dateStr, folderName, 'xCor', 'DoA_ImQ.mat');
        if ~exist(doaFile,'file'); continue; end

        % First time seeing this canonical mouse → initialise entry
        if ~isKey(diskRegistry, canonID)
            entry.group     = groupName;
            entry.dateFiles = containers.Map('KeyType','char','ValueType','char');
            diskRegistry(canonID) = entry;
        end

        entry = diskRegistry(canonID);

        if ~strcmp(entry.group, groupName)
            warning('Mouse %d: group mismatch (%s vs %s) on %s — keeping %s.', ...
                canonID, entry.group, groupName, dateStr, entry.group);
        end

        if isKey(entry.dateFiles, dateStr)
            warning('Mouse %d, %s: duplicate folder — keeping first.', canonID, dateStr);
        else
            entry.dateFiles(dateStr) = doaFile;
            fprintf('  Mouse %3d | %-15s | %s | %s\n', ...
                canonID, groupName, dateStr, folderName);
        end

        diskRegistry(canonID) = entry;
    end
end

allCanonIDs = cell2mat(keys(diskRegistry));
fprintf('\nTotal unique mice registered: %d\n\n', numel(allCanonIDs));

%% ============================================================
%  DATA COLLECTION
%  mouse → TP → dates → load medP2P → extract 10 features
%  Multiple dates in same TP are averaged into one value.
%% ============================================================
fprintf('=== Collecting data ===\n');

featData = cell(nFeat,1);
featTP   = cell(nFeat,1);
featGrp  = cell(nFeat,1);
for f = 1:nFeat
    featData{f} = [];
    featTP{f}   = {};
    featGrp{f}  = {};
end

for ii = 1:numel(allCanonIDs)
    canonID = allCanonIDs(ii);
    entry   = diskRegistry(canonID);
    grpName = entry.group;

    for t = 1:numel(TimePoints)
        tpDates = TimePoints{t};
        tpVals  = nan(nFeat, numel(tpDates));   % features × dates

        for d = 1:numel(tpDates)
            currDate = tpDates{d};
            if ~isKey(entry.dateFiles, currDate); continue; end

            doaFile = entry.dateFiles(currDate);
            try
                D = load(doaFile);              % load ALL variables
            catch ME
                warning('Cannot load %s: %s', doaFile, ME.message);
                continue;
            end

            if ~isfield(D, statChoice)
                warning('%s not found in %s', statChoice, doaFile);
                continue;
            end

            statMat = D.(statChoice);   % e.g. D.medP2P or D.meanP2P

            fprintf('  Mouse %3d | %s | %s  [size %s: %s]\n', ...
                canonID, TimeLabels{t}, currDate, ...
                statChoice, num2str(size(statMat)));

            for f = 1:nFeat
                tpVals(f,d) = doaGetVal(statMat, ...
                                        featPush(f), featFall(f), ...
                                        featFI(f),   featROI(f));
            end
        end % date loop

        tpMean = mean(tpVals, 2, 'omitnan');    % nFeat × 1
        if all(isnan(tpMean)); continue; end

        for f = 1:nFeat
            if ~isnan(tpMean(f))
                featData{f}(end+1) = tpMean(f);      %#ok<AGROW>
                featTP{f}{end+1}   = TimeLabels{t};  %#ok<AGROW>
                featGrp{f}{end+1}  = grpName;        %#ok<AGROW>
            end
        end
    end % TP loop
end % mouse loop

fprintf('\nData collection complete.\n\n');

%% ============================================================
%  PLOTTING
%% ============================================================
fprintf('=== Generating plots ===\n');

roiLabels = {'roiHMI (tumor contour)', 'roiBKD rect', 'roiINC rect'};
figAll    = figure('Color','w','Position',[50 50 2500 950]);

for f = 1:nFeat

    if isempty(featData{f})
        fprintf('  No data for: %s — skipping.\n', featName{f});
        continue;
    end

    dat  = featData{f};
    grp  = featGrp{f};
    cats = categorical(featTP{f}, TimeLabels, 'Ordinal', true);
    x    = double(cats);

    %% ---- Subplot in combined figure ------------------------
    figure(figAll);
    ax1 = subplot(2,5,f);
    hold(ax1,'on');

    for gg = 1:numel(GroupNames)
        idxG = strcmp(grp, GroupNames{gg});
        if any(idxG)
            boxchart(ax1, x(idxG)+offsets(gg), dat(idxG), ...
                'BoxWidth',         0.2, ...
                'BoxFaceColor',     colors(gg,:), ...
                'WhiskerLineColor', colors(gg,:), ...
                'MarkerStyle',      '.');
        end
    end

    title(ax1, strrep(featName{f},'_',' '), 'FontSize',11,'FontWeight','bold');
    grid(ax1,'on'); box(ax1,'on');
    xlim(ax1,[0.5 numel(TimeLabels)+0.5]);
    xticks(ax1, 1:numel(DayLabels));
    xticklabels(ax1, string(DayLabels));
    ax1.FontSize=10; ax1.FontWeight='bold';
    ax1.XColor='k';  ax1.YColor='k'; ax1.LineWidth=1.2;
    if mod(f,5)==1
        ylabel(ax1,[statChoice ' DoA – PD'],'FontSize',11,'FontWeight','bold');
    end

    %% ---- Individual figure ---------------------------------
    figSingle = figure('Color','w','Position',[200 150 900 650]);
    ax2 = axes('Parent',figSingle);
    hold(ax2,'on');

    for gg = 1:numel(GroupNames)
        idxG = strcmp(grp, GroupNames{gg});
        if any(idxG)
            boxchart(ax2, x(idxG)+offsets(gg), dat(idxG), ...
                'BoxWidth',         0.2, ...
                'BoxFaceColor',     colors(gg,:), ...
                'WhiskerLineColor', colors(gg,:), ...
                'MarkerStyle',      '.');
        end
    end

    if featPush(f)==0
        pushStr = 'Average Push';
    elseif featPush(f)==featFall(f)
        pushStr = sprintf('Push %d', featPush(f));
    else
        pushStr = sprintf('Push %d  (fallback → Push %d)', featPush(f), featFall(f));
    end
    title(ax2, sprintf('DoA – PD  |  %s  |  %s', pushStr, roiLabels{featROI(f)}), ...
        'FontSize',14,'FontWeight','bold');

    box(ax2,'on'); grid(ax2,'off');
    xlim(ax2,[0.5 numel(TimeLabels)+0.5]);
    xticks(ax2, 1:numel(DayLabels));
    xticklabels(ax2, string(DayLabels));
    ax2.FontSize=14; ax2.FontWeight='bold';
    ax2.XColor='k';  ax2.YColor='k'; ax2.LineWidth=1.5;

    xlabel(ax2,'Days after 4T1 cell implantation','FontSize',15,'FontWeight','bold');
    ylabel(ax2,[statChoice ' DoA – Peak Displacement'],'FontSize',15,'FontWeight','bold');

    hL(1)=plot(ax2,nan,nan,'s','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none','MarkerSize',10);
    hL(2)=plot(ax2,nan,nan,'s','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none','MarkerSize',10);
    hL(3)=plot(ax2,nan,nan,'s','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none','MarkerSize',10);
    lgd=legend(ax2, hL, GroupNames, 'Location','best');
    lgd.FontSize=13; lgd.FontWeight='bold';

    saveas(figSingle, fullfile(outDir,[featName{f} '.png']));
    saveas(figSingle, fullfile(outDir,[featName{f} '.fig']));
    fprintf('  Saved: %s\n', featName{f});

end % feature loop

%% ---- Combined figure ---------------------------------------
figure(figAll);
hL(1)=plot(nan,nan,'s','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none','MarkerSize',10);
hL(2)=plot(nan,nan,'s','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none','MarkerSize',10);
hL(3)=plot(nan,nan,'s','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none','MarkerSize',10);
lgd=legend(hL, GroupNames, 'Position',[0.93 0.44 0.05 0.12]);
lgd.FontSize=13; lgd.FontWeight='bold';
sgtitle(['DoA – Peak Displacement  (All Features) — ' statChoice],'FontSize',17,'FontWeight','bold');

saveas(figAll, fullfile(outDir,'DoA_PD_AllFeatures_Combined.png'));
saveas(figAll, fullfile(outDir,'DoA_PD_AllFeatures_Combined.fig'));

fprintf('\nDone. All plots saved in:\n  %s\n', outDir);
%% ============================================================
%  DIFFERENCE PLOTS  (TPX – TP1)
%% ============================================================
fprintf('=== Generating difference plots (TPX – TP1) ===\n');

outDirDiff = fullfile(outDir, 'Difference_vs_TP1');
if ~exist(outDirDiff,'dir'); mkdir(outDirDiff); end

figAllDiff = figure('Color','w','Position',[50 50 2500 950]);

for f = 1:nFeat

    if isempty(featData{f})
        fprintf('  No data for: %s — skipping.\n', featName{f});
        continue;
    end

    dat  = featData{f};
    grp  = featGrp{f};
    cats = categorical(featTP{f}, TimeLabels, 'Ordinal', true);
    x    = double(cats);

    %% ---- Compute per-mouse TP1 baseline and subtract ----------
    % We need to re-derive per-mouse values. Instead, work from the
    % already-collected flat arrays: find each unique mouse's TP1
    % value per group and subtract.
    %
    % Strategy: iterate over GroupNames, then per mouse (identified
    % by unique (grp, x) combinations that share TP1).
    % Since featData is flat (no mouse ID), we rebuild from diskRegistry.

    datDiff  = [];
    grpDiff  = {};
    xDiff    = [];

    for gg = 1:numel(GroupNames)
        gName = GroupNames{gg};

        % Find all mice in this group
        for ii = 1:numel(allCanonIDs)
            canonID = allCanonIDs(ii);
            entry   = diskRegistry(canonID);
            if ~strcmp(entry.group, gName); continue; end

            % Collect per-TP values for this mouse
            mouseVals = nan(1, numel(TimePoints));
            for t = 1:numel(TimePoints)
                tpDates = TimePoints{t};
                tpV     = nan(1, numel(tpDates));

                for d = 1:numel(tpDates)
                    currDate = tpDates{d};
                    if ~isKey(entry.dateFiles, currDate); continue; end
                    doaFile = entry.dateFiles(currDate);
                    try
                        D = load(doaFile);
                    catch
                        continue;
                    end
                    if ~isfield(D, statChoice); continue; end
                    tpV(d) = doaGetVal(D.(statChoice), ...
                                       featPush(f), featFall(f), ...
                                       featFI(f),   featROI(f));
                end
                mouseVals(t) = mean(tpV, 'omitnan');
            end

            % Need TP1 baseline
            baseline = mouseVals(1);
            if isnan(baseline); continue; end

            for t = 1:numel(TimePoints)
                if isnan(mouseVals(t)); continue; end
                datDiff(end+1) = (mouseVals(t) - baseline)/baseline; %#ok<AGROW>
                grpDiff{end+1} = gName;                   %#ok<AGROW>
                xDiff(end+1)   = t;                       %#ok<AGROW>
            end
        end
    end

    if isempty(datDiff)
        fprintf('  No difference data for: %s — skipping.\n', featName{f});
        continue;
    end

    %% ---- Subplot in combined diff figure --------------------
    figure(figAllDiff);
    ax1 = subplot(2,5,f);
    hold(ax1,'on');
    yline(ax1, 0, 'k--', 'LineWidth', 1.2);

    for gg = 1:numel(GroupNames)
        idxG = strcmp(grpDiff, GroupNames{gg});
        if any(idxG)
            boxchart(ax1, xDiff(idxG)+offsets(gg), datDiff(idxG), ...
                'BoxWidth',         0.2, ...
                'BoxFaceColor',     colors(gg,:), ...
                'WhiskerLineColor', colors(gg,:), ...
                'MarkerStyle',      '.');
        end
    end

    title(ax1, strrep(featName{f},'_',' '), 'FontSize',11,'FontWeight','bold');
    grid(ax1,'on'); box(ax1,'on');
    xlim(ax1,[0.5 numel(TimeLabels)+0.5]);
    xticks(ax1, 1:numel(DayLabels));
    xticklabels(ax1, string(DayLabels));
    ax1.FontSize=10; ax1.FontWeight='bold';
    ax1.XColor='k';  ax1.YColor='k'; ax1.LineWidth=1.2;
    if mod(f,5)==1
        ylabel(ax1,['\Delta' statChoice ' (TPX – TP1)'],'FontSize',11,'FontWeight','bold');
    end

    %% ---- Individual diff figure ----------------------------
    figSingle = figure('Color','w','Position',[200 150 900 650]);
    ax2 = axes('Parent', figSingle);
    hold(ax2,'on');
    yline(ax2, 0, 'k--', 'LineWidth', 1.5, 'Label','Baseline (TP1)','FontSize',12);

    for gg = 1:numel(GroupNames)
        idxG = strcmp(grpDiff, GroupNames{gg});
        if any(idxG)
            boxchart(ax2, xDiff(idxG)+offsets(gg), datDiff(idxG), ...
                'BoxWidth',         0.2, ...
                'BoxFaceColor',     colors(gg,:), ...
                'WhiskerLineColor', colors(gg,:), ...
                'MarkerStyle',      '.');
        end
    end

    if featPush(f)==0
        pushStr = 'Average Push';
    elseif featPush(f)==featFall(f)
        pushStr = sprintf('Push %d', featPush(f));
    else
        pushStr = sprintf('Push %d  (fallback → Push %d)', featPush(f), featFall(f));
    end
    title(ax2, sprintf('DoA – PD  \x0394(TPX–TP1)  |  %s  |  %s', ...
        pushStr, roiLabels{featROI(f)}), ...
        'FontSize',14,'FontWeight','bold');

    box(ax2,'on'); grid(ax2,'off');
    xlim(ax2,[0.5 numel(TimeLabels)+0.5]);
    xticks(ax2, 1:numel(DayLabels));
    xticklabels(ax2, string(DayLabels));
    ax2.FontSize=14; ax2.FontWeight='bold';
    ax2.XColor='k';  ax2.YColor='k'; ax2.LineWidth=1.5;

    xlabel(ax2,'Days after 4T1 cell implantation','FontSize',15,'FontWeight','bold');
    ylabel(ax2,['\Delta ' statChoice ' – Peak Displacement (TPX – TP1)'], ...
        'FontSize',15,'FontWeight','bold');

    hL(1)=plot(ax2,nan,nan,'s','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none','MarkerSize',10);
    hL(2)=plot(ax2,nan,nan,'s','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none','MarkerSize',10);
    hL(3)=plot(ax2,nan,nan,'s','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none','MarkerSize',10);
    lgd = legend(ax2, hL, GroupNames, 'Location','best');
    lgd.FontSize=13; lgd.FontWeight='bold';

    saveas(figSingle, fullfile(outDirDiff, ['Diff_' featName{f} '.png']));
    saveas(figSingle, fullfile(outDirDiff, ['Diff_' featName{f} '.fig']));
    close(figSingle);
    fprintf('  Saved diff: %s\n', featName{f});

end % feature loop (diff)

%% ---- Combined diff figure legend & save --------------------
figure(figAllDiff);
hL(1)=plot(nan,nan,'s','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none','MarkerSize',10);
hL(2)=plot(nan,nan,'s','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none','MarkerSize',10);
hL(3)=plot(nan,nan,'s','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none','MarkerSize',10);
lgd=legend(hL, GroupNames, 'Position',[0.93 0.44 0.05 0.12]);
lgd.FontSize=13; lgd.FontWeight='bold';
sgtitle(['\DeltaDoA – Peak Displacement  (TPX – TP1) — ' statChoice], ...
    'FontSize',17,'FontWeight','bold');

saveas(figAllDiff, fullfile(outDirDiff,'Diff_DoA_PD_AllFeatures_Combined.png'));
saveas(figAllDiff, fullfile(outDirDiff,'Diff_DoA_PD_AllFeatures_Combined.fig'));

fprintf('\nDone. Difference plots saved in:\n  %s\n', outDirDiff);

%% ============================================================
%  PER-MOUSE TRAJECTORY PLOTS  (one figure per group)
%  Rows: HMI | INC_rect   Cols: Push1 | Push2 | Push3 | Push4 | AvgPush
%  Each line = one mouse, colored uniquely, X marks missing pushes
%% ============================================================
fprintf('=== Generating per-mouse trajectory plots ===\n');

outDirMouse = fullfile(outDir, 'PerMouse_Trajectories');
if ~exist(outDirMouse,'dir'); mkdir(outDirMouse); end

% ---- Feature layout for this figure -------------------------
% Row 1: HMI  (featROI==1), Row 2: INC_rect (featROI==3)
% Cols:  Push1 Push2 Push3 Push4 AvgPush
pushCols  = [1 2 3 4 0];          % pushReq values
pushLabels = {'Push 1','Push 2','Push 3','Push 4','Avg Push'};
roiRows    = [1 3];               % roiIdx values
roiRowLabels = {'HMI (tumor contour)', 'INC rect'};
nCols = numel(pushCols);
nRows = numel(roiRows);

% Fallback rule: Push3/4 → Push2 if unavailable (same as FeatDef)
fallbackMap = containers.Map([1 2 3 4 0],[1 2 2 2 0]);

% Colour palette — up to 15 mice, cycling if more
basePalette = lines(15);

xTick_pos  = 1:numel(TimeLabels);
xTick_days = DayLabels;

for gg = 1:numel(GroupNames)
    gName = GroupNames{gg};

    % Collect mice in this group (sorted for reproducibility)
    grpMouseIDs = sort( allCanonIDs( arrayfun(@(id) ...
        strcmp(diskRegistry(id).group, gName), allCanonIDs) ) );

    nMice = numel(grpMouseIDs);
    if nMice == 0
        fprintf('  No mice in group %s — skipping.\n', gName);
        continue;
    end

    % Assign a unique colour per mouse
    mouseColors = basePalette(mod((1:nMice)-1, size(basePalette,1))+1, :);

    figG = figure('Color','w', ...
        'Position',[30 30 320*nCols+120  320*nRows+100], ...
        'Name', sprintf('Per-Mouse Trajectories — %s', gName));

    axG = gobjects(nRows, nCols);   % pre-allocate axes handles

    for rr = 1:nRows
        roiIdx = roiRows(rr);

        for cc = 1:nCols
            pushReq  = pushCols(cc);
            fallback = fallbackMap(pushReq);

            ax = subplot(nRows, nCols, (rr-1)*nCols + cc);
            axG(rr,cc) = ax;
            hold(ax,'on');

            legendEntries  = {};
            legendHandles  = gobjects(0);
            missingMarkers = gobjects(0);   % collect 'X' handles for legend

            for mm = 1:nMice
                canonID = grpMouseIDs(mm);
                entry   = diskRegistry(canonID);
                col     = mouseColors(mm,:);

                yVals      = nan(1, numel(TimePoints));
                isFallback = false(1, numel(TimePoints));

                for t = 1:numel(TimePoints)
                    tpDates = TimePoints{t};
                    tpV     = nan(1, numel(tpDates));

                    for d = 1:numel(tpDates)
                        currDate = tpDates{d};
                        if ~isKey(entry.dateFiles, currDate); continue; end
                        doaFile = entry.dateFiles(currDate);
                        try
                            D = load(doaFile);
                        catch
                            continue;
                        end
                        if ~isfield(D, statChoice); continue; end

                        statMat = D.(statChoice);
                        nPush   = size(statMat,1) - 1;

                        % Detect fallback use
                        if pushReq > 0 && pushReq > nPush
                            isFallback(t) = true;
                        end

                        tpV(d) = doaGetVal(statMat, pushReq, fallback, 1, roiIdx);
                    end
                    yVals(t) = mean(tpV, 'omitnan');
                end

                % Skip mouse entirely if no valid data
                if all(isnan(yVals)); continue; end

                % ---- Plot connecting line (gaps where NaN) ----
                hLine = plot(ax, xTick_pos, yVals, '-o', ...
                    'Color',           col, ...
                    'LineWidth',       1.8, ...
                    'MarkerSize',      6, ...
                    'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', 'none', ...
                    'DisplayName',     sprintf('Mouse %d', canonID));

                legendHandles(end+1) = hLine; %#ok<AGROW>
                legendEntries{end+1} = sprintf('Mouse %d', canonID); %#ok<AGROW>

                % ---- Mark fallback TPs with X ----------------
                fbIdx = find(isFallback & ~isnan(yVals));
                if ~isempty(fbIdx)
                    hX = plot(ax, xTick_pos(fbIdx), yVals(fbIdx), 'x', ...
                        'Color',      col, ...
                        'LineWidth',  2.5, ...
                        'MarkerSize', 14, ...
                        'HandleVisibility','off');   % don't clutter legend
                end

                % ---- Mark missing TPs with X in black --------
                missIdx = find(isnan(yVals));
                if ~isempty(missIdx)
                    % Place X at bottom of current y-axis (will rescale)
                    % Store x positions; we'll draw after ylim is set
                    % Use a temporary NaN-safe y = 0 placeholder
                    for mk = 1:numel(missIdx)
                        missingMarkers(end+1) = plot(ax, ...
                            xTick_pos(missIdx(mk)), NaN, 'x', ...
                            'Color',     col*0.5, ...   % darker shade
                            'LineWidth', 2, ...
                            'MarkerSize',12, ...
                            'HandleVisibility','off');   %#ok<AGROW>
                    end
                end

            end % mouse loop

            % ---- Axes formatting ----------------------------
            xlim(ax, [0.5, numel(TimeLabels)+0.5]);
            xticks(ax, xTick_pos);
            xticklabels(ax, string(xTick_days));
            ax.FontSize   = 11;
            ax.FontWeight = 'bold';
            ax.XColor = 'k'; ax.YColor = 'k';
            ax.LineWidth  = 1.3;
            grid(ax,'on'); box(ax,'on');

            % Move missing-data X markers to bottom of y-axis
            yl = ylim(ax);
            yBot = yl(1) - 0.08*(yl(2)-yl(1));   % just below axis bottom
            for hh = 1:numel(missingMarkers)
                if isvalid(missingMarkers(hh))
                    missingMarkers(hh).YData = yBot;
                end
            end

            % Column title (Push label) on top row only
            if rr == 1
                title(ax, pushLabels{cc}, 'FontSize',13,'FontWeight','bold');
            end

            % Row label (ROI) on left column only
            if cc == 1
                ylabel(ax, {roiRowLabels{rr}; [statChoice ' DoA–PD']}, ...
                    'FontSize',12,'FontWeight','bold');
            end

            % X-axis label on bottom row only
            if rr == nRows
                xlabel(ax,'Days after 4T1 implantation','FontSize',11,'FontWeight','bold');
            end

            % Legend only on the last column to save space
            if cc == nCols && ~isempty(legendHandles)
                lgd = legend(ax, legendHandles, legendEntries, ...
                    'Location','eastoutside','FontSize',10,'FontWeight','bold');
            end

        end % col loop
    end % row loop

    % ---- Fallback note & supertitle -------------------------
    % Add a small annotation explaining the X markers
    annotation(figG,'textbox',[0 0 1 0.03], ...
        'String', ...
        ['Coloured  ×  = push unavailable, fallback used   |   ' ...
         'Dark  ×  at axis bottom = no data for that TP'], ...
        'EdgeColor','none','HorizontalAlignment','center', ...
        'FontSize',10,'FontAngle','italic');

    sgtitle(figG, sprintf('Per-Mouse DoA–PD Trajectories  |  %s  |  %s', ...
        gName, statChoice), 'FontSize',16,'FontWeight','bold');

    % ---- Save -----------------------------------------------
    saveas(figG, fullfile(outDirMouse, sprintf('PerMouse_%s.png', gName)));
    saveas(figG, fullfile(outDirMouse, sprintf('PerMouse_%s.fig', gName)));
    fprintf('  Saved: PerMouse_%s\n', gName);

end % group loop

fprintf('\nPer-mouse trajectory plots saved in:\n  %s\n', outDirMouse);

%% ============================================================
%  PUSH OPTIMIZER  — Q3-based scoring
%  (replace medP2P median → 75th percentile for trend matching)
%% ============================================================

outDirOpt = fullfile(outDir, 'OptimizedPush_Q3');
if ~exist(outDirOpt,'dir'); mkdir(outDirOpt); end

%% ---- Config ------------------------------------------------
statChoices  = {'medP2P','meanP2P'};
roiList      = [1, 3];
roiListNames = {'HMI','INC_rect'};
pushOptions  = [1 2 3 4 0];
pushOptNames = {'Push1','Push2','Push3','Push4','AvgPush'};
nTP          = numel(TimeLabels);
nPushOpt     = numel(pushOptions);
fallbackMap  = containers.Map([1 2 3 4 0],[1 2 2 2 0]);

% Later TPs weighted more heavily
monoWeights = 1:nTP-1;   % [1 2 3 4 5]
sepWeights  = 1:nTP;     % [1 2 3 4 5 6]

lambda_mono = 0.6;
lambda_sep  = 0.4;

sepDirection = {'Control','Chemotherapy','Immunotherapy'};


%% ---- Build all 5^6 combos ----------------------------------
fprintf('Building push combination table (Q3 scoring)...\n');
[g1,g2,g3,g4,g5,g6] = ndgrid(1:nPushOpt,1:nPushOpt,1:nPushOpt,...
                               1:nPushOpt,1:nPushOpt,1:nPushOpt);
allCombos = [g1(:) g2(:) g3(:) g4(:) g5(:) g6(:)];
nCombos   = size(allCombos,1);
fprintf('Total combinations: %d\n\n', nCombos);

%% ============================================================
%  MAIN OPTIMIZATION LOOP
%% ============================================================
OptResults = struct();

for si = 1:numel(statChoices)
    sc = statChoices{si};

    for ri = 1:numel(roiList)
        roiIdx = roiList(ri);
        roiNm  = roiListNames{ri};

        fprintf('=== %s | %s ===\n', sc, roiNm);

        %% -- Preload all mouse data ---------------------------
        allMouseData  = cell(numel(GroupNames),1);
        allMouseIDsG  = cell(numel(GroupNames),1);

        for gi = 1:numel(GroupNames)
            gName   = GroupNames{gi};
            grpMIDs = sort(allCanonIDs(arrayfun(@(id) ...
                strcmp(diskRegistry(id).group,gName), allCanonIDs)));
            nM = numel(grpMIDs);
            allMouseIDsG{gi}  = grpMIDs;
            allMouseData{gi}  = nan(nM, nTP, nPushOpt);

            for pi = 1:nPushOpt
                pReq = pushOptions(pi);
                pFal = fallbackMap(pReq);
                [mVals,~] = collectMouseTP(grpMIDs, diskRegistry, ...
                    TimePoints, sc, pReq, pFal, 1, roiIdx);
                allMouseData{gi}(:,:,pi) = mVals;
            end
            fprintf('  Group %-15s: %d mice\n', gName, nM);
        end

        %% -- Precompute Q3 for all groups, all pushes, all TPs
        %  grpQ3All{g}(pi, t) = Q3 of group g, push option pi, TP t
        grpQ3All = cell(numel(GroupNames),1);
        for gi = 1:numel(GroupNames)
            grpQ3All{gi} = nan(nPushOpt, nTP);
            for pi = 1:nPushOpt
                for t = 1:nTP
                    vals = allMouseData{gi}(:,t,pi);
                    vals = vals(~isnan(vals));
                    if ~isempty(vals)
                        grpQ3All{gi}(pi,t) = prctile(vals, 75);
                    end
                end
            end
        end

        %% -- Per-group optimization ---------------------------
        for gi = 1:numel(GroupNames)
            gName = GroupNames{gi};
            fprintf('  Optimizing: %s | %s | %s\n', sc, roiNm, gName);

            bestScore     = -inf;
            bestCombo     = ones(1,nTP);
            bestIsMono    = false;
            bestIsSep     = false;
            bestMonoScore = 0;
            bestSepScore  = 0;
            foundStrict   = false;

            % ---- Strict pass ---------------------------------
            for ci = 1:nCombos
                pVec = allCombos(ci,:);

                % Build Q3 trajectory for all groups using this combo
                grpQ3TP = cell(numel(GroupNames),1);
                for g = 1:numel(GroupNames)
                    grpQ3TP{g} = arrayfun(@(t) grpQ3All{g}(pVec(t),t), 1:nTP);
                end

                [sc_val, isMono, isSep, mSc, sSc] = scorePushComboQ3( ...
                    grpQ3TP, monoWeights, sepWeights, ...
                    lambda_mono, lambda_sep, nTP, GroupNames, sepDirection);

                if isMono && sc_val > bestScore
                    bestScore     = sc_val;
                    bestCombo     = pVec;
                    bestIsMono    = isMono;
                    bestIsSep     = isSep;
                    bestMonoScore = mSc;
                    bestSepScore  = sSc;
                    foundStrict   = true;
                end
            end

            % ---- Soft fallback -------------------------------
            if ~foundStrict
                fprintf('    !! Strict not found → SOFT fallback\n');
                bestScore = -inf;
                for ci = 1:nCombos
                    pVec = allCombos(ci,:);
                    grpQ3TP = cell(numel(GroupNames),1);
                    for g = 1:numel(GroupNames)
                        grpQ3TP{g} = arrayfun(@(t) grpQ3All{g}(pVec(t),t), 1:nTP);
                    end
                    [sc_val, isMono, isSep, mSc, sSc] = scorePushComboQ3( ...
                        grpQ3TP, monoWeights, sepWeights, ...
                        lambda_mono, lambda_sep, nTP, GroupNames, sepDirection);
                    if sc_val > bestScore
                        bestScore     = sc_val;
                        bestCombo     = pVec;
                        bestIsMono    = isMono;
                        bestIsSep     = isSep;
                        bestMonoScore = mSc;
                        bestSepScore  = sSc;
                    end
                end
            end

            %% -- Extract optimized per-mouse values -----------
            nM = size(allMouseData{gi},1);
            optMouseVals = nan(nM, nTP);
            for t = 1:nTP
                optMouseVals(:,t) = allMouseData{gi}(:,t,bestCombo(t));
            end

            %% -- Store ----------------------------------------
            res.bestCombo     = bestCombo;
            res.pushNames     = pushOptNames(bestCombo);
            res.foundStrict   = foundStrict;
            res.score         = bestScore;
            res.isMono        = bestIsMono;
            res.isSep         = bestIsSep;
            res.monoScore     = bestMonoScore;
            res.sepScore      = bestSepScore;
            res.mouseIDs      = allMouseIDsG{gi};
            res.mouseData     = allMouseData{gi};
            res.grpQ3All      = grpQ3All{gi};     % nPushOpt×nTP
            res.optMouseVals  = optMouseVals;
            % Store Q3 trajectory of the winning combo
            res.optQ3 = arrayfun(@(t) grpQ3All{gi}(bestCombo(t),t), 1:nTP);

            OptResults(si,ri,gi).res = res;

            fprintf('    Combo  : %s\n',    strjoin(pushOptNames(bestCombo),' | '));
            fprintf('    Q3 traj: %s\n',    num2str(res.optQ3,'%.3f  '));
            fprintf('    Strict : %d | Sep: %d | Score: %.4f\n', ...
                    foundStrict, bestIsSep, bestScore);
            fprintf('    Mono=%.3f  Sep=%.3f\n\n', bestMonoScore, bestSepScore);
        end
    end
end

%% ============================================================
%  SUMMARY TABLE
%% ============================================================
fprintf('\n%s\n', repmat('=',1,90));
fprintf('Q3-BASED OPTIMIZATION SUMMARY\n');
fprintf('%s\n', repmat('=',1,90));
fprintf('%-10s %-10s %-15s %-8s %-8s %-8s %-7s  Push Selection\n', ...
    'Stat','ROI','Group','Strict','Mono','Sep','Score');
fprintf('%s\n', repmat('-',1,110));
for si = 1:numel(statChoices)
    for ri = 1:numel(roiList)
        for gi = 1:numel(GroupNames)
            res = OptResults(si,ri,gi).res;
            fprintf('%-10s %-10s %-15s %-8d %-8d %-8d %-7.3f  %s\n', ...
                statChoices{si}, roiListNames{ri}, GroupNames{gi}, ...
                res.foundStrict, res.isMono, res.isSep, res.score, ...
                strjoin(res.pushNames,' | '));
        end
        fprintf('%s\n', repmat('-',1,110));
    end
end

%% ============================================================
%  PLOT 1 — Per-group boxplots with Q3 trend line
%% ============================================================
fprintf('\n=== Generating Q3-optimized boxplots ===\n');

for si = 1:numel(statChoices)
    sc = statChoices{si};
    for ri = 1:numel(roiList)
        roiNm = roiListNames{ri};

        figBox = figure('Color','w','Position',[50 50 1600 560]);
        sgtitle(sprintf('Q3-Optimized Push — %s | %s', sc, roiNm), ...
            'FontSize',16,'FontWeight','bold');

        for gi = 1:numel(GroupNames)
            res = OptResults(si,ri,gi).res;
            col = colors(gi,:);

            ax = subplot(1,3,gi);
            hold(ax,'on');

            % Boxchart
            datVec=[]; xVec=[];
            for t = 1:nTP
                vals = res.optMouseVals(:,t);
                ok   = ~isnan(vals);
                datVec=[datVec; vals(ok)];           %#ok<AGROW>
                xVec  =[xVec;  repmat(t,sum(ok),1)]; %#ok<AGROW>
            end
            if ~isempty(datVec)
                boxchart(ax, xVec, datVec, ...
                    'BoxWidth',0.45,'BoxFaceColor',col, ...
                    'WhiskerLineColor',col,'MarkerStyle','.','MarkerColor',col);
            end

            % Q3 trend line (the optimized signal)
            plot(ax, 1:nTP, res.optQ3, '-^', ...
                'Color',           col*0.55, ...
                'LineWidth',       2.5, ...
                'MarkerSize',      9, ...
                'MarkerFaceColor', col*0.55, ...
                'DisplayName',     'Q3 trend');

            % Median line (reference, thinner)
            medLine = arrayfun(@(t) median(res.optMouseVals(:,t),'omitnan'),1:nTP);
            plot(ax, 1:nTP, medLine, '--o', ...
                'Color',           col*0.75, ...
                'LineWidth',       1.3, ...
                'MarkerSize',      6, ...
                'MarkerFaceColor', col*0.75, ...
                'DisplayName',     'Median');

            legend(ax,'Q3 trend','Median','Location','northwest', ...
                'FontSize',9,'FontWeight','bold');

            % Push annotations just inside bottom of axes
            yl = ylim(ax);
            yAnn = yl(1) + 0.04*(yl(2)-yl(1));
            for t = 1:nTP
                pN = strrep(res.pushNames{t},'Push','P');
                pN = strrep(pN,'AvgPush','Av');
                text(ax,t,yAnn,pN,'HorizontalAlignment','center', ...
                    'FontSize',8,'FontWeight','bold','Color',col*0.65);
            end

            if res.foundStrict
                badge='STRICT ✓'; bcol=[0 0.52 0];
            else
                badge='SOFT ~';   bcol=[0.72 0.37 0];
            end
            title(ax,sprintf('%s\n(%s)  Score=%.3f',GroupNames{gi},badge,res.score), ...
                'FontSize',12,'FontWeight','bold','Color',bcol);

            xlim(ax,[0.5 nTP+0.5]); xticks(ax,1:nTP);
            xticklabels(ax,string(DayLabels));
            xlabel(ax,'Days after 4T1 implantation','FontSize',11,'FontWeight','bold');
            if gi==1
                ylabel(ax,[sc ' DoA–PD'],'FontSize',11,'FontWeight','bold');
            end
            ax.FontSize=11; ax.FontWeight='bold';
            ax.XColor='k'; ax.YColor='k'; ax.LineWidth=1.3;
            grid(ax,'on'); box(ax,'on');
        end

        saveas(figBox,fullfile(outDirOpt,sprintf('Q3_OptBox_%s_%s.png',sc,roiNm)));
        saveas(figBox,fullfile(outDirOpt,sprintf('Q3_OptBox_%s_%s.fig',sc,roiNm)));
        fprintf('  Saved: Q3_OptBox_%s_%s\n',sc,roiNm);
    end
end

%% ============================================================
%  PLOT 2 — All groups overlaid, Q3 trend lines prominent
%% ============================================================
fprintf('\n=== Generating Q3 combined overlay ===\n');

for si = 1:numel(statChoices)
    sc = statChoices{si};
    for ri = 1:numel(roiList)
        roiNm = roiListNames{ri};

        figComb = figure('Color','w','Position',[50 50 980 680]);
        ax = axes('Parent',figComb);
        hold(ax,'on');

        offC = [-0.28, 0, 0.28];
        hBox = gobjects(numel(GroupNames),1);
        hQ3  = gobjects(numel(GroupNames),1);

        for gi = 1:numel(GroupNames)
            res = OptResults(si,ri,gi).res;
            col = colors(gi,:);

            % Boxchart (offset per group)
            datVec=[]; xVec=[];
            for t = 1:nTP
                vals=res.optMouseVals(:,t);
                ok=~isnan(vals);
                datVec=[datVec; vals(ok)];              %#ok<AGROW>
                xVec  =[xVec;  repmat(t,sum(ok),1)];   %#ok<AGROW>
            end
            if ~isempty(datVec)
                boxchart(ax, xVec+offC(gi), datVec, ...
                    'BoxWidth',0.22,'BoxFaceColor',col, ...
                    'WhiskerLineColor',col,'MarkerStyle','.','MarkerColor',col);
            end

            % Q3 trend — thick, prominent
            hQ3(gi) = plot(ax, (1:nTP)+offC(gi), res.optQ3, '-^', ...
                'Color',           col*0.55, ...
                'LineWidth',       3.0, ...
                'MarkerSize',      10, ...
                'MarkerFaceColor', col*0.55, ...
                'DisplayName',     sprintf('%s Q3',GroupNames{gi}));

            % Dummy handle for box legend entry
            hBox(gi) = plot(ax,nan,nan,'s', ...
                'MarkerFaceColor',col,'MarkerEdgeColor','none', ...
                'MarkerSize',11,'DisplayName',GroupNames{gi});
        end

        % Push annotation rows at bottom
        yl  = ylim(ax);
        yA  = yl(1) - [0.10 0.17 0.24]*(yl(2)-yl(1));
        ax.YLim(1) = yl(1) - 0.28*(yl(2)-yl(1));
        for gi = 1:numel(GroupNames)
            res = OptResults(si,ri,gi).res;
            for t = 1:nTP
                pN = strrep(res.pushNames{t},'Push','P');
                pN = strrep(pN,'AvgPush','Av');
                text(ax, t+offC(gi), yA(gi), pN, ...
                    'HorizontalAlignment','center','FontSize',7.5, ...
                    'FontWeight','bold','Color',colors(gi,:)*0.7);
            end
            % Row label
            text(ax, 0.3, yA(gi), GroupNames{gi}(1:min(4,end)), ...
                'HorizontalAlignment','right','FontSize',7.5, ...
                'FontWeight','bold','Color',colors(gi,:)*0.7);
        end

        % Legend: boxes + Q3 lines
        legend(ax,[hBox; hQ3], ...
            [GroupNames, strcat(GroupNames,' Q3')], ...
            'Location','northwest','FontSize',10,'FontWeight','bold', ...
            'NumColumns',2);

        title(ax, sprintf('Q3-Optimized Push — %s | %s\n▲ = Q3 trend  |  box = distribution', ...
            sc,roiNm),'FontSize',13,'FontWeight','bold');
        xlim(ax,[0.5 nTP+0.5]); xticks(ax,1:nTP);
        xticklabels(ax,string(DayLabels));
        xlabel(ax,'Days after 4T1 implantation','FontSize',13,'FontWeight','bold');
        ylabel(ax,[sc ' DoA–PD'],'FontSize',13,'FontWeight','bold');
        ax.FontSize=12; ax.FontWeight='bold';
        ax.XColor='k'; ax.YColor='k'; ax.LineWidth=1.4;
        grid(ax,'on'); box(ax,'on');

        saveas(figComb,fullfile(outDirOpt,sprintf('Q3_OptCombined_%s_%s.png',sc,roiNm)));
        saveas(figComb,fullfile(outDirOpt,sprintf('Q3_OptCombined_%s_%s.fig',sc,roiNm)));
        fprintf('  Saved: Q3_OptCombined_%s_%s\n',sc,roiNm);
    end
end

%% ============================================================
%  PLOT 3 — Per-mouse trajectories with Q3 overlay
%% ============================================================
fprintf('\n=== Generating Q3 per-mouse trajectory plots ===\n');

basePalette = lines(15);

for si = 1:numel(statChoices)
    sc = statChoices{si};
    for ri = 1:numel(roiList)
        roiNm = roiListNames{ri};
        for gi = 1:numel(GroupNames)
            res   = OptResults(si,ri,gi).res;
            gName = GroupNames{gi};
            nM    = numel(res.mouseIDs);
            mCols = basePalette(mod((1:nM)-1,size(basePalette,1))+1,:);

            figT = figure('Color','w','Position',[60 60 870 520]);
            ax   = axes('Parent',figT);
            hold(ax,'on');

            hLines = gobjects(nM,1);
            for mm = 1:nM
                yV = res.optMouseVals(mm,:);
                hLines(mm) = plot(ax,1:nTP,yV,'-o', ...
                    'Color',mCols(mm,:),'LineWidth',1.8,'MarkerSize',6, ...
                    'MarkerFaceColor',mCols(mm,:),'MarkerEdgeColor','none', ...
                    'DisplayName',sprintf('Mouse %d',res.mouseIDs(mm)));

                missIdx = find(isnan(yV));
                if ~isempty(missIdx)
                    yl  = ylim(ax);
                    yBt = yl(1)-0.05*(yl(2)-yl(1));
                    plot(ax,missIdx,repmat(yBt,1,numel(missIdx)),'x', ...
                        'Color',mCols(mm,:)*0.6,'LineWidth',2.5,'MarkerSize',12, ...
                        'HandleVisibility','off');
                end
            end

            % Q3 trend — thick prominent line
            hQ3line = plot(ax,1:nTP,res.optQ3,'-^k', ...
                'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','k', ...
                'DisplayName','Group Q3');

            % Median reference
            medT = arrayfun(@(t)median(res.optMouseVals(:,t),'omitnan'),1:nTP);
            hMed = plot(ax,1:nTP,medT,'--ok', ...
                'LineWidth',1.5,'MarkerSize',7,'MarkerFaceColor','w', ...
                'DisplayName','Group Median');

            % Push labels at top
            yl2  = ylim(ax);
            yTop = yl2(2)+0.02*(yl2(2)-yl2(1));
            ax.YLim(2) = yl2(2)+0.12*(yl2(2)-yl2(1));
            for t = 1:nTP
                pN = strrep(res.pushNames{t},'Push','P');
                pN = strrep(pN,'AvgPush','Av');
                text(ax,t,yTop,pN,'HorizontalAlignment','center', ...
                    'FontSize',9,'FontWeight','bold','Color',[0.3 0.3 0.3]);
            end

            if res.foundStrict
                modeStr='STRICT'; mc=[0 0.5 0];
            else
                modeStr='SOFT';   mc=[0.68 0.33 0];
            end
            title(ax,sprintf('%s | %s | %s  [%s  Q3-Score=%.3f]', ...
                sc,roiNm,gName,modeStr,res.score), ...
                'FontSize',12,'FontWeight','bold','Color',mc);

            xlim(ax,[0.5 nTP+0.5]); xticks(ax,1:nTP);
            xticklabels(ax,string(DayLabels));
            xlabel(ax,'Days after 4T1 implantation','FontSize',12,'FontWeight','bold');
            ylabel(ax,[sc ' DoA–PD'],'FontSize',12,'FontWeight','bold');
            ax.FontSize=11; ax.FontWeight='bold';
            ax.XColor='k'; ax.YColor='k'; ax.LineWidth=1.3;
            grid(ax,'on'); box(ax,'on');

            legend(ax,[hLines;hQ3line;hMed],'Location','eastoutside', ...
                'FontSize',10,'FontWeight','bold');

            fName = sprintf('Q3_OptTraj_%s_%s_%s',sc,roiNm,gName);
            saveas(figT,fullfile(outDirOpt,[fName '.png']));
            saveas(figT,fullfile(outDirOpt,[fName '.fig']));
            fprintf('  Saved: %s\n',fName);
        end
    end
end

fprintf('\n=== Q3 Push Optimizer complete.\n  Outputs in: %s\n',outDirOpt);
%% ============================================================
%  LOCAL FUNCTION  — must be at the END of the script file
%% ============================================================
function val = doaGetVal(medP2P, pushReq, fallback, fI, roiIdx)
% Extract one scalar from medP2P(pushIdx, fI, roiIdx).
% medP2P size: (nPush+1) x 2 x 3
%   rows 1..nPush = individual pushes
%   last row      = average across all pushes
    nPush = size(medP2P,1) - 1;
    if pushReq == 0
        pIdx = size(medP2P,1);          % average row
    elseif pushReq <= nPush
        pIdx = pushReq;                 % requested push available
    else
        pIdx = fallback;                % push not available → fallback
        fprintf('        Push%d unavailable (nPush=%d) → using Push%d\n', ...
                pushReq, nPush, fallback);
    end
    val = medP2P(pIdx, fI, roiIdx);
end
%% ---- Helper: collect per-mouse per-TP ----------------------
function [mouseVals, mouseIDs] = collectMouseTP(allCanonIDs, diskRegistry, ...
        TimePoints, statChoice, pushReq, fallbackVal, fI, roiIdx)
    mouseIDs  = [];
    mouseVals = [];
    for ii = 1:numel(allCanonIDs)
        cID   = allCanonIDs(ii);
        entry = diskRegistry(cID);
        rowV  = nan(1, numel(TimePoints));
        for t = 1:numel(TimePoints)
            tpV = nan(1, numel(TimePoints{t}));
            for d = 1:numel(TimePoints{t})
                cd = TimePoints{t}{d};
                if ~isKey(entry.dateFiles, cd); continue; end
                try D = load(entry.dateFiles(cd)); catch; continue; end
                if ~isfield(D, statChoice); continue; end
                tpV(d) = doaGetVal(D.(statChoice), pushReq, fallbackVal, fI, roiIdx);
            end
            rowV(t) = mean(tpV,'omitnan');
        end
        mouseIDs(end+1,1)  = cID;       %#ok<AGROW>
        mouseVals(end+1,:) = rowV;      %#ok<AGROW>
    end
end

%% ---- Helper: Q3-based score --------------------------------
%  grpQ3TP{g}(t) = 75th percentile of group g at TP t
%  using the selected push for that TP
function [score, isMono, isSep, monoScore, sepScore] = ...
        scorePushComboQ3(grpQ3TP, monoWeights, sepWeights, ...
                         lambda_mono, lambda_sep, nTP, GroupNames, direction)

    % ---- Monotonicity (Q3 must increase each step) ----------
    monoScore = 0;
    isMono    = true;
    for g = 1:numel(grpQ3TP)
        q3 = grpQ3TP{g};
        for t = 1:nTP-1
            v1 = q3(t); v2 = q3(t+1);
            if isnan(v1) || isnan(v2); continue; end
            if v2 > v1
                monoScore = monoScore + monoWeights(t);
            else
                isMono = false;
            end
        end
    end
    totalMonoW = sum(monoWeights) * numel(grpQ3TP);
    monoScore  = monoScore / max(totalMonoW, 1);

    % ---- Separation (Q3_Control > Q3_Chemo > Q3_Immuno) ----
    sepScore = 0;
    isSep    = true;
    nPairs   = numel(direction) - 1;
    for t = 1:nTP
        for p = 1:nPairs
            g1 = find(strcmp(GroupNames, direction{p}),   1);
            g2 = find(strcmp(GroupNames, direction{p+1}), 1);
            if isempty(g1) || isempty(g2); continue; end
            q1 = grpQ3TP{g1}(t);
            q2 = grpQ3TP{g2}(t);
            if isnan(q1) || isnan(q2); continue; end
            if q1 > q2
                sepScore = sepScore + sepWeights(t);
            else
                isSep = false;
            end
        end
    end
    totalSepW = sum(sepWeights) * nPairs;
    sepScore  = sepScore / max(totalSepW, 1);

    score = lambda_mono*monoScore + lambda_sep*sepScore;
end
