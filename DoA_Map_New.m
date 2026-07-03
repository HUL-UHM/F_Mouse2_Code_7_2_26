%% Build DoA_map_fromSaved_interp.mat three ways: Push1, Push2, and Push-Mean
% Walks: <ROOT>\StepXX\AcqYY_AngZZZ...\Process_Fun(or _PI)\xCorL*.mat
% Saves (three times, same filename in separate folders):
%   <ROOT>\Analysis_DoA_AllSteps_P1\DoA_map_fromSaved_interp.mat
%   <ROOT>\Analysis_DoA_AllSteps_P2\DoA_map_fromSaved_interp.mat
%   <ROOT>\Analysis_DoA_AllSteps_Pmean\DoA_map_fromSaved_interp.mat
%
% Each MAT contains (orientation consistent with your prior workflow):
%   DoA_map     [nSteps x nDepth]   % rows: Step01..StepNN, cols: depth_grid
%   depth_grid  [1 x nDepth] (mm)
%   step_labels {1 x nSteps} of 'Step##'
%   meta        struct with details (time, root, push_mode, agg, etc.)

clear; clc;

%% ===== USER SETTINGS =====
ROOT           = 'E:\2026.01.22_20\Phantom_2mm\ARFI\Fun\T_76_Fs_2\';
XCORL_PATTERN  = 'xCorL*.mat';     % your displacement file pattern
AGG            = 'mean';         % lateral aggregation across columns: 'median'|'mean'|'max'
LAT_WINDOW_MM  = [];               % [] = use all laterals; e.g., 2 for ±2 mm around 0 (requires 'lat' in file)
DEPTH_LIMS     = [];               % [] = use overlap automatically; or e.g., [8 30]
COMMON_DY_FALLBACK = 0.1;          % fallback depth step (mm) for common grid if needed

% Build all three outputs in one go:
PUSH_MODES = {'push1','push2','mean'};     % exact modes requested

% ===== Helper: natural sort Step## list =====
stepDirs = dir(fullfile(ROOT,'Step*'));
stepDirs = stepDirs([stepDirs.isdir]);
if isempty(stepDirs)
    error('No Step* folders under %s', ROOT);
end
nums = nan(numel(stepDirs),1);
for i = 1:numel(stepDirs)
    tk = regexp(stepDirs(i).name,'Step(\d+)','tokens','once');
    if ~isempty(tk), nums(i) = str2double(tk{1}); end
end
[~,ord] = sort(nums);
stepDirs = stepDirs(ord);

fprintf('Root: %s\n', ROOT);
fprintf('Found %d steps.\n', numel(stepDirs));

% ===== Main loop over push modes =====
for pm = 1:numel(PUSH_MODES)
    PUSH_MODE = lower(PUSH_MODES{pm});  % 'push1' | 'push2' | 'mean'
    tag = upper(strrep(PUSH_MODE,'mean','Pmean'));  % P1,P2,Pmean

    fprintf('\n==============================\n');
    fprintf('Building DoA map for: %s\n', tag);
    fprintf('==============================\n');

    ALL_DoA   = {};
    ALL_depth = {};
    ALL_step  = {};

    for s = 1:numel(stepDirs)
        stepPath = fullfile(stepDirs(s).folder, stepDirs(s).name);
        fprintf('\n== %s ==\n', stepDirs(s).name);

        % List acquisitions in natural order by Acq##
        acqDirs = dir(fullfile(stepPath, 'Acq*_Ang*'));
        acqDirs = acqDirs([acqDirs.isdir]);
        if isempty(acqDirs)
            fprintf('  (no acquisitions)\n'); continue;
        end
        anum = nan(numel(acqDirs),1);
        for i=1:numel(acqDirs)
            t = regexp(acqDirs(i).name,'Acq(\d+)','tokens','once');
            if ~isempty(t), anum(i) = str2double(t{1}); end
        end
        [~,ao] = sort(anum);
        acqDirs = acqDirs(ao);

        angList = [];
        depthVec = [];
        PeakDA   = [];
        focalDepth = [];

        for a = 1:numel(acqDirs)
            acqPath = fullfile(acqDirs(a).folder, acqDirs(a).name);

            % Parse angle from name (Ang###)
            tok = regexp(acqDirs(a).name,'Ang(\d+)','tokens','once');
            if isempty(tok), continue; end
            ang = str2double(tok{1});

            % Find a displacement file (prefer Process_Fun, then PI, then fallback recursive)
            mfiles = dir(fullfile(acqPath,'Process_Fun',XCORL_PATTERN));
            if isempty(mfiles), mfiles = dir(fullfile(acqPath,'Process_PI',XCORL_PATTERN)); end
            if isempty(mfiles), mfiles = dir(fullfile(acqPath,'**',XCORL_PATTERN)); end
            if isempty(mfiles)
                fprintf('  (no %s under %s)\n', XCORL_PATTERN, acqPath); 
                continue;
            end
            [~, idxLatest] = max([mfiles.datenum]);
            fmat = fullfile(mfiles(idxLatest).folder, mfiles(idxLatest).name);
            S = load(fmat);
            if ~isfield(S,'Peak') || ~isfield(S,'axial')
                fprintf('  Missing Peak/axial in %s (skipped)\n', fmat);
                continue;
            end

            Peak  = S.Peak;           % could be [ax x 2], [ax x lat], or [ax x lat x push]
            axial = S.axial(:);       % [ax x 1]
            lat   = []; if isfield(S,'lat') && ~isempty(S.lat), lat = S.lat(:); end
            if isempty(depthVec), depthVec = axial; end
            if isfield(S,'aimgIni') && isfield(S.aimgIni,'focalDepthmm')
                focalDepth = S.aimgIni.focalDepthmm;
            end

            % ---- Normalize Peak to [ax x lat x push] ----
            if ndims(Peak) == 2
                [axR, c] = size(Peak);
                if c==2 && (isempty(lat) || numel(lat)<=2)
                    % [ax x push] -> [ax x 1 lat x 2]
                    Peak = reshape(Peak,[axR,1,2]);
                    if isempty(lat), lat = (1:1)'; end
                else
                    % [ax x lat] -> [ax x lat x 1]
                    Peak = reshape(Peak,[axR,c,1]);
                    if isempty(lat), lat = (1:c)'; end
                end
            elseif ndims(Peak) ~= 3
                fprintf('  Unexpected Peak dims in %s\n', fmat); 
                continue;
            end
            nPush = size(Peak,3);

            % ---- Select/Combine pushes per requested mode ----
            switch PUSH_MODE
                case 'push1'
                    if nPush < 1, fprintf('  no push1 in %s\n', fmat); continue; end
                    Peak2D = Peak(:,:,1);     % [ax x lat]
                case 'push2'
                    if nPush < 2, fprintf('  no push2 in %s\n', fmat); continue; end
                    Peak2D = Peak(:,:,2);
                case 'mean'
                    if nPush >= 2
                        Peak2D = mean(Peak, 3, 'omitnan');
                    else
                        Peak2D = Peak(:,:,1);
                    end
                otherwise
                    error('Unknown PUSH_MODE: %s', PUSH_MODE);
            end

            % ---- Lateral window (optional) ----
            latMask = true(1, size(Peak2D,2));
            if ~isempty(LAT_WINDOW_MM) && ~isempty(lat) && numel(lat)==size(Peak2D,2)
                latMask = (lat >= -LAT_WINDOW_MM) & (lat <= LAT_WINDOW_MM);
                if ~any(latMask), latMask = true(size(latMask)); end
            end

            % ---- Aggregate across lateral (AGG) -> one column per angle ----
            P_use = Peak2D(:, latMask);
            switch lower(AGG)
                case 'median', pCol = nanmedian(P_use, 2);
                case 'mean',   pCol = nanmean(P_use, 2);
                case 'max',    pCol = nanmax(P_use, [], 2);
                otherwise,     pCol = nanmedian(P_use, 2);
            end

            angList(end+1,1) = ang;              %#ok<AGROW>
            PeakDA(:, end+1) = pCol;             %#ok<AGROW>
            fprintf('  + angle %g°\n', ang);
        end

        if isempty(PeakDA)
            fprintf('  (no usable acquisitions in %s)\n', stepPath);
            continue;
        end

        % Sort columns by angle
        [angList, ord] = sort(angList);
        PeakDA = PeakDA(:, ord);
        % ---- DoA vs depth: FIXED ANGLES ----
        %idx090 = find(angList == 90, 1);
        idx000 = find(angList == 0,  1);
        
%         if isempty(idx090) || isempty(idx000)
%             error('Required angles Ang090 or Ang000 not found in %s.', stepDirs(s).name);
%         end
        % ---- DoA vs depth: per depth max/min across angles ----
        %maxP = PeakDA(:, idx090);   % Ang090 nanmax(PeakDA, [], 2);
        minP = PeakDA(:, idx000);   % Ang000 nanmin(PeakDA, [], 2);
        %DoA  = maxP ./ max(minP, eps);
        DoA = minP;
        % Collect for all-steps common grid
        ALL_DoA{end+1,1}  = DoA;           %#ok<AGROW>
        ALL_depth{end+1,1}= depthVec;      %#ok<AGROW>
        ALL_step{end+1,1} = stepDirs(s).name; %#ok<AGROW>
        fprintf('  -> collected %s\n', stepDirs(s).name);
    end

    % ===== Build a common depth grid & map across steps =====
    if isempty(ALL_DoA)
        warning('No DoA curves collected for %s. Skipped.', tag);
        continue;
    end

    % Overlap limits
    mins = cellfun(@(v) min(v), ALL_depth);
    maxs = cellfun(@(v) max(v), ALL_depth);
    yMin = max(mins); yMax = min(maxs);
    if ~isempty(DEPTH_LIMS)
        yMin = max(yMin, DEPTH_LIMS(1));
        yMax = min(yMax, DEPTH_LIMS(end));
    end
    if ~(yMax > yMin)
        error('No overlapping depth region across steps for %s.', tag);
    end

    % Choose a reasonable dy from native spacings
    dAll = cellfun(@(v) median(diff(v)), ALL_depth);
    dy   = median(dAll(~isnan(dAll) & dAll>0));
    if isempty(dy) || ~isfinite(dy) || dy<=0
        dy = COMMON_DY_FALLBACK;
    end
    depth_grid = yMin:dy:yMax;        % 1 x nDepth (row)

    nSteps = numel(ALL_DoA);
    Z = nan(numel(depth_grid), nSteps);  % [nDepth x nSteps]

    for k = 1:nSteps
        z_k = ALL_depth{k}(:);
        d_k = ALL_DoA{k}(:);
        % ensure ascending
        if z_k(1) > z_k(end)
            z_k = flip(z_k); d_k = flip(d_k);
        end
        Z(:,k) = interp1(z_k, d_k, depth_grid, 'linear', NaN);
    end

    % Orient to [nSteps x nDepth] to match your previous MAT (rows=steps)
    DoA_map = Z.';                  % [nSteps x nDepth]
    step_labels = ALL_step(:)';     % 1 x nSteps cell row
    depth_grid = depth_grid(:)';    % 1 x nDepth row

    % ===== Save to push-specific folder with the same filename =====
    outDir = fullfile(ROOT, ['Analysis_DoA_AllSteps_' tag]);  % e.g., ...\Analysis_DoA_AllSteps_P1
    if ~exist(outDir,'dir'), mkdir(outDir); end

    meta.root        = ROOT;
    meta.push_mode   = PUSH_MODE;   % 'push1'|'push2'|'mean'
    meta.agg         = AGG;
    meta.lat_window  = LAT_WINDOW_MM;
    meta.depth_lims  = DEPTH_LIMS;
    meta.generated   = datestr(now, 30);
    meta.note        = 'DoA_map is [nSteps x nDepth]; depth_grid is row vector in mm; step_labels are StepXX strings.';

    save(fullfile(outDir, 'DoA_map_fromSaved_interp.mat'), ...
         'DoA_map', 'depth_grid', 'step_labels', 'meta', '-v7.3');

    fprintf('Saved: %s\n', fullfile(outDir, 'DoA_map_fromSaved_interp.mat'));
end

fprintf('\nAll three versions completed (P1, P2, Pmedian).\n');

%% Export DoA (all steps) to Excel:  Depth | Step01 | Step02 | ... 
clear; clc;

matFile = fullfile('E:\2026.01.22_20\Phantom_2mm\ARFI\Fun\T_76_Fs_2\', ...
                   'Analysis_DoA_AllSteps_PMEAN', 'DoA_map_fromSaved_interp.mat');   % <-- path to your MAT
outXLS  = fullfile(fileparts(matFile), 'DoA_map_fromSaved_interp.xlsx');           % output Excel

S = load(matFile);    % expects fields: depth_grid, DoA_map, step_labels (if available)

% --- depth vector (ensure column) ---
if isfield(S,'depth_grid')
    depth = S.depth_grid(:);
else
    error('depth_grid not found in %s', matFile);
end

% --- DoA map (robust handling of shapes/cells) ---
if isfield(S,'DoA_map')
    M = S.DoA_map;
    if iscell(M)
        % Cell array: each cell is a step vector
        nSteps = numel(M);
        Z = nan(numel(depth), nSteps);
        for k = 1:nSteps
            v = M{k}(:);
            if numel(v) == numel(depth)
                Z(:,k) = v;
            else
                % interpolate to depth_grid if needed
                warning('Step %d length mismatch -> interpolating to depth_grid', k);
                % assume the k-th cell also has a corresponding depth vector S.depth_cell{k}
                % if not available, we try to linearly stretch to match (last resort)
                if isfield(S,'depth_cells') && numel(S.depth_cells) >= k
                    Z(:,k) = interp1(S.depth_cells{k}(:), v, depth, 'linear', NaN);
                else
                    Z(:,k) = interp1(linspace(depth(1),depth(end),numel(v))', v, depth, 'linear', NaN);
                end
            end
        end
    else
        % Numeric matrix/vector
        M = double(M);
        if isvector(M)
            % Only one step stored: make it a single column
            Z = M(:);
            if numel(Z) ~= numel(depth)
                error('Vector DoA_map length (%d) != depth_grid length (%d).', numel(Z), numel(depth));
            end
            Z = Z(:); % [nDepth x 1]
        else
            % Matrix: decide orientation. We want [nDepth x nSteps]
            [r,c] = size(M);
            if r == numel(depth)
                Z = M;                   % already [nDepth x nSteps]
            elseif c == numel(depth)
                Z = M.';                 % transpose to [nDepth x nSteps]
            else
                error('DoA_map size [%d x %d] not compatible with depth length %d.', r, c, numel(depth));
            end
        end
    end
else
    error('DoA_map not found in %s', matFile);
end

% --- Step labels ---
if isfield(S,'step_labels') && ~isempty(S.step_labels)
    if iscell(S.step_labels)
        labels = S.step_labels(:).';   % cell row
    else
        labels = cellstr(S.step_labels); 
    end
else
    % fabricate labels if missing
    nSteps = size(Z,2);
    labels = arrayfun(@(k) sprintf('Step%02d',k), 1:nSteps, 'uni',0);
end

% --- Build table (Depth_mm + one column per step) ---
T = table(depth, 'VariableNames', {'Depth_mm'});

% If Z ended up as a vector, make it a 2-D [nDepth x 1]
if isvector(Z), Z = Z(:); end

% Ensure labels length matches number of step columns
nSteps = size(Z,2);
if numel(labels) ~= nSteps
    warning('Number of step_labels (%d) != columns in DoA_map (%d). Relabeling.', numel(labels), nSteps);
    labels = arrayfun(@(k) sprintf('Step%02d',k), 1:nSteps, 'uni',0);
end

for k = 1:nSteps
    T.(labels{k}) = Z(:,k);
end

% --- Write Excel ---
writetable(T, outXLS, 'Sheet', 'DoA_all_steps');
fprintf('Excel written: %s\n', outXLS);
%%
% ====== Post-processing & mapping of DoA (after Excel is saved) ======
% Axes:
%   x_mm  : real lateral axis in mm for the original DoA columns
%   xq_mm : lateral axis in mm for the interpolated/filtered grid
%   y (depth) already in mm from Excel
% All figures/ticks: show mm (not steps)

clear; clc;close all

% ================= INPUTS =================
%xlsPath   = 'C:\Users\huluh\Downloads\DoA_map_fromSaved_interp.xlsx';
xlsPath='E:\2026.01.22_20\Phantom_2mm\ARFI\Fun\T_76_Fs_2\Analysis_DoA_AllSteps_PMEAN\DoA_map_fromSaved_interp.xlsx'
sheetName = 'DoA_all_steps';
CLIM      = [0 3];
OUT_DIR   = fileparts(xlsPath);

MOVMEAN_WIN = 60;      % depth smoothing
dx_mm      = 2;     % lateral step per column (mm)
PX_STEP_MM = 0.1;     % interp lateral resolution (mm)
PX_MM      = 0.1;     % interp depth resolution (mm)

MEDFILT_WIN = [35 35];  % REQUIRED median filter window (pixels)

% ================= LOAD EXCEL =================
T = readtable(xlsPath,'Sheet',sheetName,'VariableNamingRule','preserve');
vars = string(T.Properties.VariableNames);
depthCol = find(vars=="Depth_mm" | vars=="Depth",1);

depth = T{:,depthCol};
Z     = T{:,setdiff(1:width(T),depthCol)};
[nDepth,nSteps] = size(Z);

if depth(1) > depth(end)
    depth = flip(depth);
    Z     = flipud(Z);
end

% ================= LATERAL AXIS =================
x_mm = linspace(0, dx_mm*nSteps, nSteps);
% ================= PRE-PROCESS =================
Z(Z > 1.5) = NaN;

% Fill NaNs along depth
for j = 1:nSteps
    v = isfinite(Z(:,j));
    if nnz(v) >= 3
        Z(:,j) = interp1(depth(v), Z(v,j), depth, 'linear', NaN);
    end
end

% Depth smoothing
Z_smooth = Z;
for j = 1:nSteps
    Z_smooth(:,j) = movmean(Z(:,j), MOVMEAN_WIN, 'omitnan','Endpoints','shrink');
    %Z_smooth(:,j) = movmedian(Z_smooth(:,j), MOVMEAN_WIN, 'omitnan','Endpoints','shrink');
    %Z_smooth(:,j) = movmedian(Z_smooth(:,j), MOVMEAN_WIN, 'omitnan','Endpoints','shrink');
end

% ================= RAW PLOT =================
figure('Color','w');
imagesc(x_mm, depth, flipud(Z_smooth));
set(gca,'YDir','normal'); axis tight;
xlabel('Lateral (mm)'); ylabel('Axial depth (mm)');
title('DoA — movmean smoothed');
colorbar; caxis(CLIM);

saveas(gcf, fullfile(OUT_DIR,'DoA_RAW_mm.png'));

% ================= INTERP2 =================
xq_mm = 0:PX_STEP_MM:(dx_mm*nSteps);
yq_mm = min(depth):PX_MM:max(depth);

[X,Y]   = meshgrid(x_mm, depth);
[Xq,Yq] = meshgrid(xq_mm, yq_mm);

Zq = interp2(X, Y, Z_smooth, Xq, Yq, 'makima');

figure('Color','w');
imagesc(xq_mm, yq_mm, Zq);
set(gca,'YDir','normal'); axis image;
xlabel('Lateral (mm)'); ylabel('Axial depth (mm)');
title('DoA — interp2 (makima)');
c= colorbar; %caxis(CLIM);
%c.Position=[0.9 2 3 5.5];    

saveas(gcf, fullfile(OUT_DIR,'DoA_interp2_mm.png'));

% ================= REQUIRED MEDFILT2 =================
K_clean = Zq;

% Fill NaNs before medfilt2 (required)
nanMask = isnan(K_clean);
K_clean(nanMask) = median(K_clean(:),'omitnan');

% REQUIRED METHOD
K_processed = medfilt2(K_clean, MEDFILT_WIN);

% Restore NaNs
K_processed(nanMask) = NaN;

Zq_med = K_processed;
Z_plot = Zq_med;
Zq_med_flip = flipud(Zq_med);
% % % Mask middle range
% mask_mid = (Zq_med_flip > 1.5 & Zq_med_flip < 4.0);
% Zq_med_flip(mask_mid) = NaN;

figure('Color','w');
imagesc(xq_mm, yq_mm, Zq_med);
set(gca,'YDir','reverse');   % << ONLY change is here
axis image;
xlabel('Lateral (mm)');
ylabel('Axial depth (mm)');
title(sprintf('DoA Map', MEDFILT_WIN));
colorbar; 
colormap(jet);
%caxis([0.8 5.4]);
cb.Label.String = 'DoA Ratio';
              % show from 0 to -12
uiopen('"E:\2026.01.22_20\untitled.fig"',1)

%set(gca,'XDir','reverse');   % flip x-axis direction
%xlim([-12 0]);
%ylim([15 50]);
%saveas(gcf, fullfile(OUT_DIR,'DoA_medfilt2_mm.png'));

% ================= SAVE FINAL MAP =================


meta.generated_on = datestr(now,'yyyy-mm-dd HH:MM:SS');
meta.filter       = 'medfilt2';
meta.window       = MEDFILT_WIN;
meta.dx_mm        = dx_mm;

save(fullfile(OUT_DIR,'DoA_medfilt2_final.mat'), ...
    'Zq_med','Zq_med_flip','xq_mm','yq_mm','x_mm','depth', ...
    'MEDFILT_WIN','CLIM','meta','-v7.3');

% ================= INPUT-STYLE SAVE =================
DoA_map    = Z_smooth;
Zq         = Zq_med;
x          = x_mm;
xq         = xq_mm;
yq         = yq_mm;
depth_grid = depth(:);

save(fullfile(OUT_DIR,'DoA_input_style_final.mat'), ...
    'DoA_map','Zq','x','xq','yq','depth_grid','CLIM','-v7.3');

disp('✅ Finished: medfilt2-based DoA processing');
%%
clear; clc; close all;

% ================= USER INPUTS =================
xlsPath   = 'E:\2026.01.22_20\Phantom_2mm\ARFI\Fun\T_76_Fs_2\Analysis_DoA_AllSteps_PMEAN\DoA_map_fromSaved_interp.xlsx';
sheetName = 'DoA_All_Steps';
bmode_fig = 'E:\2026.01.22_20\untitled.fig';
xShift_mm = 0;   % desired +x shift in mm

CLIM = [0.5 1.1];

dx_mm = 2;        % lateral step per column (mm)
PX_STEP_MM = 0.1;   % interp lateral resolution (mm)
PX_MM = 0.1;        % interp depth resolution (mm)

MOVMEAN_WIN = 60;
MEDFILT_WIN = [15 15];

OUT_DIR = fileparts(xlsPath);

% ================= LOAD DoA DATA =================
T = readtable(xlsPath,'Sheet',sheetName,'VariableNamingRule','preserve');
vars = string(T.Properties.VariableNames);
depthCol = find(vars=="Depth_mm" | vars=="Depth",1);

depth = T{:,depthCol};
Z = T{:,setdiff(1:width(T),depthCol)};
[nDepth,nSteps] = size(Z);

% Ensure depth increasing downward
if depth(1) > depth(end)
    depth = flip(depth);
    Z = flipud(Z);
end

x_mm = linspace(0, dx_mm*nSteps, nSteps);

% Remove outliers
Z(Z > 1.1) = NaN;

% Fill NaNs along depth
for j = 1:nSteps
    v = isfinite(Z(:,j));
    if nnz(v) >= 3
        Z(:,j) = interp1(depth(v), Z(v,j), depth, 'linear', NaN);
    end
end

% Depth smoothing
for j = 1:nSteps
    Z(:,j) = movmean(Z(:,j), MOVMEAN_WIN, 'omitnan','Endpoints','shrink');
end

% ================= INTERPOLATE DoA =================
xq_mm = 0:PX_STEP_MM:(dx_mm*nSteps);
yq_mm = min(depth):PX_MM:max(depth);

[X,Y]   = meshgrid(x_mm, depth);
[Xq,Yq] = meshgrid(xq_mm, yq_mm);

Zq = interp2(X, Y, Z, Xq, Yq, 'makima');

% Median filter (required)
nanMask = isnan(Zq);
Zq(nanMask) = median(Zq(:),'omitnan');
Zq = medfilt2(Zq, MEDFILT_WIN);
Zq(nanMask) = NaN;

Zq_med = Zq;

% ================= LOAD B-MODE FROM FIG =================
hFig = openfig(bmode_fig,'invisible');
axB  = findobj(hFig,'Type','axes');
imgB = findobj(axB,'Type','image');

Bmode = imgB.CData;
xB = imgB.XData;
yB = imgB.YData;

% Expand X/Y vectors if needed
if numel(xB)==2
    xB = linspace(xB(1), xB(2), size(Bmode,2));
end
if numel(yB)==2
    yB = linspace(yB(1), yB(2), size(Bmode,1));
end

xlim_B = axB.XLim;
ylim_B = axB.YLim;
xdir_B = axB.XDir;
ydir_B = axB.YDir;

close(hFig);

% ================= DRAW ROI ON B-MODE =================
figure('Color','w');
imagesc(xB, yB, Bmode);
set(gca,'XDir',xdir_B,'YDir',ydir_B);
axis image;
%xlim(xlim_B); ylim(ylim_B);
colormap gray;
colorbar;
xlabel('Lateral (mm)');
ylabel('Axial depth (mm)');
title('B-mode (Draw ROI)');

disp('👉 Draw ROI on B-mode (double-click to finish)');
hROI = drawpolygon('Color','r','LineWidth',2);
roi_mm = hROI.Position;     % ROI in B-mode coordinates (mm)

% ================= B-MODE WITH ROI CONTOUR (SIGN-FLIPPED DISPLAY) =================

% Flip x only for visualization (do NOT touch raw data)
%xShift_mm = 5;   % +5 mm shift

xB_plot  = xB + xShift_mm;      % shift B-mode image
roi_Bplt = roi_mm;
roi_Bplt(:,1) = roi_Bplt(:,1) + xShift_mm;   % shift ROI consistently


figure('Color','w');
imagesc(xB_plot, yB, Bmode);
set(gca,'YDir',ydir_B,'XDir','normal');
axis image;
%xlim([0 12]);
%ylim([10 35]);

colormap gray;
colorbar;

hold on;
plot(roi_Bplt(:,1), roi_Bplt(:,2), 'm--', 'LineWidth', 4);
hold off;

xlabel('Lateral (mm)');
ylabel('Axial depth (mm)');
title('B-mode with ROI contour');

% ================= SIGN FIX FOR DoA =================
% B-mode: 0 → -12 mm
% DoA:    0 → +12 mm
%roi_DoA = roi_mm;
%roi_DoA(:,1) = roi_DoA(:,1);   % SIGN FLIP ONLY HERE

% ================= DoA WITH ROI CONTOUR =================
figure('Color','w');

Zshow = Zq_med;
Zshow(isnan(Zshow)) = -1;   % visualize NaNs
Zshow_flipped = fliplr(Zshow);   % flip image data horizontally

imagesc(xq_mm, yq_mm, Zshow_flipped);
set(gca,'YDir','reverse','XDir','normal');
axis image;
%xlim([0 12]);
ylim([15,45])

colormap(jet);
%caxis([0.8 5.4]);
cb = colorbar;
cb.Label.String = 'DoA Ratio';

hold on;
%plot(roi_DoA(:,1), roi_DoA(:,2), 'black--', 'LineWidth', 4);
plot(roi_Bplt(:,1), roi_Bplt(:,2), 'black--', 'LineWidth', 4);
hold off;

xlabel('Lateral (mm)');
ylabel('Axial depth (mm)');
title('DoA Map with ROI contour');

% % ================= SAVE ROI =================
% ROI.Bmode_mm = roi_mm;
% ROI.DoA_mm   = roi_DoA;
% ROI.generated = datestr(now);
% 
% save(fullfile(OUT_DIR,'ROI_contour_Bmode_DoA.mat'),'ROI');

disp('✅ Finished: ROI contour correctly shown on B-mode and DoA');

figure('Color','w');
imagesc(xB, yB, Bmode);
set(gca,'YDir','reverse','XDir','normal');
axis image;
colormap gray;
xlabel('Lateral (mm)');
ylabel('Axial depth (mm)');
title('Draw ROIs: 1) Tumor, 2) Background');

% ---- Tumor ROI ----
disp('👉 Draw TUMOR ROI (double-click to finish)');
hROI_tumor = drawpolygon('Color','r','LineWidth',2);
roiTumor_mm = hROI_tumor.Position;

% ---- Background ROI ----
disp('👉 Draw BACKGROUND ROI (double-click to finish)');
hROI_bkg = drawpolygon('Color','b','LineWidth',2);
roiBkg_mm = hROI_bkg.Position;
[Xq, Yq] = meshgrid(xq_mm, yq_mm);

% Shift ROIs to DoA frame
roiTumor_DoA = roiTumor_mm;
roiTumor_DoA(:,1) = roiTumor_DoA(:,1) + xShift_mm;

roiBkg_DoA = roiBkg_mm;
roiBkg_DoA(:,1) = roiBkg_DoA(:,1) + xShift_mm;

% Masks
maskTumor = inpolygon(Xq, Yq, roiTumor_DoA(:,1), roiTumor_DoA(:,2));
maskBkg   = inpolygon(Xq, Yq, roiBkg_DoA(:,1),   roiBkg_DoA(:,2));
Z_doA = nan(size(Zq_med));   % everything transparent by default

% Tumor DoA
Z_doA(maskTumor) = Zq_med(maskTumor);

% Background DoA
Z_doA(maskBkg)   = Zq_med(maskBkg);
figure;

arfiBmodeOverlay( ...
    xq_mm, yq_mm, ...
    Z_doA, ...              % tumor + background DoA
    [0.7 1.2], ...
    xB + xShift_mm, yB, ...
    Bmode, ...
    0.4, ...
    'ARFI DoA (Tumor + Background) over B-mode');

hold on;
contour(xq_mm, yq_mm, maskTumor, [1 1], 'k--', 'LineWidth', 5);
contour(xq_mm, yq_mm, maskBkg,   [1 1], 'c--', 'LineWidth', 0);
hold off;

axis image;
xlabel('Lateral (mm)');
ylabel('Axial depth (mm)');
ylim([15,45]);
xlim([0,23]);
title('CIRS Phantom SMR Prediction Map overlaid on B-mode')