
%% ============================
% INPUT
% =============================
% PeakCom_1 assumed size: (Axial x Angle x Lat)
% axial is your actual axial vector

% Example:
% load('your_file.mat');   % must contain PeakCom_1 and axial

%% ============================
% CHECK SIZE
% =============================
[nAx, nAng, nLat] = size(PeakCom_1);

AxialVec = axial(:);

if length(AxialVec) ~= nAx
    error('Length of axial vector does not match first dimension of PeakCom_1.');
end

%% ============================
% CREATE AXIS VECTORS
% =============================
LatVec   = (0:nLat-1) * 1.5;   % 0, 1.5, 3, 4.5, ...
AngleIdx = 1:nAng;             % unknown angle, use index

%% ============================
% MAKE FULL ROW-WISE TABLE
% =============================
[Ax, Ang, Lat] = ndgrid(AxialVec, AngleIdx, LatVec);

T = table(Ax(:), Ang(:), Lat(:), PeakCom_1(:), ...
    'VariableNames', {'Axial','Angle','Lat','PD'});

%% ============================
% SELECT NEAREST AXIAL SLICE TO 25
% =============================
targetAxial = 16;

[~, idxAx] = min(abs(AxialVec - targetAxial));
exactAxial = AxialVec(idxAx);

fprintf('Requested axial = %.4f\n', targetAxial);
fprintf('Using nearest available axial = %.4f\n', exactAxial);

Tsub = T(T.Axial == exactAxial, :);

%% ============================
% EXPORT
% =============================
outFile = sprintf('Axial_%.4f.csv', exactAxial);
writetable(Tsub, outFile);

fprintf('Saved: %s\n', outFile);
fprintf('Number of rows saved: %d\n', height(Tsub));

%% ============================
% CHECK
% =============================
disp(unique(Tsub.Axial));