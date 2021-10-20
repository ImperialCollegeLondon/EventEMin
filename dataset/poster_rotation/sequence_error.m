clear, clc, close all;

%% read ground-truth vars
disp('please select the file containing the ground-truth (it should be the imu.txt file)');
[file, folder] = uigetfile({'*.txt';'*.*'}, 'File Selector');
fimu = fullfile(folder, file);

A = readmatrix(fimu);
Aw = [A(:, 1) rad2deg(A(:, 5)) rad2deg(A(:, 6)) rad2deg(A(:, 7))];

%% read estimated vars
disp('please select the file containing the estimates');
[file, folder] = uigetfile({'*.txt';'*.*'}, 'File Selector');
festimates = fullfile(folder, file);

B = readmatrix(festimates);
delay = 2.4e-3;
Bw = [(B(:, 1) - delay) rad2deg(B(:, 2)) -rad2deg(B(:, 3)) rad2deg(B(:, 4))];
maxT = 59.75;
Aw = Aw(Aw(:, 1) <= maxT, :);
Bw = Bw(Bw(:, 1) <= maxT, :);

%% error computation
Awi = interp1(Aw(:, 1), Aw(:, 2:4), Bw(:, 1), 'spline');
erw = Awi - Bw(:, 2:4);

% error statistics
erwAbsMean = mean(abs(erw));
fprintf('\nmean (ex, ey, ez): %f %f %f\n', erwAbsMean);
fprintf('mean (ew): %f\n', mean(erwAbsMean));
erwStd = std(erw);
fprintf('std (ex, ey, ez): %f %f %f\n', erwStd);
fprintf('std (ew): %f\n', mean(erwStd));
erwRms = rms(erw);
erwMeanRms = mean(erwRms);
fprintf('rms (ex, ey, ez): %f %f %f\n', erwRms);
fprintf('rms (ew, percentage): %f %f\n', erwMeanRms, ...
    100 * erwMeanRms/max(abs(Aw(:, 2:4)), [], 'all'));
