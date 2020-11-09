clear, clc, close all;

%% read ground-truth vars
disp('please select the file containing the ground-truth (it should be the imu.txt file)');
[file, folder]=uigetfile;
fimu=fullfile(folder, file);

fid=fopen(fimu, 'r');
Aimu=fscanf(fid, '%f %f %f %f %f %f %f', [7 Inf])';
fclose(fid);
Aw=[Aimu(:, 1) rad2deg(Aimu(:, 5)) rad2deg(Aimu(:, 6)) rad2deg(Aimu(:, 7))];

%% read estimated vars
disp('please select the file containing the estimates');
[file, folder]=uigetfile;
festimates=fullfile(folder, file);

fid=fopen(festimates, 'r');
B=fscanf(fid, '%f %f %f %f', [4 Inf])';
fclose(fid);
Bw=[B(:, 1) rad2deg(B(:, 2)) -rad2deg(B(:, 3)) rad2deg(B(:, 4))];
maxT=59.75;
Aw=Aw(Aw(:, 1)<=maxT, :);
Bw=Bw(Bw(:, 1)<=maxT, :);

%% error computation
nparams=3;
Awi=zeros(size(Bw, 1), nparams);
for p=1:nparams
    pp=p+1;
    Awi(:, p)=interp1(Aw(:, 1), Aw(:, pp), Bw(:, 1), 'spline');
end

er=Awi-Bw(:, 2:4);

% error statistics
fprintf('\nmean (ex, ey, ez): %f %f %f %f\n', mean(abs(er)));
fprintf('mean (ew): %f\n', mean(mean(abs(er))));
fprintf('std (ex, ey, ez): %f %f %f\n', std(er));
fprintf('std (ew): %f\n', mean(std(er)));
fprintf('rms (ex, ey, ez): %f %f %f\n', rms(er));
fprintf('rms (ew, percentage): %f %f\n', mean(rms(er)), 100*mean(rms(er))/max(abs(Aw(:, 2:4)), [], 'all'));
