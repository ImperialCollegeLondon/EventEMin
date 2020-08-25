clear, clc, close all;

addpath('pyPlotCMap/');
cMap=getPyPlot_cMap('copper');

%% read events
fname='events.txt';
fid=fopen(fname, 'r');
A=fscanf(fid, '%f %f %f %f %d', [5 Inf])';
fclose(fid);

nEvents=size(A, 1);
K=[226.38018519795807 0 173.6470807871759; 0 226.15002947047415 133.73271487507847; 0 0 1];
for k=1:nEvents
  A(k, 2)=A(k, 4)*(A(k, 2)-K(1, 3))/K(1, 1);
  A(k, 3)=A(k, 4)*(K(2, 3)-A(k, 3))/K(2, 2);
end

%% read events modelled
fname=uigetfile('*.txt');
fid=fopen(fname, 'r');
Am=fscanf(fid, '%f %f %f %f %d', [5 Inf])';
fclose(fid);

%% plot events
figure;
scatter3(A(:, 2), A(:, 4), A(:, 3), 1, A(:, 1));
colormap(cMap);
c=colorbar('southoutside');
c.Label.String='time (s)';
axis tight;
axis off;
view([-0.1 -1 0]);

%% plot events modelled
figure;
scatter3(Am(:, 2), Am(:, 4), Am(:, 3), 1, Am(:, 1));
colormap(cMap);
c=colorbar('southoutside');
c.Label.String='time (s)';
axis tight;
axis off;
view([-0.1 -1 0]);
