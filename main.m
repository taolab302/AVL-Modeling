% add base into MATLAB path
addpath( genpath('base') );

%% Figure 4a
PlotWholeCellCurrent(650);

%% Figure 4b
PlotTailPeakIV(650);

%% Figure 4c & S3a
PlotWholeCellVoltage('wt', 7000);

%% Figure 4d & S3b
PlotWholeCellVoltage('exp-2(lf)', 7000);

%% Figure 4e & S3c
PlotWholeCellVoltage('exp-2(gf)', 7000);

%% Figure S3d
PlotWholeCellVoltage('shl-1', 8000);

%% Figure S3e
PlotWholeCellVoltage('exp-2;shl-1', 8000);
