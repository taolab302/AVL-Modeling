% add base into MATLAB path
addpath( genpath('base') );

%% Figure 4A
PlotWholeCellCurrent(650);

%% Figure 4B
PlotTailPeakIV(650);

%% Figure 4C & S3A
PlotWholeCellVoltage('wt', 7000);

%% Figure 4D & S3B
PlotWholeCellVoltage('exp-2(lf)', 7000);

%% Figure 4E & S3C
PlotWholeCellVoltage('exp-2(gf)', 7000);

%% Figure S3D
PlotWholeCellVoltage('shl-1', 7000);

%% Figure S3E
PlotWholeCellVoltage('exp-2;shl-1', 7000);
