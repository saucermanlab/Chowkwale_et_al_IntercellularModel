%% Sets up data and parameters needed for the model
% Run this once while setting up the model
% Mukti Chowkwale, last edited 1/10/2022

dataFile = 'Datasets.xlsx';
simTime = 0:0.001:30*24; % 30 days, in hours

calibrationName = loadCalibrationData(dataFile);
disp('Created MATLAB data file of calibration data.');

validationName = loadValidationData(dataFile);
disp('Created MATLAB data file of validation data.');