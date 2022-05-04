function [ calibrationFile ] = loadCalibrationData(dataFile)
    [num, txt, raw] = xlsread(dataFile, 'Data');
    c_data = {};
    for i=1:length(raw)
        time = 24 * str2num(raw{i, 6})'; % time in days
        data = str2num(raw{i, 5})';
        errors = str2num(raw{i, 7})';
        
        if strcmp(raw{i, 1}, 'calibration') && strcmp(raw{i, 2}, 'lTGFB')
            c_data{1} = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'calibration') && strcmp(raw{i, 2}, 'TNFa') 
            c_data{2} = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'calibration') && strcmp(raw{i, 2}, 'IL1')
            c_data{3} = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'calibration') && strcmp(raw{i, 2}, 'IL6')
            c_data{4} = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'calibration') && strcmp(raw{i, 2}, 'MMP-9')
            c_data{5} = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'calibration') && strcmp(raw{i, 2}, 'GM-CSF')
            c_data{6} = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'calibration') && strcmp(raw{i, 2}, 'Collagen')    
            c_data{7} = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'calibration') && strcmp(raw{i, 2}, 'Macrophages')
            c_data{8} = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'calibration') && strcmp(raw{i, 2}, 'Fibroblasts')
            c_data{9} = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'calibration') && strcmp(raw{i, 2}, 'TIMP-1')
            c_data{10} = {time, data, errors};
        end
    end
    
    calibrationFile = 'calibrationData.mat';
    save(calibrationFile, 'c_data');
end