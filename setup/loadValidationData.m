% Validation data points
function [validationName] = loadValidationData(dataFile)

    [num, txt, raw] = xlsread(dataFile, 'Data');
    for i=1:length(raw)
        time = 24 * str2num(raw{i, 6})'; % time in days
        data = str2num(raw{i, 5})';
        errors = str2num(raw{i, 7})';
        
        if strcmp(raw{i, 1}, 'validation') && strcmp(raw{i, 2}, 'lTGFB')
            ltgfb = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'validation') && strcmp(raw{i, 2}, 'TNFa')   
            tnfa = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'validation') && strcmp(raw{i, 2}, 'IL1')
            il1 = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'validation') && strcmp(raw{i, 2}, 'IL6')
            il6 = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'validation') && strcmp(raw{i, 2}, 'MMP-9')
            mmp9 = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'validation') && strcmp(raw{i, 2}, 'GM-CSF')
            gmcsf = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'validation') && strcmp(raw{i, 2}, 'Collagen')
            collagen = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'validation') && strcmp(raw{i, 2}, 'Macrophages')
            macro = {time, data, errors};
        elseif strcmp(raw{i, 1}, 'validation') && strcmp(raw{i, 2}, 'Fibroblasts')
            fibro = {time, data, errors};
        end
    end
    
    validationName = 'validationData.mat';
    validation = {ltgfb; tnfa; il1; []; mmp9; gmcsf; collagen; macro; fibro}; % [] placeholder
    save(validationName, 'validation');
end