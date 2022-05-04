%% Validation function
function [percentAgree, resultChart] = Validation(filename, threshold)
    % Inputs:
    % filename - validation excel file
    % steadyState - steadyState.mat for input conditions
    % baseline - baseline.mat for comparison
    % Outputs:
    % percentAgree - validation percentage
    % resultChart - summary of validation information
    
    [~, txt, raw] = xlsread(filename);

    ids = txt(2:end, 1);
    perturbations = txt(2:end, 2); % Perturbations
    perturbCodes = txt(2:end, 3); % Code for perturbations
    output = raw(2:end, 5); % Output indices
    validTimes = cell2mat(raw(2:end, 6)); % timepoint (in days)
    measurement = txt(2:end, 7); % Increase/decrease/no change

    % Output variables
    percentChangeAct = zeros(1, length(output));
    percentChangeActStr = cell(1, length(output));
    prediction = cell(1, length(output));
    numMatching = 0;
    match = cell(1, length(output));
       
    % Simulation parameters
    dt = 0.001; % step size
    
    for i=1:length(ids)
        disp(i);
        % Baseline
        [ parameters, constants, receptors, knockouts ] = loadParameters();
        
        infarct_size = 1; codeSnip = '';
        [ y0, constants, codeSnip ] = initParams(infarct_size, constants, codeSnip);
        time = 0:dt:(24*validTimes(i));
        [~, y] = ode45(@(t, y) modelEquations(t, y, parameters, ...
            constants, receptors, knockouts, codeSnip),...
            time, y0);
        control = y;
        
        % Perturbation 
        [ parameters, constants, receptors, knockouts ] = loadParameters();
        infarct_size = 1; codeSnip = '';
        [ y0, constants, codeSnip ] = initParams(infarct_size, constants, codeSnip);
        
        eval(perturbCodes{i});
        time = 0:dt:(24*validTimes(i));
        [~, y] = ode45(@(t, y) modelEquations(t, y, parameters, ...
            constants, receptors, knockouts, codeSnip),...
            time, y0);

        % Compare outputs
        controlAUC = trapz(control(:, output{i}));
        perturbAUC = trapz(y(:, output{i}));
        percentChangeAct(i) = perturbAUC/controlAUC;
        if percentChangeAct(i) >= (1+threshold)
            prediction{i} = 'Increase';
        elseif percentChangeAct(i) <= (1-threshold)
            prediction{i} = 'Decrease';
        else
            prediction{i} = 'No change';
        end
        
        percentChangeActStr{i} = num2str(percentChangeAct(i));
        if isequal(prediction{i}, measurement{i})
            numMatching = numMatching + 1; match{i} = 1;
        else
            match{i} = 0;
        end
    end
    
    percentAgree = numMatching/length(ids) * 100;
    resultChart = {perturbations, output, measurement, prediction', ...
        percentChangeActStr', match'};
    resultChart = horzcat(resultChart{:});
    header = { 'perturbations' , 'output', 'measurement', 'prediction',...
        'predictedChange', 'match'};
    resultChart = vertcat(header, resultChart);
%     xlswrite('Output.xlsx', resultChart);
end