close all; 
clearvars; clc;

%% Model
[ parameters, constants, receptors, knockouts ] = loadParameters();

dt = 0.05; % step size

% Calibration data from literature
time = 0:dt:(24*30);
calibData = load('setup/calibrationData.mat');
calibData = calibData.c_data;
[c_ltgfb, c_tnfa, c_il1, c_il6, c_mmp9, c_gmcsf, c_collagen, c_macro, ...
    c_fibro] = calibData{:};  

validData = load('setup/validationData.mat');
validData = validData.validation;
[v_ltgfb, v_tnfa, v_il1, v_il6, v_mmp9, v_gmcsf, v_collagen, v_macro, ...
    v_fibro] = validData{:}; 

infarct_size = 1;
codeSnip = '';
[ y0, constants, codeSnip ] = initParams(infarct_size, constants, codeSnip);
[t, y] = ode45(@(t, y) modelEquations(t, y, parameters, constants, ...
    receptors, knockouts, codeSnip), time, y0);
normalY = y;
normMax = max(y, [], 1); %ones(1, 41); % max of each curve

%% Output
% Figure 2A
figure;
%subplot(3, 3, 1);
plot(time, normalY(:, 1)./normMax(1), 'LineWidth', 1.5); hold on; 
errorbar(c_macro{1}, c_macro{2}./max(c_macro{2}), c_macro{3}./...
    max(c_macro{2}), 'rx', 'MarkerSize', 8, 'LineWidth', 1);
errorbar(v_macro{1}, v_macro{2}./max(v_macro{2}), v_macro{3}./...
    max(v_macro{2}), 'ko', 'MarkerSize', 8, 'LineWidth', 1); hold off;
yticks(0:0.2:1.4);
imageLabels(length(time), 'Macrophages', 'Normalized Cell Count', 1);

figure;
%subplot(3, 3, 2);
plot(time, normalY(:, 2)./normMax(2), 'LineWidth', 1.5); hold on; 
errorbar(c_fibro{1}, c_fibro{2}./max(c_fibro{2}), c_fibro{3}./...
    max(c_fibro{2}), 'rx', 'MarkerSize', 8, 'LineWidth', 1);
errorbar(v_fibro{1}, v_fibro{2}./max(v_fibro{2}), v_fibro{3}./...
    max(v_fibro{2}), 'ko', 'MarkerSize', 8, 'LineWidth', 1); hold off;
imageLabels(length(time), 'Fibroblasts', 'Normalized Cell Count', 1);

figure;
%subplot(3, 3, 3);
plot(time, normalY(:, 3)./normMax(3), 'LineWidth', 1.5); hold on; 
errorbar(c_collagen{1}, c_collagen{2}./max(c_collagen{2}), c_collagen{3}./...
    max(c_collagen{2}), 'rx', 'MarkerSize', 8, 'LineWidth', 1); 
errorbar(v_collagen{1}, v_collagen{2}./max(v_collagen{2}), v_collagen{3}./...
    max(v_collagen{2}), 'ko', 'MarkerSize', 8, 'LineWidth', 1); hold off;
imageLabels(length(time), 'Collagen', 'Normalized Concentration', 1);

figure;
%subplot(3, 3, 4);
plot(time, normalY(:, 4)./normMax(4), 'LineWidth', 1.5); hold on; 
errorbar(c_il1{1}, c_il1{2}./max(c_il1{2}), c_il1{3}./max(c_il1{2}), ...
    'rx', 'MarkerSize', 8, 'LineWidth', 1); 
errorbar(v_il1{1}, v_il1{2}./max(v_il1{2}), v_il1{3}./max(v_il1{2}), ...
    'ko', 'MarkerSize', 8, 'LineWidth', 1); hold off;
imageLabels(length(time), 'Secreted IL-1', 'Normalized Concentration', 1);

figure;
%subplot(3, 3, 5);
plot(time, normalY(:, 5)./normMax(5), 'LineWidth', 1.5); hold on; 
errorbar(c_gmcsf{1}, c_gmcsf{2}./max(c_gmcsf{2}), c_gmcsf{3}./...
    max(c_gmcsf{2}), 'rx', 'MarkerSize', 8, 'LineWidth', 1); 
errorbar(v_gmcsf{1}, v_gmcsf{2}./max(v_gmcsf{2}), v_gmcsf{3}./...
    max(v_gmcsf{2}), 'ko', 'MarkerSize', 8, 'LineWidth', 1); hold off;
imageLabels(length(time), 'Secreted GM-CSF', 'Normalized Concentration', 1);

figure;
%subplot(3, 3, 6);
plot(time, normalY(:, 6)./normMax(6), 'LineWidth', 1.5); hold on; 
errorbar(c_ltgfb{1}, c_ltgfb{2}./max(c_ltgfb{2}), c_ltgfb{3}./...
    max(c_ltgfb{2}), 'rx', 'MarkerSize', 8, 'LineWidth', 1); 
errorbar(v_ltgfb{1}, v_ltgfb{2}./max(v_ltgfb{2}), v_ltgfb{3}./...
    max(v_ltgfb{2}), 'ko', 'MarkerSize', 8, 'LineWidth', 1); hold off;
imageLabels(length(time), 'Secreted latent TGFB', 'Normalized Concentration', 1);

figure;
%subplot(3, 3, 7);
plot(time, normalY(:, 8)./normMax(8), 'LineWidth', 1.5); hold on; 
errorbar(c_mmp9{1}, c_mmp9{2}./max(c_mmp9{2}), c_mmp9{3}./max(c_mmp9{2}), ...
    'rx', 'MarkerSize', 8, 'LineWidth', 1); 
errorbar(v_mmp9{1}, v_mmp9{2}./max(v_mmp9{2}), v_mmp9{3}./max(v_mmp9{2}), ...
    'ko', 'MarkerSize', 8, 'LineWidth', 1); hold off;
imageLabels(length(time), 'Secreted MMP-9', 'Normalized Concentration', 1);

figure;
%subplot(3, 3, 8);
plot(time, normalY(:, 34)./normMax(34), 'LineWidth', 1.5); hold on; 
errorbar(c_tnfa{1}, c_tnfa{2}./max(c_tnfa{2}), c_tnfa{3}./max(c_tnfa{2}), ...
    'rx', 'MarkerSize', 8, 'LineWidth', 1); 
errorbar(v_tnfa{1}, v_tnfa{2}./max(v_tnfa{2}), v_tnfa{3}./max(v_tnfa{2}), ...
    'ko', 'MarkerSize', 8, 'LineWidth', 1); hold off;
imageLabels(length(time), 'Secreted TNFa', 'Normalized Concentration', 1);

% Figure 2B
map = [247,252,240
224,243,219
204,235,197
168,221,181
123,204,196
78,179,211
43,140,190
8,104,172
8,64,129]./255;

figure;
maxY = max(y, [], 1);
indOrder = [9, 11, 12, 1, 19, 20, 2, 4, 34, 5, 6, 7, 8, 33, 3, 10, 40, ...
    42, 15:18, 35, 21:23, 36, 24, 37, 41, 26:28, 38, 30:32, 39];
imagesc((y(:, indOrder)./maxY(indOrder))'); colormap(map);
xticklabels({}); yticklabels({});
yticks(1:length(indOrder)); 
yticklabels({'Cardiomyocytes', 'Neutrophils', 'Monocytes', 'Macrophages', ...
    'Macrophage differentiation', 'Macrophage proliferation', ...
    'Fibroblasts', 'IL-1', 'TNFa', 'GM-CSF', 'Latent TGF\beta', 'TGF\beta', ...
    'MMP-9', 'Collagen secretion', 'Collagen', 'Cardiomyocyte debris', ...
    'Neutrophil debris', 'Ingested debris', 'Fibroblast IL-1', 'Fibroblast GM-CSF', ...
    'Fibroblast TGF\beta', 'Fibroblast MMP-9', 'Fibroblast TNFa', ...
    'Macrophage TGF\beta', 'Macrophage IL-1', 'Macrophage MMP-9', 'Macrophage TNFa', ...
    'Cardiomyocyte IL-1', 'Cardiomyocyte TNFa', 'Cardiomyocyte TGF\beta', ...
    'Neutrophil IL-1', 'Neutrophil GM-CSF', 'Neutrophil MMP-9', ...
    'Neutrophil TNFa', 'Monocyte IL-1', 'Monocyte GM-CSF', 'Monocyte MMP-9', ...
    'Monocyte TNFa'});
xticks(0:5*24/dt:length(y)); xticklabels(0:5:30);
c = colorbar; c.Location = 'southoutside';

%% Figure 3
[percentAgree, resultChart] = validation('Validation.xlsx', 0.05);

%% Validation loop for Supplemental figure 1
thresholds = 10.^[-3:0.25:0]; 
numNo = zeros(1, length(thresholds));
percent = zeros(1, length(thresholds));
for i=1:length(thresholds)
    disp(i);
    [percentAgree, resultChart] = validation('Validation.xlsx', thresholds(i));
    percent(1, i) = percentAgree;
    numNo(1, i) = sum(cell2mat(strfind(resultChart(:, 4), 'No change')));
end

figure; hold on;
plot(log10(thresholds), percent, '-ko', 'LineWidth', 1.5);
yline(70, '--b', 'LineWidth', 1);
xticklabels([]); yticklabels([]);
figure; plot(thresholds, numNo, '-ko', 'LineWidth', 1.5);

%% Functions
function [] = imageLabels(xmax, titleName, yName, legFlag)
    title(titleName); xticks(0:120:xmax);
    xticklabels(0:5:30); xlabel('Time (days)');
    ylabel(yName);
    if legFlag == 1
        legend('Simulation', 'Calibration', 'Validation');
    elseif legFlag == 2
        legend('Simulation', 'Experiment');
    end
end