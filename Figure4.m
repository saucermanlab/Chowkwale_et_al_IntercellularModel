tic
close all; 
clearvars; clc;

%% Figure 4b
% Model perturbations
p = {''}; l = {'Baseline'};

% Model simulations
dt = 0.05; % step size
time = 0:dt:(24*30);
[ parameters, constants, receptors, knockouts ] = loadParameters();
data = zeros(length(p), length(time), 43);

for i=1:length(p)
    disp(i);
    [ parameters, constants, receptors, knockouts ] = loadParameters();
    [ y0, constants, codeSnip ] = initParams(1, constants, p{i});
    [t, y] = ode45(@(t, y) modelEquations(t, y, parameters, constants, ...
        receptors, knockouts, codeSnip), time, y0);
    normalY = y;
    data(i, :, :) = y;
end

% Output [4, 15, 22, 24, 26, 30]
norm = squeeze(data(1, 1:2400, [4, 15, 22, 24, 26, 30]));
cellData = squeeze(trapz(norm, 1))'; % days 0-5
norm = squeeze(data(1, 2401:4800, [4, 15, 22, 24, 26, 30]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 5-10
norm = squeeze(data(1, 4801:7200, [4, 15, 22, 24, 26, 30]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 10-15
norm = squeeze(data(1, 7201:9600, [4, 15, 22, 24, 26, 30]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 15-20
norm = squeeze(data(1, 9601:12000, [4, 15, 22, 24, 26, 30]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 20-25
norm = squeeze(data(1, 12001:14401, [4, 15, 22, 24, 26, 30]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 25-30
cellData = cellData./cellData(1, :); %total, fibro, macro, cm, n, mo

figure;
b = bar(cellData(2:6, :)', 'stacked'); ylim([0 1.2]);
b(5).FaceColor = [127,205,187]./255;
b(4).FaceColor = [65,182,196]./255;
b(3).FaceColor = [29,145,192]./255;
b(2).FaceColor = [34,94,168]./255;
b(1).FaceColor = [37,52,148]./255;
legend({'Fibroblasts', 'Macrophages', 'Cardiomyocytes', 'Neutrophils', 'Monocytes'});
title('Fraction of infarct IL-1 secreted by cells');
xticklabels({'Days 0-5', 'Days 5-10', 'Days 10-15', 'Days 15-20', ...
    'Days 20-25', 'Days 25-30'}); xtickangle(45);

totalAUC = cellData(1, :);
totalAUC = totalAUC.*100./sum(totalAUC);

%% Figure 4a
[ parameters, constants, receptors, knockouts ] = loadParameters();
infSize = [0:0.05:1]; %[0:0.002:0.01, 0.02:0.02:0.1, 0.2:0.2:1];

dt = 0.05; % step size
time = 0:dt:(24*30);
data = zeros(length(infSize), length(time), 43);

for j=1:length(infSize)
    disp(j);
    infarct_size = infSize(j);
    [ parameters, constants, receptors, knockouts ] = loadParameters();
    [ y0, constants, codeSnip ] = initParams(infarct_size, constants, '');
    [t, y] = ode45(@(t, y) modelEquations(t, y, parameters, constants, ...
        receptors, knockouts, codeSnip), time, y0);
    normalY = y;
    data(j, :, :) = y; 
end

% Output
ind = [1, 6, 11, 16, 21]; 
lineStyle = {'-', '--', ':', '-.', '-'};
opIndex = 4; opData = data(:, :, opIndex);
l2 = "IS = " + split(num2str(infSize));
figure; hold on;
for i=1:size(ind, 2)
    if i == 5 color = 'k'; 
    else color = '#737373'; end
    plot(opData(ind(i), 1:6721), lineStyle{i}, 'Color', color, 'LineWidth', 1.5);
end
xticks(0:960:6720); xticklabels(split(num2str(0:2:14)));
yticks(0:360:1800); yticklabels(split(num2str(0:360:1800)));
title('IL-1'); xlabel('Time'); ylabel('Concentration'); legend(l2([1, 6, 11, 16, 21]));

peakVals = max(opData, [], 2);
figure; plot(infSize, peakVals, '-ko', 'LineWidth', 1.5);
yticks(0:360:1800); yticklabels(split(num2str(0:360:1800)));
xlabel('Infarct size'); ylabel('Peak values');
title('IL-1 peaks vs infarct size');

% Figure 4c
opIndex = 11; opData = data(:, :, opIndex);
l2 = "IS = " + split(num2str(infSize));

figure; hold on;
for i=1:size(ind, 2)
    if i == 5 color = 'k'; 
    else color = '#737373'; end
    plot(opData(ind(i), 1:6721), lineStyle{i}, 'Color', color, 'LineWidth', 1.5);
end
xticks(0:960:6720); xticklabels(split(num2str(0:2:14)));
title('Neutrophils'); xlabel('Time'); ylabel('Cell Counts'); legend(l2([1, 6, 11, 16, 21]));

peakVals = max(opData, [], 2);
figure; plot(infSize, peakVals, '-ko', 'LineWidth', 1.5);
% yticks(0:360:1800); yticklabels(split(num2str(0:360:1800)));
xlabel('IS'); ylabel('Peak values');
title('Neutrophils peaks vs infarct size');

%% Figure 4d
p = {'', 'parameters(22) = 50*parameters(22);', 'parameters(6) = 25*parameters(6);'}; 
l = {'Baseline', 'IL-1 deg', 'Decreased Ns'};

% Model simulations
dt = 0.05; % step size
time = 0:dt:(24*30);
[ parameters, constants, receptors, knockouts ] = loadParameters();
data = zeros(length(p), length(time), 43);

for i=1:length(p)
    disp(i);
    [ parameters, constants, receptors, knockouts ] = loadParameters();
    [ y0, constants, codeSnip ] = initParams(1, constants, p{i});
    [t, y] = ode45(@(t, y) modelEquations(t, y, parameters, constants, ...
        receptors, knockouts, codeSnip), time, y0);
    normalY = y;
    data(i, :, :) = y;
end

% Output
opIndex = [11, 4];
opData = data(:, :, opIndex);
labels = {'Neutrophils', 'IL-1'};
dMax = squeeze(max(opData, [], 2)); dMax = dMax./dMax(1, :);

figure; 
lCat = categorical(labels); lCat = reordercats(lCat, labels);
b1 = bar(lCat, dMax', 1, 'FaceColor', 'flat'); ylim([0 1.2]);
b1(1).FaceColor = '#1d91c0';
b1(2).FaceColor = '#225ea8';
b1(3).FaceColor = '#081d58';
title('Neutrophil-IL1 feedback loop'); legend(l);
ylabel('Peak values relative to baseline'); xlabel('Output measured');

%%
toc