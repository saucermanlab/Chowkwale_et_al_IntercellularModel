tic
close all; 
clearvars; clc;

%% Figure 5a
% Model perturbations
p1 = {'parameters(22) = 0.2*parameters(22);' 'parameters(22) = 0.5*parameters(22);',...
    '', 'parameters(22) = 1.5*parameters(22);', 'parameters(22) = 2*parameters(22);'};
l1 = {'0.2*IL-1 degradation', '0.5*IL-1 degradation', 'Baseline', ...
    '1.5*IL-1 degradation', '2*IL-1 degradation'};
p = p1'; l = l1';

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
lineStyle = {'-', '--', '-', ':', '-.'};
normIL1 = squeeze(data(:, :, 4))';
figure; hold on;
for i=1:size(normIL1, 2)
    if i == 3 color = 'k'; 
    else color = '#737373'; end
    plot(normIL1(:, i), lineStyle{i}, 'Color', color, 'LineWidth', 1.5);
    xticks(0:2400:14400); xticklabels(split(num2str(0:5:30)));
end
title('IL-1'); ylabel('Concentration'); xlabel('Time'); legend(l);

% T50 
opIndex = [4]; opData = data(:, :, opIndex); labels = {'IL-1'};
dMax = findT50(opData, dt, 0.5);

figure; 
lCat = categorical(labels); lCat = reordercats(lCat, labels);
b1 = bar(lCat, dMax', 'FaceColor', 'flat'); legend(l); 
b1(5).FaceColor = [127,205,187]./255;
b1(4).FaceColor = [65,182,196]./255;
b1(3).FaceColor = [29,145,192]./255;
b1(2).FaceColor = [34,94,168]./255;
b1(1).FaceColor = [37,52,148]./255;
title('T50'); ylabel('T_{50} in days'); xlabel('Output measured');

%% Figure 5b
p1 = {'parameters(6) = 0.2*parameters(6);' 'parameters(6) = 0.5*parameters(6);',...
    '', 'parameters(6) = 1.5*parameters(6);', 'parameters(6) = 2*parameters(6);'};
l1 = {'0.2*Neutrophil rem', '0.5*Neutrophil rem', 'Baseline', ...
    '1.5*Neutrophil rem', '2*Neutrophil rem'};
p = p1'; l = l1';

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
lineStyle = {'-', '--', '-', ':', '-.'};
normIL1 = squeeze(data(:, :, 4))';
figure; hold on;
for i=1:size(normIL1, 2)
    if i == 3 color = 'k'; 
    else color = '#737373'; end
    plot(normIL1(:, i), lineStyle{i}, 'Color', color, 'LineWidth', 1.5);
    xticks(0:2400:14400); xticklabels(split(num2str(0:5:30)));
end
% yticks(0:10:50); yticklabels(split(num2str(0:10:50)));
title('IL-1'); ylabel('Concentration'); xlabel('Time'); legend(l);

% T50 
opIndex = [4]; opData = data(:, :, opIndex); labels = {'IL-1'};
dMax = findT50(opData, dt, 0.5);

figure; 
lCat = categorical(labels); lCat = reordercats(lCat, labels);
b1 = bar(lCat, dMax', 'FaceColor', 'flat'); legend(l); 
b1(5).FaceColor = [127,205,187]./255;
b1(4).FaceColor = [65,182,196]./255;
b1(3).FaceColor = [29,145,192]./255;
b1(2).FaceColor = [34,94,168]./255;
b1(1).FaceColor = [37,52,148]./255;
title('T50'); ylabel('T_{50} in days'); xlabel('Output measured');

%% Figure 5c
p1 = {'parameters(37) = 0.05*parameters(37);' 'parameters(37) = 0.5*parameters(37);',...
    '', 'parameters(37) = 5*parameters(37);', 'parameters(37) = 50*parameters(37);'};
l1 = {'0.05*MMP-9 deg', '0.5*MMP-9 deg', 'Baseline', '5*MMP-9 deg', '50*MMP-9 deg'};
p = p1'; l = l1';

% Model simulations
dt = 0.05; % step size
time = 0:dt:(24*60);
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
lineStyle = {'-', '--', '-', ':', '-.'};
normIL1 = squeeze(data(:, :, 4))';
figure; hold on;
for i=1:size(normIL1, 2)
    if i == 3 color = 'k'; 
    else color = '#737373'; end
    plot(normIL1(:, i), lineStyle{i}, 'Color', color, 'LineWidth', 1.5);
    xticks(0:4800:28800); xticklabels(split(num2str(0:10:60)));
end
yticks(0:500:2000); yticklabels(split(num2str(0:500:2000)));
title('IL-1'); ylabel('Concentration'); xlabel('Time'); legend(l);

% T90 
opIndex = [4]; opData = data(:, :, opIndex); labels = {'IL-1'};
dMax = findT50(opData, dt, 0.1);

figure; 
lCat = categorical(labels); lCat = reordercats(lCat, labels);
b1 = bar(lCat, dMax', 'FaceColor', 'flat'); legend(l); 
b1(5).FaceColor = [127,205,187]./255;
b1(4).FaceColor = [65,182,196]./255;
b1(3).FaceColor = [29,145,192]./255;
b1(2).FaceColor = [34,94,168]./255;
b1(1).FaceColor = [37,52,148]./255;
title('T90'); ylabel('T_{90} in days'); xlabel('Output measured');

%% Figure 5d
p = {'', 'parameters(32) = 50*parameters(32);'};
l = {'Baseline', 'TGFB inhibition'};

dt = 0.05; % step size
time = 0:dt:(24*60);
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
normIL1 = squeeze(data(:, :, 4))';
figure; hold on;
for i=1:size(normIL1, 2)
    if i == 1 color = 'k'; 
    else color = '#737373'; end
    plot(normIL1(:, i), 'Color', color, 'LineWidth', 1.5);
    xticks(0:4800:28800); xticklabels(split(num2str(0:10:60)));
end
yticks(0:300:1800); yticklabels(split(num2str(0:300:1800)));
title('IL-1'); ylabel('Concentration'); xlabel('Time'); legend(l);

% T90 
opIndex = [4]; opData = data(:, :, opIndex); labels = {'IL-1'};
dMax = findT50(opData, dt, 0.1);

figure; 
lCat = categorical(labels); lCat = reordercats(lCat, labels);
b1 = bar(lCat, dMax', 'FaceColor', 'flat'); legend(l); ylim([0 50]);
yticks(0:10:50);
b1(2).FaceColor = [34,94,168]./255;
b1(1).FaceColor = [37,52,148]./255;
title('T90'); ylabel('T_{90} in days'); xlabel('Output measured');

% Output
norm = squeeze(data(2, 1:2400, [4, 15, 22, 24, 26, 30]));
cellData = squeeze(trapz(norm, 1))'; % days 0-5
norm = squeeze(data(2, 2401:4800, [4, 15, 22, 24, 26, 30]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 5-10
norm = squeeze(data(2, 4801:7200, [4, 15, 22, 24, 26, 30]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 10-15
norm = squeeze(data(2, 7201:9600, [4, 15, 22, 24, 26, 30]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 15-20
norm = squeeze(data(2, 9601:12000, [4, 15, 22, 24, 26, 30]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 20-25
norm = squeeze(data(2, 12001:14401, [4, 15, 22, 24, 26, 30]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 25-30

totalAUC = cellData(1, :);
totalAUC = totalAUC.*100./sum(totalAUC);
cellData = cellData./cellData(1, :); %total, fibro, macro, cm, n, mo

figure;
b = bar(cellData(2:6, :)', 'stacked'); ylim([0 1.2]);
yticks(0:0.2:1.2);
b(5).FaceColor = [127,205,187]./255;
b(4).FaceColor = [65,182,196]./255;
b(3).FaceColor = [29,145,192]./255;
b(2).FaceColor = [34,94,168]./255;
b(1).FaceColor = [37,52,148]./255;
legend({'Fibroblasts', 'Macrophages', 'Cardiomyocytes', 'Neutrophils', 'Monocytes'});
title('Fraction of infarct IL-1 secreted by cells');
xticklabels({'Days 0-5', 'Days 5-10', 'Days 10-15', 'Days 15-20', ...
    'Days 20-25', 'Days 25-30'}); xtickangle(45);

%% Functions
function [ t50 ] = findT50(data, dt, fraction)
    t50 = zeros(size(data, 1), size(data, 3));
    [maxVal, maxTime] = max(data, [], 2);
    maxVal = squeeze(maxVal); maxTime = squeeze(maxTime);
    halfVal = fraction*maxVal;
    for i=1:size(data, 1)
        for j=1:size(data, 3)
            t = find(data(i, :, j) < halfVal(i, j), length(data), 'last');
            t50(i, j) = t(find(t > maxTime(i, j), 1))*dt/24; % in hours
        end
    end
end