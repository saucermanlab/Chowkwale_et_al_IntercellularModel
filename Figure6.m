tic
close all; 
clearvars; clc;

%% Figure 6a-d
% Model perturbations
p1 = {'', 'parameters(22) = 50*parameters(22);', 'parameters(43) = 50*parameters(43);' ...
    'parameters(26) = 50*parameters(26);'};
l1 = {'Baseline', 'IL-1 inhibition', 'TNFa inhibition', 'GM-CSF inhibition'};
p = p1';
l = l1';

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

%% Output
% Figure 6a
norm = squeeze(data(1, 1:2400, [6, 17, 21, 41]));
cellData = squeeze(trapz(norm, 1))'; % days 0-5
norm = squeeze(data(2, 2401:4800, [6, 17, 21, 41]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 5-10
norm = squeeze(data(2, 4801:7200, [6, 17, 21, 41]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 10-15
norm = squeeze(data(2, 7201:9600, [6, 17, 21, 41]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 15-20
norm = squeeze(data(2, 9601:12000, [6, 17, 21, 41]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 20-25
norm = squeeze(data(2, 12001:14401, [6, 17, 21, 41]));
cellData = [cellData squeeze(trapz(norm, 1))']; % days 25-30

totalAUC = cellData(1, :);
totalAUC = totalAUC.*100./sum(totalAUC);
cellData = cellData./sum(cellData(2:4, :)); %total, fibro, macro, cm

figure;
b = bar(cellData(2:4, :)', 'stacked'); ylim([0 1.2]);
yticks(0:0.2:1.2);
b(3).FaceColor = [127,205,187]./255;
b(2).FaceColor = [29,145,192]./255;
b(1).FaceColor = [37,52,148]./255;
% b(3).FaceColor = [29,145,192]./255;
% b(2).FaceColor = [34,94,168]./255;
% b(1).FaceColor = [37,52,148]./255;
legend({'Fibroblasts', 'Macrophages', 'Cardiomyocytes'});
title('Fraction of infarct TGFB secreted by cells');
xticklabels({'Days 0-5', 'Days 5-10', 'Days 10-15', 'Days 15-20', ...
    'Days 20-25', 'Days 25-30'}); xtickangle(45);

% Figure 6b
% opIndex = [6, 7, 1, 9, 21, 41]; opData = data(:, :, opIndex);
% labels = {'Latent TGFB', 'TGFB', 'Macrophages', ...
%     'Cardiomyocytes', 'Macrophage TGFB', 'CM TGFB'};
opIndex = [6, 7, 17, 21, 41]; opData = data(:, :, opIndex);
labels = {'Latent TGFB', 'TGFB', 'Fibroblast TGFB', 'Macrophage TGFB', 'CM TGFB'};
dMax = squeeze(trapz(opData, 2)); dMax = dMax./dMax(1,:);

figure; 
lCat = categorical(labels); lCat = reordercats(lCat, labels); 
b1 = bar(lCat, dMax, 'FaceColor', 'flat'); legend(l);
b1(4).FaceColor = [199,233,180]./255;
b1(3).FaceColor = [127,205,187]./255;
b1(2).FaceColor = [29,145,192]./255;
b1(1).FaceColor = [37,52,148]./255;
title('Effect of inflammatory cytokines on TGFB');
ylabel('AUC relative to baseline'); xlabel('Output measured');

% Figure 6c
opIndex = [1, 2, 10, 40, 42]; opData = data(:, :, opIndex);
dMax = squeeze(max(opData, [], 2)); 
dMax = dMax./dMax(1,:);
labels = {'Macrophages', 'Fibroblasts', 'Cardiomyocyte debris', 'Neutrophil debris', 'Ingested debris'};

figure; 
lCat = categorical(labels); lCat = reordercats(lCat, labels); 
b1 = bar(lCat, dMax', 'FaceColor', 'flat'); ylim([0 1.6]);
yticks(0:0.4:1.6);
legend(l);
b1(4).FaceColor = [199,233,180]./255;
b1(3).FaceColor = [127,205,187]./255;
b1(2).FaceColor = [29,145,192]./255;
b1(1).FaceColor = [37,52,148]./255;
title('Effect of inflammatory cytokines on phagocytosis');
ylabel('Peak values relative to baseline'); xlabel('Output measured');

%% Figure 6d
% Model perturbations
p1 = {'', 'mp2 = 0; mp1 = 0;', 'mp1 = 0; mp5 = 1;', 'mp2 = 0;', 'mp3 = 0;', 'mp4 = 0;'}; pInd = 1:length(p1);
l1 = {'Baseline', 'Phagocytosis inhibition', 'TNFa-mediated phagocytosis inhibition', ...
    'TGFB-mediated phagocytosis inhibition', 'Cardiomyocyte phagocytosis inhibition', ...
    'Neutrophil phagocytosis inhibition'};
p = p1';
l = l1';

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
opIndex = [6, 7, 21]; opData = data(:, :, opIndex);
labels = {'Latent TGFB', 'TGFB', 'Macrophage TGFB'};
dMax = squeeze(trapz(opData, 2)); dMax = dMax./dMax(1,:);

figure; 
lCat = categorical(labels); lCat = reordercats(lCat, labels); 
b1 = bar(lCat, dMax, 'FaceColor', 'flat'); legend(l);
b1(1).CData = [34,94,168]./255; b1(2).CData = [8,29,88]./255;
b1(5).FaceColor = [127,205,187]./255;
b1(4).FaceColor = [65,182,196]./255;
b1(3).FaceColor = [29,145,192]./255;
b1(2).FaceColor = [34,94,168]./255;
b1(1).FaceColor = [37,52,148]./255;
title('Effect of phagocytosis inhibition on TGFB');
ylabel('AUC relative to baseline'); xlabel('Output measured');

%%
toc