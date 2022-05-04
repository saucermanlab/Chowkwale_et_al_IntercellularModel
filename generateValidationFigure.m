%function [] = generateValidationFigure(exp, sim, perturb, output)
% exp = experimental data
% sim = simulation data

filename = 'Validation.xlsx';
[~, txt, raw] = xlsread(filename);

% 7 - measurement, 8 - prediction
% 2 - perturbations, 4 - output

perturb = txt(2:end, 2); % Perturbations
output = txt(2:end, 4); % Output indices
exp = txt(2:end, 7); % Increase/decrease/no change
sim = txt(2:end, 8); % Increase/decrease/no change

imData = zeros(length(exp), 2);
for i=1:length(exp)
    if strcmp(exp{i}, 'Increase')
        imData(i, 2) = 1;
    elseif strcmp(exp{i}, 'Decrease')
        imData(i, 2) = -1;
    else
        imData(i, 2) = 0;
    end

    if strcmp(sim{i}, "'Increase'")
        imData(i, 1) = 1;
    elseif strcmp(sim{i}, "'Decrease'")
        imData(i, 1) = -1;
    else
        imData(i, 1) = 0;
    end
end

figureRange = 84:86;
figure(1);
%     cmaprange = 0.5:0.01:1;
%     blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
%     myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
myrgbcmap = [27,7,147; 247,247,247; 103,0,31]./255;
colormap(myrgbcmap);
maxVal = 0.25;
caxis([-maxVal, maxVal]);
im = imagesc(imData(figureRange, :), [-maxVal,maxVal]);
set(gca, 'YTick', 1:length(perturb(figureRange)), 'YTickLabel',...
    perturb(figureRange), 'FontSize', 12);
yyaxis right
im = imagesc(imData(figureRange, :),[-maxVal,maxVal]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Simulation', 'Experiment'}, ...
    'FontSize', 12);
set(gca, 'YTick', 1:length(output(figureRange)), 'YTickLabel', ...
    output(figureRange), 'YColor', [0, 0, 0], 'FontSize', 12);

%     image = getframe(gcf);
%     imwrite(image.cdata, myrgbcmap, 'Output.png');
%    % hold on;
%[rows, columns, numChan] = size(im);
%end