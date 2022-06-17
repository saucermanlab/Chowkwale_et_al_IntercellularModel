close all; 
clearvars; clc;

%% Figure 7A and 7C
p = {''};
[ parameters, constants, receptors, knockouts ] = loadParameters();
infSize = [0:0.01:1];

dt = 0.05; % step size
time = 0:dt:(24*30);
data = zeros(length(p)*length(infSize), length(time), 43);
crowdingData = zeros(length(p)*length(infSize), length(time), 8);

for j=1:length(infSize)
    disp(j);
    infarct_size = infSize(j);
    [ parameters, constants, receptors, knockouts ] = loadParameters();
    [ y0, constants, codeSnip ] = initParams(infarct_size, constants, '');
    [t, y, crowding] = ode45(@(t, y) modelEquations(t, y, parameters, constants, ...
        receptors, knockouts, codeSnip), time, y0);
    normalY = y;
    data(j, :, :) = y; 
    
    for tstep=1:length(t)
        [~, x(tstep, :)] = modelEquations(t(tstep), y(tstep, :), ...
            parameters, constants, receptors, knockouts, codeSnip);
    end
    crowdingData(j, :, :) = x;
end

% Figure 7A
opIndex = 3; opData = data(:, :, opIndex);
peakVals = max(opData, [], 2);
figure; plot(infSize(1:5:end), peakVals(1:5:end), '-ko', 'LineWidth', 1.5); hold on;
xlabel('Infarct size'); ylabel('Peak values');
title('Collagen peaks vs infarct size');

% Fit to Hill function
% max(peakVals); = 45.0527
% k = find(peakVals >= (0.5*45.0527), 1); ec50 = infSize(k); % = 0.73
hill_func = fittype('45.0527*((a^n)/((0.73^n) + (a^n)))', 'dependent', {'e'}, ...
    'independent', {'a'}, 'coefficients', {'n'});
soln = fit(infSize', peakVals, hill_func);
plot(infSize, soln(infSize), 'r', 'LineWidth', 1.5); hold off;
% n = 9.569 Positively cooperative binding
yticks(0:10:50); yticklabels(split(num2str(0:10:50)));
xlabel('Infarct size'); ylabel('Peak values');
title('Collagen peaks vs infarct size');

% Figure 7C
opIndex = 2; opData = data(:, :, opIndex);
peakValsFib = max(opData, [], 2);
figure; plot(infSize(1:5:end), peakValsFib(1:5:end)./max(peakValsFib(1:5:end)), '-ko', 'LineWidth', 1.5); hold on;
yticks(0:0.2:1); yticklabels(split(num2str(0:0.2:1)));

opIndex = 7; opData = data(:, :, opIndex);
peakValsFib = max(opData, [], 2);
plot(infSize(1:5:end), peakValsFib(1:5:end)./max(peakValsFib(1:5:end)), '--ko', 'LineWidth', 1.5);
yticks(0:0.2:1); yticklabels(split(num2str(0:0.2:1)));

opIndex = 33; opData = data(:, :, opIndex);
peakValsFib = max(opData, [], 2);
plot(infSize(1:5:end), peakValsFib(1:5:end)./max(peakValsFib(1:5:end)), ':ko', 'LineWidth', 1.5); hold off;
yticks(0:0.2:1); yticklabels(split(num2str(0:0.2:1)));
legend({'Fibroblasts', 'TGFB', 'Collagen secretion'}); 
xlabel('Infarct size'); ylabel('Peak values');
title('Normalized peak values vs infarct size');

%% Figure 7B
p1 = {'', 'parameters(22) = 50*parameters(22);', 'parameters(43) = 50*parameters(43);' ...
    'parameters(26) = 50*parameters(26);'};
l1 = {'Baseline', 'IL-1 inhibition', 'TNFa inhibition', 'GM-CSF inhibition'};
p = p1'; l = l1';
[ parameters, constants, receptors, knockouts ] = loadParameters();
infSize = [0:0.05:1];

dt = 0.05; % step size
time = 0:dt:(24*30);
data = zeros(length(p)*length(infSize), length(time), 43);

for j=1:length(infSize)
    for i=1:length(p)
        disp(j); disp(i);
        infarct_size = infSize(j);
        [ parameters, constants, receptors, knockouts ] = loadParameters();
        [ y0, constants, codeSnip ] = initParams(infarct_size, constants, p{i});
        [t, y] = ode45(@(t, y) modelEquations(t, y, parameters, constants, ...
            receptors, knockouts, codeSnip), time, y0);
        normalY = y;
        data(length(p)*(j-1) + i, :, :) = y; 
    end
end

% Figure 7B 
opIndex = 3; opData = data(:, :, opIndex);
peakVals = max(opData, [], 2);
figure; plot(infSize, peakVals(1:4:end), '-ko', 'LineWidth', 1.5); hold on;
plot(infSize, peakVals(2:4:end), '-o', infSize, peakVals(3:4:end), '--o', ...
    infSize, peakVals(4:4:end), ':o', 'Color', '#737373', 'LineWidth', 1.5); hold off;
yticks(0:10:50); yticklabels(split(num2str(0:10:50)));
xlabel('Infarct size'); ylabel('Peak values');
title('Collagen peaks vs infarct size'); legend(l);

%% Figure 7D - Bifurcation diagram with TGFB as parameter
t = [0:1:20]; f_max = 7000; lam = 0.15; d = 0.0035; sat = 3;
inf_size = 1;
fixed_points = zeros(2, length(t));
stability = zeros(2, length(t));

for i=1:length(t)
    % equation --> dFdt = lam*tgfb*F*(1 - F/(f_max*inf_size)) - dF
    % (lam*tgfb - d)F - (lam*tgfb/(f_max*inf_size))*(F^2)
    tgfb = t(i)./(t(i)+sat); 
    eq = [-(lam*tgfb/(f_max*inf_size)), (lam*tgfb - d), 0];
    %eq(isinf(eq)) = 0;
    r = roots(eq);
    fixed_points(:, i) = r';

    % find stability of each point, 1 if stable and 0 if unstable
    % derivative
    syms f
    for i=1:length(t)
        tgfb = t(i)./(t(i)+sat);
        dfdt = lam*tgfb*f*(1 - (f/(inf_size*f_max))) - d*f;
        g = diff(dfdt, f);
        stability(1, i) = vpa(subs(g,f,fixed_points(1, i)));
        stability(2, i) = vpa(subs(g,f,fixed_points(2, i)));
    end
end

% negative = stable, positive = unstable
stability(stability > 0) = 0;
stability(stability < 0) = 1;

figure; hold on;
for i=1:length(t)
    for j=1:2
        p = plot(t(i), fixed_points(j, i), '-ko');
        if stability(j, i) == 1 % stable
            p.MarkerFaceColor = [0 0 0];
            p.MarkerSize = 7;
            p.MarkerEdgeColor = [0 0 0];
        end
    end
end
plot(t, fixed_points(1, :), '-k', 'LineWidth', 1.5);
plot(t, fixed_points(2, :), '-k', 'LineWidth', 1.5);
title('Transcritical bifurcation');

%% Figure 7D - Inset
t = [0:0.01:0.1, 0.2:0.1:1]; f_max = 7000; lam = 0.15; d = 0.0035; sat = 3;
inf_size = 1;
fixed_points = zeros(2, length(t));
stability = zeros(2, length(t));

for i=1:length(t)
    % equation --> dFdt = lam*tgfb*F*(1 - F/(f_max*inf_size)) - dF
    % (lam*tgfb - d)F - (lam*tgfb/(f_max*inf_size))*(F^2)
    tgfb = t(i)./(t(i)+sat); 
    eq = [-(lam*tgfb/(f_max*inf_size)), (lam*tgfb - d), 0];
    %eq(isinf(eq)) = 0;
    r = roots(eq);
    fixed_points(:, i) = r';

    % find stability of each point, 1 if stable and 0 if unstable
    % derivative
    syms f
    for i=1:length(t)
        tgfb = t(i)./(t(i)+sat);
        dfdt = lam*tgfb*f*(1 - (f/(inf_size*f_max))) - d*f;
        g = diff(dfdt, f);
        stability(1, i) = vpa(subs(g,f,fixed_points(1, i)));
        stability(2, i) = vpa(subs(g,f,fixed_points(2, i)));
    end
end

stability(stability > 0) = 0;
stability(stability < 0) = 1;

figure; hold on;
for i=1:length(t)
    for j=1:2
        p = plot(t(i), fixed_points(j, i), '-ko');
        if stability(j, i) == 1
            p.MarkerFaceColor = [0 0 0];
            p.MarkerSize = 7;
            p.MarkerEdgeColor = [0 0 0];
        end
    end
end
plot(t, fixed_points(1, :), '-k', 'LineWidth', 1.5);
plot(t, fixed_points(2, :), '-k', 'LineWidth', 1.5);
title('Transcritical bifurcation');

%% Figure 7E Bifurcation diagram - removing F
t = [0:1:20]; f_max = 7000; lam = 0.15; d = 0.0035; sat = 3;
inf_size = 1;
fixed_points = zeros(1, length(t));
stability = zeros(1, length(t));

for i=1:length(t)
    % equation --> dFdt = lam*tgfb*F*(1 - F/(f_max*inf_size)) - dF
    % (lam*tgfb - d)F - (lam*tgfb/(f_max*inf_size))*(F^2)
    tgfb = t(i)./(t(i)+sat); 
    eq = [-((lam*tgfb/(f_max*inf_size)) + d), (lam*tgfb)];
    eq(isinf(eq)) = 0;
    r = roots(eq);
    fixed_points(:, i) = r';

    % find stability of each point, 1 if stable and 0 if unstable
    % derivative
    syms f
    for i=1:length(t)
        tgfb = t(i)./(t(i)+sat);
        dfdt = lam*tgfb*(1 - (f/(inf_size*f_max))) - d*f;
        g = diff(dfdt);
        stability(1, i) = vpa(subs(g,f,fixed_points(1, i)));
    end
end

stability(stability > 0) = 0;
stability(stability < 0) = 1;

figure; hold on;
for i=1:length(t)
    p = plot(t(i), fixed_points(1, i), '-ko');
    if stability(1, i) == 1
        p.MarkerFaceColor = [0 0 0];
        p.MarkerSize = 6;
        p.MarkerEdgeColor = [0 0 0];
    end
end
plot(t, fixed_points(1, :), '-k', 'LineWidth', 1.5);
xticks(0:4:20);
title('Transcritical bifurcation');

%% Figure 7F Bifurcation diagram - removing crowding
t = [0:1:20]; f_max = 7000; lam = 0.15; d = 0.0035; sat = 3;
inf_size = 1;
fixed_points = zeros(1, length(t));
stability = zeros(1, length(t));

for i=1:length(t)
    % equation --> dFdt = lam*tgfb*F*(1 - F/(f_max*inf_size)) - dF
    % (lam*tgfb - d)F - (lam*tgfb/(f_max*inf_size))*(F^2)
    tgfb = t(i)./(t(i)+sat); 
    eq = [lam*tgfb - d, 0];
    eq(isinf(eq)) = 0;
    r = roots(eq);
    fixed_points(:, i) = r';

    % find stability of each point, 1 if stable and 0 if unstable
    % derivative
    syms f
    for i=1:length(t)
        tgfb = t(i)./(t(i)+sat);
        dfdt = lam*tgfb*f - d*f;
        g = diff(dfdt);
        stability(1, i) = vpa(subs(g,f,fixed_points(1, i)));
    end
end

stability(stability > 0) = 0;
stability(stability < 0) = 1;

figure; hold on;
for i=1:length(t)
    p = plot(t(i), fixed_points(1, i), '-ko');
    if stability(1, i) == 1
        p.MarkerFaceColor = [0 0 0];
        p.MarkerSize = 6;
        p.MarkerEdgeColor = [0 0 0];
    end
end
plot(t, fixed_points(1, :), '-k', 'LineWidth', 1.5);
xticks(0:4:20); yticks(-1:0.5:1);
title('Transcritical bifurcation');