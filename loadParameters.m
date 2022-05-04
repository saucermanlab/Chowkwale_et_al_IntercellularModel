%% Parameters
% Last updated: 10/25/21
function [ parameters, constants, receptors, knockouts ] = loadParameters()
    
    parameters = []; 
    
    % Macrophage parameters
    parameters(1) = 0.0592; % differentiation rate from monocytes 
    parameters(2) = 2.3e-04; % macrophage proliferation rate
    parameters(3) = 0.0175; % removal rate of macrophages 
    parameters(4) = 0.5*0.03; % macrophage phagocytosis rate
    
    % Neutrophil parameters
    parameters(5) = 0.095; % % neutrophil recruitment rate
    parameters(6) = 0.525; % neutrophil removal rate
    parameters(7) = 0.0022; % neutrophil phagocytosis rate
    
    % Monocyte parameters
    parameters(8) = 0.0025; % % monocyte recruitment rate
    parameters(9) = 0.01; % monocytes removal rate
    parameters(10) = 0.00001; % monocyte proliferation in blood
   
    % Other cell parameters
    parameters(11) = 0.065; 0.15; % proliferation rate of fibroblasts
    parameters(12) = 0.025; % death rate of cardiomyocytes 
    parameters(13) = 4*0.25; % debris conversion

    % Collagen parameters
    parameters(14) = 0.095; % maturation rate of collagen
    parameters(15) = 0.0025; % collagen secretion by fibroblasts 
    parameters(16) = 4.85e-03; % collagen degradation by MMP-9
   
    % IL-1 parameters
    parameters(17) = 7.55e-5; % IL-1 secretion by cardiomyocytes
    parameters(18) = 0.06; % IL-1 secretion rate by neutrophils
    parameters(19) = 0.95854e-2; % IL-1 secretion rate by monocytes
    parameters(20) = 1.554e-2; % IL-1 secretion rate by macrophages
    parameters(21) = 0.225e-2; % IL-1 secretion rate by fibroblasts
    parameters(22) = 0.4086; % Degradation rate of IL-1
    
    % GM-CSF parameters
    parameters(23) = 0.96e-2; % GM-CSF secretion by neutrophils
    parameters(24) = 1.95e-3; % GM-CSF secretion by monocytes
    parameters(25) = 0.524e-3; % GM-CSF secretion by fibroblasts
    parameters(26) = 3.15; % degradation rate of GM-CSF
     
    % Latent TGFB parameters
    parameters(27) = 0.01e-5; % latent TGFB secretion by cardiomyocytes
    parameters(28) = 0.64285e-2; % latent TGFB secretion by macrophages
    parameters(29) = 2.95e-4; % latent TGFB secretion by fibroblasts
    parameters(30) = 0.0075; % degradation rate of latent TGFB
    
    % TGFB parameters
    parameters(31) = 1e-2; % activation rate of TGFB
    parameters(32) = 0.0769; % degradation rate of TGFB 
    
    % MMP-9 parameters
    parameters(33) = 7.35e-4; % MMP-9 secretion by neutrophils
    parameters(34) = 0.000065; % MMP-9 secretion by monocytes
    parameters(35) = 0.00003; % MMP-9 secretion by macrophages
    parameters(36) = 0.000065; % MMP-9 secretion by fibroblasts
    parameters(37) = 0.292; % degradation rate of MMP-9
    
    % TNFa parameters
    parameters(38) = 0.05e-5; % TNFa secretion rate by cardiomyocytes
    parameters(39) = 0.0012; % TNFa secretion rate by neutrophils
    parameters(40) = 1.5854e-4; % TNFa secretion rate by monocytes
    parameters(41) = 3.514e-4; % TNFa secretion rate by macrophages
    parameters(42) = 0.95e-4; % TNFa secretion rate by fibroblasts
    parameters(43) = 0.4786; % degradation rate of TNFa
    
    parameters(44) = 0.0015; % removal rate for fibroblasts
        
    % Constants
    constants = [];
    constants(1) = 100; % initial fibroblast density
    constants(2) = 7000; % maximum fibroblast density
    constants(3) = 5*2; % saturation constant for TNFa, 
    constants(4) = 10; % saturation constant for GM-CSF
    constants(5) = 3; % saturation constant for TGFB
    constants(6) = 10; % saturation constant for MMP-9
    constants(7) = 400; % saturation constant for IL-1
    constants(8) = 106601; % maximum macrophage density
    constants(9) = 5e2; % maximum collagen density
    constants(10) = 1.0714e+07; % maximum cardiomyocytes
    constants(11) = 3e4; % maximum neutrophils
    constants(12) = 3e4; % maximum monocytes
    constants(13) = 1; % infarct
    constants(14) = 1; % infarct size

    % Receptors
    receptors = [];
    receptors(1) = 1; % GM-CSF receptors
    receptors(2) = 1; % IL-1 receptors
    receptors(3) = 1; % TGFb receptors
    receptors(4) = 1; % TNFa receptors

    % Knockouts for secreted factors
    knockouts = ones(1, 6);
end