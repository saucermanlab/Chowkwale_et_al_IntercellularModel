%% Digitization of calibration and validation datasets
% Using webplotdigitizer

%% Calibration lTGFB - PMID 15537506, Figure 3A
lTGFB_calib = [2.193, 2.645, 3.075, 3.440, 5.268, 5.569, 2.193, 2.559, ....
    2.086, 2.344];
data1 = lTGFB_calib(1:2:end); % 2.1930, 3.0750, 5.2680, 2.1930, 2.0860
error1 = lTGFB_calib(2:2:end) - lTGFB_calib(1:2:end); % 0.452, 0.365, 0.301, 0.366, 0.258

%% Calibration TNFa - PMID 12457250, Figure 8
tnfa_calib = [6.911, 8.788, 10.580, 15.955, 17.576, 19.283, 11.433, ...
    12.713, 14.163, 2.389, 3.754, 5.119];
data2 = tnfa_calib(2:3:end); % 8.7880, 17.5760, 12.7130, 3.7540
error2_l = tnfa_calib(2:3:end) - tnfa_calib(1:3:end);
error2_u = tnfa_calib(3:3:end) - tnfa_calib(2:3:end); % 1.8770, 1.6210, 1.2800, 1.3650

%% Calibration IL-1 - PMID 14742270, Figure 9B
il1_calib = [0.009, 0.012, 0.058, 0.075, 0.193, 0.250, 0.030, 0.033, ...
    0.009, 0.011, 0.004, 0.006];
data3 = il1_calib(1:2:end); % 0.0090, 0.0580, 0.1930, 0.0300, 0.0090, 0.0040
error3 = il1_calib(2:2:end) - il1_calib(1:2:end); % 0.0030, 0.0170, 0.0570, 0.0030, 0.0020, 0.0020
error3 = error3./data3(1); 
data3 = data3./data3(1);

%% Calibration MMP-9 - PMID 18022590, Figure 3
sham = [0.1569]; shamError = [0.2167];
mmp9_calib = [1.905, 2.156, 2.156, 2.387, 0.213, 0.254, 0.175, 0.209];
data = mmp9_calib(1:2:end); % 1.9050, 2.1560, 0.2130, 0.1750
error = mmp9_calib(2:2:end) - mmp9_calib(1:2:end); % 0.2510, 0.2310, 0.0410, 0.0340
data = data./sham; 
error = error./sham; 

%% Calibration GM-CSF - PMID 17208206, Figure 5
gmcsf_calib = [63.771, 86.244, 30.249, 43.921, 5.800, 6.568, 2.800, ...
    4.605, 1.376, 1.998, 1.017, 1.036];
data4 = gmcsf_calib(1:2:end); % 63.771, 30.249, 5.800, 2.800, 1.376, 1.017
error4 = gmcsf_calib(2:2:end) - gmcsf_calib(1:2:end); % 22.473, 13.672, 0.768, 1.805, 0.622, 0.019

%% Calibration Collagen - PMID 14729404, Figure 3
% Sham, sham error, data, data error
coll_calib = [7.676, 7.949, 8.058, 8.384, 7.676, 8.493, 7.186, 7.785, ...
    6.914, 7.186, 11.488 13.992, 7.676, 8.711, 20.798, 23.248, 7.949, ...
    8.548, 29.836, 35.771];
sham5 = coll_calib(1:4:end); % 7.6760, 7.6760, 6.9140, 7.6760, 7.9490
shamError5 = coll_calib(2:4:end) - coll_calib(1:4:end); % 0.2730, 0.8170, 0.2720, 1.0350, 0.5990
data5 = coll_calib(3:4:end); % 8.0580, 7.1860, 11.4880, 20.7980, 29.8360
error5 = coll_calib(4:4:end) - coll_calib(3:4:end); % 0.3260, 0.5990, 2.5040, 2.4500, 5.9350
error5 = error5./sham5(1); data5 = data5./sham5(1);

%% Calibration Macrophages - PMID 23644221, Figure 7C
% [sham, sham error, data(1, 3, 5, 7), error]
macro_calib = [44109.589, 52876.712, 70684.931, 81095.890, ...
    106027.397, 121095.890, 93972.602, 107123.287, 59726.027, 73972.602];
sham6 = macro_calib(1); shamError6 = macro_calib(2) - macro_calib(1);
data6 = macro_calib(3:2:end); % 1.0e+05 * [0.7068, 1.0603, 0.9397, 0.5973]
error6 = macro_calib(4:2:end) - macro_calib(3:2:end); % 1.0e+04 * [1.0411, 1.5068, 1.3151, 1.4247]
% sham 4.4110e+04 +/- 8.7671e+03
data6 = data6./sham6; error6 = error6./sham6; shamError6 = shamError6./sham6;

%% Calibration fibroblasts - PMID 12481929, Figure 1
sham7 = [83.478];
fibro_calib = [76.521, 97.391, 79.999, 97.391, 354.782, 386.086, 1036.521, ...
    1088.695, 1276.521, 1353.043, 1784.347, 1933.913, 2000, 2180.869];
data7 = fibro_calib(1:2:end); % 1.0e+03 * [0.0765, 0.0800, 0.3548, 1.0365, 1.2765, 1.7843, 2.0000]
error7 = fibro_calib(2:2:end) - fibro_calib(1:2:end); % 20.8700, 17.3920, 31.3040, 52.1740, 76.5220, 149.5660, 180.8690
data7 = data7./sham7; error7 = error7./sham7; 

%% Validation lTGFB - PMID 22260784, Figure 1D
sham8 = 1;
ltgfb_valid = [3.716, 3.984, 1.535, 1.669, 1.652, 1.820];
data8 = ltgfb_valid(1:2:end); % 3.7160, 1.5350, 1.6520
error8 = ltgfb_valid(2:2:end) - ltgfb_valid(1:2:end); % 0.2680, 0.1340, 0.1680

%% Validation TNFa - PMID 14742270, Figure 9D
sham9 = [0.00476];
tnfa_valid = [0.00924, 0.01065, 0.01460, 0.01707, 0.01141, 0.013125, ...
    0.00532, 0.00641, 0.00342, 0.00407];
data9 = tnfa_valid(1:2:end); % 0.0092, 0.0146, 0.0114, 0.0053, 0.0034
error9 = tnfa_valid(2:2:end) - tnfa_valid(1:2:end); % 0.0014, 0.0025, 0.0017, 0.0011, 0.0006
data9 = data9./sham9; error9 = error9./sham9;

%% Validation IL-1 - PMID 11691538, Figure 3
% sham, sham error, data, data error
il1_valid = [0.0456, 0.0534, 0.0762, 0.0912, 0.0462, 0.0540, 0.1519, ...
    0.1831, 0.0504, 0.0606, 0.2252, 0.2714, 0.0426, 0.0570, 0.2120, ...
    0.2516, 0.0396, 0.0546, 0.0828, 0.1027, 0.0288, 0.0360, 0.0576, ...
    0.0732, 0.0270, 0.0372, 0.0516, 0.0642, 0.0306, 0.0372, 0.0042, 0.0072];
sham10 = il1_valid(1:4:end); % 0.0456, 0.0462, 0.0504, 0.0426, 0.0396, 0.0288, 0.0270, 0.0306
shamError10 = il1_valid(2:4:end) - il1_valid(1:4:end); % 0.0078, 0.0078, 0.0102, 0.0144, 0.0150, 0.0072, 0.0102, 0.0066
data10 = il1_valid(3:4:end); % 0.0762, 0.1519, 0.2252, 0.2120, 0.0828, 0.0576, 0.0516, 0.0042
error10 = il1_valid(4:4:end) - il1_valid(3:4:end); % 0.0150, 0.0312, 0.0462, 0.0396, 0.0199, 0.0156, 0.0126, 0.0030
data10 = data10./sham10(1); error10 = error10./sham10(1);

%% Validation MMP-9 - PMID 10502816, Figure 1
sham11 = [0.1534]; shamError11 = [0.2727];
mmp9_valid = [2.3352, 2.7556, 0.6534, 1.2045, 0.4261, 0.7670];
data11 = mmp9_valid(1:2:end); % 2.3352, 0.6534, 0.4261
error11 = mmp9_valid(2:2:end) - mmp9_valid(1:2:end); % 0.4204, 0.5511, 0.3409
data11 = data11./sham11; error11 = error11./sham11;

%% Validation GM-CSF - PMID 28978634, Figure 2
sham12 = [2.4960]; shamError12 = [3.7441-2.4960]./sham12;
gmcsf_valid = [88.6115, 103.5881, 16.5366, 19.9687, 1.2480, 1.8720];
data12 = gmcsf_valid(1:2:end); % 88.6115, 16.5366, 1.2480
error12 = gmcsf_valid(2:2:end) - gmcsf_valid(1:2:end); % 14.9766, 3.4321, 0.6240
data12 = data12./sham12; error12 = error12./sham12;

%% Validation Collagen - PMID 26080361, Figure 5
coll1 = [0.0018, 0.0027, 0, 0, 0, 0, 0.0677, 0.0807, 0.1424, 0.1762];
coll3 = [0.3783, 0.4114, 0.1981, 0.2312, 0.7987, 1.0480, 1.1291, 1.2372, ...
    1.2222, 1.4594];
coll_valid = coll1 + coll3;
data13 = coll_valid(1:2:end); % 0.3801, 0.1981, 0.7987, 1.1968, 1.3646
error13 = coll_valid(2:2:end) - coll_valid(1:2:end); % 0.0340, 0.0331, 0.2493, 0.1211, 0.2710
error13 = error13./data13(1); data13 = data13./data13(1); 

%% Validation Macrophages - PMID 30538339, Figure 4B
macro_valid = [100, 136.054, 4897.959, 5873.015, 3356.009, 3922.902, ...
    2380.952, 2721.088, 45.351, 158.730];
data14 = macro_valid(1:2:end); % 1.0e+03 * 0.1000, 4.8980, 3.3560, 2.3810, 0.0454
error14 = macro_valid(2:2:end) - macro_valid(1:2:end); % 36.0540, 975.0560, 566.8930, 340.1360, 113.3790
error14 = error14./data14(1); data14 = data14./data14(1);

%% Validation Fibroblasts - PMID 29664017, Figure 1J
fibro_valid = [722.807, 785.964, 940.350, 1024.561, 2898.245, 2996.491, ...
    2863.157, 2961.403, 2743.859, 2926.315];
data15 = fibro_valid(1:2:end); % 1.0e+03 * 0.7228, 0.9404, 2.8982, 2.8632, 2.7439
error15 = fibro_valid(2:2:end) - fibro_valid(1:2:end); % 63.1570, 84.2110, 98.2460, 98.2460, 182.4560
error15 = error15./data15(1); data15 = data15./data15(1); 
