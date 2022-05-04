%% ODE Function
% Last updated: 5/4/22
function [ dydt, crowding ] = modelEquations(t, y, parameters, constants, ...
    receptors, knockouts, codeSnip)
    dydt = zeros(42, 1);

    % Parameters for perturbations
    [n1, n2, n3, n4, n5] = deal(1); [n6, n7] = deal(0);
    [nt1, nt2, nt4, nt6, nt7] = deal(1); 
    [mp1, mp2, mp3, mp4] = deal(1); mp5 = 0; ci1 = 1; ci2 = 0;
    if isempty(codeSnip) == 0
        eval(codeSnip);
    end
   
    % Saturation constants
    k_tnfa = constants(3); k_tgfb = constants(5); k_il1 = constants(7); 
    k_gmcsf = constants(4); k_mmp = constants(6);
    
    % Hill function terms
    h_tnfa = receptors(4)*y(34)/(y(34) + k_tnfa); h_tnfa_i = k_tnfa/(receptors(4)*y(34) + k_tnfa);
    h_tgfb = receptors(3)*y(7)/(y(7) + k_tgfb); h_tgfb_i = k_tgfb/(receptors(3)*y(7) + k_tgfb);
    h_il1 = receptors(2)*y(4)/(y(4) + k_il1); 
    h_gmcsf = receptors(1)*y(5)/(y(5) + k_gmcsf); 
    h_mmp = y(8)/(y(8) + k_mmp); h_mmp_i = k_mmp/(y(8) + k_mmp);

    % Crowding effects
    kCrowd = (y(9) + y(10))/constants(10) + (y(11) + y(40))/constants(11) + ...
        y(12)/constants(12) + y(1)/constants(8) + y(3)/constants(9)...
            + (y(2)/constants(2)) + 0.17;
    crowdEff = kCrowd;
    
    % Cell behavior terms
    mac_diff = (h_gmcsf + h_il1)*parameters(1)*y(12);
    mac_prolif = (h_gmcsf*(h_tgfb + h_tnfa) + h_tgfb)*parameters(2)*y(1)*(1-crowdEff);
    mac_rem = (h_tnfa_i + h_mmp_i)*parameters(3);
    phago_mac = parameters(4)*(mp2*h_tgfb + mp1*h_tnfa_i + mp5)*h_mmp*y(1); 
    n_inf = (n1*h_gmcsf + n2*h_il1 + n3*h_tgfb_i + n6)*parameters(5)*y(13)*(1 - crowdEff);
    n_rem = (n4*h_tnfa_i + n5*h_mmp_i + n7)*parameters(6)*y(11);
    phago_n = parameters(7)*y(11);
    mo_inf = (h_gmcsf + h_il1)*parameters(8)*y(14)*(1 - crowdEff);
    mo_rem = (h_tnfa_i + h_mmp_i)*parameters(9)*y(12);
    mo_blood_fb = parameters(10)*h_il1*y(14);
    fib_prolif = h_tgfb*parameters(11)*y(2)*(1-crowdEff);
    fib_death = parameters(44)*y(2);
    cm_death = parameters(12)*(h_tnfa + h_il1 + h_tgfb_i)*y(9);
    coll_mat = parameters(14)*h_mmp_i*y(33);
    pro_coll = h_tgfb*parameters(15)*y(2)*(ci1*(1/(y(3)+1)) + ci2);
    
    % IL-1 secretion terms
    il1_cm = parameters(17)*y(9);
    il1_n = parameters(18)*y(11);
    il1_mo = parameters(19)*y(12);
    il1_mac = h_tgfb_i*parameters(20)*y(1);
    il1_fib = parameters(21)*y(2);
    il1_deg = parameters(22);
    
    % GM-CSF secretion terms
    gm_n = parameters(23)*y(11);
    gm_mo = parameters(24)*y(12);
    gm_fib = parameters(25)*y(2);
    gm_deg = parameters(26);
    
    % TGFB secretion terms
    tgfb_cm = parameters(27)*y(9);
    tgfb_mac = h_tgfb*(y(42)/(5e4 + y(42)))*parameters(28)*y(1);
    tgfb_fib = h_tgfb*parameters(29)*y(2);
    tgfb_deg = parameters(30);
    
    % MMP-9 secretion terms
    mmp_n = parameters(33)*y(11);
    mmp_mo = parameters(34)*y(12);
    mmp_mac = parameters(35)*y(1);
    mmp_fib = (h_tgfb)*parameters(36)*y(2);
    mmp_deg = parameters(37);
    
    % TNFa secretion terms
    tnfa_cm = nt7*parameters(38)*y(9);
    tnfa_n = nt1*parameters(39)*y(11);
    tnfa_mo = nt2*parameters(40)*y(12);
    tnfa_m = h_tgfb_i*nt4*parameters(41)*y(1);
    tnfa_fib = nt6*h_tgfb_i*parameters(42)*y(2);
    tnfa_deg = parameters(43);
    
    % DEs
    macrophages = mac_diff + mac_prolif - mac_rem*y(1); % Macrophages
    fibroblasts = fib_prolif - fib_death; % Fibroblasts
    collagen = coll_mat - parameters(16)*y(3)*y(8); % Mature collagen
    il_1 = knockouts(1)*(il1_cm + il1_n + il1_mo + il1_mac + il1_fib - il1_deg*y(4)); % IL-1
    gm_csf = knockouts(2)*(gm_n + gm_mo + gm_fib - gm_deg*y(5)); % GM-CSF
    latent_tgfb = knockouts(3)*(tgfb_cm + tgfb_mac + tgfb_fib - tgfb_deg*y(6)); % Latent TGFB
    tgfb = knockouts(4)*(parameters(31)*(h_mmp+0.05)*y(6)...
        - parameters(32)*y(7)); % TGFB
    mmp_9 = knockouts(5)*(mmp_n + mmp_mo + mmp_mac + mmp_fib - mmp_deg*y(8)); % MMP-9
    cardiomyo = -cm_death; % Cardiomyocytes
    cardiomyo_deb = parameters(13)*cm_death - mp3*phago_mac*y(10) - 2*phago_n*y(10); % Debris of dead cardiomyocytes
    inf_neutro = n_inf - n_rem; % infarct neutrophils
    inf_mono = mo_inf - mo_rem; % infarct monocytes 
    blood_neutro = -n_inf; % blood neutrophils
    blood_mono = -mo_inf + mo_blood_fb; % blood monocytes
    tnfa = knockouts(6)*(tnfa_cm + tnfa_n + tnfa_mo + tnfa_m + ...
        tnfa_fib - tnfa_deg*y(34)); % TNFa 
    neutro_deb = parameters(13)*n_rem - mp4*phago_mac*y(40) - phago_n*y(40); % Neutrophil debris
    ingested_deb = mp4*phago_mac*y(40) + mp3*phago_mac*y(10) - 0.25*y(42); % Debris ingested by macrophages
    infiltrating_mono = -mo_inf; % infiltrating monocytes
    
    % Fibroblast behavior
    il1_fibroblasts = knockouts(1)*(il1_fib - il1_deg*y(15)); % IL-1 secretion
    gmcsf_fibroblasts = knockouts(2)*(gm_fib - gm_deg*y(16)); % GM-CSF secretion
    tgfb_fibroblasts = knockouts(3)*(tgfb_fib - tgfb_deg*y(17)); % TGFB secretion
    mmp9_fibroblasts = knockouts(5)*(mmp_fib - mmp_deg*y(18)); % MMP-9 secretion
    tnfa_fibroblasts = knockouts(6)*(tnfa_fib - tnfa_deg*y(35)); % TNFa secretion
 
    % Macrophage behavior
    diff_macrophages = mac_diff - mac_rem*y(19); % Macrophage differentiation
    prolif_macrophages = mac_prolif - mac_rem*y(20); % Macrophage proliferation
    tgfb_macrophages = knockouts(3)*(tgfb_mac - tgfb_deg*y(21)); % TGFB secretion 
    il1_macrophages = knockouts(1)*(il1_mac - il1_deg*y(22)); % IL-1 secretion
    mmp9_macrophages = knockouts(5)*(mmp_mac - mmp_deg*y(23)); % MMP-9 secretion
    tnfa_macrophages = knockouts(6)*(tnfa_m - tnfa_deg*y(36)); % TNFa secretion

    % Cardiomyocyte behavior
    il1_cardiomyo = knockouts(1)*(il1_cm - il1_deg*y(24)); % IL-1 secretion
    tnfa_cardiomyo = knockouts(6)*(tnfa_cm - tnfa_deg*y(37)); % TNFa secretion
    tgfb_cardiomyo = knockouts(3)*(tgfb_cm - tgfb_deg*y(41)); % TGFB secretion
    
    % Neutrophil behavior
    il1_neutro = knockouts(1)*(il1_n - il1_deg*y(26)); % IL-1 secretion
    gmcsf_neutro = knockouts(2)*(gm_n - gm_deg*y(27)); % GM-CSF secretion
    mmp9_neutro = knockouts(5)*(mmp_n - mmp_deg*y(28)); % MMP-9 secretion
    tnfa_neutro = knockouts(6)*(tnfa_n - tnfa_deg*y(38)); % TNFa secretion
    
    % Monocyte behavior
    il1_mono = knockouts(1)*(il1_mo - il1_deg*y(30)); % IL-1 secretion
    gmcsf_mono = knockouts(2)*(gm_mo - gm_deg*y(31)); % GM-CSF secretion
    mmp9_mono = knockouts(5)*(mmp_mo - mmp_deg*y(32)); % MMP-9 secretion
    tnfa_mono = knockouts(6)*(tnfa_mo - tnfa_deg*y(39)); % TNFa secretion
    
    % Collagen
    secreted_collagen = pro_coll - coll_mat; % Collagen secretion by fibroblasts
    
    dydt = [macrophages; fibroblasts; collagen; il_1; gm_csf; latent_tgfb; ...
        tgfb; mmp_9; cardiomyo; cardiomyo_deb; inf_neutro; inf_mono; blood_neutro; ...
        blood_mono; il1_fibroblasts; gmcsf_fibroblasts; tgfb_fibroblasts; ...
        mmp9_fibroblasts; diff_macrophages; prolif_macrophages; tgfb_macrophages; ...
        il1_macrophages; mmp9_macrophages; il1_cardiomyo; 0; il1_neutro; ...
        gmcsf_neutro; mmp9_neutro; 0; il1_mono; gmcsf_mono; mmp9_mono; ...
        secreted_collagen; tnfa; tnfa_fibroblasts; tnfa_macrophages; ...
        tnfa_cardiomyo; tnfa_neutro; tnfa_mono; neutro_deb; tgfb_cardiomyo; ...
        ingested_deb; infiltrating_mono];
    
    % Return crowding effect variables
    crowding = [ crowdEff, (cardiomyo + cardiomyo_deb)/constants(10), ...
        (inf_neutro + neutro_deb)/constants(11), ...
        inf_mono/constants(12), macrophages/constants(8), collagen/constants(9), ...
        (fibroblasts/constants(2)), 0.17];
    
    % For the perturbations
    if isempty(codeSnip) == 0
        eval(codeSnip);
    end
end