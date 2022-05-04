function [ y0, constants, codeSnip ] = initParams(infarct_size, constants, codeSnip)
    y0 = zeros(1, 43); 
    constants(14) = infarct_size;
    y0(2) = floor(constants(1)*infarct_size);
    y0(9) = 6e6*infarct_size;
    y0(6) = 40*infarct_size;
    y0(3) = 15*infarct_size; 
    y0(13) = 1e5*infarct_size;
    y0(14) = 1e5*infarct_size;
    y0(43) = 1e5*infarct_size; 
    codeSnip = codeSnip;
    
    constants([2, 8:12]) = constants([2, 8:12])*infarct_size;
end