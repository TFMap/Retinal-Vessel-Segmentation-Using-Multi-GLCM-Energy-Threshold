function [IQR_Diff_Var_AVG  IQR_Sum_AVG_AVG IQR_Sum_Entropy_AVG IQR_Contrast_AVG IQR_Diff_Entropy_AVG  IQR_Correl_2_AVG IQR_Sum_Var_AVG IQR_Entropy_AVG]=glcm_all_Average_IQR_Call(IG2)
% HARALICK Fast Calculation of Haralick Features
%   IN:   glcm = Co-Occurrence Matrix     
%   OUT:  F = Feature Vector   
%
%   Stefan Winzeck 2012   
%   winzeck@hm.edu
% 
%   Feature Calculation according to:
%   [1] R. Haralick: 'Textural Feature for Image Classification' (1979)
%   [2] E. Miyamoto: 'Fast Calculation of Haralick Texture Features' 
% 
% MISSING:   f14  [1]

%% ALLOCATION
clc;
[ Mean_Contrast_val_Zero IQR_Diff_Var_val_ZERO_Degree IQR_Sum_AVG_val_ZERO_Degree IQR_Sum_Entropy_val_ZERO_Degree IQR_Contrast_val_ZERO_Degree IQR_Diff_Entropy_val_ZERO_Degree IQR_Correl_2_val_ZERO_Degree IQR_Sum_Var_val_ZERO_Degree IQR_Entropy_val_ZERO_Degree ]=glcm_Diff_all_Zero_Degree_IQR_Call(IG2);

[Mean_Contrast_val_45 IQR_Diff_Var_val_45_Degree IQR_Sum_AVG_val_45_Degree IQR_Sum_Entropy_val_45_Degree IQR_Contrast_val_45_Degree IQR_Diff_Entropy_val_45_Degree IQR_Correl_2_val_45_Degree IQR_Sum_Var_val_45_Degree IQR_Entropy_val_45_Degree ]=glcm_Diff_all_45_Degree_IQR_Call(IG2);

[Mean_Contrast_val_90 IQR_Diff_Var_val_90_Degree IQR_Sum_AVG_val_90_Degree IQR_Sum_Entropy_val_90_Degree IQR_Contrast_val_90_Degree IQR_Diff_Entropy_val_90_Degree IQR_Correl_2_val_90_Degree IQR_Sum_Var_val_90_Degree IQR_Entropy_val_90_Degree ]=glcm_Diff_all_90_Degree_IQR_Call(IG2);

[Mean_Contrast_val_135 IQR_Diff_Var_val_135_Degree IQR_Sum_AVG_val_135_Degree IQR_Sum_Entropy_val_135_Degree IQR_Contrast_val_135_Degree IQR_Diff_Entropy_val_135_Degree IQR_Correl_2_val_135_Degree IQR_Sum_Var_val_135_Degree IQR_Entropy_val_135_Degree ]=glcm_Diff_all_135_Degree_IQR_Call(IG2);


    IQR_Diff_Var_AVG =  (IQR_Diff_Var_val_ZERO_Degree + IQR_Diff_Var_val_45_Degree + IQR_Diff_Var_val_90_Degree + IQR_Diff_Var_val_135_Degree)/4;
    
    IQR_Sum_AVG_AVG = (IQR_Sum_AVG_val_ZERO_Degree + IQR_Sum_AVG_val_45_Degree + IQR_Sum_AVG_val_90_Degree + IQR_Sum_AVG_val_135_Degree)/4;
    
    IQR_Sum_Entropy_AVG = (IQR_Sum_Entropy_val_ZERO_Degree + IQR_Sum_Entropy_val_45_Degree + IQR_Sum_Entropy_val_90_Degree + IQR_Sum_Entropy_val_135_Degree)/4;

    IQR_Contrast_AVG = (IQR_Contrast_val_ZERO_Degree + IQR_Contrast_val_45_Degree + IQR_Contrast_val_90_Degree + IQR_Contrast_val_135_Degree )/4;
    
    IQR_Diff_Entropy_AVG = (IQR_Diff_Entropy_val_ZERO_Degree + IQR_Diff_Entropy_val_45_Degree + IQR_Diff_Entropy_val_90_Degree + IQR_Diff_Entropy_val_135_Degree )/4;

    IQR_Correl_2_AVG = (IQR_Correl_2_val_ZERO_Degree + IQR_Correl_2_val_45_Degree + IQR_Correl_2_val_90_Degree + IQR_Correl_2_val_135_Degree )/4;

    IQR_Sum_Var_AVG = (IQR_Sum_Var_val_ZERO_Degree + IQR_Sum_Var_val_45_Degree + IQR_Sum_Var_val_90_Degree + IQR_Sum_Var_val_135_Degree )/4;

    IQR_Entropy_AVG = (IQR_Entropy_val_ZERO_Degree + IQR_Entropy_val_45_Degree + IQR_Entropy_val_90_Degree + IQR_Entropy_val_135_Degree  )/4;

    