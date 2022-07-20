function [Sum_AVG_AVG Diff_Var_AVG Sum_Entropy_AVG Contrast_AVG Diff_Entropy_AVG  Correl_2_AVG Sum_Var_AVG Entropy_AVG]=glcm_all_Average_Range_Call(IG2)
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
[ Mean_Contrast_val_Zero Diff_Var_val_ZERO_Degree Sum_AVG_val_ZERO_Degree Sum_Entropy_val_ZERO_Degree Contrast_val_ZERO_Degree Diff_Entropy_val_ZERO_Degree Correl_2_val_ZERO_Degree Sum_Var_val_ZERO_Degree Entropy_val_ZERO_Degree ]=glcm_Diff_all_Zero_Degree_Range_Call(IG2);

[Mean_Contrast_val_45 Diff_Var_val_45_Degree Sum_AVG_val_45_Degree Sum_Entropy_val_45_Degree Contrast_val_45_Degree Diff_Entropy_val_45_Degree Correl_2_val_45_Degree Sum_Var_val_45_Degree Entropy_val_45_Degree ]=glcm_Diff_all_45_Degree_Range_Call(IG2);

[Mean_Contrast_val_90 Diff_Var_val_90_Degree Sum_AVG_val_90_Degree Sum_Entropy_val_90_Degree Contrast_val_90_Degree Diff_Entropy_val_90_Degree Correl_2_val_90_Degree Sum_Var_val_90_Degree Entropy_val_90_Degree ]=glcm_Diff_all_90_Degree_Range_Call(IG2);

[Mean_Contrast_val_135 Diff_Var_val_135_Degree Sum_AVG_val_135_Degree Sum_Entropy_val_135_Degree Contrast_val_135_Degree Diff_Entropy_val_135_Degree Correl_2_val_135_Degree Sum_Var_val_135_Degree Entropy_val_135_Degree ]=glcm_Diff_all_135_Degree_Range_Call(IG2);


    Diff_Var_AVG =  (Diff_Var_val_ZERO_Degree + Diff_Var_val_45_Degree + Diff_Var_val_90_Degree + Diff_Var_val_135_Degree)/4;
    
    Sum_AVG_AVG = (Sum_AVG_val_ZERO_Degree + Sum_AVG_val_45_Degree + Sum_AVG_val_90_Degree + Sum_AVG_val_135_Degree)/4;
    
    Sum_Entropy_AVG = (Sum_Entropy_val_ZERO_Degree + Sum_Entropy_val_45_Degree + Sum_Entropy_val_90_Degree + Sum_Entropy_val_135_Degree)/4;

    Contrast_AVG = (Contrast_val_ZERO_Degree + Contrast_val_45_Degree + Contrast_val_90_Degree + Contrast_val_135_Degree )/4;
    
    Diff_Entropy_AVG = (Diff_Entropy_val_ZERO_Degree + Diff_Entropy_val_45_Degree + Diff_Entropy_val_90_Degree + Diff_Entropy_val_135_Degree )/4;

    Correl_2_AVG = (Correl_2_val_ZERO_Degree + Correl_2_val_45_Degree + Correl_2_val_90_Degree + Correl_2_val_135_Degree )/4;

    Sum_Var_AVG = (Sum_Var_val_ZERO_Degree + Sum_Var_val_45_Degree + Sum_Var_val_90_Degree + Sum_Var_val_135_Degree )/4;

    Entropy_AVG = (Entropy_val_ZERO_Degree + Entropy_val_45_Degree + Entropy_val_90_Degree + Entropy_val_135_Degree  )/4;

    