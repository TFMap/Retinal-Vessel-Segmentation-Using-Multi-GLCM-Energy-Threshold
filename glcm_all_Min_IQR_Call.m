function [Min_IQR_Sum_AVG Min_IQR_Diff_Var Min_IQR_Entropy Min_IQR_Diff_Entropy Min_IQR_Sum_Var Min_IQR_Sum_Entropy Min_IQR_Correl_2 Min_IQR_Contrast ]=glcm_all_Min_IQR_Call(IG2)
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
IQR_Min_Sum_AVG_val = zeros(4,1);
IQR_Min_Diff_Var_val= zeros(4,1);
IQR_Min_Entropy_val= zeros(4,1);
IQR_Min_Diff_Entropy_val= zeros(4,1);
IQR_Min_Sum_Var_val = zeros(4,1);
IQR_Min_Sum_Entropy_val = zeros(4,1);
IQR_Min_Correl_2_val = zeros(4,1);
IQR_Min_Contrast_val = zeros(4,1);

[ Mean_Contrast_val_Zero IQR_Diff_Var_val_ZERO_Degree IQR_Sum_AVG_val_ZERO_Degree IQR_Sum_Entropy_val_ZERO_Degree IQR_Contrast_val_ZERO_Degree IQR_Diff_Entropy_val_ZERO_Degree IQR_Correl_2_val_ZERO_Degree IQR_Sum_Var_val_ZERO_Degree IQR_Entropy_val_ZERO_Degree ]=glcm_Diff_all_Zero_Degree_IQR_Call(IG2);

IQR_Min_Sum_AVG_val(1)= IQR_Sum_AVG_val_ZERO_Degree;
IQR_Min_Diff_Var_val(1)= IQR_Diff_Var_val_ZERO_Degree;
IQR_Min_Entropy_val(1)= IQR_Entropy_val_ZERO_Degree;
IQR_Min_Diff_Entropy_val(1)= IQR_Diff_Entropy_val_ZERO_Degree;
IQR_Min_Sum_Var_val(1) = IQR_Sum_Var_val_ZERO_Degree;
IQR_Min_Sum_Entropy_val(1) = IQR_Sum_Entropy_val_ZERO_Degree;
IQR_Min_Correl_2_val(1) = IQR_Correl_2_val_ZERO_Degree;
IQR_Min_Contrast_val(1) = IQR_Contrast_val_ZERO_Degree;

[Mean_Contrast_val_45 IQR_Diff_Var_val_45_Degree IQR_Sum_AVG_val_45_Degree IQR_Sum_Entropy_val_45_Degree IQR_Contrast_val_45_Degree IQR_Diff_Entropy_val_45_Degree IQR_Correl_2_val_45_Degree IQR_Sum_Var_val_45_Degree IQR_Entropy_val_45_Degree ]=glcm_Diff_all_45_Degree_IQR_Call(IG2);

IQR_Min_Sum_AVG_val(2)= IQR_Sum_AVG_val_45_Degree;
IQR_Min_Diff_Var_val(2)= IQR_Diff_Var_val_45_Degree;
IQR_Min_Entropy_val(2)= IQR_Entropy_val_45_Degree;
IQR_Min_Diff_Entropy_val(2)= IQR_Diff_Entropy_val_45_Degree;
IQR_Min_Sum_Var_val(2) = IQR_Sum_Var_val_45_Degree;
IQR_Min_Sum_Entropy_val(2) = IQR_Sum_Entropy_val_45_Degree;
IQR_Min_Correl_2_val(2) = IQR_Correl_2_val_45_Degree;
IQR_Min_Contrast_val(2) = IQR_Contrast_val_45_Degree;

[Mean_Contrast_val_90 IQR_Diff_Var_val_90_Degree IQR_Sum_AVG_val_90_Degree IQR_Sum_Entropy_val_90_Degree IQR_Contrast_val_90_Degree IQR_Diff_Entropy_val_90_Degree IQR_Correl_2_val_90_Degree IQR_Sum_Var_val_90_Degree IQR_Entropy_val_90_Degree ]=glcm_Diff_all_90_Degree_IQR_Call(IG2);

IQR_Min_Sum_AVG_val(3)= IQR_Sum_AVG_val_90_Degree;
IQR_Min_Diff_Var_val(3)= IQR_Diff_Var_val_90_Degree;
IQR_Min_Entropy_val(3)= IQR_Entropy_val_90_Degree;
IQR_Min_Diff_Entropy_val(3)= IQR_Diff_Entropy_val_90_Degree;
IQR_Min_Sum_Var_val(3) = IQR_Sum_Var_val_90_Degree;
IQR_Min_Sum_Entropy_val(3) = IQR_Sum_Entropy_val_90_Degree;
IQR_Min_Correl_2_val(3) = IQR_Correl_2_val_90_Degree;
IQR_Min_Contrast_val(3) = IQR_Contrast_val_90_Degree;

[Mean_Contrast_val_135 IQR_Diff_Var_val_135_Degree IQR_Sum_AVG_val_135_Degree IQR_Sum_Entropy_val_135_Degree IQR_Contrast_val_135_Degree IQR_Diff_Entropy_val_135_Degree IQR_Correl_2_val_135_Degree IQR_Sum_Var_val_135_Degree IQR_Entropy_val_135_Degree ]=glcm_Diff_all_135_Degree_IQR_Call(IG2);

IQR_Min_Sum_AVG_val(4)= IQR_Sum_AVG_val_135_Degree;
IQR_Min_Diff_Var_val(4)= IQR_Diff_Var_val_135_Degree;
IQR_Min_Entropy_val(4)= IQR_Entropy_val_135_Degree;
IQR_Min_Diff_Entropy_val(4)= IQR_Diff_Entropy_val_135_Degree;
IQR_Min_Sum_Var_val(4) = IQR_Sum_Var_val_135_Degree;
IQR_Min_Sum_Entropy_val(4) = IQR_Sum_Entropy_val_135_Degree;
IQR_Min_Correl_2_val(4) = IQR_Correl_2_val_135_Degree;
IQR_Min_Contrast_val(4) = IQR_Contrast_val_135_Degree;


Min_IQR_Sum_AVG = min(IQR_Min_Sum_AVG_val);
Min_IQR_Diff_Var = min( IQR_Min_Diff_Var_val);
Min_IQR_Entropy = min(IQR_Min_Entropy_val);
Min_IQR_Diff_Entropy = min(IQR_Min_Diff_Entropy_val);
Min_IQR_Sum_Var = min(IQR_Min_Sum_Var_val);
Min_IQR_Sum_Entropy = min(IQR_Min_Sum_Entropy_val);
Min_IQR_Correl_2 = min(IQR_Min_Correl_2_val);
Min_IQR_Contrast = min(IQR_Min_Contrast_val);


    