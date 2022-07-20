function [Max_Mean_Contrast  Min_Mean_Contrast AVG_Mean_Contrast]=glcm_all_All_Mean_Contrast_Call(IG2)
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
All_Mean_Contrast_val = zeros(4,1);

[ Mean_Contrast_val_Zero ]=glcm_Diff_all_Zero_Degree_Range_Call(IG2);


All_Mean_Contrast_val(1) = mean(Mean_Contrast_val_Zero);

[Mean_Contrast_val_45 ]=glcm_Diff_all_45_Degree_Range_Call(IG2);

All_Mean_Contrast_val(2) = mean(Mean_Contrast_val_45);

[Mean_Contrast_val_90 ]=glcm_Diff_all_90_Degree_Range_Call(IG2);

All_Mean_Contrast_val(3) = mean(Mean_Contrast_val_90);

[Mean_Contrast_val_135 ]=glcm_Diff_all_135_Degree_Range_Call(IG2);

All_Mean_Contrast_val(4) = mean(Mean_Contrast_val_135);


Max_Mean_Contrast = max(All_Mean_Contrast_val);
Min_Mean_Contrast = min(All_Mean_Contrast_val);
AVG_Mean_Contrast = mean(All_Mean_Contrast_val);


    