function [Max_Sum_AVG Max_Diff_Var Max_Entropy Max_Diff_Entropy Max_Sum_Var Max_Sum_Entropy Max_Correl_2 Max_Contrast ]=glcm_all_Max_Range_Call(IG2)
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
Max_Sum_AVG_val = zeros(4,1);
Max_Diff_Var_val= zeros(4,1);
Max_Entropy_val= zeros(4,1);
Max_Diff_Entropy_val= zeros(4,1);
Max_Sum_Var_val = zeros(4,1);
Max_Sum_Entropy_val = zeros(4,1);
Max_Correl_2_val = zeros(4,1);
Max_Contrast_val = zeros(4,1);

[ Mean_Contrast_val_Zero Diff_Var_val_ZERO_Degree Sum_AVG_val_ZERO_Degree Sum_Entropy_val_ZERO_Degree Contrast_val_ZERO_Degree Diff_Entropy_val_ZERO_Degree Correl_2_val_ZERO_Degree Sum_Var_val_ZERO_Degree Entropy_val_ZERO_Degree ]=glcm_Diff_all_Zero_Degree_Range_Call(IG2);

Max_Sum_AVG_val(1)= Sum_AVG_val_ZERO_Degree;
Max_Diff_Var_val(1)= Diff_Var_val_ZERO_Degree;
Max_Entropy_val(1)= Entropy_val_ZERO_Degree;
Max_Diff_Entropy_val(1)= Diff_Entropy_val_ZERO_Degree;
Max_Sum_Var_val(1) = Sum_Var_val_ZERO_Degree;
Max_Sum_Entropy_val(1) = Sum_Entropy_val_ZERO_Degree;
Max_Correl_2_val(1) = Correl_2_val_ZERO_Degree;
Max_Contrast_val(1) = Contrast_val_ZERO_Degree;

[Mean_Contrast_val_45 Diff_Var_val_45_Degree Sum_AVG_val_45_Degree Sum_Entropy_val_45_Degree Contrast_val_45_Degree Diff_Entropy_val_45_Degree Correl_2_val_45_Degree Sum_Var_val_45_Degree Entropy_val_45_Degree ]=glcm_Diff_all_45_Degree_Range_Call(IG2);

Max_Sum_AVG_val(2)= Sum_AVG_val_45_Degree;
Max_Diff_Var_val(2)= Diff_Var_val_45_Degree;
Max_Entropy_val(2)= Entropy_val_45_Degree;
Max_Diff_Entropy_val(2)= Diff_Entropy_val_45_Degree;
Max_Sum_Var_val(2) = Sum_Var_val_45_Degree;
Max_Sum_Entropy_val(2) = Sum_Entropy_val_45_Degree;
Max_Correl_2_val(2) = Correl_2_val_45_Degree;
Max_Contrast_val(2) = Contrast_val_45_Degree;

[Mean_Contrast_val_90 Diff_Var_val_90_Degree Sum_AVG_val_90_Degree Sum_Entropy_val_90_Degree Contrast_val_90_Degree Diff_Entropy_val_90_Degree Correl_2_val_90_Degree Sum_Var_val_90_Degree Entropy_val_90_Degree ]=glcm_Diff_all_90_Degree_Range_Call(IG2);

Max_Sum_AVG_val(3)= Sum_AVG_val_90_Degree;
Max_Diff_Var_val(3)= Diff_Var_val_90_Degree;
Max_Entropy_val(3)= Entropy_val_90_Degree;
Max_Diff_Entropy_val(3)= Diff_Entropy_val_90_Degree;
Max_Sum_Var_val(3) = Sum_Var_val_90_Degree;
Max_Sum_Entropy_val(3) = Sum_Entropy_val_90_Degree;
Max_Correl_2_val(3) = Correl_2_val_90_Degree;
Max_Contrast_val(3) = Contrast_val_90_Degree;

[Mean_Contrast_val_135 Diff_Var_val_135_Degree Sum_AVG_val_135_Degree Sum_Entropy_val_135_Degree Contrast_val_135_Degree Diff_Entropy_val_135_Degree Correl_2_val_135_Degree Sum_Var_val_135_Degree Entropy_val_135_Degree ]=glcm_Diff_all_135_Degree_Range_Call(IG2);

Max_Sum_AVG_val(4)= Sum_AVG_val_135_Degree;
Max_Diff_Var_val(4)= Diff_Var_val_135_Degree;
Max_Entropy_val(4)= Entropy_val_135_Degree;
Max_Diff_Entropy_val(4)= Diff_Entropy_val_135_Degree;
Max_Sum_Var_val(4) = Sum_Var_val_135_Degree;
Max_Sum_Entropy_val(4) = Sum_Entropy_val_135_Degree;
Max_Correl_2_val(4) = Correl_2_val_135_Degree;
Max_Contrast_val(4) = Contrast_val_135_Degree;


Max_Sum_AVG = max(Max_Sum_AVG_val);
Max_Diff_Var = max( Max_Diff_Var_val);
Max_Entropy = max(Max_Entropy_val);
Max_Diff_Entropy = max(Max_Diff_Entropy_val);
Max_Sum_Var = max(Max_Sum_Var_val);
Max_Sum_Entropy = max(Max_Sum_Entropy_val);
Max_Correl_2 = max(Max_Correl_2_val);
Max_Contrast = max(Max_Contrast_val);


    