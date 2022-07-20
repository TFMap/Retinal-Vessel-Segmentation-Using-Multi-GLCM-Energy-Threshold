function [IQR_Energy_AVG IQR_IDM_AVG]=glcmEnergy_AND_IDM_Average_IQR_Call(IG2)
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
[Energy_Zero_Degree IDM_Zero_Degree]=glcmEnergy_AND_IDM_Zero_Degree_IQR_Call(IG2);
[Energy_45_Degree IDM_45_Degree]=glcmEnergy_AND_IDM_45_Degree_IQR_Call(IG2);
[Energy_90_Degree IDM_90_Degree]=glcmEnergy_AND_IDM_90_Degree_IQR_Call(IG2);
[Energy_135_Degree IDM_135_Degree]=glcmEnergy_AND_IDM_135_Degree_IQR_Call(IG2);

    IQR_Energy_AVG =  (Energy_Zero_Degree + Energy_45_Degree + Energy_90_Degree + Energy_135_Degree)/4;
    IQR_IDM_AVG =  (IDM_Zero_Degree + IDM_45_Degree + IDM_90_Degree + IDM_135_Degree)/4;