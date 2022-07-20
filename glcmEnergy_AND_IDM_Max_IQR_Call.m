function [IQR_Max_Energy IQR_Max_IDM]=glcmEnergy_AND_IDM_Max_IQR_Call(IG2)
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
Energy_Max_Val=zeros(1,4);
IDM_Max_Val=zeros(1,4);
[Energy_Zero_Degree IDM_Zero_Degree]=glcmEnergy_AND_IDM_Zero_Degree_IQR_Call(IG2);
Energy_Max_Val(1)= Energy_Zero_Degree;
IDM_Max_Val(1)= IDM_Zero_Degree;
[Energy_45_Degree IDM_45_Degree]=glcmEnergy_AND_IDM_45_Degree_IQR_Call(IG2);
Energy_Max_Val(2)= Energy_45_Degree;
IDM_Max_Val(2)= IDM_45_Degree;
[Energy_90_Degree IDM_90_Degree]=glcmEnergy_AND_IDM_90_Degree_IQR_Call(IG2);
Energy_Max_Val(3)= Energy_90_Degree;
IDM_Max_Val(3)= IDM_90_Degree;
[Energy_135_Degree IDM_135_Degree]=glcmEnergy_AND_IDM_135_Degree_IQR_Call(IG2);
Energy_Max_Val(4)= Energy_135_Degree;
IDM_Max_Val(4)= IDM_135_Degree;

  IQR_Max_Energy =  max(Energy_Max_Val);
  IQR_Max_IDM =  max(IDM_Max_Val);
  %[IQR_Energy_AVG IQR_IDM_AVG]=glcmEnergy_AND_IDM_Average_IQR_Call(IG2);
  %[IQR_Max_Energy IQR_Max_IDM]=glcmEnergy_AND_IDM_Max_IQR_Call(IG2)