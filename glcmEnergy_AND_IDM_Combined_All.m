%function [Energy_Zero_Degree IDM_Zero_Degree]=glcmEnergy_AND_IDM_Combined_All(IG2)
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
Energy_val = zeros(4,1);
IDM_val = zeros(4,1);
fx5=strcat('AVERAGE OF ENERGY AND LOCAL HOMOGENEITY DIFFERENCES FOR ALL DIRECTIONS', '.txt');
fName=fx5; %# A file name
fid = fopen(fName,'w');            %# Open the file

fx8=('AVERAGE OF ENERGY DIFFERENCES AND LOCAL HOMOGENEITY DIFFERENCES FOR ALL DIRECTIONS'); %ImagesTestSegREAL_AVERAGE_COMBINEMORPHOLOGY3
fprintf(fid,'__________________________________________________________________________________________\r\n');
fprintf(fid,'%s\r\n',fx8);  %# Print the string
fprintf(fid,'__________________________________________________________________________________________\r\n');
for kat=1:20
% original file name
  if (kat<=9)
    str = int2str(kat); % 01_test.tif
    fi=strcat('0',str, '_test.tif');
    u_filename = fi;
  else
    str = int2str(kat); %  01_test.tif
    fi=strcat( str, '_test.tif');
    u_filename = fi;
  end;

  
% read files 
F = double(imread(u_filename)) / 255;

% figure, imshow([u_GT, u_bw]) 
H=F(:,:,1);
IG2=F(:,:,2);
J=F(:,:,3);



[Energy_Zero_Degree IDM_Zero_Degree]=glcmEnergy_AND_IDM_Zero_Degree_Call(IG2);
[Energy_45_Degree IDM_45_Degree]=glcmEnergy_AND_IDM_45_Degree_Call(IG2);
[Energy_90_Degree IDM_90_Degree]=glcmEnergy_AND_IDM_90_Degree_Call(IG2);
[Energy_135_Degree IDM_135_Degree]=glcmEnergy_AND_IDM_135_Degree_Call(IG2);

    Energy_AVG =  (Energy_Zero_Degree + Energy_45_Degree + Energy_90_Degree + Energy_135_Degree)/4;
    IDM_AVG =  (IDM_Zero_Degree + IDM_45_Degree + IDM_90_Degree + IDM_135_Degree)/4;
    str = int2str(kat);
    fx2=strcat('MEASUREMENT OF ENERGY DIFFERENCE AND LOCAL HOMOGENIETY DIFFERENCE ALL DIRECTIONS OF IMAGE: ',str);
    fprintf(fid,'____________________________________________________________________________________________\r\n');
    fprintf(fid,'%s\r\n',fx2);  %# Print the string
    fprintf(fid,'____________________________________________________________________________________________ \r\n');
    fprintf(fid,'Energy Difference of Zero Degree = %9.5f \r\n',Energy_Zero_Degree);
    fprintf(fid,'Local Homogeneity Difference of Zero Degree = %9.5f \r\n', IDM_Zero_Degree);
    fprintf(fid,' ........................................................................... \r\n');
    fprintf(fid,'Energy Difference of 45 Degree = %9.5f \r\n',Energy_45_Degree);
    fprintf(fid,'Local Homogeneity Difference of 45 Degree = %9.5f \r\n', IDM_45_Degree);
    fprintf(fid,' ........................................................................... \r\n');
    fprintf(fid,'Energy Difference of 90 Degree = %9.5f \r\n',Energy_90_Degree);
    fprintf(fid,'Local Homogeneity Difference of 90 Degree = %9.5f \r\n', IDM_90_Degree);
    fprintf(fid,' ........................................................................... \r\n');
    fprintf(fid,'Energy Difference of 135 Degree = %9.5f \r\n',Energy_135_Degree);
    fprintf(fid,'Local Homogeneity Difference of 135 Degree = %9.5f \r\n', IDM_135_Degree);
    fprintf(fid,' ........................................................................... \r\n');
    fprintf(fid,' ........................................................................... \r\n');
    fprintf(fid,' ........................................................................... \r\n');
    fprintf(fid,'Average of Total Energy Difference of All Directions = %9.5f \r\n',Energy_AVG);
    fprintf(fid,'Average of Total Local Homogeneity Difference of All Directions = %9.5f \r\n', IDM_AVG);
    
    fprintf(fid,' \r\n');
    fprintf(fid,'\r\n');
    fprintf(fid,' \r\n');
 
%-------------------------------------------------------------------------%
    

end

fclose(fid);
    