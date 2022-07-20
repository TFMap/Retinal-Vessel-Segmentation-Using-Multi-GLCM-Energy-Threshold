%function [ F ] = haralick( glcm )
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
fx5=strcat('GLCM_PROPERTIES_OF_RETINAL_IMAGES_GREENCHANNEL453333jj', '.txt');
fName=fx5; %# A file name
fid = fopen(fName,'w');            %# Open the file

fx8=('GLCM_PROPERTIES_OF_RETINAL_IMAGES_GREENCHANNEL45333JJJ'); %ImagesTestSegREAL_AVERAGE_COMBINEMORPHOLOGY3
fprintf(fid,'___________________________\r\n');
fprintf(fid,'%s\r\n',fx8);  %# Print the string
fprintf(fid,'___________________________\r\n');
for kat=1:20
% original file name
  if (kat<=9)
    str = int2str(kat);
    fi=strcat('0',str, '_test.tif');
    u_filename = fi;
  else
    str = int2str(kat);
    fi=strcat(str, '_test.tif');
    u_filename = fi;
  end;
  
% read files 
F = double(imread(u_filename)) / 255;

% figure, imshow([u_GT, u_bw]) 
H=F(:,:,1);
I=F(:,:,2);
J=F(:,:,3);
%GLCM2 = graycomatrix(I,'Offset',[2 0;0 2]);
%str = int2str(k);
%fi=strcat('COMPUTATION FOR IMAGE:',str);
disp(['COMPUTATION FOR IMAGE: ', num2str(kat)]);
%fprintf('COMPUTATION FOR IMAGE:\n',kat);
fprintf('These offsets define pixel relationships of direction 0 degree and distances 1-4'); 
GLCM2 = graycomatrix(I,'Offset', [ 0 1; 0 2; 0 3; 0 4]);
stats = graycoprops(GLCM2,{'Contrast','Homogeneity','Correlation','Energy'})
GLCM2 = graycomatrix(I,'Offset', [0 1]);
%disp (' ENTROPY: ');
I1= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [0 2]);
%disp (' ENTROPY: ');
I2= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [0 3]);
%disp (' ENTROPY: ');
I3= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [0 4]);
%disp (' ENTROPY: ');
I4= entropy(GLCM2);
fprintf('         Entropy: [');
fprintf('%6.4f ', I1, I2,I3, I4);
fprintf(']\n');
fprintf('\n');
fprintf('\n');
fprintf('These offsets define pixel relationships of direction 45 degree and distances 1-4');
GLCM2 = graycomatrix(I,'Offset', [-1 1; -2 2; -3 3; -4 4]);
stats = graycoprops(GLCM2,{'Contrast','Homogeneity','Correlation','Energy'})
GLCM2 = graycomatrix(I,'Offset', [-1 1]);
%disp (' ENTROPY: ');
I1= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [-2 2]);
%disp (' ENTROPY: ');
I2= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [-3 3]);
%disp (' ENTROPY: ');
I3= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [-4 4]);
%disp (' ENTROPY: ');
I4= entropy(GLCM2);
fprintf('         Entropy: [');
fprintf('%6.4f ', I1, I2,I3, I4);
fprintf(']\n');
fprintf('\n');
fprintf('\n');
fprintf('These offsets define pixel relationships of direction 90 degree and distances 1-4');
GLCM2 = graycomatrix(I,'Offset', [-1 0; -2 0; -3 0; -4 0]);
stats = graycoprops(GLCM2,{'Contrast','Homogeneity','Correlation','Energy'})
GLCM2 = graycomatrix(I,'Offset', [-1 0]);
%disp (' ENTROPY: ');
I1= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [-2 0]);
%disp (' ENTROPY: ');
I2= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [-3 0]);
%disp (' ENTROPY: ');
I3= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [-4 0]);
%disp (' ENTROPY: ');
I4= entropy(GLCM2);
fprintf('         Entropy: [');
fprintf('%6.4f ', I1, I2,I3, I4);
fprintf(']\n');
fprintf('\n');
fprintf('\n');
fprintf('These offsets define pixel relationships of direction 135 degree and distances 1-4');
GLCM2 = graycomatrix(I,'Offset', [-1 -1; -2 -2; -3 -3; -4 -4]);
stats = graycoprops(GLCM2,{'Contrast','Homogeneity','Correlation','Energy'})
GLCM2 = graycomatrix(I,'Offset', [-1 -1]);
%disp (' ENTROPY: ');
I1= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [-2 -2]);
%disp (' ENTROPY: ');
I2= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [-3 -3]);
%disp (' ENTROPY: ');
I3= entropy(GLCM2);
GLCM2 = graycomatrix(I,'Offset', [-4 -4]);
%disp (' ENTROPY: ');
I4= entropy(GLCM2);
fprintf('         Entropy: [');
fprintf('%6.4f ', I1, I2,I3, I4);
fprintf(']\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('\n');
end
fclose(fid); 



