%function [Energy_135_Degree IDM_135_Degree]=glcmEnergy_AND_IDM_135_Degree(IG2)
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
fx5=strcat('GLCM_PROPERTIES FOR DIRECTION ONE HUNDRED AND THIRTY FIVE DEGREE IN ALL DISTANCES', '.txt');
fName=fx5; %# A file name
fid = fopen(fName,'w');            %# Open the file

fx8=('GLCM_PROPERTIES FOR DIRECTION 135 DEGREE IN ALL DISTANCES'); %ImagesTestSegREAL_AVERAGE_COMBINEMORPHOLOGY3
fprintf(fid,'___________________________\r\n');
fprintf(fid,'%s\r\n',fx8);  %# Print the string
fprintf(fid,'___________________________\r\n');
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
I=F(:,:,2);
J=F(:,:,3);
f1=0;
f2=0;
f3=0;
f4=0;
f5=0;
f6=0;
f7=0;
f8=0;
f9=0;
f10=0;
f11=0;
f12=0;
f13=0;
% [ 0 1; 0 2; 0 3; 0 4]);
glcm = graycomatrix(I,'Offset',[-1 -1]);
%I = (rgb2gray( imread('rock sample2.jpg')));
%glcm = graycomatrix(I,'Offset',[0 1]);
%glcm = [0 1 2 3;1 1 2 3;1 0 2 0;0 0 0 3];
%'Offset', [-1 1; -2 2; -3 3; -4 4]
% 'Offset', [-1 0; -2 0; -3 0; -4 0]);
% 'Offset', [-1 -1; -2 -2; -3 -3; -4 -4])
S=size(glcm,1);

f_2=zeros(S);
f_3=zeros(S);
f_4=zeros(S);
f_5=zeros(S);
f_6=zeros(1,2*S);
f_7=zeros(1,2*S);
f_8=zeros(1,2*S);
f_9=zeros(S);
f_11=zeros(1,S);

pxy=zeros(1,2*S);
px_y=zeros(1,S);

HX_=zeros(1,S);
HY_=zeros(1,S);
HXY_1=zeros(S);
HXY_2=zeros(S);

%% CALCULATION
% Normalization
M = glcm/sum(glcm(:));

% Energy
f_1 = M.^2;
f1 = sum(f_1(:));
Energy = f1;
Energy_val(1)=Energy;
%-------------------------------------------------------------------------%

u = mean2(M);
py = sum(M,1);
px = sum(M,2);

for i=1:S
    for j=1:S
       
        f_3(i,j) = i*j*M(i,j);
        f_4(i,j) = (i-u)^2*M(i,j);
        f_5(i,j) = M(i,j)/(1+(i-j)^2);
        f_9(i,j) = M(i,j)*log(M(i,j)+eps);
    
        pxy(i+j) = pxy(i+j)+M(i,j);
        px_y(abs(i-j)+1) = px_y(abs(i-j)+1)+M(i,j);
             
        HX_(i)= px(i)*log(px(i)+eps);
        HY_(j)= py(j)*log(py(j)+eps);
        HXY_1(i,j) = M(i,j)*log(px(i)*py(j)+eps);
        HXY_2(i,j) = px(i)*py(j)*log(px(i)*py(j)+eps);
    end
end


% Correlation
ux = mean(px); sx=std(px);
uy = mean(py); sy=std(py);
f3 =(sum(f_3(:))-(ux*uy))/(sx*sy);
Correlation = f3;

% Sum of Variances
 f4 = sum(f_4(:));
 Sum_of_Variances = f4;

% Inverse Difference Moment
 f5 = sum(f_5(:));
 Inverse_Difference_Moment = f5;
 IDM_val(1)= Inverse_Difference_Moment;

% Entropy
 f9 = -sum(f_9(:));
 Entropy = f9;

% Information Measures of Correlation 1&2
HX = -sum(HX_);
HY = -sum(HY_);
HXY = Entropy;
HXY1 = -sum(HXY_1(:));
HXY2 = -sum(HXY_2(:));

 f12 = (HXY-HXY1)/max([HX, HY]);
 Correlation_1 = f12;
 f13 = (1 - exp((-2)*(HXY2 - HXY)))^0.5;
 Correlation_2 = f13;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_6(i) = i*pxy(i);
    f_8(i) = pxy(i)*log(pxy(i)+eps);
end

% Sum Average
%f_6(1) = [];       % not necessary f_6(1) is zero anyway
 f6 = sum(f_6);
 Sum_Average =f6;

% Sum Entropy
%f_8(1)=[];         % not necessary f_8(1) is zero anyway
f8 = -sum(f_8);
Sum_Entropy = f8;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_7(i)=(i-f8)^2*pxy(i);
end

% Sum Variance
%f_7(1) = [];       % not necessary f_7(1) is zero anyway
 f7 = sum(f_7);
 Sum_Variance = f7;

% Difference Variance
 f10 = var(px_y);
 Difference_Variance = f10;

%-------------------------------------------------------------------------%

for k=1:S
    f_2(k) = (k-1)^2*px_y(k);
    f_11(k) = px_y(k)*log(px_y(k)+eps);
end


% Contrast
f2 = sum(f_2(:));
Contrast =f2;

% Difference Entropy
 f11 = -sum(f_11);
 Difference_Entropy = f11;
 
    str = int2str(kat);
    fx2=strcat('MEASUREMENT OF DIRECTION 135 DEGREE IN DISTANCE ONE OF IMAGE: ',str);
    fprintf(fid,'_____________________________ \r\n');
    fprintf(fid,'%s\r\n',fx2);  %# Print the string
    fprintf(fid,'_____________________________ \r\n');
    fprintf(fid,'Energy = %9.5f \r\n',Energy);
    fprintf(fid,'Correlation = %9.5f \r\n', Correlation);
    fprintf(fid,'Sum_of_Variances  = %9.5f \r\n',Sum_of_Variances );
    fprintf(fid,'Inverse_Difference_Moment = %9.5f \r\n' , Inverse_Difference_Moment);
    fprintf(fid,'Entropy = %9.5f \r\n', Entropy);
    fprintf(fid,'Correlation_1 = %9.5f \r\n', Correlation_1);
    fprintf(fid,'Correlation_2 = %9.5f \r\n', Correlation_2);
    fprintf(fid,'Sum_Average = %9.5f \r\n', Sum_Average);
    fprintf(fid,'Sum_Entropy = %9.5f \r\n', Sum_Entropy);
    fprintf(fid,'Sum_Variance = %9.5f \r\n', Sum_Variance);
    fprintf(fid,'Difference_Variance = %9.5f \r\n', Difference_Variance);
    fprintf(fid,'Difference_Entropy = %9.5f \r\n', Difference_Entropy);
    fprintf(fid,'Contrast = %9.5f \r\n', Contrast);
    fprintf(fid,' \r\n');
    fprintf(fid,'\r\n');
    fprintf(fid,' \r\n');
 
%-------------------------------------------------------------------------%
    
%F = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13];



glcm = graycomatrix(I,'Offset',[-2 -2]);
%I = (rgb2gray( imread('rock sample2.jpg')));
%glcm = graycomatrix(I,'Offset',[0 1]);
%glcm = [0 1 2 3;1 1 2 3;1 0 2 0;0 0 0 3];
S=size(glcm,1);

f_2=zeros(S);
f_3=zeros(S);
f_4=zeros(S);
f_5=zeros(S);
f_6=zeros(1,2*S);
f_7=zeros(1,2*S);
f_8=zeros(1,2*S);
f_9=zeros(S);
f_11=zeros(1,S);

pxy=zeros(1,2*S);
px_y=zeros(1,S);

HX_=zeros(1,S);
HY_=zeros(1,S);
HXY_1=zeros(S);
HXY_2=zeros(S);

%% CALCULATION
% Normalization
M = glcm/sum(glcm(:));

% Energy
f_1 = M.^2;
f1 = sum(f_1(:));
Energy = f1;
Energy_val(2) = Energy;
%-------------------------------------------------------------------------%

u = mean2(M);
py = sum(M,1);
px = sum(M,2);

for i=1:S
    for j=1:S
       
        f_3(i,j) = i*j*M(i,j);
        f_4(i,j) = (i-u)^2*M(i,j);
        f_5(i,j) = M(i,j)/(1+(i-j)^2);
        f_9(i,j) = M(i,j)*log(M(i,j)+eps);
    
        pxy(i+j) = pxy(i+j)+M(i,j);
        px_y(abs(i-j)+1) = px_y(abs(i-j)+1)+M(i,j);
             
        HX_(i)= px(i)*log(px(i)+eps);
        HY_(j)= py(j)*log(py(j)+eps);
        HXY_1(i,j) = M(i,j)*log(px(i)*py(j)+eps);
        HXY_2(i,j) = px(i)*py(j)*log(px(i)*py(j)+eps);
    end
end


% Correlation
ux = mean(px); sx=std(px);
uy = mean(py); sy=std(py);
f3 =(sum(f_3(:))-(ux*uy))/(sx*sy);
Correlation2 = f3;
Correlation = Correlation - Correlation2;

% Sum of Variances
 f4 = sum(f_4(:));
 Sum_of_Variances2 = f4;
 Sum_of_Variances = Sum_of_Variances - Sum_of_Variances2;
% Inverse Difference Moment
 f5 = sum(f_5(:));
 Inverse_Difference_Moment = f5;
 IDM_val(2)= Inverse_Difference_Moment;
 %Inverse_Difference_Moment = Inverse_Difference_Moment - Inverse_Difference_Moment2;

% Entropy
 f9 = -sum(f_9(:));
 Entropy2 = f9;
 Entropy = Entropy - Entropy2;

% Information Measures of Correlation 1&2
HX = -sum(HX_);
HY = -sum(HY_);
HXY = Entropy;
HXY1 = -sum(HXY_1(:));
HXY2 = -sum(HXY_2(:));

 f12 = (HXY-HXY1)/max([HX, HY]);
 Correlation_12 = f12;
 Correlation_1 = Correlation_1 - Correlation_12;
 f13 = (1 - exp((-2)*(HXY2 - HXY)))^0.5;
 Correlation_22 = f13;
 Correlation_2 = Correlation_2- Correlation_22;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_6(i) = i*pxy(i);
    f_8(i) = pxy(i)*log(pxy(i)+eps);
end

% Sum Average
%f_6(1) = [];       % not necessary f_6(1) is zero anyway
 f6 = sum(f_6);
 Sum_Average2 =f6;
 Sum_Average = Sum_Average - Sum_Average2;

% Sum Entropy
%f_8(1)=[];         % not necessary f_8(1) is zero anyway
f8 = -sum(f_8);
Sum_Entropy2 = f8;
Sum_Entropy = Sum_Entropy - Sum_Entropy2;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_7(i)=(i-f8)^2*pxy(i);
end

% Sum Variance
%f_7(1) = [];       % not necessary f_7(1) is zero anyway
 f7 = sum(f_7);
 Sum_Variance2 = f7;
 Sum_Variance = Sum_Variance - Sum_Variance2;

% Difference Variance
 f10 = var(px_y);
 Difference_Variance2 = f10;
 Difference_Variance = Difference_Variance - Difference_Variance2;

%-------------------------------------------------------------------------%

for k=1:S
    f_2(k) = (k-1)^2*px_y(k);
    f_11(k) = px_y(k)*log(px_y(k)+eps);
end


% Contrast
f2 = sum(f_2(:));
Contrast2 =f2;
Contrast = Contrast - Contrast2;

% Difference Entropy
 f11 = -sum(f_11);
 Difference_Entropy2 = f11;
 Difference_Entropy = Difference_Entropy - Difference_Entropy2;
 
    str = int2str(kat);
    fx2=strcat('MEASUREMENT OF DIRECTION 135 DEGREE IN DISTANCE TWO OF IMAGE: ',str);
    fprintf(fid,'_____________________________ \r\n');
    fprintf(fid,'%s\r\n',fx2);  %# Print the string
    fprintf(fid,'_____________________________ \r\n');
    fprintf(fid,'Energy = %9.5f \r\n',Energy);
    fprintf(fid,'Correlation = %9.5f \r\n', Correlation);
    fprintf(fid,'Sum_of_Variances  = %9.5f \r\n',Sum_of_Variances );
    fprintf(fid,'Inverse_Difference_Moment = %9.5f \r\n' , Inverse_Difference_Moment);
    fprintf(fid,'Entropy = %9.5f \r\n', Entropy);
    fprintf(fid,'Correlation_1 = %9.5f \r\n', Correlation_1);
    fprintf(fid,'Correlation_2 = %9.5f \r\n', Correlation_2);
    fprintf(fid,'Sum_Average = %9.5f \r\n', Sum_Average);
    fprintf(fid,'Sum_Entropy = %9.5f \r\n', Sum_Entropy);
    fprintf(fid,'Sum_Variance = %9.5f \r\n', Sum_Variance);
    fprintf(fid,'Difference_Variance = %9.5f \r\n', Difference_Variance);
    fprintf(fid,'Difference_Entropy = %9.5f \r\n', Difference_Entropy);
    fprintf(fid,'Contrast = %9.5f \r\n', Contrast);
    fprintf(fid,' \r\n');
    fprintf(fid,'\r\n');
    fprintf(fid,' \r\n');
 
%-------------------------------------------------------------------------%
    


glcm = graycomatrix(I,'Offset',[-3 -3]);
%I = (rgb2gray( imread('rock sample2.jpg')));
%glcm = graycomatrix(I,'Offset',[0 1]);
%glcm = [0 1 2 3;1 1 2 3;1 0 2 0;0 0 0 3];
S=size(glcm,1);

f_2=zeros(S);
f_3=zeros(S);
f_4=zeros(S);
f_5=zeros(S);
f_6=zeros(1,2*S);
f_7=zeros(1,2*S);
f_8=zeros(1,2*S);
f_9=zeros(S);
f_11=zeros(1,S);

pxy=zeros(1,2*S);
px_y=zeros(1,S);

HX_=zeros(1,S);
HY_=zeros(1,S);
HXY_1=zeros(S);
HXY_2=zeros(S);

%% CALCULATION
% Normalization
M = glcm/sum(glcm(:));

% Energy
f_1 = M.^2;
f1 = sum(f_1(:));
Energy = f1;
Energy_val(3) = Energy;
%-------------------------------------------------------------------------%

u = mean2(M);
py = sum(M,1);
px = sum(M,2);

for i=1:S
    for j=1:S
       
        f_3(i,j) = i*j*M(i,j);
        f_4(i,j) = (i-u)^2*M(i,j);
        f_5(i,j) = M(i,j)/(1+(i-j)^2);
        f_9(i,j) = M(i,j)*log(M(i,j)+eps);
    
        pxy(i+j) = pxy(i+j)+M(i,j);
        px_y(abs(i-j)+1) = px_y(abs(i-j)+1)+M(i,j);
             
        HX_(i)= px(i)*log(px(i)+eps);
        HY_(j)= py(j)*log(py(j)+eps);
        HXY_1(i,j) = M(i,j)*log(px(i)*py(j)+eps);
        HXY_2(i,j) = px(i)*py(j)*log(px(i)*py(j)+eps);
    end
end


% Correlation
ux = mean(px); sx=std(px);
uy = mean(py); sy=std(py);
f3 =(sum(f_3(:))-(ux*uy))/(sx*sy);
Correlation2 = f3;
Correlation = Correlation - Correlation2;

% Sum of Variances
 f4 = sum(f_4(:));
 Sum_of_Variances2 = f4;
 Sum_of_Variances = Sum_of_Variances - Sum_of_Variances2;
% Inverse Difference Moment
 f5 = sum(f_5(:));
 Inverse_Difference_Moment = f5;
 IDM_val(3)= Inverse_Difference_Moment;
 %Inverse_Difference_Moment = Inverse_Difference_Moment - Inverse_Difference_Moment2;

% Entropy
 f9 = -sum(f_9(:));
 Entropy2 = f9;
 Entropy = Entropy - Entropy2;

% Information Measures of Correlation 1&2
HX = -sum(HX_);
HY = -sum(HY_);
HXY = Entropy;
HXY1 = -sum(HXY_1(:));
HXY2 = -sum(HXY_2(:));

 f12 = (HXY-HXY1)/max([HX, HY]);
 Correlation_12 = f12;
 Correlation_1 = Correlation_1 - Correlation_12;
 f13 = (1 - exp((-2)*(HXY2 - HXY)))^0.5;
 Correlation_22 = f13;
 Correlation_2 = Correlation_2- Correlation_22;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_6(i) = i*pxy(i);
    f_8(i) = pxy(i)*log(pxy(i)+eps);
end

% Sum Average
%f_6(1) = [];       % not necessary f_6(1) is zero anyway
 f6 = sum(f_6);
 Sum_Average2 =f6;
 Sum_Average = Sum_Average - Sum_Average2;

% Sum Entropy
%f_8(1)=[];         % not necessary f_8(1) is zero anyway
f8 = -sum(f_8);
Sum_Entropy2 = f8;
Sum_Entropy = Sum_Entropy - Sum_Entropy2;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_7(i)=(i-f8)^2*pxy(i);
end

% Sum Variance
%f_7(1) = [];       % not necessary f_7(1) is zero anyway
 f7 = sum(f_7);
 Sum_Variance2 = f7;
 Sum_Variance = Sum_Variance - Sum_Variance2;

% Difference Variance
 f10 = var(px_y);
 Difference_Variance2 = f10;
 Difference_Variance = Difference_Variance - Difference_Variance2;

%-------------------------------------------------------------------------%

for k=1:S
    f_2(k) = (k-1)^2*px_y(k);
    f_11(k) = px_y(k)*log(px_y(k)+eps);
end


% Contrast
f2 = sum(f_2(:));
Contrast2 =f2;
Contrast = Contrast - Contrast2;

% Difference Entropy
 f11 = -sum(f_11);
 Difference_Entropy2 = f11;
 Difference_Entropy = Difference_Entropy - Difference_Entropy2;
 
     str = int2str(kat);
    fx2=strcat('MEASUREMENT OF DIRECTION 135 DEGREE IN DISTANCE THREE OF IMAGE: ',str);
    fprintf(fid,'_____________________________ \r\n');
    fprintf(fid,'%s\r\n',fx2);  %# Print the string
    fprintf(fid,'_____________________________ \r\n');
    fprintf(fid,'Energy = %9.5f \r\n',Energy);
    fprintf(fid,'Correlation = %9.5f \r\n', Correlation);
    fprintf(fid,'Sum_of_Variances  = %9.5f \r\n',Sum_of_Variances );
    fprintf(fid,'Inverse_Difference_Moment = %9.5f \r\n' , Inverse_Difference_Moment);
    fprintf(fid,'Entropy = %9.5f \r\n', Entropy);
    fprintf(fid,'Correlation_1 = %9.5f \r\n', Correlation_1);
    fprintf(fid,'Correlation_2 = %9.5f \r\n', Correlation_2);
    fprintf(fid,'Sum_Average = %9.5f \r\n', Sum_Average);
    fprintf(fid,'Sum_Entropy = %9.5f \r\n', Sum_Entropy);
    fprintf(fid,'Sum_Variance = %9.5f \r\n', Sum_Variance);
    fprintf(fid,'Difference_Variance = %9.5f \r\n', Difference_Variance);
    fprintf(fid,'Difference_Entropy = %9.5f \r\n', Difference_Entropy);
    fprintf(fid,'Contrast = %9.5f \r\n', Contrast);
    fprintf(fid,' \r\n');
    fprintf(fid,'\r\n');
    fprintf(fid,' \r\n');
 
%-------------------------------------------------------------------------%
    
glcm = graycomatrix(I,'Offset',[-4 -4]);
%I = (rgb2gray( imread('rock sample2.jpg')));
%glcm = graycomatrix(I,'Offset',[0 1]);
%glcm = [0 1 2 3;1 1 2 3;1 0 2 0;0 0 0 3];
S=size(glcm,1);

f_2=zeros(S);
f_3=zeros(S);
f_4=zeros(S);
f_5=zeros(S);
f_6=zeros(1,2*S);
f_7=zeros(1,2*S);
f_8=zeros(1,2*S);
f_9=zeros(S);
f_11=zeros(1,S);

pxy=zeros(1,2*S);
px_y=zeros(1,S);

HX_=zeros(1,S);
HY_=zeros(1,S);
HXY_1=zeros(S);
HXY_2=zeros(S);

%% CALCULATION
% Normalization
M = glcm/sum(glcm(:));

% Energy
f_1 = M.^2;
f1 = sum(f_1(:));
Energy = f1;
Energy_val(4) = Energy;
%-------------------------------------------------------------------------%

u = mean2(M);
py = sum(M,1);
px = sum(M,2);

for i=1:S
    for j=1:S
       
        f_3(i,j) = i*j*M(i,j);
        f_4(i,j) = (i-u)^2*M(i,j);
        f_5(i,j) = M(i,j)/(1+(i-j)^2);
        f_9(i,j) = M(i,j)*log(M(i,j)+eps);
    
        pxy(i+j) = pxy(i+j)+M(i,j);
        px_y(abs(i-j)+1) = px_y(abs(i-j)+1)+M(i,j);
             
        HX_(i)= px(i)*log(px(i)+eps);
        HY_(j)= py(j)*log(py(j)+eps);
        HXY_1(i,j) = M(i,j)*log(px(i)*py(j)+eps);
        HXY_2(i,j) = px(i)*py(j)*log(px(i)*py(j)+eps);
    end
end


% Correlation
ux = mean(px); sx=std(px);
uy = mean(py); sy=std(py);
f3 =(sum(f_3(:))-(ux*uy))/(sx*sy);
Correlation2 = f3;
Correlation = Correlation - Correlation2;

% Sum of Variances
 f4 = sum(f_4(:));
 Sum_of_Variances2 = f4;
 Sum_of_Variances = Sum_of_Variances - Sum_of_Variances2;
% Inverse Difference Moment
 f5 = sum(f_5(:));
 Inverse_Difference_Moment = f5;
 IDM_val(4)= Inverse_Difference_Moment;
 %Inverse_Difference_Moment = Inverse_Difference_Moment - Inverse_Difference_Moment2;

% Entropy
 f9 = -sum(f_9(:));
 Entropy2 = f9;
 Entropy = Entropy - Entropy2;

% Information Measures of Correlation 1&2
HX = -sum(HX_);
HY = -sum(HY_);
HXY = Entropy;
HXY1 = -sum(HXY_1(:));
HXY2 = -sum(HXY_2(:));

 f12 = (HXY-HXY1)/max([HX, HY]);
 Correlation_12 = f12;
 Correlation_1 = Correlation_1 - Correlation_12;
 f13 = (1 - exp((-2)*(HXY2 - HXY)))^0.5;
 Correlation_22 = f13;
 Correlation_2 = Correlation_2- Correlation_22;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_6(i) = i*pxy(i);
    f_8(i) = pxy(i)*log(pxy(i)+eps);
end

% Sum Average
%f_6(1) = [];       % not necessary f_6(1) is zero anyway
 f6 = sum(f_6);
 Sum_Average2 =f6;
 Sum_Average = Sum_Average - Sum_Average2;

% Sum Entropy
%f_8(1)=[];         % not necessary f_8(1) is zero anyway
f8 = -sum(f_8);
Sum_Entropy2 = f8;
Sum_Entropy = Sum_Entropy - Sum_Entropy2;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_7(i)=(i-f8)^2*pxy(i);
end

% Sum Variance
%f_7(1) = [];       % not necessary f_7(1) is zero anyway
 f7 = sum(f_7);
 Sum_Variance2 = f7;
 Sum_Variance = Sum_Variance - Sum_Variance2;

% Difference Variance
 f10 = var(px_y);
 Difference_Variance2 = f10;
 Difference_Variance = Difference_Variance - Difference_Variance2;

%-------------------------------------------------------------------------%

for k=1:S
    f_2(k) = (k-1)^2*px_y(k);
    f_11(k) = px_y(k)*log(px_y(k)+eps);
end


% Contrast
f2 = sum(f_2(:));
Contrast2 =f2;
Contrast = Contrast - Contrast2;

% Difference Entropy
 f11 = -sum(f_11);
 Difference_Entropy2 = f11;
 Difference_Entropy = Difference_Entropy - Difference_Entropy2;
 
     str = int2str(kat);
    fx2=strcat('MEASUREMENT OF DIRECTION 135 DEGREE IN DISTANCE FOUR OF IMAGE: ',str);
    fprintf(fid,'_____________________________ \r\n');
    fprintf(fid,'%s\r\n',fx2);  %# Print the string
    fprintf(fid,'_____________________________ \r\n');
    fprintf(fid,'Energy = %9.5f \r\n',Energy);
    fprintf(fid,'Correlation = %9.5f \r\n', Correlation);
    fprintf(fid,'Sum_of_Variances  = %9.5f \r\n',Sum_of_Variances );
    fprintf(fid,'Inverse_Difference_Moment = %9.5f \r\n' , Inverse_Difference_Moment);
    fprintf(fid,'Entropy = %9.5f \r\n', Entropy);
    fprintf(fid,'Correlation_1 = %9.5f \r\n', Correlation_1);
    fprintf(fid,'Correlation_2 = %9.5f \r\n', Correlation_2);
    fprintf(fid,'Sum_Average = %9.5f \r\n', Sum_Average);
    fprintf(fid,'Sum_Entropy = %9.5f \r\n', Sum_Entropy);
    fprintf(fid,'Sum_Variance = %9.5f \r\n', Sum_Variance);
    fprintf(fid,'Difference_Variance = %9.5f \r\n', Difference_Variance);
    fprintf(fid,'Difference_Entropy = %9.5f \r\n', Difference_Entropy);
    fprintf(fid,'Contrast = %9.5f \r\n', Contrast);
    fprintf(fid,' \r\n');
    fprintf(fid,'\r\n');
    fprintf(fid,' \r\n');
 
 
 
 
    
 Energy_135_Degree= max(Energy_val)-min(Energy_val); 
 IDM_135_Degree= max(IDM_val)-min(IDM_val); 
    fprintf(fid,'......................................................... \r\n');
    fprintf(fid,'.........................................................\r\n');
    fprintf(fid,'......................................................... \r\n');
    fprintf(fid,'......................................................... \r\n');
    str55 = int2str(kat);
    fx2=strcat('SUMMARY OF DIFFERENCE FOR DIRECTION 135 DEGREE OF IMAGE ',str55);
    fprintf(fid,'_____________________________ \r\n');
    fprintf(fid,'%s\r\n',fx2);  %# Print the string
    fprintf(fid,'_____________________________ \r\n');
    fprintf(fid,'ENERGY/UNIFORMITY DIFFRENCE = %9.5f \r\n',Energy_135_Degree);
    fprintf(fid,'Inverse Difference Moment/ Local Homogeneity Difference = %9.5f \r\n' , IDM_135_Degree);
    
    fprintf(fid,' \r\n');
    fprintf(fid,'\r\n');
    fprintf(fid,' \r\n');
    fprintf(fid,' \r\n');
    fprintf(fid,'\r\n');
    fprintf(fid,' \r\n');

end

fclose(fid);
    