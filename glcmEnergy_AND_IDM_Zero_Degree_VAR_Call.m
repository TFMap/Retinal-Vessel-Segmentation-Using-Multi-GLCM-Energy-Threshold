function [Energy_IQR IDM_IQR]=glcmEnergy_AND_IDM_Zero_Degree_VAR_Call(IG2)
% HARALICK Fast Calculation of Haralick Features

 %Energy_Var= var(Energy_Var_All(:)); 
 %IDM_Var= var(IDM_Var_All(:)); 

 %Energy_Var_Min = min(Energy_Var_All);
 %Energy_Var_Max = max(Energy_Var_All);
 %Energy_Var_Mean = mean(Energy_Var_All);
 
 %IDM_Var_Min = min(IDM_Var_All);
 %IDM_Var_Max = max(IDM_Var_All);
 %IDM_Var_Mean = mean(IDM_Var_All);
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
I=IG2;
 Energy_Var_All= zeros(16,1);
 IDM_Var_All= zeros(16,1);
Energy_val = zeros(4,1);
IDM_val = zeros(4,1);
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
glcm = graycomatrix(I,'Offset',[0 1]);
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
Energy_val(1)=Energy;
Energy_Var_All(1)= Energy_val(1);
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
 IDM_Var_All(1)= IDM_val(1);

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
 
   
%-------------------------------------------------------------------------%
    
%F = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13];



glcm = graycomatrix(I,'Offset',[0 2]);
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
Energy_Var_All(2)= Energy_val(2);
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
 IDM_Var_All(2)= IDM_val(2);
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
 
   
%-------------------------------------------------------------------------%
    


glcm = graycomatrix(I,'Offset',[0 3]);
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
Energy_Var_All(3)= Energy_val(3);
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
 IDM_Var_All(3)= IDM_val(3);
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
 
    
%-------------------------------------------------------------------------%
    
glcm = graycomatrix(I,'Offset',[0 4]);
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
Energy_Var_All(4)= Energy_val(4);
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
 IDM_Var_All(4)= IDM_val(4);
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
 
     
    
 %Energy_Var_All(1)= var(Energy_val(:)); 
 %IDM_Var_All(1)= var(IDM_val(:));                      %  var(u(:))
 
 
 
 
 
 
 
 
 
 
 
 Energy_val45 = zeros(4,1);
IDM_val45 = zeros(4,1);
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
glcm = graycomatrix(I,'Offset',[-1 1]);
%I = (rgb2gray( imread('rock sample2.jpg')));
%glcm = graycomatrix(I,'Offset',[0 1]);
%glcm = [0 1 2 3;1 1 2 3;1 0 2 0;0 0 0 3];
%'Offset', [-1 1; -2 2; -3 3; -4 4]
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
Energy_val45(1)=Energy;
Energy_Var_All(5)= Energy_val45(1);
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
 IDM_val45(1)= Inverse_Difference_Moment;
 IDM_Var_All(5)= IDM_val45(1);
 

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
 
    
%-------------------------------------------------------------------------%
    
%F = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13];



glcm = graycomatrix(I,'Offset',[-2 2]);
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
Energy_val45(2) = Energy;
Energy_Var_All(6)= Energy_val45(2);
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
 IDM_val45(2)= Inverse_Difference_Moment;
 IDM_Var_All(6)= IDM_val45(2);
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
 
    
%-------------------------------------------------------------------------%
    


glcm = graycomatrix(I,'Offset',[-3 3]);
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
Energy_val45(3) = Energy;
Energy_Var_All(7)= Energy_val45(3);
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
 IDM_val45(3)= Inverse_Difference_Moment;
 IDM_Var_All(7)= IDM_val45(3);
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
 
    
%-------------------------------------------------------------------------%
    
glcm = graycomatrix(I,'Offset',[-4 4]);
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
Energy_val45(4) = Energy;
Energy_Var_All(8)= Energy_val45(4);
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
 IDM_val45(4)= Inverse_Difference_Moment;
 IDM_Var_All(8)= IDM_val45(4);
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
 
 
  %Energy_Var_All(2)= var(Energy_val45(:)); 
  %IDM_Var_All(2)= var(IDM_val45(:));   
 
 
 
 
 
 
 
 
 
 
 Energy_val90 = zeros(4,1);
IDM_val90 = zeros(4,1);

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
glcm = graycomatrix(I,'Offset',[-1 0]);
%I = (rgb2gray( imread('rock sample2.jpg')));
%glcm = graycomatrix(I,'Offset',[0 1]);
%glcm = [0 1 2 3;1 1 2 3;1 0 2 0;0 0 0 3];
%'Offset', [-1 1; -2 2; -3 3; -4 4]
% 'Offset', [-1 0; -2 0; -3 0; -4 0]);
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
Energy_val90(1)=Energy;
Energy_Var_All(9)= Energy_val90(1);
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
 IDM_val90(1)= Inverse_Difference_Moment;
 IDM_Var_All(9)= IDM_val90(1);

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
 

%-------------------------------------------------------------------------%
    
%F = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13];



glcm = graycomatrix(I,'Offset',[-2 0]);
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
Energy_val90(2) = Energy;
Energy_Var_All(10)= Energy_val90(2);
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
 IDM_val90(2)= Inverse_Difference_Moment;
 IDM_Var_All(10)= IDM_val90(2);
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
 

%-------------------------------------------------------------------------%
    


glcm = graycomatrix(I,'Offset',[-3 0]);
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
Energy_val90(3) = Energy;
Energy_Var_All(11)= Energy_val90(3);
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
 IDM_val90(3)= Inverse_Difference_Moment;
 IDM_Var_All(11)= IDM_val90(3);
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
 

%-------------------------------------------------------------------------%
    
glcm = graycomatrix(I,'Offset',[-4 0]);
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
Energy_val90(4) = Energy;
Energy_Var_All(12)= Energy_val90(4);
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
 IDM_val90(4)= Inverse_Difference_Moment;
 IDM_Var_All(12)= IDM_val90(4);
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
 
 
 %Energy_Var_All(3)= var(Energy_val90(:)); 
 %IDM_Var_All(3)= var(IDM_val90(:)); 
  
 
 
 
 
 
 
 
 
 
 
 
 
 Energy_val135 = zeros(4,1);
IDM_val135 = zeros(4,1);

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
Energy_val135(1)=Energy;
Energy_Var_All(13)= Energy_val135(1);
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
 IDM_val135(1)= Inverse_Difference_Moment;
 IDM_Var_All(13)= IDM_val135(1);

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
Energy_val135(2) = Energy;
Energy_Var_All(14)= Energy_val135(2);
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
 IDM_val135(2)= Inverse_Difference_Moment;
 IDM_Var_All(14)= IDM_val135(2);
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
Energy_val135(3) = Energy;
Energy_Var_All(15)= Energy_val135(3);
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
 IDM_val135(3)= Inverse_Difference_Moment;
 IDM_Var_All(15)= IDM_val135(3);
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
Energy_val135(4) = Energy;
Energy_Var_All(16)= Energy_val135(4);
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
 IDM_val135(4)= Inverse_Difference_Moment;
 IDM_Var_All(16)= IDM_val135(4);
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
 
  %Energy_Var_All= zeros(16,1);
  %IDM_Var_All= zeros(16,1);
  
 Energy_IQR= range(Energy_Var_All)/2; 
 IDM_IQR= range(IDM_Var_All)/2; 
 
% Energy_Var_Min = min(Energy_Var_All);
% Energy_Var_Max = max(Energy_Var_All);
% Energy_Var_Mean = mean(Energy_Var_All);
 
 %IDM_Var_Min = min(IDM_Var_All);
 %IDM_Var_Max = max(IDM_Var_All);
 %IDM_Var_Mean = mean(IDM_Var_All);
 
 
    