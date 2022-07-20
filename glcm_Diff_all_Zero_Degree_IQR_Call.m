function [ Mean_Contrast_val_Zero IQR_Diff_Var_val_ZERO_Degree IQR_Sum_AVG_val_ZERO_Degree IQR_Sum_Entropy_val_ZERO_Degree IQR_Contrast_val_ZERO_Degree IQR_Diff_Entropy_val_ZERO_Degree IQR_Correl_2_val_ZERO_Degree IQR_Sum_Var_val_ZERO_Degree IQR_Entropy_val_ZERO_Degree ]=glcm_Diff_all_Zero_Degree_IQR_Call(IG2)
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
I=IG2;
Sum_AVG_val = zeros(4,1);
Diff_Var_val= zeros(4,1);
Entropy_val= zeros(4,1);
Diff_Entropy_val= zeros(4,1);
Sum_Var_val = zeros(4,1);
Sum_Entropy_val = zeros(4,1);
Correl_2_val = zeros(4,1);
Contrast_val = zeros(4,1);
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
 Sum_Var_val(1)=  Sum_of_Variances;
 %Sum_Var_val(1)=Sum_of_Variances  

% Inverse Difference Moment
 f5 = sum(f_5(:));
 Inverse_Difference_Moment = f5;
 %IDM_val(1)= Inverse_Difference_Moment;

% Entropy
 f9 = -sum(f_9(:));
 Entropy = f9;
 Entropy_val(1)= Entropy;

% Information Measures of Correlation 1&2
HX = -sum(HX_);
HY = -sum(HY_);
HXY = Entropy;
HXY1 = -sum(HXY_1(:));
HXY2 = -sum(HXY_2(:));

 f12 = (HXY-HXY1)/max([HX, HY]);
 Correlation_1 = f12;
 Correlation_val(1)=Correlation_1 ;
 f13 = (1 - exp((-2)*(HXY2 - HXY)))^0.5;
 Correlation_2 = f13;
 Correl_2_val(1)= Correlation_2;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_6(i) = i*pxy(i);
    f_8(i) = pxy(i)*log(pxy(i)+eps);
end

% Sum Average
%f_6(1) = [];       % not necessary f_6(1) is zero anyway
 f6 = sum(f_6);
 Sum_Average =f6;
 Sum_AVG_val(1)= Sum_Average;

% Sum Entropy
%f_8(1)=[];         % not necessary f_8(1) is zero anyway
f8 = -sum(f_8);
Sum_Entropy = f8;
Sum_Entropy_val(1) = Sum_Entropy;


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
 Diff_Var_val(1)= Difference_Variance;

%-------------------------------------------------------------------------%

for k=1:S
    f_2(k) = (k-1)^2*px_y(k);
    f_11(k) = px_y(k)*log(px_y(k)+eps);
end


% Contrast
f2 = sum(f_2(:));
Contrast =f2;
Contrast_val(1) = Contrast;

% Difference Entropy
 f11 = -sum(f_11);
 Difference_Entropy = f11;
 Diff_Entropy_val(1) = Difference_Entropy;
 
   
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
 Sum_Var_val(2)=  Sum_of_Variances2;
 %Sum_of_Variances = Sum_of_Variances - Sum_of_Variances2;
% Inverse Difference Moment
 f5 = sum(f_5(:));
 Inverse_Difference_Moment = f5;
 IDM_val(2)= Inverse_Difference_Moment;
 %Inverse_Difference_Moment = Inverse_Difference_Moment - Inverse_Difference_Moment2;

% Entropy
 f9 = -sum(f_9(:));
 Entropy2 = f9;
 Entropy_val(2)= Entropy2;
 %Entropy = Entropy - Entropy2;

% Information Measures of Correlation 1&2
HX = -sum(HX_);
HY = -sum(HY_);
HXY = Entropy;
HXY1 = -sum(HXY_1(:));
HXY2 = -sum(HXY_2(:));

 f12 = (HXY-HXY1)/max([HX, HY]);
 Correlation_12 = f12;
 Correlation_val(2)= Correlation_12;
 f13 = (1 - exp((-2)*(HXY2 - HXY)))^0.5;
 Correlation_22 = f13;
 Correl_2_val(2)= Correlation_22;
 %Correlation_2 = Correlation_2- Correlation_22;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_6(i) = i*pxy(i);
    f_8(i) = pxy(i)*log(pxy(i)+eps);
end

% Sum Average
%f_6(1) = [];       % not necessary f_6(1) is zero anyway
 f6 = sum(f_6);
 Sum_Average2 =f6;
 Sum_AVG_val(2)= Sum_Average2;

% Sum Entropy
%f_8(1)=[];         % not necessary f_8(1) is zero anyway
f8 = -sum(f_8);
Sum_Entropy2 = f8;
Sum_Entropy_val(2) = Sum_Entropy2;



%-------------------------------------------------------------------------%

for i=2:2*S
    f_7(i)=(i-f8)^2*pxy(i);
end

% Sum Variance
%f_7(1) = [];       % not necessary f_7(1) is zero anyway
 f7 = sum(f_7);
 Sum_Variance2 = f7;
 

% Difference Variance
 f10 = var(px_y);
 Difference_Variance2 = f10;
 %Difference_Variance = Difference_Variance - Difference_Variance2;
 Diff_Var_val(2)= Difference_Variance2;

%-------------------------------------------------------------------------%

for k=1:S
    f_2(k) = (k-1)^2*px_y(k);
    f_11(k) = px_y(k)*log(px_y(k)+eps);
end


% Contrast
f2 = sum(f_2(:));
Contrast2 =f2;
Contrast_val(2) = Contrast2;
%Contrast = Contrast - Contrast2;

% Difference Entropy
 f11 = -sum(f_11);
 Difference_Entropy2 = f11;
 Diff_Entropy_val(2) = Difference_Entropy2;
 %Difference_Entropy = Difference_Entropy - Difference_Entropy2;
 
   
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
 Sum_of_Variances3 = f4;
 Sum_Var_val(3)=  Sum_of_Variances3;
 %Sum_of_Variances = Sum_of_Variances - Sum_of_Variances2;
% Inverse Difference Moment
 f5 = sum(f_5(:));
 Inverse_Difference_Moment = f5;
 IDM_val(3)= Inverse_Difference_Moment;
 %Inverse_Difference_Moment = Inverse_Difference_Moment - Inverse_Difference_Moment2;

% Entropy
 f9 = -sum(f_9(:));
 Entropy3 = f9;
 Entropy_val(3)= Entropy3;
 %Entropy = Entropy - Entropy2;

% Information Measures of Correlation 1&2
HX = -sum(HX_);
HY = -sum(HY_);
HXY = Entropy;
HXY1 = -sum(HXY_1(:));
HXY2 = -sum(HXY_2(:));

 f12 = (HXY-HXY1)/max([HX, HY]);
 Correlation_13 = f12;
 Correlation_val(3)= Correlation_13;
 f13 = (1 - exp((-2)*(HXY2 - HXY)))^0.5;
 Correlation_23 = f13;
 Correl_2_val(3)= Correlation_23;
 %Correlation_2 = Correlation_2- Correlation_22;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_6(i) = i*pxy(i);
    f_8(i) = pxy(i)*log(pxy(i)+eps);
end

% Sum Average
%f_6(1) = [];       % not necessary f_6(1) is zero anyway
 f6 = sum(f_6);
 Sum_Average3 =f6;
 Sum_AVG_val(3)= Sum_Average3;
  
% Sum Entropy
%f_8(1)=[];         % not necessary f_8(1) is zero anyway
f8 = -sum(f_8);
Sum_Entropy3 = f8;
Sum_Entropy_val(3) = Sum_Entropy3;
%Sum_Entropy = Sum_Entropy - Sum_Entropy2;


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
 Difference_Variance3 = f10;
 Diff_Var_val(3)= Difference_Variance3;

%-------------------------------------------------------------------------%

for k=1:S
    f_2(k) = (k-1)^2*px_y(k);
    f_11(k) = px_y(k)*log(px_y(k)+eps);
end


% Contrast
f2 = sum(f_2(:));
Contrast3 =f2;
Contrast_val(3) = Contrast3;
%Contrast = Contrast - Contrast2;

% Difference Entropy
 f11 = -sum(f_11);
 Difference_Entropy3 = f11;
 Diff_Entropy_val(3) = Difference_Entropy3;
 %Difference_Entropy = Difference_Entropy - Difference_Entropy2;
 
    
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
 Sum_of_Variances4 = f4;
 Sum_Var_val(4)=  Sum_of_Variances4;
 %Sum_of_Variances = Sum_of_Variances - Sum_of_Variances2;
 
% Inverse Difference Moment
 f5 = sum(f_5(:));
 Inverse_Difference_Moment = f5;
 IDM_val(4)= Inverse_Difference_Moment;
 %Inverse_Difference_Moment = Inverse_Difference_Moment - Inverse_Difference_Moment2;

% Entropy
 f9 = -sum(f_9(:));
 Entropy4 = f9;
 Entropy_val(4)= Entropy4;
 %Entropy = Entropy - Entropy2;

% Information Measures of Correlation 1&2
HX = -sum(HX_);
HY = -sum(HY_);
HXY = Entropy;
HXY1 = -sum(HXY_1(:));
HXY2 = -sum(HXY_2(:));

 f12 = (HXY-HXY1)/max([HX, HY]);
 Correlation_14 = f12;
 Correlation_val(4)= Correlation_14;
 f13 = (1 - exp((-2)*(HXY2 - HXY)))^0.5;
 Correlation_24 = f13;
 Correl_2_val(4)= Correlation_24;
 %Correlation_2 = Correlation_2- Correlation_22;


%-------------------------------------------------------------------------%

for i=2:2*S
    f_6(i) = i*pxy(i);
    f_8(i) = pxy(i)*log(pxy(i)+eps);
end

% Sum Average
%f_6(1) = [];       % not necessary f_6(1) is zero anyway
 f6 = sum(f_6);
 Sum_Average4 =f6;
 Sum_AVG_val(4)= Sum_Average4;
  
% Sum Entropy
%f_8(1)=[];         % not necessary f_8(1) is zero anyway
f8 = -sum(f_8);
Sum_Entropy4 = f8;
Sum_Entropy_val(4) = Sum_Entropy4;

%Sum_Entropy = Sum_Entropy - Sum_Entropy2;


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
 Difference_Variance4 = f10;
 Diff_Var_val(4)= Difference_Variance4;

%-------------------------------------------------------------------------%

for k=1:S
    f_2(k) = (k-1)^2*px_y(k);
    f_11(k) = px_y(k)*log(px_y(k)+eps);
end


% Contrast
f2 = sum(f_2(:));
Contrast4 =f2;
Contrast_val(4) = Contrast4;
%Contrast = Contrast - Contrast2;

% Difference Entropy
 f11 = -sum(f_11);
 Difference_Entropy4 = f11;
 Diff_Entropy_val(4) = Difference_Entropy4;
 %Difference_Entropy = Difference_Entropy - Difference_Entropy2;
 
     
    
 %Energy_Zero_Degree= max(Energy_val)-min(Energy_val); 
 %IDM_Zero_Degree= max(IDM_val)-min(IDM_val); 
 %Correl_1_Degree=max(Correlation_val)-min(Correlation_val)
 IQR_Diff_Var_val_ZERO_Degree = iqr(Diff_Var_val);
 IQR_Sum_AVG_val_ZERO_Degree = abs(iqr(Sum_AVG_val));
 IQR_Sum_Entropy_val_ZERO_Degree = abs(iqr(Sum_Entropy_val));
 IQR_Contrast_val_ZERO_Degree = abs(iqr(Contrast_val));
 IQR_Diff_Entropy_val_ZERO_Degree = abs(iqr(Diff_Entropy_val));
 IQR_Correl_2_val_ZERO_Degree = abs(iqr(Correl_2_val));
 IQR_Sum_Var_val_ZERO_Degree = abs(iqr(Sum_Var_val));
 IQR_Entropy_val_ZERO_Degree = abs(iqr(Entropy_val));
 Mean_Contrast_val_Zero=Contrast_val/4;