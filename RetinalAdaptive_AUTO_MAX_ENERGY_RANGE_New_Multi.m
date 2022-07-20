%  Retinal vessel segmentation using Multi-GLCM Energy info-Based thresholding
%  Approach
clc;
mkdir ('RetinalAdaptive_AUTO_MAX_ENERGY_RANGE_New__Multi');
SegmentedImage = 'SegmentedImages';
%sdirectory = 'imagesTrain';
sdirectory = 'ImagesTest';
bmpFiles = dir([sdirectory '/*.tif']);
for k = 1:length(bmpFiles)
  filename = [sdirectory '/' bmpFiles(k).name];
  str = int2str(k);
 fi=strcat(str, '.tif');
  fk=strcat(str, 'BW.tif');
 BW10 = imread(fk);
  ImageInput = imread(filename);
IG2 = rgb2gray(ImageInput);
IG2=im2double(IG2);
[Max_Energy Max_IDM]=glcmEnergy_AND_IDM_Max_Call(IG2);
[IQR_Max_Energy IQR_Max_IDM]=glcmEnergy_AND_IDM_Max_IQR_Call(IG2);
I = imresize(ImageInput, 1.0);
  IG = I(:,:,2);
  IG=imfilter(IG,fspecial('unsharp'),'replicate');
  IG=imfilter(IG,fspecial('average'),'replicate');
   IG = imadjust(IG); 
  IG1=im2double(IG);
  m1= min(median(IG2));
  level = graythresh(IG1);
  IM=IG1;
  ws=13;
  ws2=17;
  C=IQR_Max_Energy;
  C2= Max_Energy;
  tm=1;
if (tm~=0 && tm~=1)
    error('tm must be 0 or 1.');
end

IM=mat2gray(IM);

if tm==0
    mIM=imfilter(IM,fspecial('average',ws),'replicate');
else
    mIM=medfilt2(IM,[ws ws]);
    mIM2=medfilt2(IM,[ws2 ws2]);
end
sIM=mIM-IM;
sIM2=mIM2-IM;
bw=im2bw(sIM,C);
bw2=im2bw(sIM2,C2);

bw=imcomplement(bw);
bw2=imcomplement(bw2);

BW2 = bwareaopen(bw,18,4);
BW22 = bwareaopen(bw2,18,4);

K=~BW2;
K2=~BW22;

BW5 = bwareaopen(K,90,4);
BW52 = bwareaopen(K2,150,8);


 BW7 = BW5 - BW10;
 BW72 = BW52 - BW10;
 
 
 BW77 = BW72 + BW7;
 
 
 BW77= medfilt2(BW77, [2 2]);
 
 figure, imshow(BW77);  title('Best Result Map Quantization Without Masking');
 
 
 imwrite(BW77,['RetinalAdaptive_AUTO_MAX_ENERGY_RANGE_New__Multi/' fi]);
end
