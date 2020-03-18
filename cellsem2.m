%2D cell segmentation for confocal images

workingDir='C:\Users\David\Desktop\confocalSYT121\';
cd(workingDir)

data=bfopen('C:\Users\David\Desktop\confocalSYT121/Image1_slide2_slice9.OIB'); %by adding the square braquets, we make it into a single piece of text. If we didn't do that, the comas separating the three pieces of text would make the different texts appear as arguments from the imread function 

series1 = data{1, 1}; %3x2 matrix containing the pixel data for the three channels

bluechannel = series1{1, 1}; %data from plane 1 (which I think coresponds to te blue channel)
greenchannel = series1{2, 1}; %data from plane 2 (which I think corresponds to the green channel)
redchannel = series1{3, 1}; %data from plane 3 (which I think corresponds to the red channel)
farredchannel = series1{4, 1};

series1_colorMaps = data{1, 3}; %data 1,3 conains color lookup tables for each plane

figure('Name','Four channels');
subplot(2,2,1);
if (isempty(series1_colorMaps{1}))
  colormap(gray);
else
  colormap(series1_colorMaps{1}(1,:));
end
imagesc(bluechannel);
imshow(bluechannel, []);
title('blue channel','Color','blue')

subplot(2,2,2);
if (isempty(series1_colorMaps{2}))
  colormap(gray);
else
  colormap(series1_colorMaps{2}(1,:));
end
imagesc(greenchannel);
imshow(greenchannel, []);
title('green channel','Color','green')

subplot(2,2,3);
if (isempty(series1_colorMaps{3}))
  colormap(gray);
else
  colormap(series1_colorMaps{3}(1,:));
end
imagesc(redchannel);
imshow(redchannel, []);
title('red channel','Color','red');


subplot(2,2,4);
if (isempty(series1_colorMaps{4}))
  colormap(gray);
else
  colormap(series1_colorMaps{4}(1,:));
end
imagesc(farredchannel);
imshow(farredchannel, []);
title('far red channel','Color','magenta');

%% segmentation for the blue channel

%I-Global binarization

%1) median filter
% Prior to the application of any binarization algorithm, the image is
%processed with a median filter to reduce noise while trying to preserve edges

medianblue = medfilt2(bluechannel); %aply the default media filter, with siez [3,3]

figure('Name','Effect of median filter');
subplot(1,2,1);
imagesc(bluechannel);
imshow(bluechannel, []);
title('original image');

subplot(1,2,2);
imagesc(medianblue);
imshow(medianblue, []);
title('after median filter')

%2)CLAHE algorithm
% and its contrast is enhanced with the CLAHE algorithm

figure('Name','Effect of CLAHE algorithm');
subplot(1,2,1);
imagesc(bluechannel);
imshow(bluechannel, []);
title('original image');

claheblue = adapthisteq(medianblue,'ClipLimit',0.0005,'NBins',4000);
subplot(1,2,2);
imagesc(claheblue);
imshow(claheblue, []);
title('after clahe')


%3) erosion/dilation

% use Otsu’s method (Otsu, 1979) to obtain a first
% binary image, which is subject to the typical morphological trans-
% formations (erosion/ dilation) 

mask = [0 1 0; 1 1 1; 0 1 0]; %erosion dilation mask
erodila = imopen(claheblue, mask); %imopen first erodes ad then dilates

figure('Name','Effect of erosion/dilation');
subplot(1,2,1);
imagesc(bluechannel);
imshow(bluechannel, []);
title('original image');

subplot(1,2,2);
imagesc(erodila);
imshow(erodila, []);
title('after ersion/dilation');


%4) binarization
level = graythresh(erodila);
binarized = imbinarize(erodila,level);

figure('Name','Effect of binarization');
subplot(1,2,1);
imagesc(bluechannel);
imshow(bluechannel, []);
title('original image');

subplot(1,2,2);
imagesc(binarized);
imshow(binarized, []);
title('after binarization');



% % cleaning process aiming to remove small structures that come from the
% % noise still present in the image.
% areas = regionprops(BW,'area');
% %set limit at 20 fo it to be a cell
% objects = areas.Area>20;

figure;
subplot(2,3,1)
imagesc(bluechannel);
imshow(bluechannel, []);
title('original image');

subplot(2,3,2);
imagesc(medianblue);
imshow(medianblue, []);
title('after median filter');

subplot(2,3,3);
imagesc(claheblue);
imshow(claheblue, []);
title('after clahe');

subplot(2,3,4)
imagesc(erodila);
imshow(erodila, []);
title('eroded/dilated');

subplot(2,3,5);
imagesc(binarized);
imshow(binarized, []);
title('binarized');

%% 
gaussian = imgaussfilt(double(binarized),'FilterSize',[3,3]);
figure;
imagesc(gaussian);
imshow(gaussian, []);
title('Gaussian filter')

erodila2 = imopen(gaussian,mask);
figure;
imagesc(erodila2);
imshow(erodila2, []);

%it has a lot of processes, so we need to filter or clean by size...
title('erosion/dilation 2');
