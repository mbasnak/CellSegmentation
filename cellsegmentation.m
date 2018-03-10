% %2D cell segmentation for confocal images
% 
workingDir='C:\Users\melanie\Documents\Doctorado\rotations\Chris\lab meeting\mb11\slice3\';
cd(workingDir)

data=bfopen('C:\Users\melanie\Documents\Doctorado\rotations\Chris\lab meeting\mb11\slice3/maxprojstack2.tif'); %by adding the square braquets, we make it into a single piece of text. If we didn't do that, the comas separating the three pieces of text would make the different texts appear as arguments from the imread function 

series1 = data{1, 1}; %3x2 matrix containing the pixel data for the three channels

farred = mat2gray(series1{3, 1}); %data from plane 1 (which I think coresponds to te blue channel)
green = mat2gray(series1{2, 1}); %data from plane 2 (which I think corresponds to the green channel)
red = mat2gray(series1{1, 1}); %data from plane 3 (which I think corresponds to the red channel)

series1_colorMaps = data{1, 3}; %data 1,3 conains color lookup tables for each plane

figure('Name','Three channels');
subplot(1,4,1);
imshow(farred, []);
title('farred channel','Color','blue')

subplot(1,4,2);
imshow(green, []);
title('green channel','Color','green')

subplot(1,4,3);
imshow(red, []);
title('red channel','Color','red');

channels=zeros(512,512,3);
channels(:,:,1)=red;
channels(:,:,2)=green;
channels(:,:,3)=farred;
subplot(1,4,4);
imshow(channels);
title('the three channels together');

%% If they are in tif format

% clear all;close all;
% 
% workingDir='C:\Users\David\Desktop\segmentation\';
% 
% blue=mat2gray(imread('test_C0.TIF'));
% figure
% subplot(1,4,1);
% imshow(blue);
% title('blue channel');
% 
% green=mat2gray(imread('test_C1.TIF'));
% subplot(1,4,2);
% imshow(green);
% title('green channel');
% 
% red=mat2gray(imread('test_C2.TIF'));
% subplot(1,4,3);
% imshow(red)
% title('red channel');

% channels=zeros(512,512,3);
% channels(:,:,1)=red;
% channels(:,:,2)=green;
% channels(:,:,3)=blue;
% subplot(1,4,4);
% imshow(channels);
% title('the three channels together');

%% 
%intro to this animal
%name of animal? regions injected with each color?

inputTitle = char('Name of animal & region imaged');
prompt = {['Name of animal/region imaged ' ]};
answer = inputdlg(prompt, inputTitle, 1, {' '});
animal=regexp(answer,' ','split');

inputTitle = char('Injection sites');
prompt = {'Red injection ','Green injection ','Blue injection '};
injections = inputdlg(prompt, inputTitle, 1, {'','',''});


%% segmentation for the green channel


%I-Global binarization

%1) median filter
% Prior to the application of any binarization algorithm, the image is
%processed with a median filter to reduce noise while trying to preserve edges

figure('Name','Global binarization');
subplot(2,3,1);
imshow(green); title('original image');

mediangreen = medfilt2(green); %apply the default median filter, with size [3,3]
subplot(2,3,2);
imshow(mediangreen);
title('after median filter')

%2)CLAHE algorithm
% and its contrast is enhanced with the CLAHE algorithm

clahegreen = adapthisteq(mediangreen,'ClipLimit',0.00005,'NBins',6000);
subplot(2,3,3);
imshow(clahegreen); title('after clahe');
%the cliplimit needs to be very low for it to work. I don't understand what
%it's doing.

%3) binarization
% The graythresh function uses Otsu's method,
% which chooses the threshold to minimize the intraclass variance of the black and white pixels
level = graythresh(clahegreen);
binarized = imbinarize(clahegreen,level);

subplot(2,3,4);
imshow(binarized); title('after binarization');


%4) erosion/dilation

se = strel('disk',2); %erosion/dilation circular mask
erodila = imopen(binarized, se);
subplot(2,3,5);
imshow(erodila); title('after ersion/dilation');


%% 
%cell count

%to identify the cells as objects using a logical array...
[labeledImage,N] = bwlabel(erodila);
figure, imshow(labeledImage);

%area of cells
D = regionprops(labeledImage,'Area','Perimeter');
myArea = [D.Area];
figure, hist(myArea) %areas close to 0 can be dots that are not grains or things that are on the
%boundaries

%defining a threshold to get rid of small bobs that aren't cells
thresh = median(myArea) - std(myArea);

%to get rid of the dots...
for iCell = 1:length(myArea)
    if myArea(iCell) < thresh
        labeledImage(labeledImage==iCell)=0;
    end
end
labeledImage = bwlabel(labeledImage);

%to get rid of things that are on the borders...
labeledImage = imclearborder(labeledImage);
figure, imshow(labeledImage);

%% watershed or separating clusters of cells

% bwdist, which computes the distance transform. The distance transform of a binary image is
% the distance from every pixel to the nearest nonzero-valued pixel.

D = -bwdist(~labeledImage);
D(~labeledImage) = -Inf;
L = watershed(D);
figure, subplot(1,2,1);
imshow(green); title('original image');
subplot(1,2,2);
imshow(label2rgb(L,'jet','w')); title('watersheded image');

[cells,N] = bwlabel(L);
celnum = max(cells(:));

%the small blobs ae cells in a different plane. That could be taken into
%account in a 3D cell count.
%% 

%add something to be able to add manually the layers of cortex to the
%images.

%For layer II 

figure,
imshow(green)
layerII = roipoly(green);


figure,
subplot(1,3,1);imshow(green)
greenII = green;
greenII(~layerII) = 0;
subplot(1,3,2);imshow(greenII)
LII = L;
LII(~layerII) = 0;
subplot(1,3,3);imshow(label2rgb(LII,'jet','w'));

[cellsII,N] = bwlabel(LII);
cellnumII = max(cellsII(:))

%% 
%II-local binarization

% attempts to improve the
% accuracy of the binarization process by using the binary image ob-
% tained in the first step to identify individual binary blobs.

% CC = bwconncomp(BW);
% s = regionprops(BW,'centroid');
% 
% %Each of these blobs is then used to compute its bounding box in the origi-
% % nal image.
% boundingbox = regionprops(BW,'BoundingBox');
% cero=zeros(512,512)
% bb=cero(boundingbox.BoundingBox(1))
% imshow(boundingbox)


%2) gaussian filter

% These bounding boxes are then subject to a Gaussian
% low-pass filter,

% gaussian = imgaussfilt(double(binarized),'FilterSize',[3,3]);
% figure;
% imagesc(gaussian);
% imshow(gaussian, []);
% title('Gaussian filter')

%3) Niblack's binarization
% and binarized again but, in this case, with Ni-
% black’s algorithm

%filter k constant -0.2
% niblack=zeros(512,512);
% for i=1:512
%     for j=1:512
% if gaussian(i,j) > abs((mean(mean(gaussian)) - 0.2 * std2(gaussian)))
% niblack(i,j) = 1; 
% else 
% niblack(i,j) = 0; 
% end
%     end
% end
% 
% ones=gaussian(:) > abs((mean(mean(gaussian)) - 0.2 * std2(gaussian)));
% niblackb=zeros(512,512);
% niblackb=double(ones);
% b = reshape(niblackb, 512, 512);
% figure;
% imagesc(b);
% imshow(b, []);
% title('Niblack')


% %niblack with the predetermined function
% nib = niblack(gaussian);
%it asks me for positive inteers or logicals


% %erosion dilation mask
% erodered2 = imerode(b,mask);
% imagesc(erodered2);
% dilated2 = imdilate(erodered2,mask);
% figure=imshow(dilated2, []);
% 
% 
% final result
% figure;
% imshowpair(redchannel,dilated2,'montage')



