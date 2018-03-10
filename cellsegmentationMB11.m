% %2D cell segmentation for confocal images

close all; clear all;

workingDir='C:\Users\melanie\Documents\Doctorado\rotations\Chris\lab meeting\anatomy\mb11\confocal\slice3/'; %define working directory
cd(workingDir) %set working directory

data=bfopen('C:\Users\melanie\Documents\Doctorado\rotations\Chris\lab meeting\anatomy\mb11\confocal\slice3/maxprojstack2.tif'); %open the image 

series1 = data{1, 1}; %3x2 matrix containing the pixel data for the three channels

farred = mat2gray(series1{3, 1}); %data from far red channel
green = mat2gray(series1{2, 1}); %data from green channel
red = mat2gray(series1{1, 1}); %data from red channel

series1_colorMaps = data{1, 3}; %data 1,3 conains color lookup tables for each plane

figure('Name','Three channels'); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,4,1);
imshow(farred, []);
title('Far red channel','Color','blue')

subplot(1,4,2);
imshow(green, []);
title('Green channel','Color','green')

subplot(1,4,3);
imshow(red, []);
title('Red channel','Color','red');

channels=zeros(512,512,3);
channels(:,:,1)=red;
channels(:,:,2)=green;
channels(:,:,3)=farred;
subplot(1,4,4);
imshow(channels);
title('The three channels together');

% figuresdir = 'C:\Users\melanie\Documents\Doctorado\rotations\Chris\lab meeting\anatomy\mb11\confocal\figures\'; %general folder for every figure
% 
% figureName1 = 'DifferentChannels.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName1)); %save the figure


%% segmentation for the farred channel: how many neurons are there?

%I-Global binarization

%1) median filter
% Prior to the application of any binarization algorithm, the image is
%processed with a median filter to reduce noise while trying to preserve edges

figure('Name','Global binarization'); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1);
imshow(farred); title('Original image');
medianfarred = medfilt2(farred); %apply the default median filter, with size [3,3]
subplot(2,3,2);
imshow(medianfarred);
title('After median filter')

%2)CLAHE algorithm
% Contrast is enhanced with the CLAHE algorithm

clahefarred = adapthisteq(medianfarred,'ClipLimit',0.00005,'NBins',6000);
subplot(2,3,3);
imshow(clahefarred); title('After CLAHE');
%the cliplimit needs to be very low for it to work. I don't understand what
%it's doing.

%3) binarization
[frpixelCount, frgrayLevels] = imhist(farred); 
subplot(2,3,4), bar(frpixelCount); %plot histogram of grayscale values
title('Grayscale values distribution');
% The graythresh function uses Otsu's method,
% which chooses the threshold to minimize the intraclass variance of the black and white pixels
level = graythresh(clahefarred);
binarized = imbinarize(clahefarred,level);
subplot(2,3,5);
imshow(binarized); title('After binarization');


%4) erosion/dilation

se = strel('disk',2); %erosion/dilation circular mask
erodila = imopen(binarized, se);
subplot(2,3,6);
imshow(erodila); title('After erosion/dilation');

% figureName2 = 'FarRedSegmentation1.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName2)); %save the figure
%% cell count

%to identify the cells as objects using a logical array...
[labeledImage,N] = bwlabel(erodila);
figure, set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1),imshow(labeledImage); title('Binarized image');

%area of cells
D = regionprops(labeledImage,'Area','Perimeter');
myArea = [D.Area];
subplot(2,2,2), hist(myArea);
title('Areas histogram');
%areas close to 0 can be dots that are not grains or things that are on the
%boundaries

%defining a threshold to get rid of small bobs that aren't cells
thresh = mean(myArea) - std(myArea);

%to get rid of the dots...
for iCell = 1:length(myArea)
    if myArea(iCell) < thresh
        labeledImage(labeledImage==iCell)=0;
    end
end
labeledImage = bwlabel(labeledImage);

%to get rid of things that are on the borders...
labeledImage = imclearborder(labeledImage);
subplot(2,2,3), imshow(labeledImage);
title('Without noise and borders');


%watershed or separating clusters of cells

% bwdist, which computes the distance transform. The distance transform of a binary image is
% the distance from every pixel to the nearest nonzero-valued pixel.
D = -bwdist(~labeledImage);
D(~labeledImage) = -Inf;
L = watershed(D); %compute the watershed
[cells,N] = bwlabel(L);
celnum = N;

subplot(2,2,4);
imshow(label2rgb(L,'jet','w')); title('Watersheded image');
%annotation('textbox', [0.47, 0.8, 0, 0], 'string', {'Number of cells: ';celnum},'FontWeight','bold')

% figureName3 = 'FarRedSegmentation2.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName3)); %save the figure

%% find centroids and add missing cells

cent = regionprops(L,'centroid'); %find the centroids of the segmented neurons
centroids = cat(1, cent.Centroid);
%Display original image and superimpose centroids.
figure,imshow(farred);
title('Centroids of the neurons found');
hold(imgca,'on')
plot(imgca,centroids(:,1), centroids(:,2), 'r*')
hold(imgca,'off')
%manually add centroids of cells missed by the segmentation
[x, y] = getpts;
addedcent = [x,y]; %store the new centroids

%display the original image and the image with both centroids added
figure,
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1),imshow(farred); %show original figure
title('Original image');
subplot(1,2,2),imshow(farred);
title('Image with centroids');
hold(imgca,'on')
plot(imgca,centroids(:,1), centroids(:,2), 'r*') %add the segmentation centroids
plot(imgca,addedcent(:,1), addedcent(:,2), 'b*') %add the manually added centroids
hold(imgca,'off')
totalcell = N + size(x,1);
annotation('textbox', [0.49, 0.8, 0, 0], 'string', {'Total number of neurons: ';totalcell},'FontWeight','bold','FontSize',11);

% figureName4 = 'FarRedCentroids.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName4)); %save the figure

%add something to delete centroids as well
%% segmentation of the other channels

close all;
%green channel

%I-Global binarization
%1) median filter
figure('Name','Global binarization'); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1);
imshow(green); title('Original image');

mediangreen = medfilt2(green); %apply the default median filter, with size [3,3]
subplot(2,3,2);
imshow(mediangreen);
title('After median filter')

%2)CLAHE algorithm
% and its contrast is enhanced with the CLAHE algorithm
clahegreen = adapthisteq(mediangreen,'ClipLimit',0.00005,'NBins',6000);
subplot(2,3,3);
imshow(clahegreen); title('After CLAHE');
%the cliplimit needs to be very low for it to work. I don't understand what
%it's doing.

%3) binarization
[pixelCount, grayLevels] = imhist(green); 
subplot(2,3,4), bar(pixelCount); %plot histogram of grayscale values
title('Grayscale values distribution');
%glevel = graythresh(clahegreen); %this method didn't work so well here
glevel = 0.25;
gbinarized = imbinarize(clahegreen,glevel);
subplot(2,3,5);
imshow(gbinarized); title('After binarization');


%4) erosion/dilation
se = strel('disk',2); %erosion/dilation circular mask
gerodila = imopen(gbinarized, se);
subplot(2,3,6);
imshow(gerodila); title('After erosion/dilation');

% figureName5 = 'GreenSegmentation1.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName5)); %save the figure


%% cell count

%to identify the cells as objects using a logical array...
[glabeledImage,gN] = bwlabel(gerodila);
figure, set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1), imshow(glabeledImage);
title('Binarized image');

%area of cells
gD = regionprops(glabeledImage,'Area','Perimeter');
gmyArea = [gD.Area];

%defining an area threshold to get rid of small bobs that aren't cells
gthresh = mean(gmyArea) - std(gmyArea);
%to get rid of the dots...
for iCell = 1:length(gmyArea)
    if gmyArea(iCell) < gthresh
        glabeledImage(glabeledImage==iCell)=0;
    end
end
subplot(2,2,2), hist(gmyArea)
%line([gthresh gthresh], [0 17]);
title('Areas histogram');

glabeledImage = bwlabel(glabeledImage);

%to get rid of things that are on the borders...
glabeledImage = imclearborder(glabeledImage);
subplot(2,2,3), imshow(glabeledImage);
title('Without noise and borders');

%watershed
gD = -bwdist(~glabeledImage);
gD(~glabeledImage) = -Inf;
gL = watershed(gD);
[gcells,gN] = bwlabel(gL);
gcelnum = gN;
subplot(2,2,4);
imshow(label2rgb(gL,'jet','w')); title('Watersheded image');

% figureName6 = 'GreenSegmentation2.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName6)); %save the figure

%% find centroids and add missing cells

gcent = regionprops(gL,'centroid'); %find the centroids of the segmented neurons
gcentroids = cat(1, gcent.Centroid);
%Display original image and superimpose centroids.
figure,imshow(green);
title('Centroids of the neurons found');
hold(imgca,'on')
plot(imgca,gcentroids(:,1), gcentroids(:,2), 'r*')
hold(imgca,'off')
%manually add centroids of cells missed by the segmentation
[gx, gy] = getpts;
gaddedcent = [gx,gy]; %store the new centroids

%display the original image and the image with both centroids added + the
%far red channel centroids
figure,
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1),imshow(green); %show original figure
title('Original image');
subplot(1,2,2),imshow(green);
title('Image with centroids');
hold(imgca,'on')
plot(imgca,gcentroids(:,1), gcentroids(:,2), 'r*') %add the segmentation centroids
plot(imgca,gaddedcent(:,1), gaddedcent(:,2), 'b*') %add the manually added centroids
allgreencent = [gcentroids;gaddedcent];
hold(imgca,'off')
gcell = gN + size(gx,1);
annotation('textbox', [0.49, 0.8, 0, 0], 'string', {'Number of green labeled neurons: ';gcell},'FontWeight','bold');

% figureName7 = 'GreenCentroids.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName7)); %save the figure

%% segmentation of the red channel

close all;
%I-Global binarization
%1) median filter
figure('Name','Global binarization'); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1);
imshow(red); title('Original image');

medianred = medfilt2(red); %apply the default median filter, with size [3,3]
subplot(2,3,2);
imshow(medianred);
title('After median filter')

%2)CLAHE algorithm
% and its contrast is enhanced with the CLAHE algorithm
clahered = adapthisteq(medianred,'ClipLimit',0.00005,'NBins',6000);
subplot(2,3,3);
imshow(clahered); title('After CLAHE');
%the cliplimit needs to be very low for it to work. I don't understand what
%it's doing.

%3) binarization
[rpixelCount, rgrayLevels] = imhist(red); 
subplot(2,3,4), bar(rpixelCount); %plot histogram of grayscale values
title('Grayscale values distribution');
%rlevel = graythresh(clahered); 
rlevel = 0.5;
rbinarized = imbinarize(clahered,rlevel);
subplot(2,3,5);
imshow(rbinarized); title('After binarization');

%4) erosion/dilation
se = strel('disk',2); %erosion/dilation circular mask
rerodila = imopen(rbinarized, se);
subplot(2,3,6);
imshow(rerodila); title('After erosion/dilation');

% figureName8 = 'RedSegmentation1.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName8)); %save the figure

%% cell count

%to identify the cells as objects using a logical array...
[rlabeledImage,rN] = bwlabel(rerodila);
figure, set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1), imshow(rlabeledImage);
title('Binarized image');

%area of cells
rD = regionprops(rlabeledImage,'Area','Perimeter');
rmyArea = [rD.Area];

%defining an area threshold to get rid of small bobs that aren't cells
rthresh = mean(rmyArea) - std(rmyArea);
%to get rid of the dots...
for iCell = 1:length(rmyArea)
    if rmyArea(iCell) < rthresh
        rlabeledImage(rlabeledImage==iCell)=0;
    end
end
subplot(2,2,2), hist(rmyArea)
%line([rthresh rthresh], [0 12]);
title('Areas histogram');

rlabeledImage = bwlabel(rlabeledImage);

%to get rid of things that are on the borders...
rlabeledImage = imclearborder(rlabeledImage);
subplot(2,2,3), imshow(rlabeledImage);
title('Without noise and borders');

%watershed
rD = -bwdist(~rlabeledImage);
rD(~rlabeledImage) = -Inf;
rL = watershed(rD);
[rcells,rN] = bwlabel(rL);
rcelnum = rN;
subplot(2,2,4);
imshow(label2rgb(rL,'jet','w')); title('Watersheded image');

% figureName9 = 'RedSegmentation2.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName9)); %save the figure

%% %% find centroids and add missing cells

rcent = regionprops(rL,'centroid'); %find the centroids of the segmented neurons
rcentroids = cat(1, rcent.Centroid);
%Display original image and superimpose centroids.
figure,imshow(red);
title('Centroids of the neurons found');
hold(imgca,'on')
plot(imgca,rcentroids(:,1), rcentroids(:,2), 'r*')
hold(imgca,'off')
%manually add centroids of cells missed by the segmentation
[rx, ry] = getpts;
raddedcent = [rx,ry]; %store the new centroids

%display the original image and the image with both centroids added + the
%far red channel centroids
figure,
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1),imshow(red); %show original figure
title('Original image');
subplot(1,2,2),imshow(red);
title('Image with centroids');
hold(imgca,'on')
plot(imgca,rcentroids(:,1), rcentroids(:,2), 'r*') %add the segmentation centroids
plot(imgca,raddedcent(:,1), raddedcent(:,2), 'b*') %add the manually added centroids
allredcent = [rcentroids;raddedcent];
hold(imgca,'off')
rcell = rN + size(rx,1);
annotation('textbox', [0.49, 0.8, 0, 0], 'string', {'Number of green labeled neurons: ';rcell},'FontWeight','bold');

% figureName10 = 'RedCentroids.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName10)); %save the figure

%% see where labeling falls and its properties
close all;

figure,set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
imshow(farred);
hold(imgca,'on')
plot(imgca,allredcent(:,1), allredcent(:,2), 'r*','MarkerSize',8) 
plot(imgca,allgreencent(:,1), allgreencent(:,2), 'g*','MarkerSize',8) 
%plot(imgca,allbluecent(:,1), allbluecent(:,2), 'b*','MarkerSize',8)
hold(imgca,'off')
annotation('textbox', [0.2, 0.8, 0, 0], 'string', {'Proportion of SC projecting neurons: ';(100*gcell/totalcell)},'FontWeight','bold','Color','green');
annotation('textbox', [0.74, 0.8, 0, 0], 'string', {'Proportion of ORBvl projecting neurons: ';(100*rcell/totalcell)},'FontWeight','bold','Color','red');
%annotation('textbox', [0.74, 0.5, 0, 0], 'string', {'Proportion of ORBm projecting neurons: ';(100*bcell/totalcell)},'FontWeight','bold','Color','blue');
title('Proportion of projecting neurons');

 
% figureName11 = 'ProportionOfNeurons.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName11)); %save the figure
%% 
%plot cell distribution
%invert y pixels for the centroids
ygreencent = 512-abs(allgreencent(:,2));
yredcent = 512-abs(allredcent(:,2));
%ybluecent = 512-abs(allbluecent(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%matlab vs fiji segmentation

figure, set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
plot(allredcent(:,1), yredcent, 'ro','MarkerSize',5) 
hold on
plot(allgreencent(:,1), ygreencent, 'g*','MarkerSize',5)
xlim([0 512]);ylim([0 512]);
title('Cell count using MATLAB code')

red = table2array(stats);
green = table2array(statsgreen);

subplot(1,2,2);
plot(red(:,12),(512-abs(red(:,13))),'ro');
hold on
plot(green(:,12),(512-abs(green(:,13))),'g*');
xlim([0 512]); ylim([0 512]);
title('Cell count using 3D object counter in fiji');

figuresdir = 'C:\Users\melanie\Documents\Doctorado\rotations\Chris\lab meeting\anatomy\mb11\confocal\figures\'; %general folder for every figure
figureNamex = 'CellCountComparison.png'; %name the figure
saveas(gcf, strcat(figuresdir,figureNamex)); %save the figure


%overlay
figure, set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
plot(allredcent(:,1), yredcent, 'ro','MarkerSize',5) 
hold on
plot(red(:,12),(512-abs(red(:,13))),'bo');
xlim([0 512]); ylim([0 512]);
title('Overlayed red channel');

subplot(1,2,2)
plot(allgreencent(:,1), ygreencent, 'go','MarkerSize',5) 
hold on
plot(green(:,12),(512-abs(green(:,13))),'bo');
xlim([0 512]); ylim([0 512]);
title('Overlayed green channel');

figureNamez = 'OverlayedCellCounts.png'; %name the figure
saveas(gcf, strcat(figuresdir,figureNamez)); %save the figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%histograms of cell distribution
figure,set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
histogram(ygreencent,'FaceColor','green','BinWidth',25,'Orientation','horizontal')
hold on
histogram(yredcent,'FaceColor','red','BinWidth',25,'Orientation','horizontal')
%histogram(ybluecent,'FaceColor','blue','BinWidth',25,'Orientation','horizontal')
title('Cell distribution'); ylim([0 512]);
xlabel('Number of cells');ylabel('Pixel');

figureName12 = 'CellDistribution.png'; %name the figure
saveas(gcf, strcat(figuresdir,figureName12)); %save the figure

%separated
figure,set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
histogram(ygreencent,'FaceColor','green','BinWidth',25,'Orientation','horizontal')
ylim([0 512]);xlim([0 25]);
title('SC projecting neurons distribution');
subplot(2,1,2)
histogram(yredcent,'FaceColor','red','BinWidth',25,'Orientation','horizontal')
ylim([0 512]);xlim([0 25]);
title('ORBvl projecting neurons distribution');

% figureName13 = 'CellDistribution2.png';
% saveas(gcf, strcat(figuresdir,figureName13)); %save the figure

%% 
%cytoarchitecture: distribution of al neurons
allcent = [centroids;addedcent];
ycent = 512-abs(allcent(:,2));

figure,set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1), imshow(farred); title('Original image');
subplot(1,2,2)
histogram(ycent,'BinWidth',25,'Orientation','horizontal')
title('Distribution of all the neurons');
xlabel('Number of cells');ylabel('Pixel'); ylim([0 525]);

% figureName14 = 'Cytoarchitecture.png'; %name the figure
% saveas(gcf, strcat(figuresdir,figureName14)); %save the figure
