% %2D cell segmentation for confocal images
% 
workingDir='C:\Users\David\Desktop\segmentation\';
cd(workingDir)

data=bfopen('C:\Users\David\Desktop\segmentation/Image2_z1_slide2_slice9.OIB'); %by adding the square braquets, we make it into a single piece of text. If we didn't do that, the comas separating the three pieces of text would make the different texts appear as arguments from the imread function 

series1 = data{1, 1}; %3x2 matrix containing the pixel data for the three channels

blue = mat2gray(series1{1, 1}); %data from plane 1 (which I think coresponds to te blue channel)
green = mat2gray(series1{2, 1}); %data from plane 2 (which I think corresponds to the green channel)
red = mat2gray(series1{3, 1}); %data from plane 3 (which I think corresponds to the red channel)


if size(series1,1)==4
   farred = mat2gray(series1{4, 1});
end

figure('Name','Three channels');
subplot(1,4,1);
imshow(blue, []);
title('blue channel','Color','blue')

subplot(1,4,2);
imshow(green, []);
title('green channel','Color','green')

subplot(1,4,3);
imshow(red, []);
title('red channel','Color','red');

channels=zeros(512,512,3);
channels(:,:,1)=red;
channels(:,:,2)=green;
channels(:,:,3)=blue;
subplot(1,4,4);
imshow(channels);
title('the three channels together');

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

%% 

%I-Global binarization

%1) median filter
% Prior to the application of any binarization algorithm, the image is
%processed with a median filter to reduce noise while trying to preserve edges

names={'red','green','blue'}

for i=1:3

figure('Name',['Global binarization ',names{i}]);
subplot(2,3,1);
imshow(channels(:,:,i)); title(['original image ',names{i}]);

medianf{i} = medfilt2(channels(:,:,i)); %apply the default median filter, with size [3,3]
subplot(2,3,2);
imshow(medianf{i});
title(['after median filter ',names{i}])

%2)CLAHE algorithm
% and its contrast is enhanced with the CLAHE algorithm

clahef{i} = adapthisteq(medianf{i},'ClipLimit',0.00005,'NBins',6000);
subplot(2,3,3);
imshow(clahef{i}); title(['after clahe ',names{i}]);
%the cliplimit needs to be very low for it to work. I don't understand what
%it's doing.

%3) binarization
% The graythresh function uses Otsu's method,
% which chooses the threshold to minimize the intraclass variance of the black and white pixels
level{i} = graythresh(clahef{i});
binarized{i} = imbinarize(clahef{i},level{i});
subplot(2,3,4);
imshow(binarized{i}); title(['after binarization ',names{i}]);


%4) erosion/dilation

se = strel('disk',2); %erosion/dilation circular mask
erodila{i} = imopen(binarized{i}, se);
subplot(2,3,5);
imshow(erodila{i}); title(['after ersion/dilation ',names{i}]);


end

%this is working very bad for images with a lot of processes.

%% watershed or separating clusters of cells

% bwdist, which computes the distance transform. The distance transform of a binary image is
% the distance from every pixel to the nearest nonzero-valued pixel.

for i=1:3
D{i} = -bwdist(~erodila{i});
D{i}(~erodila{i}) = -Inf;
L{i} = watershed(D{i});
figure, subplot(1,2,1);
imshow(green); title(['original image ',names{i}]);
subplot(1,2,2);
imshow(label2rgb(L{i},'jet','w')); title(['watersheded image ',names{i}]);

% [cells,N] = bwlabel(L);
% celnum = max(cells(:));
end

%the small blobs ae cells in a different plane. That could be taken into
%account in a 3D cell count.
%% 



%cell count

%to identify the cells as objects using a logical array...

for i=1:3
[labeledImage{i},N{i}] = bwlabel(erodila{i});
figure,
subplot(1,3,1);
imshow(labeledImage{i});  title(['objects ',names{i}]);

%area of cells
D{i} = regionprops(labeledImage{i},'Area','Perimeter');
myArea{i} = [D{i}.Area];
subplot(1,3,2);
hist(myArea{i}); title(['object areas ',names{i}]); %areas close to 0 can be dots that are not grains or things that are on the
%boundaries

%defining a threshold to get rid of small bobs that aren't cells
thresh{i} = median(myArea{i}) - std(myArea{i});

% %to get rid of the dots...
%     for iCell = 1:length(myArea)
%         if myArea(iCell) < thresh
%          labeledImage(labeledImage==iCell)=0;
%         end
%     end
% labeledImage = bwlabel(labeledImage);
% 
% %to get rid of things that are on the borders...
% labeledImage = imclearborder(labeledImage);
% figure, imshow(labeledImage);

end


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