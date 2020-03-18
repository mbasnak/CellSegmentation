clear all;close all;

blue = double(tiffRead('/Users/shihyitseng/Documents/MATLAB/confocal/test_blue6'));
green = double(tiffRead('/Users/shihyitseng/Documents/MATLAB/confocal/test_green6'));
red = double(tiffRead('/Users/shihyitseng/Documents/MATLAB/confocal/test_red6'));


%% Blue channel processing
figure(1), 
subplot(221);
imshow(blue(20:492,20:492)/prctile(blue(:),99));
ln_blue=imgaussfilt(localnormalize(blue,10,10));
subplot(222);
imshow(ln_blue(20:492,20:492)/prctile(ln_blue(:),99));

T_blue = adaptthresh(ln_blue,0.4);
I_blue = imbinarize(ln_blue, T_blue);

%I_blue = imbinarize(ln_blue, 0.3);
I_blue = I_blue(20:492,20:492);
subplot(223);
imshow(I_blue);

se = strel('disk',2);
J_blue = imerode(I_blue,se);
se2 = strel('disk',1);
J2_blue = imdilate(J_blue, se2);
%J2_blue = imreconstruct(J_blue,I_blue);
% subplot(224);
% imshow(J2_blue/prctile(J2_blue(:),99));

%% Indentify blue cells
% bOverlay = zeros(473,473,3);
% bOverlay(:,:,3) = blue(20:492,20:492)/prctile(blue(:),99);
% bOverlay(:,:,1) = J2_blue*0.3/prctile(J2_blue(:),99);
% figure(2); imshow(bOverlay);

conn = 8;
CC_blue = bwconncomp(J2_blue,conn);

S_blue = regionprops(CC_blue, 'Area');

L_blue = labelmatrix(CC_blue);

mean_blue = median(struct2array(S_blue));
std_blue = std(struct2array(S_blue));

thresh = mean_blue - std_blue;

for ii = 1:length(S_blue)
    if S_blue(ii).Area < thresh
        J2_blue(find(L_blue == ii)) = 0;
    end
end


subplot(224);
imshow(J2_blue/prctile(J2_blue(:),99));

bOverlay = zeros(473,473,3);
bOverlay(:,:,3) = blue(20:492,20:492)/prctile(blue(:),99);
bOverlay(:,:,1) = J2_blue*0.5/prctile(J2_blue(:),99);
figure(2); imshow(bOverlay);

conn = 8;
CC_blue = bwconncomp(J2_blue,conn);

S_blue = regionprops(CC_blue, 'Area');

L_blue = labelmatrix(CC_blue);

% figure,
% imshow(label2rgb(L_blue));

%% Red channel processing
figure(3); 
subplot(221);
imshow(red(20:492,20:492)/prctile(red(:),99));
ln_red=imgaussfilt(localnormalize(red,10,100),0.8);
% ln_red=localnormalize(red,10,100);
subplot(222); 
imshow(ln_red(20:492,20:492)/prctile(ln_red(:),99));

I_red = imbinarize(ln_red, 0.4);
I_red = I_red(20:492,20:492);
subplot(223);
imshow(I_red);

se = strel('disk',2);
J_red = imerode(I_red,se);
se2 = strel('disk',3);
J2_red = imdilate(J_red, se2);

% subplot(224);
% imshow(J2_red/prctile(J2_red(:),99));



%% Indentify red cells
% bOverlay2 = zeros(473,473,3);
% bOverlay2(:,:,1) = red(20:492,20:492)/prctile(red(:),99);
% bOverlay2(:,:,2) = J2_red*0.5/prctile(J2_red(:),99);
% figure(4); imshow(bOverlay2);

conn = 8;
CC_red = bwconncomp(J2_red,conn);

S_red = regionprops(CC_red,'Area');

L_red = labelmatrix(CC_red);

thresh2 = 60;
for ii = 1:length(S_red)
    if S_red(ii).Area < thresh2
        J2_red(find(L_red == ii)) = 0;
    end
end

subplot(224);
imshow(J2_red/prctile(J2_red(:),99));


bOverlay2 = zeros(473,473,3);
bOverlay2(:,:,1) = red(20:492,20:492)/prctile(red(:),99);
bOverlay2(:,:,2) = J2_red*0.3/prctile(J2_red(:),99);
figure(4); imshow(bOverlay2);

conn = 8;
CC_red = bwconncomp(J2_red,conn);

S_red = regionprops(CC_red,'Area');

L_red = labelmatrix(CC_red);

%% Green channel processing
figure(5),
subplot(221);
imshow(green(20:492,20:492)/prctile(green(:),99));
ln_green=adapthisteq(localnormalize(green,10,100));
subplot(222);
imshow(ln_green(20:492,20:492)/prctile(ln_green(:),99));

I_green = imbinarize(ln_green, 0.4);
I_green = I_green(20:492,20:492);
subplot(223);
imshow(I_green);

se = strel('disk',3);
J_green = imerode(I_green,se);
se2 = strel('disk',3);
J2_green = imdilate(J_green, se2);
% subplot(224);
% imshow(J2_green/prctile(J2_green(:),99));


%% Indentify green cells
% bOverlay3 = zeros(473,473,3);
% bOverlay3(:,:,2) = green(20:492,20:492)/prctile(green(:),99);
% bOverlay3(:,:,1) = J2_green*0.8/prctile(J2_green(:),99);
% figure(6); imshow(bOverlay3);

conn = 8;
CC_green = bwconncomp(J2_green,conn);

S_green = regionprops(CC_green,'Area');

L_green = labelmatrix(CC_green);


thresh3 = 30;
for ii = 1:length(S_green)
    if S_green(ii).Area < thresh3
        J2_green(find(L_green == ii)) = 0;
    end
end

subplot(224);
imshow(J2_green/prctile(J2_green(:),99));


bOverlay3 = zeros(473,473,3);
bOverlay3(:,:,2) = green(20:492,20:492)/prctile(green(:),99);
bOverlay3(:,:,1) = J2_green*0.8/prctile(J2_green(:),99);
figure(6); imshow(bOverlay3);

conn = 8;
CC_green = bwconncomp(J2_green,conn);

S_green = regionprops(CC_green,'Area');

L_green = labelmatrix(CC_green);

%% Composite image
green2 = double(J2_green);
blue2 = double(J2_blue);
red2 = double(J2_red);

comp2 = cat(3, red2, green2, blue2);
figure(7),
subplot(122);
imshow(comp2);

comp = cat(3,red/prctile(red(:),99),...
        green/prctile(green(:),99),blue/prctile(blue(:),99));

subplot(121);imshow(comp(20:492,20:492,:));

numBlue = sum(struct2array(S_blue))/median(struct2array(S_blue));

fractionRed = length(S_red)/numBlue
fractionGreen = length(S_green)/numBlue