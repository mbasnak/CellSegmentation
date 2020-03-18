%for the 3D segmentation

workingDir='C:\Users\David\Desktop\segmentation\';
cd(workingDir)

data=bfopen('C:\Users\David\Desktop\segmentation/stack.OIB');


for i=1:51

series{i} = data{1, 1}(i+(2*(i-1)):i+(2*(i-1))+2);

bluechannel(i) = series{1,i}(1,1);
greenchannel(i) = series{1,i}(1,2);
redchannel(i) = series{1,i}(1,3);

end
%% 


for i=1:51
seriescolorMaps{i} = data{1, 3}(1,i);

figure;
subplot(1,3,1);
% if (isempty(seriescolorMaps{1,i}(1)))
%   colormap(gray);
% else
%   colormap(seriescolorMaps{1,i}(1,:));
% end
imagesc(bluechannel{1,i});
imshow(bluechannel{1,i}, []);
title('blue channel','Color','blue')
% 
% subplot(1,3,2);
% if (isempty(series1_colorMaps{2}))
%   colormap(gray);
% else
%   colormap(series1_colorMaps{2}(1,:));
% end
% imagesc(greenchannel);
% imshow(greenchannel, []);
% title('green channel','Color','green')
% 
% subplot(1,3,3);
% if (isempty(series1_colorMaps{3}))
%   colormap(gray);
% else
%   colormap(series1_colorMaps{3}(1,:));
% end
% imagesc(redchannel);
% imshow(redchannel, []);
% title('red channel','Color','red');

end


[imageStack]=imreadBF('C:\Users\David\Desktop\segmentation/stack.OIB',[1:51],:,RFP)