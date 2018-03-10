% plots of the statistics obtained with imageJ

red = table2array(stats);
green = table2array(statsgreen);

figure,
plot(red(:,12),(512-abs(red(:,13))),'ro');
hold on
plot(green(:,12),(512-abs(green(:,13))),'g*');
xlim([0 512]); ylim([0 512]);
title('Cell count using 3D object counter in fiji');

figuresdir = 'C:\Users\melanie\Documents\Doctorado\rotations\Chris\lab meeting\anatomy\mb11\confocal\figures\'; %general folder for every figure
figureName1 = 'FijiCellCount.png'; %name the figure
saveas(gcf, strcat(figuresdir,figureName1)); %save the figure


figure,
scatter3(red(:,12),red(:,13),red(:,14),'ro');
hold on
scatter3(green(:,12),green(:,13),green(:,14),'g*');
xlim([0 512]); ylim([0 512]); zlim([0 50])
xticks([0 512]); yticks([0 512]); zticks([0 50]); 
% xticklabels({'x = Lateral','x = Medial'});
% yticklabels({'y = Dorsal','y = Ventral'});
% zticklabels({'z = Posterior','z = Anterior'});



%for making a "histogram" that shows according to the size of the dots an
%estimate of how many green labelled and red labelled neurons there are as a
%function of y and z...


% figure,
% plot(red(:,3),red(:,2),'ro');
% hold on
% plot(green(:,3),green(:,2),'g*');
% xlim([0 51]); ylim([0 512]);

%maybe using hist3






