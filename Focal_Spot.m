clear all;
close all;

%%%%%%%%%%------Open and display the RCF scan-----%%%%%%%%%%
[nom,chemin]=uigetfile('*.tif','Fichier scan RCF (16 bits, .tif)');
fic_im=strcat(chemin,nom);
scan=imread(fic_im);
h=figure;
imagesc(scan)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
colormap(gca,'jet')
colorbar
%set(gca,'ColorScale','log')
title(colorbar, 'Dose (Gy)')
axis equal

resolution = inputdlg('Calibration (in µm/px)','Resolution (dpi)',1,{'0.176'});
resolution = str2double(char(resolution));

energy = inputdlg('Mean laser energy (in J)','Laser energy',1,{'10.9'});
energy = str2double(char(energy));

%%%%%%%%%%-----------Background substraction----------%%%%%%%%%%
h=warndlg('Click to delimit a region without signal (background subtraction)');
waitfor(h)
[mask1,x1,y1]=roipoly;
scan=double(scan);
imagebackground = scan.*mask1;
imagebackground(imagebackground==0)= NaN;
background=mean(mean(imagebackground(:,:),'omitnan'),'omitnan');

scan=scan-background;
scan(scan<0)=0;

imagesc(scan)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
colormap(gca,'jet')
colorbar
%set(gca,'ColorScale','log')
title(colorbar, 'Dose (Gy)')
axis equal

%%%%%%%%%%------------Zoom on the focal spot-----------%%%%%%%%%%
Message = warndlg('Click on the top left corner and bottom right corner of the ROI');
waitfor(Message)                                                           %Wait for dialog box to close
[x1,y1]=ginput(2);
x1(1)=round(x1(1));x1(2)=round(x1(2));y1(1)=round(y1(1));y1(2)=round(y1(2));

ROI=scan(y1(1):y1(2),x1(1):x1(2));

% close;
% figure;
% imagesc(ROI)
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% colormap(gca,'jet')
% colorbar
% %set(gca,'ColorScale','log')
% title(colorbar, 'Dose (Gy)')
% axis equal



%%%%%%%%%%------Calculation of the intensity values-----%%%%%%%%%%
Intensity=ROI.*energy/sum(sum(ROI))/24e-15/(resolution/10000*resolution/10000);

min_axis=min(min(Intensity));
max_axis=max(max(Intensity));

imagesc(Intensity);
caxis([min_axis max_axis])
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
colormap(gca,'hot');
colorbar;
%set(gca,'ColorScale','log')
title(colorbar, 'Intensity (W/cm²)')
  
axis on
axis equal