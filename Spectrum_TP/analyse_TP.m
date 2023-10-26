clear all
close all

% data_place='C:\Users\lelievre.r\Documents\Manips\2023-04 APOLLON\Qualification 4PW\TP';
% dircur=pwd
% cd(data_place);
% [nom,chemin]=uigetfile('*.tif','');
% nom_seul=nom(1:max(size(nom))-4);
% 
% data_placenoise='C:\Users\lelievre.r\Documents\Manips\2023-04 APOLLON\Qualification 4PW\TP';
% dircurnoise=pwd
% cd(data_place);
% [nomnoise,cheminnoise]=uigetfile('*.tif','');
% radeyenoise = imread(strcat(cheminnoise,nomnoise));

h2=1;
taille_pixel=48*1e-6;
%contrast = 0.1;
contrast = 2;

% radeyeraw = imread(strcat(chemin,nom));
% radeyenoise = imread(strcat(cheminnoise,nomnoise));
% 
% radeye = radeyeraw-radeyenoise;
%radeye = radeyeraw;

%figure(h1);imagesc(radeye)
%caxis([0 contrast]);
%colorbar;

% radeye1 = radeye(1:1000,1025:1536);%Découpe de la partie 1
%figure(h2+1);imagesc(radeye1)
%caxis([0 contrast]);
%colorbar;

% radeye2 = radeye(1:1000,1537:end);%Découpe de la partie 2
% radeye2 = imrotate(radeye2, 180,'bilinear','loose');
%figure(h2+2);imagesc(radeye2)
%caxis([0 contrast]);
%colorbar;

% newradeye = zeros(2016,512);%Assemblage des parties 1 et 2      %2016 = 2000 + 16 pixel de gap entre les 2 CMOS
% newradeye(1:1000,:) = radeye1;
% newradeye(1017:end,:) = radeye2;            %Enlever la partie 2 pour ne pas ajouter
%1001 + gap de 16 pixel entre les 2 CMOS    %les valeurs des pixels du radeye au valeur de pixels de l'IP (pour éviter
                                            %une surestimation du nb de particules)

%newradeye = newradeye/10000;               %Diminution des valeurs sur l'image Radeye car initialement bien plus élevées que celles de l'IP

%figure(h2+3);imagesc(newradeye)
%caxis([0 contrast]);
%colorbar;

%Ouverture du scan IP
data_place='C:\Users\lelievre.r\Documents\Manips\2023-04 APOLLON\Qualification 4PW\TP';
dircur=pwd
cd(data_place);
[nom,chemin]=uigetfile('*.inf','Fichier Thomson Parabola IP(16 bits)'); % Dialog box that allows user to specify an .inf file
fic_im=strcat(chemin,nom);
nom_seul=nom(1:max(size(nom))-4);
file_content=textread(fic_im,'%s','delimiter','\r');
cd(dircur);
taille_pixel=str2num(file_content{4})*1e-6;
bits=str2num(file_content{6});
precision=strcat('uint',file_content{6});
nl=str2num(file_content{7});
nc=str2num(file_content{8});
sensibilite=str2num(file_content{9});
latitude=str2num(file_content{10});

image1 = multibandread(strcat(chemin,nom_seul,'.img'),[nc, nl, 1],precision, 0, 'bsq', 'ieee-be');

resolution=taille_pixel*1e6;
image1=(resolution/100)^2*4000/sensibilite*10.^(latitude*(image1/(2^bits)-1/2));
%image1 = imrotate(image1, 180,'bilinear','loose');
%image1 = flipdim(image1, 2);
figure(h2+4);imagesc(image1)
caxis([0 contrast]);
colorbar;

%IP1 = IP(206:1205,370:869); 
%%IP1 = imresize(IP1, 1);
%figure(h2+4);imagesc(IP1)
%caxis([0 contrast]);
%colorbar;

%L'image recomposée du radeye dans une matrice plus grande
% b=zeros(3000,1500);
% newradeye = imresize(newradeye, 0.96);                     %Redimensionnement de l'image du Radeye pour s'adapter à la résolution de l'IP (50µm/px contre 48µm/px pour le Radeye)
% b(501:2436,501:992) = newradeye;                           %2436 = 2420 + 16 pixel de gap entre les 2 CMOS
% figure(h2+5);imagesc(b)
% caxis([0 contrast]);
% colorbar;

%Fit de l'image rognée de l'IP sur l'image du radeye
%b(1106:2105,494:993)= b(1106:2105,494:993)+ IP1;            %Ajout du scan de l'IP sur l'image Radeye

%Rognage des bords inutiles
% b = b(480:2450,480:1000);
% figure(h2+6);imagesc(b)
% caxis([0 contrast]);
% colorbar;

%image1=b;
h1=1;
% Rotation de l'image si nécessaire
% NB: the "normal" orientation is the one where the spectrum is horizontal and the slit is on the left of the spectrum
buttonr = questdlg('Besoin de tourner l''image ?',...
    'none','Oui','Non','Oui');
if strcmp(buttonr,'Oui')
    button='Non';
    nouvelle=image1;
    nouvelle=imrotate(nouvelle, 180,'bilinear','loose');
    figure(h1);imagesc(nouvelle)
    caxis([0 contrast])
    colorbar
    
    button = questdlg('Trace horizontale ?',...
            'none','Oui','Non','Oui');
    while strcmp(button,'Non')
        figure(h1);imagesc(nouvelle);%pixval on;colormap(red);
        caxis([0 contrast])
        h=warndlg('Click at two points along the IP to rotate it straight');
        waitfor(h)
        %rotation globale pour remettre l'image droite
        [xr,yr]=ginput(2);%x=H ; y=V
        teta=atan((yr(2)-yr(1))/(xr(2)-xr(1)))
        nouvelle=imrotate(nouvelle,teta*180/pi,'bilinear','crop');
        figure(h1);imagesc(nouvelle)
        caxis([0 contrast])
        colorbar
        button = questdlg('Content ?',...
            'none','Oui','Non','Oui');
    end
    image1=nouvelle;
    
%     h=warndlg('Click at the separation point between the part without filter and the part with filter');
%     waitfor(h)
%     [xsp,ysp]=ginput(1);
%     h=warndlg('Click to delimit the width of the nominal proton trace');
%     waitfor(h)
%     [xnw,ynw]=ginput(2);
%     nominalwidth = abs(ynw(1)-ynw(2))
%     h=warndlg('Click to delimit the width of the scattered proton trace');
%     waitfor(h)
%     [xsw,ysw]=ginput(2);
%     scatteredwidth = abs(ysw(1)-ysw(2))
%     scatterratio = scatteredwidth/nominalwidth
    
    clear nouvelle
end

    
%Prise de la coupe
h=warndlg('Click to delimit the ROI');
waitfor(h)
[mask0,x0,y0]=roipoly;
image1bis = image1.*mask0;
image1bis(image1bis==0)= NaN;

%Prendre en compte que les 3 plus grandes valeurs de chaque colonne
%[M1,I1]=max(image1bis);                     %où M la valeur et I la position du maximum dans l'image1bis (ROI)
%K1=image1bis==M1;
%imtemp=image1bis;
%imtemp(imtemp==M1)=NaN;
%[M2,I2]=max(imtemp);
%K2=imtemp==M2;
%imtemp(imtemp==M2)=NaN;
%[M3,I3]=max(imtemp);
%K3=imtemp==M3;
%imtemp(imtemp==M3)=NaN;
%K=K1+K2+K3;
%image1bis = image1bis.*K;                   %K étant un masque avec les positions des 3 plus grandes valeurs pour chaque colonne
%image1bis(image1bis==0)= NaN;

% coupe2=mean(image1bis(:,:),'omitnan');
 
% coupe2(xsp:end)=coupe2(xsp:end).*scatterratio;

coupe2=sum(image1bis(:,:),'omitnan')

coupewidth =size(coupe2,2);

   for i = 1:coupewidth
     PP=image1bis(:,i);
     PP2=isnan(PP);
     PP(PP2)=[];
     PP3(i)=size(PP,1)*taille_pixel*1e3;    %to convert in mm  
     largeur_trace=PP3';                    %transposition de ligne à colonne
   end
   
figure(h1+1);
plot(coupe2,'r');
hold on

%h=warndlg('Click at two points along the IP to where you want to take the lineout [which will be averaged, like in ImageJ]');
%waitfor(h)
%figure(h1);
%[xr,yr]=(ginput(2));%x=H(col) ; y=V(lig)
%xr=round(xr)
%yr=round(yr)
%coupe2=mean(image1(yr(1):yr(2),1:xr(2)));
%figure(h1+1);
%plot(coupe2,'r');
%hold on



%REMOVAL OF SATURATION WITH SECOND SCAN
%SEE EXCEL FILE
%spectre(:,3)=2.7675.*(spectre(:,3));
coupe=coupe2.*1;
%coupe=coupe2.*5.95;                     %scaling from 2nd scan to 1st one of shot#61 of 2022 RAL exp
%coupe=coupe2.*6.8494;                   %scaling from 2nd scan to 1st one of shot#65 of 2022 RAL exp
%coupe=coupe2.*4.0624;                   %scaling from 2nd scan to 1st one of shot#71 of 2022 RAL exp
%coupe=coupe2.*6.6;                      %scaling from 2nd scan to 1st one of shot#79 of 2022 RAL exp
%coupe=coupe2.*5.8;                      %scaling from 2nd scan to 1st one of shot#81 of 2022 RAL exp
%coupe=coupe2.*5.7;                      %scaling from 2nd scan to 1st one of shot#83 of 2022 RAL exp
%coupe=coupe2.*6;                        %scaling from 2nd scan to 1st one of shot#90 of 2022 RAL exp

%recherche de la fente et calcul largeur fente
h=warndlg('Click at the location of the slit (around)');
waitfor(h)
[xf,yf]=(ginput(1));%x=H(col) ; y=V(lig)
xf=round(xf);
demarre=xf-50;
[valeur,indice]=max(coupe(demarre:xf+100)); %bien choisir l'indice de fin
[valeur2,indice2]=find(coupe(demarre:xf+100)>=valeur/2);
bord1=min(indice2)+demarre;
bord2=max(indice2)+demarre;
milieu=round((bord1+bord2)/2);
%recherche du bruit ? droite
h=warndlg('Click to select the zone of the noise to remove on the right of the slit (around)');
waitfor(h)
[xdr,ydr]=(ginput(1));%x=H(col) ; y=V(lig)
xdr=round(xdr);
fond1=(coupe(xdr));
plot(milieu:max(size(coupe)),fond1,'green')
%recherche du bruit ? gauche
h=warndlg('Click to select the zone of the noise to remove on the left of the slit (around)');
waitfor(h)
[xg,yg]=(ginput(1));%x=H(col) ; y=V(lig)
xg=round(xg);
fond2=coupe(xg);
plot(1:milieu,fond2,'magenta')
plot(1:milieu,fond2,'magenta')
fond=(fond1+fond2)/2;
sans_bruit=coupe-fond;
[valeur2,indice2]=find(sans_bruit(demarre:xf+100)>=(valeur-fond)/2);
bord1=min(indice2)+demarre;
bord2=max(indice2)+demarre;
%centre=round((bord1+bord2)/2);
figure(h1+2);
plot(sans_bruit(bord1:bord2))
%largeur_fente = (bord2-bord1)*taille_pixel%largeur fente en m
centre=round((bord1+bord2)/2);

%Pinhole de 120 µm utilisée lors de la campagne de tirs APOLLON - F1
largeur_fente = 1.2e-4    %largeur fente en m en FWHM

%centre=333;%pixel du centre de la fente


%retirer le bruit
taille=size(coupe);
x=(1:taille(2))-centre;
button='Non';
while strcmp(button,'Non')
    figure(h1)
    %h=warndlg('Click at height on the IP to where you want to take the lineout of the noise to remove [which will be averaged, like in ImageJ]');
    %waitfor(h)
    %[xb,yb]=(ginput(1));%x=H(col) ; y=V(lig)
    %xb=round(xb);
    %yb=round(yb);
    %coupe_bruit=mean(image1(yb-20:yb+20,:));
    
    h=warndlg('Click to delimit the ROI (noise)');
    waitfor(h)
    [mask1,x1,y1]=roipoly;
    image1bis2 = image1.*mask1;
    image1bis2(image1bis2==0)= NaN;
    coupe_bruit=mean(image1bis2(:,:),'omitnan');
    
    clf(h1+1);
    figure(h1+1);
    hold on
    plot(x,coupe)
    plot(x,coupe_bruit,'red');
    plot(x,coupe-coupe_bruit,'green');
    legend('Signal','Noise','S-N','Location','BEST');        
    button = questdlg('Content ?',...
        'none','Oui','Non','Oui');
end

%Calibration en énergie, doit être organisé comme suit:
%1er col, distance depuis l'ordre 0, en cm, 2ème col, E_p en MeV
spectre(:,1)=zeros(max(taille),1);

spectre(centre+10:end,1)=x(1,centre+10:end).*taille_pixel*1e2; %en cm

%Calibration obtenue à partir du calcul de la déviation des protons en
%fonction de leur énergie, puis fit de la courbe obtenue pour obtenir la
%relation E=f(x):
spectre(centre+10:end,2)=(((1.602e-19*1.045*0.07*(0.252+0.07/2))^2)./((2*(spectre(centre+10:end,1)*0.01).^2)*1.67e-27*1.6e-19))/1e6;%E in MeV

ts=size(spectre);

%3rd col of "spectre" is the spectrum in PSL
coupe3=coupe-coupe_bruit;
%here we put to zero eveything that is at higher energy (lower px value) than the E_cutoff of the protons
%coupe3(1:500)=0;
spectre(:,3)=coupe3;%max(coupe-coupe_bruit,0);

% %compensation de l'effet de chromatisme magn?tique
% Ri=xi-centre;
% R=-1.*(x-Ri);
% R=(R/Ri)';
% spectre(:,3)=max(coupe-coupe_bruit,0);
% spectre(centre:end,3)=spectre(centre:end,3).*R(centre:end);


%Calcul du FIT: PSL / nbr protons
%using the equations of Mancic RSI 2008

calcul='oui';% calcul des FIT et angle solide
if calcul=='oui'
    deb_E=min(find(spectre(:,2)));
    Energie_p=spectre(deb_E:end,2);
    [valeur,indice]=min(abs(Energie_p-2.11));
    %(E < 2.11 MeV)
    FIT1_PSL_2Mev = 0.22039*exp(-(Energie_p(indice:end)-1.5049).^2 /(1.1842^2 ));
    % 2.11 MeV < E < 20 MeV
    FIT2_PSL_20Mev = 0.33357 * Energie_p(1:indice-1).^-0.91377;
    FIT_PSL(1:indice-1,1)=FIT2_PSL_20Mev;
    FIT_PSL(indice:max(size(Energie_p)),1)= FIT1_PSL_2Mev;
    save 'FIT_PSL' FIT_PSL
else
    load FIT_PSL
end
fin=max(max(size(spectre)));
spectre(deb_E:fin,4) =  spectre(deb_E:fin,3)./...
    FIT_PSL(:);

%spectre(deb_E:fin,4) =
%spectre(deb_E:fin,3).*(spectre(deb_E:fin,2).^1.5)/1150;    %Calibration
%Radeye

%Calcul de l'angle solide
%adapt to the particular distances of the exp!!!
%APOLLON Novembre 2022: distance TCC-front of spectro#2=495 mm
%dist TCC-Radeye=495 + 381.55 mm
%dist front-magnet=58 mm
%magnet length=50 mm
%dist magnet end-back(=IP)=32 mm
%since total length of TP=140 mm
if calcul=='oui'
    m_proton=1.672649e-27;%(kg)
    e_electron=1.602e-19;%(C)
    B_T=1.045;%B_field of the magnets (see Alice's PDF, eq 7)
    W_C=e_electron*B_T/m_proton;%qB/m (rad/s)
    L1=0.07;%(m)%distance dans le champ B
    L2=0.252;%(m)%distance depuis la fin du champ B jusqu'a l'IP
    
    Vitesse_p = sqrt((2e6*e_electron/m_proton).*Energie_p);
    Theta_rad = asin((L1*W_C)./Vitesse_p);
    L_champ = (1/W_C).*Vitesse_p.*Theta_rad+L2./cos(Theta_rad);
    
    % answer = inputdlg('entrez la distance au CC','en metre',1,{'1.255'});
    Longueure_tot=(427.5+381.55)*1e-3;%str2num(char(answer));%dist TCC-IP
    
    CC_ChampB = Longueure_tot-(L1+L2);%(m)
    L_tot = L_champ + CC_ChampB;
    
    largeur_fente=largeur_fente*1e3;%to convert in mm
    largeur_px=taille_pixel*1e3;%to convert in mm
    largeur_trace = largeur_trace(deb_E:end);
    
%     Omega_sr = 4*asin(largeur_fente/2./sqrt((largeur_fente/2)^2+...
%         (L_tot.*1000).^2)*largeur_px./2./sqrt((largeur_px/2)^2+...
%         (L_tot.*1000).^2));
    
%     Omega_sr = 4*asin(largeur_fente/2./sqrt((largeur_fente/2)^2+...
%         (L_tot.*1000).^2).*largeur_trace./2./sqrt((largeur_trace./2).^2+...
%         (L_tot.*1000).^2));
    
    Omega_sr = 6.188425e-8; %pinhole de 120µm de diamètre à 42,75cm du TCC
    
    save 'Omega_sr' Omega_sr
else
    load Omega_sr
end

%calcul du dE
%D=length from beginning of magnet to IP[m]
%r=length of magnet[m]
%so that (D-r) is the length (along the O-order axis) of the free trajectory (out of the magnets)
alpha=2*m_proton/(e_electron*B_T)^2;
r=L1;
D=L1+L2;
%below is the derivative of the calibration which gives
%the position on the IP vs E
%see "TP_calibrationProgram"
%r_g=sqrt(2.*E[J]*m_p[kg])/(Z*e[C]*Btheo[T]);
%l is the distance along the IP, from the Oth-order, at which protons hit it, depending on their energy
%l[m]=abs(sqrt(r_g.^2-r^2)-r_g-r*(D-r)./sqrt(r_g.^2-r^2));
%the detailed derivation can be found in Maxence's thesis

dldE=abs((alpha./(2.*sqrt(alpha.*(Energie_p.*1e6.*1.6e-19)-r^2))...
    -alpha./(2.*sqrt(alpha.*(Energie_p.*1e6.*1.6e-19)))+...
    alpha.*r.*(D-r)./(2.*(alpha.*(Energie_p.*1e6.*1.6e-19)-r^2).^1.5)));
%in m/J
%le pas en énergie au niveau de l'IP est le 1/(dl/dE)*1 px
%donc il faut multiplier par 1 px = taille_pixel (en m)
%et on veut le résultat en MeV, donc il faut le multiplier par 1/(1e6.*1.6e-19) 
facteur=taille_pixel/(1e6.*1.6e-19);
DE_energie=1./dldE.*facteur;



%Calcul du nbr_p/(str*Mev)

nbr_p_str_Mev = zeros(max(max(size(spectre))),1);
nbr_p_str_Mev  = spectre(deb_E:fin,4) ./ ...
    (Omega_sr.*DE_energie);

spectre(deb_E:fin,5) =  nbr_p_str_Mev  ;

%now to get in part/MeV, we multiply by the solid angle given by Mancic
E_max=37%MeV
E_norm=spectre(1:fin,2)./E_max;
demi_angle=(7.2327+114.54.*E_norm -321.16.*E_norm.^2+364.57.*E_norm.^3-165.08.*E_norm.^4);
angle_solide=2.*pi.*(1-cos(demi_angle./180.*pi));%in sr
spectre(1:fin,6) = spectre(1:fin,5).*angle_solide;

[ve,ie]=find(spectre(:,2)>0);

%calcul de l'énergie totale (en J) contenue dans les protons
%la conversion en =énergie laser -> protons doit Ãªtre autour du %
energie_p_totale=sum(abs(spectre(ve(1):end,6)).*abs(DE_energie).*Energie_p)*1e6*1.6e-19;
efficacite=energie_p_totale/60*100;%en % Ã  partir de 60 J (LULI2000)

%enregistrement
save(strcat(nom_seul,'_spectre.txt'),'spectre','-ascii','-tabs');
%1st col= cm along the IP
%2nd col = E(MeV)
%3rd col= PSL
%4th col= part
%5th col= part/MeV/sr
%6th col=part/MeV
clf(h1+1);
figure(h1+1);

%loglog(spectre(ve(1):end,2),spectre(ve(1):end,6))
%plot(spectre(ve(1):end,2),spectre(ve(1):end,4))
plot(spectre(ve(1):end,2),spectre(ve(1):end,6))
%axis([0 25  0.1 5e2])
axis([0 50  1e6 1e13])
title('Proton spectrum - Shot-60')
xlabel('Energy (MeV)')
%ylabel('Nb protons')
%ylabel('dN/(dE.dOmega) (part/(MeV.sr))')
ylabel('dN/dE (part/MeV)')
set(gca, 'XScale', 'linear')
set(gca, 'YScale', 'log')

cd(pwd);
return

%analyse
figure(h1)
imagesc(image1)
caxis([0 8])
colorbar

droit='oui';
if droit=='oui',
    h=warndlg('Click where you want the angular lineout (2 points in vertical)');
    waitfor(h)
    [xp,yp]=ginput(2);%x=H ; y=V
    coupe_angle(:,1)=image1(:,round(xp(1)))-0.5;%NOISE TO BE ADJUSTED !!!
    angle=((1:max(size(coupe_angle)))-1750).*taille_pixel;%1500 IS WHERE THE PEAK IS IN PX
    %TO BE ADJUSTED !!!
    angle=1e3*atan(angle./0.36);%en mrad
    figure
    plot(angle,coupe_angle./max(coupe_angle))
    coupe_angle(:,2)=coupe_angle(:,1);
    coupe_angle(:,1)=angle;
else
    h=warndlg('Click at the apex and at one edge of the desired angular lineout');
    waitfor(h)
    %rotation globale pour remettre l'image droite
    [xp,yp]=ginput(2);%x=H ; y=V
    a=abs((xp(1)-xp(2))/(yp(1)-yp(2))^2);
    p=1/(4*a);
    tmax=(yp(2)-yp(1))/(2*p);
    tabs=abs(tmax);
    t=-1*tabs:tabs/1e3:tabs;
    xangle=round(xp(1)+p.*t.^2);
    yangle=round(yp(1)+2.*p.*t);
    hold on
    plot(xangle,yangle,'c.')
    %largeur de la fente de la TP EAST (@0 deg)=0.0433 m
    angle_max=1e3*atan(0.0433/0.36);%en mrad
    coupe_angle(:,1)=t./tabs.*angle_max;
    for i=1:max(size(xangle))
        coupe_angle(i,2)=(image1(yangle(i),xangle(i)))-1;%0.5 is the noise (PSL) on the IP close to the signal
    end
    figure
    plot(coupe_angle(:,1),coupe_angle(:,2)./max(coupe_angle(:,2)))
end

xlabel('Angle (mrad)')
ylabel('Number of protons (A.U.)')


save(strcat(nom_seul,'_angle.txt'),'coupe_angle','-ascii','-tabs');

cd(pwd);
