% VERSION september 2014 by Anna-Marie Schroer


%########## calculation of the deposited energy for each film####
% used in main line 13


figure;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

colorbar
imagesc(a);
axis equal;

%recentrage de la zone à étudier


clear y1,clear x1
%n_filmb=num2str(n_film);
H = warndlg('Click on the top left corner and bottom right corner of your studying zone',strcat('RCF #',n_filmb));
waitfor(H) %attend la fermeture de la boite
[y1,x1]=ginput(2);%y1=row=y x1=column=x !!
x1(1)=round(x1(1));x1(2)=round(x1(2));y1(1)=round(y1(1));y1(2)=round(y1(2));

if x1(1)<1
    x1(1)=1;
end
if x1(2)>size(a,1)
    x1(2)=size(a,1);
end
if y1(1)<1
    y1(1)=1;
end
if y1(2)>size(a,2)
    y1(2)=size(a,2);
end

c=a(x1(1):x1(2),y1(1):y1(2));

%soustraction du bruit
%noise subtraction
close;
h=figure;
set(h,'Position',[100 100 1000 1000]);
imagesc(a);
axis equal

close;
figure;
imagesc(c)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
colorbar
title(colorbar, 'Transmission')
axis on
axis equal
hold on


% Test start -------------
% choose two points from the RCF layer
% rotate and move each RCF layer into the same original position

clear y1,clear x1
H = warndlg('Click on the top left corner and top right corner of your studying zone for rotation',strcat('RCF #',n_filmb));
waitfor(H) %attend la fermeture de la boite
[y1,x1]=ginput(2);%y1=row=y x1=column=x !!
x1(1)=round(x1(1));x1(2)=round(x1(2));y1(1)=round(y1(1));y1(2)=round(y1(2));

fprintf('\n the coordinate of the top left corner is: \n');
fprintf('x1(1) = %1.0f \n', x1(1));
fprintf('y1(1) = %1.0f \n', y1(1));

fprintf('\n the coordinate of the top right corner is: \n');
fprintf('x1(2) = %1.0f \n', x1(2));
fprintf('y1(2) = %1.0f \n', y1(2));

% if x1(1)<1
%     x1(1)=1;
% end
% if x1(2)>size(c,1)
%     x1(2)=size(c,1);
% end
% if y1(1)<1
%     y1(1)=1;
% end
% if y1(2)>size(c,2)
%     y1(2)=size(c,2);
% end



% rotate this image so that it is horizontal
% ang_rot = abs(atand(y1(2) - y1(1)) / (x1(2) - x1(1))); % in degrees
ang_rot = atand(abs(x1(2) - x1(1)) / abs(y1(2) - y1(1))); % in degrees
fprintf('ang_rot = %1.3f degree \n', ang_rot);
if x1(1) > x1(2)
    c = imrotate(c,-ang_rot,'nearest','crop');  % rotate the image by clockwise ang_rot degrees
else
    c = imrotate(c,ang_rot,'nearest','crop');  % rotate the image by anti-clockwise ang_rot degrees
end



close;
figure;
imagesc(c)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
colorbar
title(colorbar, 'Transmission')
axis on
axis equal
hold on


% now we use circshift to move the top left cornor to (1,1)
clear y1,clear x1
%n_filmb=num2str(n_film);
H = warndlg('Click on the top left corner of your studying zone to move coordinates',strcat('RCF #',n_filmb));
waitfor(H) %attend la fermeture de la boite
[y1,x1]=ginput(1);%y1=row=y x1=column=x !!
x1(1)=round(x1(1));
y1(1)=round(y1(1));

if x1(1)<1
    x1(1)=1;
end
if y1(1)<1
    y1(1)=1;
end

% move the left conor coordinate of this RCF layer to the original position (1,1)
delx = x1(1) - 1;
dely = y1(1) - 1;
fprintf('delx = %1.1f \n', delx);
fprintf('dely = %1.1f \n', dely);

if delx > 0
    c = circshift(c, -delx, 1); % attention to the dimension
end
if dely > 0
    c = circshift(c, -dely, 2);
end

close;
figure;
imagesc(c)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
colorbar
title(colorbar, 'Transmission')
axis on
axis equal
hold on

% Test end -------------


max_axis=max(max(c));
min_axis=min(min(c));
button_axis = questdlg('Need to change the contrast?',...
    'Contrast','Yes','No','Yes');



while strcmp(button_axis,'Yes')
    prompt = {'Maximum value','Minimum value'};
    dlg_title=['valeurs de',num2str(round(min_axis)),' a ',num2str(round(min_axis))];
    max_ax = inputdlg(prompt,dlg_title,1,...
        {num2str((round(max_axis))),num2str(round(max_axis))});
    if ~isempty(max_ax)
    min_axis = (str2double(char(max_ax{2})));
         
        max_axis = (str2double(char(max_ax{1})));
    end
    if min_axis<max_axis
         caxis([min_axis max_axis])
    else
        msgbox('MUST BE min_axis<max_axis')
        
    end

    colorbar
    title(colorbar, 'Transmission')
    axis on
    button_axis = questdlg('Again?',...
        'none','Yes','No','Yes');
end

%Radius = inputdlg('Enter the pinhole radius','radius (mm)',1,{'1'});
%Radius = str2double(char(Radius));
%Radius = floor(double(Radius./(25.4/300)))
% Radius = inputdlg('Enter a value for the radius','radius',1,{'1'});
Radius = 25;
[rows, columns] = size(c);

content='No';

while strcmp(content,'No')
    % FOR POLYGON uncomment bottom three lines
    %H = warndlg('Click to delimit the ROI','click now');
    waitfor(H) %attend la fermeture de la boite
    % [mask0,x0,y0]=roipoly; % de l'image actuelle [c] , delimitation generale
    
    % FOR CIRCLE uncomment bottom five lines
    % [y0,x0]=ginput(1); % this is to get the coordinate of the center of the circle
    
    % after get the coordinate, fix them here
    x0 = 25;
    y0 = 64;
    
    fprintf('\n the coordinate of the circle is: \n');
    fprintf('x0 = %1.0f \n', x0);
    fprintf('y0 = %1.0f \n', y0);

    angles = linspace(0, 2*pi, 10000);
    X_c = cos(angles) * Radius + y0;
    Y_c = sin(angles) * Radius + x0;
    mask0 = poly2mask(X_c, Y_c, rows, columns);

    % ... figure avec le mask
    b = c.*mask0;
    close;figure;
    imagesc(b)
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    colorbar;
    caxis([min_axis max_axis]);
    axis on
    content = questdlg('Happy?','none','Yes','No','Yes');
    if strcmp(content,'No')
        clear b
        close;
        figure;
        imagesc(c);
        caxis([min_axis max_axis]);
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        colorbar
        title(colorbar, 'Transmission');
        axis on
        axis equal
        hold on
    end
end

mesure_energie='oui';
if strcmp(mesure_energie,'oui')
    clear c
    c=b;
    clear b; clear big;
    
    %suppression de zones non désirées.
    %%%...
    % ETRANGE !!A REVOIR !
    %%%...
    button='Yes';
    while strcmp(button,'Yes')
        button = questdlg('Do you want to remove a part of the RCF?','none','Yes','No','Yes');
        if strcmp(button,'Yes')
            [bw,xp,yp]=roipoly;%de l'image actuelle=c
            clear bw
            lig=(round(min(yp)):round(min(max(yp))))';
            for i=round(max(min(xp))):round(min(max(xp))),
                col=ones(max(size(lig)),1).*i;
                in=inpolygon(lig,col,yp,xp);
                for k=1:max(size(in)),
                    if in(k)==1,%on est dans le polygone
                        xx=col(k);
                        yy=lig(k);
                        c(yy,xx)=1;
                    end
                end
            end
        end
        
        imagesc(c)
        caxis([min_axis max_axis])
        colorbar
        title(colorbar, 'Transmission')
        axis on
        axis equal
        
        pause(0.5)
    end
    
    %passage aux log.
    OD=-1.*log10(c);    % change to optical density!
    clear c,clear u
    
    imagesc(OD);
    colorbar;
    hh = colorbar;
    title(hh, 'OD')
    
    axis on
    hold on
    pause(0.5)
    


    %... Mesure de la SURFACE de la ROI et de l'angle solide (str)
       
    N_pxl = sum(sum(mask0)); % nombre total de pixel 
    surf_form = N_pxl*surf_indiv; %surface total de forme
    str_pxl = 4*asin(surf_indiv * 100 / 4 / ((distance_dbl+new_couche_active(n_film))*(distance_dbl+new_couche_active(n_film)))); 
    %Ωcarré ~ 4 Arcsin (p²/h²) expression pour un carré où p est la demi longueur d'un carré et h la hauteur
    %eexpression simplifiée pour h >> p
    %Attention surf_indiv est exprimé en cm^2 et distance_dbl en mm
    str_form = N_pxl * str_pxl;
    
    %stack_thickness % for calulating the divergence of the beam for different types of stacks
    str_mat(n_film-min_film+1)  = str_form;
    surf_mat(n_film-min_film+1) = surf_form;
    % div inutile mais bon ...
    divergence_mat(n_film-min_film+1)=acos(1 - str_form / (2 * pi));
     
      
 
    
    
    % %remplissage des zones masquées
    
    % %button = questdlg('Reconstruire le film?','reconstruction','Tout','En partie','Non','Tout');   %normally it is on
    % % pause(0.5);
    % button = 'Tout';
    % if strcmp(button,'Tout')
    %     min_zone_i=1;max_zone_i=2*rmax;
    %     min_zone_j=1;max_zone_j=2*rmax;
    %     rmax_zone=rmax;M(rmax_zone-1:rmax_zone)=0;
    % elseif strcmp(button,'En partie')
    %     %    H = warndlg('cliquez pour désigner le coin haut gauche de la zone à reconstruire puis le coin opposé','cliquez maintenant');
    %     %    waitfor(H) %attend la fermeture de la boite
    %     [y2,x2]=ginput(2);%y2=row=y x1=column=x !!
    %     min_zone_i=fix(x2(1));max_zone_i=fix(x2(2));
    %     min_zone_j=fix(y2(1));max_zone_j=fix(y2(2));
    %     %  rmax_zone1=max(fix(sqrt((x2(1)-rmax)^2+(y2(1)-rmax)^2)),fix(sqrt((x2(2)-rmax)^2+(y2(2)-rmax)^2)));
    %     %  rmax_zone2=max(fix(sqrt((x2(1)-rmax)^2+(y2(2)-rmax)^2)),fix(sqrt((x2(2)-rmax)^2+(y2(1)-rmax)^2)));
    %     %  rmax_zone=max(rmax_zone1,rmax_zone2);M(rmax_zone-2+pas_R)=0;
    % end
    % if strcmp(button,'Tout') || strcmp(button,'En partie')
    %     for i=min_zone_i:max_zone_i
    %         for j=min_zone_j:max_zone_j
    %             Rn=round(sqrt((i-rmax)^2+(j-rmax)^2));
    %             if OD(i,j)==0 && Rn <= max(size(M)) && Rn~=0
    %                 OD(i,j)=M(Rn);
    %             end
    %         end
    %     end
    % end
    
    % pause(2); close;
    % figure;
    % imagesc(OD)
    % axis equal
    % colorbar
    
    %ede=str2double(char(n_filmb)); % variable to give to stack_RCF, don't ask why the name cde...
    
    %RCF = questdlg('type de RCF?','RCF','HD','MD','MD');
    
    
    %stack_RCF   % gives the different type of RCF film for each film_layer
    RCF=RCF_only{n_film};
    
    
    %transformation des OD_HeNe en dose en utilisant une calibration puis en énergie
    dose_abs = 0;
    compteur = 0;
    somax    = 0;
    tmpLIMIT = 0;

    [NX,NY] = size(OD);

    for i=1:NX
        for j=1:NY
            
            if mask0(i,j) == 1 %Dans le polyn�me d�fini alors :

                if isinf(OD(i,j))==1
                    OD(i,j)=OD(i-1,j);
                end
                
                ax=OD(i,j);     % ax = optical density
                somax=somax+ax;
                
                if ax >1.5 
                    tmpLIMIT = tmpLIMIT +1;
                    % On peut ajouter une ligne pour saturer à OD = 1.7
                end

                if strcmp(RCF,'RCF_HD')     % all must have same byte size!
                    if ax <= 0.4
                        somme=ax*309.29+1511.5*ax^2+15722*ax^3-52556*ax^4+44378*ax^5;%HD in Gray=J/k  %%%  fit red for OD < 0.4, green OD >0.4
                    elseif ax > 0.4
                        somme = ax /8e-4 - 20; % essentially linear with dose up to 250Gy (According spec from Gafchromic)
                        %"- 2" add that the two functions overlapse at the vicinity of 0.4
                    end
                elseif strcmp(RCF,'RCF_MD')
                    if ax <= 0.4
                        somme=ax*183.65 -1506.1*ax^2+7663.5*ax^3-12588*ax^4+7450.6*ax^5;%MD in Gray=J/kg
                    elseif ax > 0.4
                        somme=ax*2.557e2 - 2.53e1;
                    end
                    %Last version uses a linear function, according to the performance data of Gafchromic, the media is essentially linear with dose up to 50 Gy
                    %somme = dose/3.5e-2; % <- I don't understand where it is coming from
                elseif strcmp(RCF,'RCF_EBT2')
                    somme=-6.838+ax*247.59 +1859*ax^2+6589*ax^3-10580*ax^4+6434.5*ax^5; % EBT2 in Gray=J/kg ------ of Marc Glesser 14/03/11
                elseif strcmp(RCF,'RCF_HD_V2')
                    Ymax = 1.4;
                    Beta = 0.75;
                    x0 = 2400;
                    somme=10^(log10(x0)-(1/Beta)*log10((Ymax/ax)-1));
%                     if ax <=0.33
%                         somme=3.5426+ax*521.881 +5250.93*ax^2-26211.8*ax^3+117880.*ax^4-167372.4*ax^5; % HDv2 in Gray=J/kg ------ from Chen et al., 2016, RSI ----- For EPSON 2450
%                     else
%                         somme=2360.06*ax-230.21; %from 0.33, we can extrapolate a linear fitting
%                     end   
                elseif strcmp(RCF,'RCF_EBT3')
                    Ymax = 1.455;
                    Beta = 0.8;
                    x0 = 38;
                    somme=10^(log10(x0)-(1/Beta)*log10((Ymax/ax)-1));
%                     if ax >0.7
%                         somme=17383.47*ax-8130.61; % This expression is already expressed in Gy, not as the following
%                     elseif ax<=0.7
%                         somme=6.326+ax*1465.94 -4962.68*ax^2+47412.7*ax^3-100782*ax^4+79637.3*ax^5; % EBT3 in Gray=J/kg ------ from Chen et al., 2016, RSI ----- For EPSON 2450               
%                         somme=somme*1e-2;%It's expressed in cGy in this fucking article! Take care of stupidities in this fucking world!!
%                     end
                end
                
                dose_abs=dose_abs+somme;
                compteur=compteur+1;
                
            end
        end
    end
    
    
    if strcmp(RCF,'RCF_MD')
        ener=dose_abs*0.908e-3*surf_indiv*30e-4;%si MD    % 1.3e-3 is the density of the active MD-55-2 layer in g/cm^3   -----    30*10^-4 is the thickness of the active MD-55-2 layer in cm
        doseAvg = ener / (surf_mat*0.908e-3* 30e-4);
        %dose_fct_angle=Y1.*18   %   why times 18? Should be wrong!
        ener_fct_angle = 0 ;% To carry out! No important, now!
    elseif strcmp(RCF,'RCF_HD')%RCF=='HD   '
        % energie = dose * density * volum 
        %Density expressed in kg.cm-2 previously value used = 1.3e-3 but the rest of the code doesn't used this value; let's coherent
        % surf_indiv, surface of one pixel in cm2 and the last is the thickness of the active layer in cm
        ener=dose_abs*0.908e-3*surf_indiv*6.5e-4;%si HD
        doseAvg = ener / (surf_mat*0.908e-3* 6.5e-4);
        dose_fct_angle=y1.*309.29+1511.5*y1.^2+15722*y1.^3-52556*y1.^4+44378*y1.^5;    % numbers of the somme above (calibration!)
        ener_fct_angle=dose_fct_angle*0.908e-3*surf_indiv*6.5e-4;        %1.3e-3 is the density of the active HD layer in g/cm^3 - 6.5e-4 is the thickness of the active HD layer in cm
    elseif strcmp(RCF,'RCF_HD_V2')
        %   elseif strcmp(RCF,'HD_V2')
        ener=dose_abs*1.2e-3*surf_indiv*12e-4;%if HD_V2
        doseAvg = ener / (surf_mat*1.2e-3* 12e-4);
        %dose_fct_angle=;    % numbers of the somme above (calibration!)
        ener_fct_angle = 0; % To carry out! No important, now!
        %ener_fct_angle=dose_fct_angle*1.3e-3*surf_indiv*8e-4;     % 1.3e-3 is the density of the active HD layer in g/cm^3 - 8e-4 is the thickness of the active HD_V2 layer in cm
    elseif strcmp(RCF,'RCF_EBT2')
        ener=dose_abs*1.2e-3*surf_indiv*30e-4;%if EBT2
        doseAvg = ener / (surf_mat*1.2e-3* 28e-4);
        %dose_fct_angle=;    % numbers of the somme above (calibration!)
        ener_fct_angle = 0;% To carry out! No important, now!
        %ener_fct_angle=dose_fct_angle*1.3e-3*surf_indiv*30e-4;     % 1.3e-3 is the density of the active HD layer in g/cm^3 - 30e-4 is the thickness of the active EBT2 layer in cm
    else strcmp(RCF,'RCF_EBT3')
        ener=dose_abs*1.2e-3*surf_indiv*28e-4;%if EBT2
        doseAvg = ener / (surf_mat*1.2e-3* 28e-4);
        
        %dose_fct_angle=;    % numbers of the somme above (calibration!)
        ener_fct_angle = 0; % To carry out! No important, now!
        %ener_fct_angle=dose_fct_angle*1.3e-3*surf_indiv*30e-4;     %
        %1.3e-3 is the density of the active HD layer in g/cm^3 - 30e-4 is the thickness of the active EBT2 layer in cm
    end
    %en Joules
    
    
    %1.3 g/cm^3 moyenne de la densite entre muscle et os
    %30 micron=epaisseur des 2 couches sensibles dans chaque film MD
   % ener
    %ener_str;
    %ener_str_IP;
    %ener_moy_pxl=ener/compteur;
    %clear OD
   % pause(1)
end
