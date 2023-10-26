% in ensemble2 line 14

function [Coef_A_B densite]=propriete(x,compteur)
% msgbox(RCF_type2)
% msgbox(num2str(compteur))
% msgbox(num2str(x))
global RCF_type
% compteur=compteur+1
% iscell(RCF_type2)
% iscellstr(RCF_type2)
% isstruct(RCF_type2)
% isglobal(RCF_type2)
% ischar(RCF_type2) 
% 
% who

%aa=char(RCF_type{compteur})
%stopping unit is here MeV/mm
   if strcmp(RCF_type{compteur},'al')
             Coef_A_B=[0.16644, -.71372];%[0.18765,-.75417];%%Aluminium
            densite=2.7019E+02;     
   elseif strcmp(RCF_type{compteur},'mylar')
             Coef_A_B=[0.25626,-0.78372];%%Mylar
            densite=1.397E+02;
   elseif strcmp(RCF_type{compteur},'ta')
             Coef_A_B=[ 0.06885341, -0.54111393];%%Tantalum
            densite=1.66E+03;
   elseif strcmp(RCF_type{compteur},'cu')
             Coef_A_B=[0.13008,-.7023];%%Cuivre
             densite=8.96E+02;
%              Coef_A_B=[0.21,-.77];%%Cuivre
%             densite=2.0E-04;
   elseif strcmp(RCF_type{compteur},'ni')
              Coef_A_B=[0.12693, -0.633];%%nickel      % <-maybe wrong, calculated with SRIM and an excel/matlab curve fit
             densite=8.91E+02;
   elseif strcmp(RCF_type{compteur},'PET')
             Coef_A_B=[0.25626,-0.78372];%%Mylar/PET
            densite=1.397E+02;
   elseif strcmp(RCF_type{compteur},'fe') 
             Coef_A_B=[0.14323,-.71568];%%Fer
            densite=7.96E+02;
   elseif strcmp(RCF_type{compteur},'pb')
             Coef_A_B=[0.07077,-.62748];%%Plomb
            densite=1.134E+03;
   elseif strcmp(RCF_type{compteur},'zr')
             Coef_A_B=[0.09514525, -0.52576844];%%Zirconium
            densite=6.49E+03;
   elseif strcmp(RCF_type{compteur},'CR39')   
             Coef_A_B=[0.26855,-.78898];%%CR39
            densite=1E+02;
   elseif strcmp(RCF_type{compteur},'mo')   
             Coef_A_B=[0.10969,-.6813];%%molybdene
            densite=1.02E+03;
   elseif strcmp(RCF_type{compteur},'PLA')
           Coef_A_B=[0.25626,-0.78372];
           densite=1.397E+02;       
%    elseif strcmp(RCF_type{compteur},'OLD') %OLD layer       %OLD comparé ac Active SRIM
%            Coef_A_B=[0.24118,-0.7221];
%            densite=1.20E+02;
%    elseif strcmp(RCF_type{compteur},'ACT') %Active layer
%            Coef_A_B=[0.2887,-0.7979];
%            densite=1.20E+02;
%    elseif strcmp(RCF_type{compteur},'ADH') %Adhesive layer
%            Coef_A_B=[0.2916,-0.7988];
%            densite=1.20E+02;
%    elseif strcmp(RCF_type{compteur},'MPoly') %Matte Polyester layer
%            Coef_A_B=[0.2709,-0.7948];
%            densite=1.35E+02;
%    elseif strcmp(RCF_type{compteur},'OPoly') %Polyester Overlaminate layer
%            Coef_A_B=[0.267,-0.7936];
%            densite=1.35E+02;
%    elseif strcmp(RCF_type{compteur},'TCoat') %Top Coat layer
%            Coef_A_B=[0.2836,-0.7963];
%            densite=1.20E+02;           
   else
            %Coef_A_B=[0.24118,-0.72241];%[0.25626,-0.78372];%%Mylar
            [Coef_A_B,densite] = densite_in_RCF(x,compteur); %other program "densite_in_RCF"
   end      

%disp(['densite ' num2str(densite) ' of x=' num2str(x) ' for compteur' num2str(compteur)])


