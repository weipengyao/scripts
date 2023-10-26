% used in calcul_profil_depo_E  - line 63 etc

function [layer1,layer2]=limite_sens2(f)
%

global RCF_type
global RCF_only
global new_couche

global ep_mylar_MD
global ep_sens_MD
global ep_epoxy_MD
global ep_mylar_milieu_MD

global ep_surf_HD
global ep_mylar_HD
global ep_sens_HD

global ep_sens_HD_V2
global ep_mylar_HD_V2   %maybe not necessary?

global ep_mylar_over_EBT2
global ep_adhesive_EBT2
global ep_surf_EBT2
global ep_sens_EBT2
global ep_mylar_under_EBT2
global ep_topcoat_EBT2

global ep_surf_EBT3
global ep_sens_EBT3

%
RCF_only;
%disp(['f ' num2str(f)])

if strcmp(RCF_type{f},'RCF_HD') %HD
    
    debut_film=new_couche(f);   %entire thickness?
    layer1(1)=ep_surf_HD+debut_film;
    layer1(2)=ep_surf_HD+debut_film + ep_sens_HD;     %where is the mylar of HD?
    layer2(1)=0;
    layer2(2)=0;
    
elseif strcmp(RCF_type{f},'RCF_MD') %MD
    debut_film=new_couche(f);
    layer1(1)=ep_mylar_MD+debut_film;
    layer1(2)=ep_mylar_MD+debut_film + ep_sens_MD;
    layer2(1)=ep_mylar_MD+debut_film + ep_sens_MD + 2*ep_epoxy_MD + ep_mylar_milieu_MD;
    layer2(2)=ep_mylar_MD+debut_film + 2*ep_sens_MD+2*ep_epoxy_MD+ep_mylar_milieu_MD;
    
elseif strcmp(RCF_type{f},'RCF_HD_V2') %HD_V2
    debut_film=new_couche(f);
    layer1(1)=debut_film;
    %layer1(2)=0;    
    layer1(2)=ep_sens_HD_V2+debut_film;     %maybe not necessary?
    layer2(1)=0;
    layer2(2)=0;
    
elseif strcmp(RCF_type{f},'RCF_EBT2') %EBT2
    debut_film=new_couche(f); 
    layer1(1)=ep_mylar_over_EBT2+ep_adhesive_EBT2+debut_film;
    layer1(2)=ep_mylar_over_EBT2+ep_adhesive_EBT2+ep_sens_EBT2+debut_film;
    layer2(1)=0;
    layer2(2)=0;
    
elseif strcmp(RCF_type{f},'RCF_EBT3') %EBT3
    debut_film=new_couche(f); 
    layer1(1)=ep_surf_EBT3+debut_film;
    layer1(2)=ep_surf_EBT3+ep_sens_EBT3+debut_film;
    layer2(1)=0;
    layer2(2)=0;
    
    
end
%disp(['limitesens ' num2str(debut_film) ' ' num2str(layer1(1)) ' ' num2str(layer1(2))])