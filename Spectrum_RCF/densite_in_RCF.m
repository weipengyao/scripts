% in propriete line 36

function [Coef_A_B res] =densite_in_RCF(x,compteur)


%pour avoir finalement S(E,x)*rho en MeV/mm, res=100*rho(g/cc)
%tout en mm !!
%global ep_Plast

global RCF_type
global epaisseur
global new_couche

global ep_mylar_MD
global ep_sens_MD
global ep_epoxy_MD
global ep_mylar_milieu_MD
global total_MD

global ep_surf_HD
global ep_mylar_HD
global ep_sens_HD
global total_HD

global ep_sens_HD_V2
global ep_mylar_HD_V2

global ep_mylar_over_EBT2
global ep_adhesive_EBT2
%global ep_surf_EBT2
global ep_sens_EBT2
%global ep_mylar_under_EBT2
%global ep_topcoat_EBT2

global ep_surf_EBT3
global ep_sens_EBT3

% it is pointing the regions with the mylar density and the sens-layer-activity

x_rel=(x-new_couche(compteur));

%disp(['x_Rel' num2str(x_rel)]);
if x_rel<0
    disp('***************************')
end

%error_str=0;

if strcmp(RCF_type(compteur),'RCF_HD') %on est dans les films HD
    total_HD=epaisseur(compteur);
    if or( x_rel<=ep_surf_HD , and( x_rel>=total_HD-ep_mylar_HD , x_rel<=total_HD )),
        res=1.3970E+02;%Mylar   density
        Coef_A_B=[0.25626,-0.78372];%%Mylar
        %txtstr='mylardensity';
    elseif and(x_rel>ep_surf_HD,...
            x_rel<=ep_surf_HD+ep_sens_HD),
        res=9.0810E+01;%sens layer  ?
        Coef_A_B=[0.25626,-0.78372];%%Active layer      
        %txtstr='sensitive layer';
    else
       disp(['*****************' num2str(x_rel)])
        %error_str=error_str+1;
       % fid = fopen(['Error' num2str(error_str) '.txt'],'wt');
        res=100;
    end
    
   % teststr=(['density in HD' RCF_type(compteur) ' x= ' num2str(x_rel) ' - ' txtstr]);
    
elseif strcmp(RCF_type(compteur),'RCF_MD')%on est dans les films MD
    total_MD=epaisseur(compteur);
    if or( x_rel<ep_mylar_MD , (and( x_rel>=total_MD-ep_mylar_MD , x_rel<total_MD)   )),
        res=1.3970E+02;%Mylar  external layers
        Coef_A_B=[0.25626,-0.78372];%%Mylar
        %txtstr='external layers';
    elseif and(x_rel>=ep_mylar_MD+ep_sens_MD+ep_epoxy_MD,...
                x_rel<ep_mylar_MD+ep_sens_MD+ep_epoxy_MD+ep_mylar_milieu_MD),
        res=1.3970E+02;%Mylar   in between the two sensitive layers
        Coef_A_B=[0.25626,-0.78372];%%Mylar
        %txtstr='between 2 sens layers';
    elseif or(and(x_rel>=ep_mylar_MD,x_rel<ep_mylar_MD+ep_sens_MD),...
            and(x_rel>=ep_mylar_MD+ep_sens_MD+2*ep_epoxy_MD+ep_mylar_milieu_MD,...
            x_rel<ep_mylar_MD+2*ep_sens_MD+2*ep_epoxy_MD+ep_mylar_milieu_MD)),
        res=9.0810E+01;%sens layer
        Coef_A_B=[0.25626,-0.78372];%%Active layer            
        %txtstr='in sensible layers';
    elseif or(and(x_rel>=ep_mylar_MD+ep_sens_MD,x_rel<ep_mylar_MD+ep_sens_MD+ep_epoxy_MD),...
            and(x_rel>=ep_mylar_MD+ep_sens_MD+ep_epoxy_MD+ep_mylar_milieu_MD,...
            x_rel<ep_mylar_MD+ep_sens_MD+2*ep_epoxy_MD+ep_mylar_milieu_MD)),
        res=1.1800E+02;%epoxy
        Coef_A_B=[0.25626,-0.78372];%%Epoxy      
        %txtstr='epoxy';
    end
   % teststr=(['density in MD' RCF_type(compteur) ' x= ' num2str(x_rel) ' - ' txtstr]);
    
elseif strcmp(RCF_type(compteur),'RCF_HD_V2')
        %total_HD_V2=epaisseur(compteur);
    if x_rel<=ep_sens_HD_V2
        res=1.20E+02;%sens layer
        Coef_A_B=[0.2887,-0.7979];
        %txtstr='in sensible layers';    
    elseif and(x_rel>ep_sens_HD_V2,x_rel<=ep_sens_HD_V2 + ep_mylar_HD_V2), %and(x<total_HD_V2),
        res=1.35E+02;%Mylar   density
        Coef_A_B=[0.267,-0.7936];
      %txtstr='external layers';
      else
        disp(['*****************' num2str(x_rel)])
        %error_str=error_str+1;
       % fid = fopen(['Error' num2str(error_str) '.txt'],'wt');
        res=100;
    end
    
elseif strcmp(RCF_type(compteur),'RCF_EBT2')
    %total_EBT2=epaisseur(compteur);    
    if or( (x_rel<ep_mylar_over_EBT2) , (x_rel>(ep_mylar_over_EBT2+ep_adhesive_EBT2+ep_sens_EBT2))) %+ep_topcoat_EBT2
        res=1.35E+02;%Mylar 
        Coef_A_B=[0.267,-0.7936];
    elseif and( (x_rel>=ep_mylar_over_EBT2) , (x_rel<(ep_mylar_over_EBT2+ep_adhesive_EBT2)))
        res=1.2E+02;%epoxy
        Coef_A_B=[0.2916,-0.7988];        
%    elseif and( (x_rel>=ep_mylar_over_EBT2+ep_adhesive_EBT2) , (x_rel<(ep_mylar_over_EBT2+ep_adhesive_EBT2))) %+ep_topcoat_EBT2
%        res=1.2E+02;%topcoat    
%        Coef_A_B=[0.2836,-0.7963];  
    elseif and( (x_rel>=(ep_mylar_over_EBT2+ep_adhesive_EBT2)) , (x_rel<=(ep_mylar_over_EBT2+ep_adhesive_EBT2+ep_sens_EBT2))) %+ep_topcoat_EBT2
        res=1.2E+02;%sens layer
        Coef_A_B=[0.2887,-0.7979];
    end 

elseif strcmp(RCF_type(compteur),'RCF_EBT3')   
    if or( (x_rel<ep_surf_EBT3) , (x_rel>(ep_surf_EBT3+ep_sens_EBT3)))
        res=1.35E+02;%surf layer (matte polyester)   
        Coef_A_B=[0.2709,-0.7948];
    elseif and( (x_rel>=ep_surf_EBT3) , (x_rel<(ep_surf_EBT3+ep_sens_EBT3)))
        res=1.2E+02;%sens layer   
        Coef_A_B=[0.2887,-0.7979];
    end 
    

    
end
%disp(txtstr)
