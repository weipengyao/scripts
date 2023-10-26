function res=ensemble2(x,E,new_couche)

global Coef_A_B
global densite
global new_couche
global RCF_type
global epaisseur_totale
%E en MeV
%x en mm

if x<0
    res=NaN;
elseif x<epaisseur_totale
    temp=[];
    temp=find((x-new_couche)>=0);
    compteur=max(size(temp)); %cherche dans quelle couche on est.
    %disp(['x ' num2str(x) 'compteur ' num2str(compteur)])
    [Coef_A_B densite]= propriete(x,compteur);
    res=-1.*Coef_A_B(1).*E.^(Coef_A_B(2)).*densite;
   %disp(res)
else
    
    res=0;
end




% for x=0:0.01:0.5;
%     temp=[];
%     temp=find((x-new_couche)>=0);
%     compteur=max(size(temp)); %cherche dans quelle couche on est.
%     disp(['x ' num2str(x) 'compteur ' num2str(compteur)])
% end
