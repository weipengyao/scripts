%VERSION MARS 2008
%######## Calcul l'écart entre les �nergies d�pos�es mesur�es et cacul�es ######

function [Eth,erreur1,max_film_couche1,couche1,Ereel]=calc_erreur(spectre,indice,energie)

global n_MD max_film RCF_counter
global in_layer
global min_film

%calcule l'�nergie cumul�e d�pos�e dans chaque film (en MeV) c'est la tendance de E_d�pos�e=f(E incident)
%qui est � comparer � l'�nergie mesur�e par dosim�trie pour voir si le spectre suppos� est correct
 for f=min_film:max_film    
      if f<= (RCF_counter - n_MD) %n_HD          % maybe here it is everytime like this, because we only have film with one active layer, and not two (like MD)  - same like in calcul_spectre
        couche1(f-min_film+1)=trapz(in_layer(1:indice,f).*spectre(1:indice,1));
      else        
         m=-(RCF_counter - n_MD) +2*f; %n_HD+(f-n_HD)*2;   
         couche1(f-min_film+1)=trapz(in_layer(1:indice,m).*spectre(1:indice,1))+... 
             trapz(in_layer(1:indice,m-1).*spectre(1:indice,1));   
         
      end
 end

km=find(couche1);
max_film_couche1=max(km);
Eth = couche1(1:max_film_couche1)'./couche1(max_film_couche1-1);
Ereel = energie(1:max_film_couche1)./energie(max_film_couche1-1)';

%Mesure de l'erreur commise avec les valeures prises
erreur1=0;
    for f=1:max_film_couche1
        erreur1=erreur1+sqrt((log(Eth(f))-log(Ereel(f)))^2);
    end
     erreur1= erreur1/max_film;
