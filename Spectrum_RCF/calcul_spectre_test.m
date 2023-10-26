%VERSION fev10
%###### On cheche le meilleur spectre par it�ration ########

% used in main - line 279

%DERNI�RE PARTIE OU ON CALCULE L'INT�GRALE DE L'�NERGIE D�POS�E DANS CHAQUE COUCHE SENSIBLE
%POUR LE COMPARER � LA MESURE D'�NERGIE (JOULES)ET FINALEMENT EN D�DUIRE LE SPECTRE INCIDENT
%EN PART/MeV
%

%global energie Ereel
%global max_film n_HD
global in_layer %e_inc

%u=find(energie);    % gives a vector with the numbers of layers (1,2,3,4....,last exposed film)
%max_film = max(u);  %number of exposed films
%max_film = max_film - 1;

%ici, il faut "deviner" le spectre incident (car cela influe sur l'int�grale de l'�nergie)
%donc il faut proc�der par essai-erreur

%cacul des bornes pour l'it�ration en Energie.

if n_MD == 0                 % maybe here it is everytime like this, because we only have film with one active layer, and not two (like MD)
    max_layer=max_film;
    kk=1;
else
    max_layer=max_film+(n_MD - RCF_counter + max_film); % n_MD - RCF_counter + max_film : nbr of MD w/ signal
    kk=2;
end
size(in_layer)

indice_deb=min(find(in_layer(:,max_layer)));


if max_layer <= size(in_layer,2)-2
    indice_fin=fix(min(find(in_layer(:,max_layer+kk))))-1;
else
    indice_fin=fix(min(find(in_layer(:,max_layer))))+15;  % ???
end

if indice_deb > indice_fin
    indice_fin=indice_deb+2;
end
%% This limits represents the limits of the interval which the Energie cut off is estimated. 
%% Smaller is the interval, less calcul is proceedded by the script

h=figure;
hold on
q=size(in_layer);
for i=1:q(2)%2*max_film,
    plot(e_inc,in_layer(:,i));%.*spectre)
end

title('RCF Pack Response Function')
ylabel('Energy Deposited (MeV/n)')
xlabel('Proton Energy (MeV)')


%interval_energie=[indice_deb,indice_fin]/10
pause(1)


%Construction du profil.


profil='Non';
profil2='Oui';
while strcmp(profil,'Non')
    close
    erreur=1000;
    profil2 = questdlg('Besoin de 2 tendances?','none','Oui','Non','Oui');
    if strcmp(profil2,'Oui')
    %if profil2=='Oui'
        % for e_bord=2:1:0.8*indice_fin
        e_bord = inputdlg('NRJ changement de tendance','e_bord',1,{'7.2'});
        e_bord=str2num(char(e_bord));
        tmp = abs(e_inc-e_bord);
        [val idx_bord] = min(tmp);
        clear tmp val

        i_test= 0;
        stoke = [];

        for indice=indice_deb:1:indice_fin
            for E0=0.01:0.01:2
                for E1=E0+0.1:0.01:2.4
                    clear spectre couche1 Eth
                    for p=1:length(e_inc)
                        spectre(p) = exp(-1.*sqrt(2.*e_inc(p)./E0))./sqrt(2.*e_inc(p).*E0) + exp(-1.*sqrt(2.*e_inc(p)./E1))./sqrt(2.*e_inc(p).*E1);
                    end
    %                alpha=sqrt(E0/E1)*exp(-1*sqrt(2*e_bord)*(1/sqrt(E1)-1/sqrt(E0)));
    %                for p=1:idx_bord
    %                    spectre(p)=alpha.*exp(-1.*sqrt(2.*e_inc(p)./E0))./sqrt(2.*e_inc(p).*E0);
    %                end
    %                for p=idx_bord:length(e_inc)
    %                    spectre(p)=exp(-1.*sqrt(2.*e_inc(p)./E1))./sqrt(2.*e_inc(p).*E1);
    %                end
                    spectre= spectre';
                  %  spectre=spectre(p)';
                    [Eth,erreur1,max_film_couche1,couche1,Ereel]=calc_dif_spec_th_reel(spectre,indice,energie);
                    
                    i_test = i_test +1;
                    stoke(i_test,1:4)=[indice,E0,E1,erreur1];

                    if  erreur1 < erreur
                        erreur=erreur1;
                        E0_b=E0; E1_b=E1; indice_b=indice; couche1_b=couche1;  Eth_b=Eth;...
                            spectre_b=spectre; e_bord_b=e_bord;  Ereel_b=Ereel;
                    end
                end
            end
        end
        %end
    else
        'bin drin'
        E0=0; E0_b=0;E1_b=0;
        for indice=indice_deb:1:indice_fin
            for E1=0.02:0.01:5    %loop to find the best fit for the proton temp (E1).
                clear spectre couche1 Eth
                for p=1:length(e_inc),
                    spectre(p)=exp(-1.*sqrt(2.*e_inc(p)./E1))./sqrt(2.*e_inc(p).*E1);
                end
                spectre=spectre';
                [Eth,erreur1,max_film_couche1,couche1,Ereel]=calc_dif_spec_th_reel(spectre,indice,energie);
                if  erreur1 < erreur   %It goes through the loop and calculates the error at each point (erreur1).
                    clear couche1_b Eth_b spectre_b
                    erreur=erreur1;    %In each loop it compares the error from that loop to the lowest error found so far (erreur).
                    E1_b=E1; indice_b=indice; couche1_b=couche1;  Eth_b=Eth; spectre_b=spectre;  Ereel_b=Ereel;
                end
            end
        end
    end
    E1=E1_b;
    E0=E0_b;
    energie_max=e_inc(indice_b);
    erreur;    %Since the energies are in 0.1 MeV increments then it is multiplied by 10. (E= min:0.1:max)
    
    Eth=Eth_b; 
    Ereel=Ereel_b;
    size_E=max(size(Ereel));
    
    save('spectre_exp.mat','spectre')
    
    h=figure;
    hold on
    plot(1:size_E,log(Ereel),'b*:')
    plot(1:size_E,log(Eth),'r*:')
    caxis('auto')
    title('bleu=log(Ereel/Ereel-max) ; Rouge=log(Eth/Eth-max)');
    
    E0=num2str(E0,'%10.2f');E1=num2str(E1,'%10.2f');energie_max=num2str(energie_max,'%10.1f');...
        erreur=num2str(erreur,'%10.3f');

    graph_name = strcat('E0=',E0,', E1=',E1,', Emax=',energie_max,'MeV, erreur='...
        ,erreur);
    xlabel(graph_name)
    ylabel('log(E dep (MeV))')
    profil = questdlg('Happy?','Non','Oui','Non','Oui');


    if strcmp(profil,'Non')
        E_min_max = questdlg('Modifier Emin et Emax?','none','Oui','Non','Oui');
        if E_min_max=='Oui'
            
            nrj_deb = inputdlg('NRJ min','interval',1,{'5'});
            tmp = abs(e_inc-str2num(char(nrj_deb)));
            [val indice_deb] = min(tmp);

            nrj_fin = inputdlg('NRJ max','interval',1,{'6'});
            tmp = abs(e_inc-str2num(char(nrj_fin)));
            [val indice_fin] = min(tmp);
            
        end
    end
end

clear f

indice=indice_b; couche1=couche1_b; spectre=spectre_b;

%    save 'Eth' e_inc -ASCII -DOUBLE -TABS

%calcule l'�nergie moyenne INCIDENTE des protons qui ont d�pos� de
%l'�nergie, ceci dans chaque couche.POUR UN SPECTRE THEORIQUE!!!

autre_e_inc=(e_inc)';
for f=min_film:max_film,
         if f<=(RCF_counter-n_MD) % maybe here it is everytime like this, because we only have film with one active layer, and not two (like MD)  - same like in line 25
            e_moy(f-min_film+1)=trapz((in_layer(1:indice,f).*spectre(1:indice)...
                ).*autre_e_inc(1:indice))./...
            trapz((in_layer(1:indice,f).*spectre(1:indice)));
         else %2 active layers
            m=-(RCF_counter - n_MD) + 2*f; %max_film - n_MD+2*f  %(f-n_HD)*2+n_HD; 
            e_moy(f-min_film+1)=trapz((in_layer(1:indice,m).*spectre(1:indice)+...
            in_layer(1:indice,m-1).*spectre(1:indice)).*autre_e_inc(1:indice))./...
            trapz((in_layer(1:indice,m).*spectre(1:indice)+...
            in_layer(1:indice,m-1).*spectre(1:indice)));
        end
end
e_moy';

pos_moy = zeros(size(e_moy));
for k = 1: length(e_moy),

    tmp = abs(e_inc-e_moy(k));
    [val idx] = min(tmp);
    pos_moy(k)=idx;
end

%enfin, voici le chiffre � multiplier par l'�nergie d�pos�e mesur�e (en J)
%pour obtenir le nbre de part/MeV
%cherche facteur entre �nergie dep et le DN/DE
%Edep*cst(Emoy)=DN/DE(Emoy)!! calcul� pour un f(E) donc une couche donn�!

cst=spectre(pos_moy)./(1e6.*1.6e-19.*(couche1(1:size_E)').*e_inc_step);
cst=cst';
DN_DE=cst.*energie(1:size_E)';

%if strcmp(type_spectre,'spectre_VS_div')~=1
%    h1=figure;
%    %semilogy(e_moy(1:max(size(e_moy))-1),DN_DE(1:max(size(e_moy))-1),'b*:')
%    semilogy(e_moy,DN_DE,'b*:')
%    label=strcat('energie (Mev)                       Tir',nbr_tir);
%    xlabel(label);
%    
%    %plot(1:max(size(DN_DE_b)),DN_DE_b,'r*:')
%    %caxis('auto')
%end
clear DN_DE1
%DN_DE1(:,1)= DN_DE';
%DN_DE1(:,2)= e_moy';


%calcul l'�nergie totale d�pos� dans le pack RCF.

if strcmp(type_spectre,'spectre_total')
    
    h1=figure;
    %semilogy(e_moy(1:max(size(e_moy))-1),DN_DE(1:max(size(e_moy))-1),'b*:')
    semilogy(e_moy,DN_DE,'b*:')
    label=strcat('Energy (MeV)                       Tir',nbr_tir);
    xlabel(label);

    DN_DE1(:,1)= DN_DE';
    DN_DE1(:,2)= e_moy';

    Edepos_Mev = 0;
    coef_dir = (log10(DN_DE(2))-log10(DN_DE(1)))/(e_moy(2)-e_moy(1));
    ord_orig = log10(DN_DE(1))-coef_dir*e_moy(1);
    e_moy(1) = 1.5; % l'int�grale doit toujours d�buter au m�me point
    DN_DE(1) = 10^(coef_dir*e_moy(1)+ ord_orig);
    %h2=figure;
    %semilogy(e_moy,DN_DE,'b*:')
    for i=1:max(size(e_moy))-1
        delta_E=e_moy(i+1)-e_moy(i);
        Edepos_Mev = Edepos_Mev + trapz(DN_DE(i:i+1).*e_moy(i:i+1))*delta_E;
    end
    Edepos_J= Edepos_Mev*1e6*(1.6e-19);
end

%calcul l'�nergie par str d�pos� dans le pack RCF.

if strcmp(type_spectre,'spectre_par_str')
    
    %Loading angle of the study area
    cd(dirfichier)
    
    DN_DE = DN_DE./str_mat';

    h1=figure;
    %semilogy(e_moy(1:max(size(e_moy))-1),DN_DE(1:max(size(e_moy))-1),'b*:')
    semilogy(e_moy,DN_DE,'b*:')
    label=strcat('Energy (MeV)                       Tir',nbr_tir);
    xlabel(label);

    DN_DE1(:,1)= DN_DE';
    DN_DE1(:,2)= e_moy';

    Edepos_Mev = 0;
    coef_dir = (log10(DN_DE(2))-log10(DN_DE(1)))/(e_moy(2)-e_moy(1));
    ord_orig = log10(DN_DE(1))-coef_dir*e_moy(1);
    e_moy(1) = 1.5; % l'int�grale doit toujours d�buter au m�me point
    DN_DE(1) = 10^(coef_dir*e_moy(1)+ ord_orig);
    %h2=figure;
    %semilogy(e_moy,DN_DE,'b*:')
    for i=1:max(size(e_moy))-1
        delta_E=e_moy(i+1)-e_moy(i);
        Edepos_Mev = Edepos_Mev + trapz(DN_DE(i:i+1).*e_moy(i:i+1))*delta_E;
    end
    Edepos_J= Edepos_Mev*1e6*(1.6e-19);
end
