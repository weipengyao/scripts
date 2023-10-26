%VERSION décembre 2009
%###### CETTE PREMIÈRE PARTIE CALCULE LE DEPOT D'ENERGIE DANS LES RCF ######

% used in main line 226



%incident energy of the proton in the RCF layers expressed in MeV
global epaisseur_totale
global next_layer
next_layer=0;
%global new_couche;

clear energie_deposee;
%energie_deposee=zeros(i:max(size(e_inc)));
pas=0.0005;
pack=0:pas:epaisseur_totale;%
h=waitbar(0,'Please wait...');
for i=1:max(size(e_inc)),
    %for each incident energy given  for each layer for different type of
    %RCF stacks.
    %total thickness is defined taking into account the Alluminium foil 
    %the evaluation stops when the proton energy comes to be zero.
    %The resolution of the method Runge kutta 4 is given by the variable
    %pas;
   
    %disp([num2str(e_inc(i)) ' / ' num2str(e_inc(end))])
    options = odeset('Events',@events2,'InitialStep',pas,'MaxStep',0.005);
    [x,E,xe,ee,ie]=ode45(@ensemble2,[0 epaisseur_totale],e_inc(i),options);
    %if 
% If the last element of calculation gives an imaginary energy  , the last
% element is neglected;
   
    u=find(imag(E)~=0);
    if isempty(u)==0,
        E1=E(1:u(1)-1);
        x1=x(1:u(1)-1);
        clear x E
        E=E1;
        x=x1;
        clear E1 x1
    end
    
    %At the end , the interpolation gives you the deposited energy solving the differential equation dE/dx (energy/thichness)  
    %from the value of  the stopping power S. 
    %et comme on cherche l'énergie déposée à chaque endroit
    %(pour reconstruire la dose dans chaque film)
    %on la calcule à partir de l'énergie de la particule
    %à chaque endroit
    
    
    energie_deposee(i,:) = interp1(x(1:end-1),diff(E)./diff(x),pack,'pchip',0);
    
    %compar(i)=sum(-1.*energie_deposee(i,:))*pas;
    %sum(-1.*energie_deposee(i,:))*pas;
    %E(end);
    %=energie initiale normalement !!
    %if abs((compar(i)-e_inc(i))/e_inc(i))>0.3
    %   break
    %end
    waitbar(i/max(size(e_inc)),h)
end
close(h)

% for i=1:size(energie_deposee,2)    
%     if sum(abs(energie_deposee(:,i)))==0
%         msgbox(['Energy not high enough: ' num2str(i) 'you need to increase the max energy'])
%         break
%     end
%     
% end

%   *******************************************************************************************************   %%
%   **********************   SECOND PART ******************************************************************   %%
%   *******************************************************************************************************   %%
%   From the calculed deposited energy (previously), the deposit energy is calculed for each sensible layer   %% 


global in_layer;
'numberlayer';
number_layer;
clear in_layer;

%multiply the value dE/dx times pas you finally can have the value of the
%enrgy for the different active layer
%donne l'énergie déposée, en MeV, dans chaque couche sensible
k=1;

for f=1:length(RCF_type)
    f;
    
    if strcmp(RCF_type{f},'RCF_HD')
        [l1,l2]=limite_sens2(f);     % in limite_sens2 are the different layers of the film
        w=find(and(pack>l1(1),pack<l1(2)));
        in_layer(:,k)=-1.*sum(energie_deposee(:,min(w):max(w)),2).*pas;
        %disp (['k=' num2str(k) ' f=' num2str(f) 'RCF_type{f}=' RCF_type{f} 'w=' num2str(w) ' l1 ' num2str(l1) ' l2 ' num2str(l2)])
        % msgbox(num2str(in_layer(:,k)))
        %pause
        k=k+1;
   elseif strcmp(RCF_type{f},'RCF_MD')%soit MD
        [l1,l2]=limite_sens2(f);
        w=find(and(pack>l1(1),pack<l1(2)));
        in_layer(:,k)=-1.*sum(energie_deposee(:,min(w):max(w)),2).*pas;
       %disp (['k=' num2str(k) ' f=' num2str(f) 'RCF_type{f}=' RCF_type{f} 'w=' num2str(w) ' l1 ' num2str(l1) ' l2 ' num2str(l2)])
        % msgbox(num2str(in_layer(:,k)))
        %pause
        w=find(and(pack>l2(1),pack<l2(2)));
        in_layer(:,k+1)=-1.*sum(energie_deposee(:,min(w):max(w)),2).*pas;
        %disp (['k=' num2str(k) ' f=' num2str(f) 'RCF_type{f}=' RCF_type{f} 'w=' num2str(w) ' l1 ' num2str(l1) ' l2 ' num2str(l2)])
        % msgbox(num2str(in_layer(:,k+1)))
        %pause
        k=k+2;
    elseif strcmp(RCF_type{f},'RCF_HD_V2')
        [l1,l2]=limite_sens2(f);
        w=find(and(pack>l1(1),pack<l1(2))); %maybe not working because l1(2)=0, one can change this in limite_sens2.m
        in_layer(:,k)=-1.*sum(energie_deposee(:,min(w):max(w)),2).*pas;
        %disp (['k=' num2str(k) ' f=' num2str(f) 'RCF_type{f}=' RCF_type{f} 'w=' num2str(w)])
        % msgbox(num2str(in_layer(:,k)))
        %pause
        k=k+1;
    elseif strcmp(RCF_type{f},'RCF_EBT2')
        [l1,l2]=limite_sens2(f);
        w=find(and(pack>l1(1),pack<l1(2)));
        in_layer(:,k)=-1.*sum(energie_deposee(:,min(w):max(w)),2).*pas;
        %disp (['k=' num2str(k) ' f=' num2str(f) 'RCF_type{f}=' RCF_only{f} 'w=' num2str(w)])
        % msgbox(num2str(in_layer(:,k)))
        %pause
        k=k+1;
    elseif strcmp(RCF_type{f},'RCF_EBT3')
        [l1,l2]=limite_sens2(f);
        w=find(and(pack>l1(1),pack<l1(2)));
        in_layer(:,k)=-1.*sum(energie_deposee(:,min(w):max(w)),2).*pas;
        k=k+1;
    end
%    disp_v{f}= (['k=' num2str(k) ' f=' num2str(f) 'RCF_type{f}=' RCF_type{f} 'w=' num2str(w) ' l1 ' num2str(l1) ' l2 ' num2str(l2)])
    
end
%
%trace les énergies déposées dans chaque couche sensible, en fonction de l'énergie incidente
%(graphe en MeV/MeV)
h1=figure;
hold on
q=size(in_layer);

%plot  Bragg peak versus energy deposited in each layer;
for i=1:q(2)
    plot(e_inc,in_layer(:,i));
end

title('RCF Pack Response Function')
ylabel('Energy Deposited (MeV/n)')
xlabel('Proton Energy (MeV)')

saveas (h1,[save_name 'RCF_Response'],'fig');
print('-dtiff',[save_name 'RCF_Response']);




