%    VERSION  Mars 2015 Patrizio Antici
%#################################     PROGRAMME PRINCIPAL       ###############################


% scanner parameters: Pos film  black & white photo 16bits
%       exposure=-1   gamma=1    highlight=255    shadow=0




%
%incident energy 
%#e_inc: defined in calc_profil_depo; 

%energy matrix used to evaluated the theoretical energy (energy deposit-
%stopping power)
%#energy defined in calc_profil_depo;

%minimal energy necessary to pass through the film
%# E_nfilm: defined in main

%energy define in dosec_ptfilm
%rate between Emin to pass through the film n and the Emin energy to pass through the
%last film 
% #En_divi_Emax: defined in dosecfilm;

%energy deposited in each film 
% #energie: defined in  dosec_ptfilm;

%the energy of the proton in the active layer
%#in_layer: defined in  calc_profil_depo;

%number of film H_D and MD defined in init_values;
% #n_HD et n_MD: 


clc;
close all;
clear all;



dircur=pwd;

global RCF_only
global RCF_type
global RCF_counter
global n_HD n_MD
global energie ener
global doseROI doseAvg
global max_film
global min_film
global a
global epaisseur_totale   %spessore -thickness
%strato -layer
global new_couche
global epaisseur


%All units in mm

% Parameters for different types of RCFs; 
global ep_PET  %thickness of PET-layer in mm
ep_PET=0.500;
global ep_ni
ep_ni=0.05;    %thickness of Ni-layer in mm

global ep_ta   %thickness of Ta-layer in mm
ep_ta=0.050;
global ep_al   %thickness of Al-layer in mm
ep_al=0.015;  
global ep_al_filter1
ep_al_filter1=0.1;
global ep_al_filter2
ep_al_filter2=0.5;
global ep_al_filter3
ep_al_filter3=1;
global ep_al_filter4
ep_al_filter4=1.5;

global ep_fe_filter   %thickness of Fe-layer in mm
ep_fe_filter=0.025;
global ep_cu
ep_cu=2.0;    %thickness of Cu-layer in mm
global ep_pb
ep_pb=1.0;    %thickness of Pb-layer in mm


% Film MD
global ep_mylar_MD
ep_mylar_MD=0.067;
%couche sensible aux protons   - sensitive layer 
global ep_sens_MD
ep_sens_MD=0.015;
%colle                      - colla / glue
global ep_epoxy_MD
ep_epoxy_MD=0.025;
%mylar
global ep_mylar_milieu_MD
ep_mylar_milieu_MD=0.025;
%thickness of MD55
global total_MD
total_MD=0.239;
%global n_MD

%film HD-810
global ep_surf_HD   %thickness surface layer in mm
ep_surf_HD=0.00075;
global ep_sens_HD   %thickness active layer
ep_sens_HD=0.0065;
global ep_mylar_HD  %thickness polyester substrate
ep_mylar_HD=0.097;
global total_HD     %total thickness
total_HD=0.10425;
 
%film HD-V2
global ep_sens_HD_V2
ep_sens_HD_V2=0.012;
global ep_mylar_HD_V2
ep_mylar_HD_V2=0.097;
global total_HD_V2
total_HD_V2=0.109;
%global n_HD_V2
%load n_HD_V2

%film EBT2 (1 sensitive layer)
global ep_mylar_over_EBT2
ep_mylar_over_EBT2=0.05;
global ep_adhesive_EBT2
ep_adhesive_EBT2=0.025;
%    global ep_topcoat_EBT2
%    ep_topcoat_EBT2=0.005;
global ep_sens_EBT2
ep_sens_EBT2=0.028;
global ep_mylar_under_EBT2
ep_mylar_under_EBT2=0.175;
global total_EBT2
total_EBT2=0.283;
%global n_EBT2;
%load n_EBT2;

%film EBT3 
global ep_surf_EBT3
ep_surf_EBT3=0.125;
global ep_sens_EBT3
ep_sens_EBT3=0.028;
global total_EBT3
total_EBT3=ep_sens_EBT3+2.*ep_surf_EBT3;
%global n_EBT
%load n_HD

%  Alu=Aluminum, HD1=HD810, HD2=HD_V2, MD1=MD55, EBT2=EBT2, PET= plastic
%  layer , PLA= 2*PET (double of thickness)- same material;
save_name='Analyse_RCF_F1_2023';
etude='Yes';
while strcmp(etude,'Yes')
    analyse = questdlg('What do you want to do?','Please select one option:','1 Define RCF pack','2 Retrieve Energy in RCF','3 Reconstruct Spectre','1 Define RCF pack');
    
  
    if strcmp(analyse,'1 Define RCF pack')
         %Default stack composition
%          defaultanswer=  {'Alu;HD1;MD1;MD1;',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','15'};

       %TITAN
          %defaultanswer=  {'Alu;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;HD2;PET;HD2;PET;HD2;PET;PET;HD2;PET;PET;HD2;PET;PET;HD2;PET;PET;HD2;PET;PET;HD2;PET;PET;HD2;PET;PET;HD2;PET;PET;HD2',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','35'};
          
       %SFA Commissioning
          %defaultanswer=  {'Alu;EB3;EB3;EB3;AL1;EB3;AL1;EB3;AL1;EB3;AL1;AL1;EB3;AL1;AL1;EB3;AL1;AL1;EB3;AL1;AL1;EB3;AL1;AL1;EB3;AL1;AL1;AL1;EB3;AL1;AL1;AL1;EB3;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;AL1;EB3;',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','35'};

       %LULI - July2020
          %defaultanswer=  {'Alu;PET;PET;PET;HD2;PET;PET;HD2;PET;PET;EB3;PET;PET;EB3;PET;PET;EB3',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','25'};
          %defaultanswer=  {'Alu;HD2;HD2;PET;HD2;PET;HD2;PET;PET;HD2;PET;PET;PET;HD2;PET;PET;PET;PET;HD2;PET;PET;PET;PET;HD2;PET;PET;PET;PET;EB3;PET;PET;PET;PET;EB3;PET;PET;PET;PET;EB3;PET;PET;PET;PET;PET;EB3;PET;PET;PET;PET;PET;EB3;PET;PET;PET;PET;PET;PET;EB3;PET;PET;PET;PET;PET;PET;PET;EB3;',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','30'};
          %defaultanswer=  {'Alu;HD2;HD2;HD2;HD2;HD2;HD2;EB3;EB3;EB3;EB3;EB3;',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','30'};    
          %defaultanswer=  {'Alu;EB3;EB3;EB3;Alu;Alu;Alu;Alu;Alu;Alu;Alu;Alu;Alu;Alu;EB3;',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','30'};
          %defaultanswer=  {'AL1;PET;PET;PET;HD2;AL1;AL1;AL1;HD2;PET;PET;PET;PET;EB3;PET;PET;PET;PET;EB3;PET;PET;PET;PET;PET;PET;EB3;PET;PET;PET;PET;PET;PET;EB3;PET;PET;PET;PET;PET;PET;EB3;AL1;EB3;AL1;EB3;',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','35'};
          %defaultanswer=  {'Tan;EB3;EB3;EB3;EB3;EB3;EB3;EB3;EB3;EB3;EB3;EB3;EB3;EB3;',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','40'}; 
           
          %defaultanswer=  {'Alu;Al1;Al1;HD2;Al1;Al1;Al1;Al1;HD2;Al1;Al2;HD2;Al1;Al1;Al1;Al2;HD2;Al3;HD2;Al1;Al3;HD2;Al1;Al1;Al3;EB3;Al1;Al1;Al1;Al3;EB3;Al4;EB3;Al1;Al4;EB3;Al1;Al1;Al4;EB3;Al1;Al1;Al1;Al1;Al4;EB3;Al3;AL3;EB3;Al1;Al3;Al3;EB3;Al1;Al1;Al1;Al3;Al3;EB3;Al1;Al1;Al1;Al3;Al3;EB3',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','85'};
          %defaultanswer=  {'Alu;Al1;Al1;HD2;Al1;Al1;Al1;Al1;HD2;Al1;Al2;HD2;Al1;Al1;Al1;Al2;HD2;Al3;HD2;Al1;Al3;HD2;Al1;Al1;Al3;EB3;Fer;EB3;Al1;Al1;Fer;EB3;Al1;Al1;Al1;Al1;Fer;EB3;Al2;Fer;EB3;Al1;Al2;Fer;EB3;Al1;Al1;Al2;Fer;EB3;Al1;Al1;Al1;Al1;Al2;Fer;EB3;Al3;Fer;EB3;Al3;Fer;EB3',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','85'};
         
       %APOLLON - DPM 2022   
           % defaultanswer=  {'Alu;EB3;EB3;EB3;AL1;EB3;AL1;EB3;AL1;EB3;AL1;AL1;EB3;AL1;AL1;EB3;AL1;AL1;EB3;AL1;AL1;EB3;AL1;AL1;EB3;AL1;AL1;AL1;EB3;AL1;AL1;AL1;EB3;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;AL1;EB3;AL1;AL1;AL1;AL1;AL1;EB3;',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','45'};

       %Qualif F1
          
          % design #2
          % defaultanswer=  {'Alu;COP;EB3;AL1;AL1;AL1;FER;EB3;AL4;EB3;PLO;EB3;PLO;EB3;PLO;EB3;PLO;EB3;PLO;EB3;PLO;EB3;PLO;EB3',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','90'};
          
          % design #3
          % defaultanswer=  {'Alu;Al1;Al1;HD2;Al1;Al1;Al1;Al1;HD2;Al1;Al2;HD2;Al1;Al1;Al1;Al2;HD2;Al3;HD2;Al1;Al3;HD2;FER;EB3;AL4;EB3;PLO;EB3;PLO;EB3;PLO;EB3;PLO;EB3;PLO;EB3;PLO;EB3;PLO;EB3',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','90'}; 
          
          % design #4
          defaultanswer=  {'Alu;Al1;Al1;HD2;Al1;Al1;HD2;Al1;Al1;Al1;Al1;HD2;Al1;Al1;Al1;Al1;HD2;Al2;HD2;Al2;HD2;Al3;EB3;Al4;EB3;Al4;EB3;Al4;EB3;Al4;EB3;Al4;EB3;Al3;EB3;Al3;EB3',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','80'};
          
          %defaultanswer=  {'Alu;Al1;Al1;HD2;Al1;Al1;HD2;Al1;Al1;Al1;Al1;HD2;Al1;Al1;Al1;Al1;HD2',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','25'};

       %RAL 2023      
          %defaultanswer=  {'Alu;FER;HD2;FER;HD2;FER;HD2;FER;FER;HD2;FER;FER;HD2;FER;FER;HD2;EB3;FER;FER;HD2;EB3;FER;FER;FER;FER;HD2;EB3;FER;FER;FER;FER;HD2;EB3;FER;FER;FER;FER;HD2;EB3',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','90'};

                   
       %Test
          %defaultanswer=  {'Alu;HD2;PET;HD2;PET;HD2',num2str(ep_al*1000),num2str(ep_PET*1000),num2str(ep_cu*1000),num2str(ep_ta*1000),'0.1','0.1','35'};

          
           Stack_info= inputdlg({'Shot composition (ALWAYS 3 DIGITS)','Thickness Alu Filter in um;','Thickness PET Filter in um:','Thickness Cu (ou Zr) Filter in um:','Thickness Ta Filter in um:','Minimum proton energy:','Step','Maximum proton energy'},'FIX Parameters',[3, 150],defaultanswer);
        stack_data=char(Stack_info(1));
        ep_al =str2double(char(Stack_info(2)))/1000;    %new max_film for while-loop!!!!
        ep_PET=str2double(char(Stack_info(3)))/1000;
        ep_cu =str2double(char(Stack_info(4)))/1000;
        ep_ta =str2double(char(Stack_info(5)))/1000;
        e_inc_start=str2double(char(Stack_info(6)));
        e_inc_step =str2double(char(Stack_info(7)));
        e_inc_stop =str2double(char(Stack_info(8)));
        
        
        e_inc=e_inc_start:e_inc_step:e_inc_stop;
      
        counter=1;    %don't take aluminum in consideration
        %initiation of the values
        RCF_counter=0;    
        number_layer=0;
        datafields = cell(size(stack_data));
        new_couche(1)=0;
        new_couche_active(1)=0;

        %check the different type of RCFs and passive material :
        %Nichel,Copper,Plastic layer 
        n_MD=0;
        n_HD=0;
        n_EBT2=0;
        n_HD2=0;
        n_EBT3=0;
        
        for i=1:4:length(stack_data)
            datafields{counter}=upper(stack_data(i:i+2));
            
            if strcmp(datafields{counter},'ALU')
                RCF_type{counter}='al';
                epaisseur(counter)=ep_al;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
            
            elseif strcmp(datafields{counter},'AL1')
                RCF_type{counter}='al';
                epaisseur(counter)=ep_al_filter1;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
                
            elseif strcmp(datafields{counter},'AL2')
                RCF_type{counter}='al';
                epaisseur(counter)=ep_al_filter2;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter); 
                
            elseif strcmp(datafields{counter},'AL3')
                RCF_type{counter}='al';
                epaisseur(counter)=ep_al_filter3;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);    
             
            elseif strcmp(datafields{counter},'AL4')
                RCF_type{counter}='al';
                epaisseur(counter)=ep_al_filter4;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
                
           elseif strcmp(datafields{counter},'FER')
                RCF_type{counter}='fe';
                epaisseur(counter)=ep_fe_filter;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
                
           elseif strcmp(datafields{counter},'PLO')
                RCF_type{counter}='pb';
                epaisseur(counter)=ep_pb;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);     
            
            elseif strcmp(datafields{counter},'NIC')
                RCF_type{counter}='ni';
                epaisseur(counter)=ep_ni;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
            
            elseif strcmp(datafields{counter},'COP')
               RCF_type{counter}='cu';
               epaisseur(counter)=ep_cu;
               new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
            
            elseif strcmp(datafields{counter},'ZIR')
               RCF_type{counter}='zr';
               epaisseur(counter)=ep_cu;
               new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
            
            elseif strcmp(datafields{counter},'TAN')
                RCF_type{counter}='ta';
                epaisseur(counter)=ep_ta;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
            
            elseif strcmp(datafields{counter},'HD1')
                n_HD=n_HD+1;
                RCF_counter=RCF_counter+1;
                RCF_type{counter}='RCF_HD';
                epaisseur(counter)=total_HD;
                RCF_only{RCF_counter}='RCF_HD';
                number_layer=number_layer+1;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
                new_couche_active(RCF_counter)=new_couche(counter);
            
            elseif strcmp(datafields{counter},'HDR')
                RCF_counter=RCF_counter+1;
                RCF_type{counter}='RCF_HDrev';
                epaisseur(counter)=total_HD;
                RCF_only{RCF_counter}='RCF_HDrev';
                number_layer=number_layer+1;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
                new_couche_active(RCF_counter)=new_couche(counter);
            
            elseif strcmp(datafields{counter},'MD1')
                n_MD=n_MD+1;
                RCF_counter=RCF_counter+1;
                RCF_type{counter}='RCF_MD';
                epaisseur(counter)=total_MD;
                RCF_only{RCF_counter}='RCF_MD';
                number_layer=number_layer+2;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
                new_couche_active(RCF_counter)=new_couche(counter);
            
            elseif strcmp(datafields{counter},'PET')
                RCF_type{counter}='PET';%PET'; % 
                epaisseur(counter)=ep_PET;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
           
           elseif strcmp(datafields{counter},'PLA')
                RCF_type{counter}='PLA';
                epaisseur(counter)=2*ep_PET;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);

            elseif strcmp(datafields{counter},'HD2')
                n_HD2 = n_HD2+1;
                RCF_counter=RCF_counter+1;
                RCF_type{counter}='RCF_HD_V2';
                epaisseur(counter)=total_HD_V2;
                RCF_only{RCF_counter}='RCF_HD_V2';
                number_layer=number_layer+1;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
                new_couche_active(RCF_counter)=new_couche(counter);

            elseif strcmp(datafields{counter},'EB2')
                n_EBT2 = n_EBT2 +1;
                RCF_counter=RCF_counter+1;
                RCF_type{counter}='RCF_EBT2';
                epaisseur(counter)=total_EBT2;
                RCF_only{RCF_counter}='RCF_EBT2';              
                number_layer=number_layer+1;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
                new_couche_active(RCF_counter)=new_couche(counter);

            elseif strcmp(datafields{counter},'EB3')
                n_EBT3 = n_EBT3 +1;
                RCF_counter=RCF_counter+1;
                RCF_type{counter}='RCF_EBT3';
                epaisseur(counter)=total_EBT3;      
                RCF_only{RCF_counter}='RCF_EBT3';  
                number_layer=number_layer+1;
                new_couche(counter+1)=new_couche(counter)+epaisseur(counter);
                new_couche_active(RCF_counter)=new_couche(counter);

            end
            
            counter=counter+1;
        end
        
                
        %Calculates total thickness of pack
        epaisseur_totale=sum(epaisseur);
        
        %PAS LA MEILLEUR METHODE
        %{
        %initiation for the different layer
        new_couche(1)=0;
        
        %epaisseur(1);
        
        for i=2:max(size(epaisseur))
            if epaisseur(i)>0 
               new_couche(i)=new_couche(i-1)+epaisseur(i-1);     % stack thickness as sum in a vector
            end
        end
%}

        
        % to display  number of layers,  different RCF types, total thickness; 

        for ttt =1:numel(epaisseur)-1
            fprintf('%0.0f: ',ttt);
            fprintf('%s = ',cell2mat(RCF_type(ttt)));
            fprintf('%0.3f mm %0.3f mm \n ',epaisseur(ttt),new_couche(ttt+1));
        end
        ttt=ttt+1;
        fprintf('%0.0f: ',ttt);
        fprintf('%s = ',cell2mat(RCF_type(ttt)));
        fprintf('%0.3f mm %0.3f mm \n ',epaisseur(ttt),new_couche(ttt)+epaisseur(ttt));
          
       
        calcul_profil_depo_E
        
        %save('profile_depo_E.mat','e_inc','in_layer')
        
        % Find minimum Energy to which RCF is sensitive to. MD have 2
        % active 
        % layers, HD only one.
        
        m=1; 
        %number of layers
        i=1;  
        
        E_nfilm= zeros;
        for f=1:RCF_counter
            try
                if strcmp(RCF_type{f},'RCF_HD')
                    E_nfilm(i)=e_inc(find(in_layer(:,m),1,'first')); % energie min for reaching film m.
                    m=m+1;i=i+1;
                elseif strcmp(RCF_type{f},'RCF_MD')     %-----------!!!!!!!!!!----------
                    E_nfilm(i)=e_inc(find(in_layer(:,m),1,'first'));
                    m=m+2;i=i+1;
                elseif strcmp(RCF_type{f},'RCF_HD_V2')
                    E_nfilm(i)=e_inc(find(in_layer(:,m),1,'first')); % energie min for reaching film m.
                    m=m+1;i=i+1;
                elseif strcmp(RCF_type{f},'RCF_EBT2')
                    E_nfilm(i)=e_inc(find(in_layer(:,m),1,'first')); % energie min pour atteindre le film m.
                    m=m+1;i=i+1;
                elseif strcmp(RCF_type{f},'RCF_EBT3')
                    E_nfilm(i)=e_inc(find(in_layer(:,m),1,'first')); % energie min pour atteindre le film m.
                    m=m+1;i=i+1;
                end
            catch ME
                fprintf(ME.identifier,'\r\n');
                fprintf('\nArret au film no %d \n',f)
                fprintf('Type du RCF (premier RCF non-affect�): %s \n',RCF_type{f})
                break
            end
        end
        
       
        max_film=RCF_counter;
        
        save(save_name)
        %pause(1)
    end
        % **********************************************************************************************************  %%
        % *************** END OF CALCULATION OF DEPOSITED ENERGY   *************************************************  %%
        % **********************************************************************************************************  %% 
        
        
        
        % *********************************************************************************************************  %%
        % ***************       CALCULATE ENERGY DEPOSITED IN RCF      ********************************************  %%
        % *********************************************************************************************************  %%
        
        
        if strcmp(analyse,'2 Retrieve Energy in RCF')

            %Exctraction form the saved data
            load(save_name)

            clear energie
            clear doseROI
            clear energie_fct_angle_str
            clear rmax_mat      
            clear divergence_mat
            clear surf_mat
            clear str_mat
            clear sommeOD
            
            resolution = inputdlg('Enter the resolution of the RCF-scan in dpi','resolution (dpi)',1,{'200'});
            resolution=str2double(char(resolution));
            
            
            %resolution factor 1px= 25.4 micrometers
            facteur_mm_px=resolution/25.4;%px/mm
            %surface resolution
            surf_indiv=(2.54/resolution)^2;% d'un pixel carre en cm^2       cm^2/px^2
            
            
            %stack_comp
            
            loading=0;
            n_film=1;
            if max_film<2
                max_film=10;
            end
            
            max_sheet_film=18;         %max_sheet_film: amount of films per scan-sheet; max_film: what ever..., it will ask you about the amount of films with a signal, this will be max_film in future
            
            current_background_type='';
            while n_film <= max_film    %max_film amount of films with signal
                
                clear En_divi_Emax
                
                if loading == 0 % very first use of main, first loading
                    loading=1;
                    n_film = inputdlg('Enter the number of the film with which you want to start','film number',1,{'1'});
                    n_film1=str2double(char(n_film));  % n_film is not a variable, it is a string.
                    if n_film1 == 1  % recharge les fichiers de sortie pour les completer
                        
                        %load E_nfilm;  ### A SAUVER DANS LE WORKSPACE
                        defaultanswer=  {'01','0',num2str(RCF_counter),'25','1'} ;                     
                        info_RCF= inputdlg({'Enter shot number','How many RCFs saturated(count only active film,no plastic)?','What is the number of the last RCFs with signal?','distance to RCF pack in (mm)','number of films per scanned sheet'},'Properties',3,defaultanswer);
                        nbr_tir=char(info_RCF(1));
                        min_film=str2double(char(info_RCF(2)))+1;
                        max_film=str2double(char(info_RCF(3)));    %new max_film for while-loop!!!!
                        distance_dbl=str2double(char(info_RCF(4)));
                        max_sheet_film=str2double(char(info_RCF(5)));
                        
                        if max_film>RCF_counter
                            msgbox('Your Pack calculated in Step 01 has less RCFs than you have to analize. Please repeat Step 1 adapting the RCF number')
                            break
                        end
                        
                        
%                       % Initialize values of Energy string
                        energie=zeros(max_film-min_film+1,1);
                        doseROI=zeros(max_film-min_film+1,1);
                        surf_mat = zeros(max_film-min_film+1,1);
                        str_mat  = zeros(max_film-min_film+1,1);
                        divergence_mat=zeros(max_film-min_film+1,1);
                        rmax_mat=zeros(max_film-min_film+1,1);
                        
                        close all;
                        
                        energie_str(1:max_film-min_film+1)=0;
                        energie_moy_pxl(1:max_film-min_film+1)=0; 
                        X=atan(((1:10000)/facteur_mm_px)/distance_dbl);
                        X1=atan(((1:500:10000)/facteur_mm_px)/distance_dbl);
                        energie_fct_angle_str(1:max_film-min_film+1,1:max(size(X1)))=0;
%                       
                        clear n_film1
                        save workspace
                        n_film=min_film;
                        %load the data file
                        try
                        if exist('dirfichier', 'file')
                           cd(dirfichier)
                        end 
                        catch
 
                            warningMessage = sprintf('Warning: file does not exist:\n%s', 'dirfichier');
                            uiwait(msgbox(warningMessage));
                        end 
                       
                       

                        [nom,chemin]=uigetfile('*.tif','Fichier scan RCF (16 bits, .tif)');
                        cd(chemin)
                        mkdir(datestr(now,29));
                        cd(datestr(now,29))
                        dirfichier=pwd;
                        cd(dircur);
                        save 'dirfichier' dirfichier
                        fic_im=strcat(chemin,nom);
                        a=imread(fic_im);
                        
                        
                        
                        
                        
                        
                    else % if you not start with n_film=1
                        load dirfichier %si n'�xiste pas ou erron� il faut le cr�er!
                        load workspace
                        n_film=n_film1;
                        cd(dirfichier)
                        [nom,chemin]=uigetfile('*.tif','Fichier scan RCF (16 bits, .tif)');    % could load whole scan-sheet
                        cd(dircur);
                        fic_im=strcat(chemin,nom);
                        a=imread(fic_im);
                        
                        
                        %----background "function"
                        h=figure;
                        %set(h,'Position',[100 100 1000 1000])
                        imagesc(a);
                        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
                        axis equal
                        
                        background2
                        %----
                        
                    end
                end
                
                
                if (n_film-1)/max_sheet_film==fix((n_film-1)/max_sheet_film) && n_film ~= 1 % charge une nouvelle image (ici max_film films par image) ATTENTION; IF MAX_FILM IS
                    cd(dirfichier)
                    [nom,chemin]=uigetfile('*.tif','Fichier scan RCF (16 bits, .tif)');
                    cd(dircur);
                    fic_im=strcat(chemin,nom);
                    nom_seul=nom(1:max(size(nom))-4);
                    a=imread(fic_im);
                    
                    colorbar
                    %----background "function"
                    h=figure;
                    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
                    %set(h,'Position',[100 100 1000 1000])
                    imagesc(a);
                    axis equal
                    
                    background2
                    %----
                end
                
                
                
                if ~strcmp(RCF_only{n_film},  current_background_type)
                    h=figure;
                    title(current_background_type);
                    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
                    imagesc(a);
                    axis equal
                    background2     % opens background function
                end
                
                
                dosec_ptfilm
                
                sommeOD(n_film-min_film+1)   = somax;
                sommeLIMIT(n_film-min_film+1)= tmpLIMIT; % OD > 1.5 aren't well calibrated(), Have a idea the nb of pixies in this case
                energie(n_film-min_film+1)   = ener;
                %doseROI(n_film-min_film+1)   = doseAvg;
                
                %TO DO ! now it is a null !
                fprintf('now play with the dN/dE/dsr');
                energie_fct_angle_str(n_film-min_film+1,:)=ener_fct_angle./(surf_indiv*100/(distance_dbl)^2);
                
                fprintf('\n******************************\n')
                fprintf('Solid Angle      = %1.6f sr\n',str_mat(n_film-min_film+1))
                fprintf('Surface          = %1.6f cm�\n',surf_mat(n_film-min_film+1))
                fprintf('Energy           = %1.6f J\n',energie(n_film-min_film+1))
                %fprintf('Dose (Gy)          = %1.6f\n',doseROI(n_film-min_film+1))
                fprintf('Saturated pixels = %d\n',tmpLIMIT)
                fprintf('******************************\n\n')
                %pause(1)
                
                n_film=n_film+1;
                
                clear n_film1
                save workspace
            end
            
            
            fichier_str(1:max_film-min_film+1,1) = min_film:max_film;%Give no of the film 
            fichier_str(1:max_film-min_film+1,2) = rmax_mat;         %Give the half angle on where the signal has been processed
            fichier_str(1:max_film-min_film+1,3) = divergence_mat;   %Give the half angle of divergence of the beam at the half maximum
            
            cd(dirfichier)
            
            save 'fichier_solid_angle.txt'  fichier_str -ASCII -TABS
            save (strcat('energie_tir',nbr_tir),'energie')
            save (strcat('str-tir',nbr_tir),'str_mat')
            
            %To save when "ener_fct_angle" is implemented
            %save (strcat('energie_fct_angle_str_tir',nbr_tir),'energie_fct_angle_str')
            %save 'indice_VS_divergence.txt'  X1 -ASCII -TABS
            
            save(save_name) %Save the workspace, if the script is stopped
            
            %Save the deposit energy (already saved) but in the data directory 
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
            print('-dtiff','RCF_Response');
            
            cd(dircur)
            save(save_name) % Save in two differents directories

        end
       %%  **************************************************************************************************************** %%
       %%  **********************************END OF THE CALCUL OF THE DEPOSIT ENERGY ON THE RCFs ************************  %%
       %%  **************************************************************************************************************** %%
       

       %%  ***************************************************************************************************************  %%
       %%  ********************* CALCUL OF DEPOSIT ENERGIE IN EACH FILM FOR A GIVEN ENERGIE OF THE PROTON ****************  %%
       %%  ***************************************************************************************************************  %%
       %%   Fit a power law which match the previous measured deposition of energy in the RCF
        
        if strcmp(analyse,'3 Reconstruct Spectre')
            
                 
            %type_spectre='spectre_VS_div';
            %  load dirfichier
            %  cd(dirfichier)
            
            %... Retrieving the data from the associated workspace
            [nomData,cheminData]=uigetfile('*.mat','WORKSPACE in the RCF folder');
            dirfichier=cheminData;
            
            fichier=strcat(cheminData,nomData);
            
            load(fichier);

            cd(dircur);
            

            type_spectre = questdlg('Which type of spectrum?','Analyse','spectre_total', ...
            'spectre_par_str','spectre_moy_pxl','spectre_total');

            clear e_moy delta_E
            
            global in_layer;
            
            %Bloque semble obsolete et sans intérêt
            %{
            if strcmp(type_spectre','spectre_VS_div')
                            nbr_tir=nomData(21:end-4);
                            DN_DE_fct_angle_str=zeros(size(energie_fct_angle_str));
                            
                            size_fct=1;
                            while energie_fct_angle_str(1,size_fct)~=0
                                energie=energie_fct_angle_str(:,size_fct)';
                                calcul_spectre % calcul le spectre en comparant un fit theorique et les points exp.
                                DN_DE_fct_angle_str(1:size(DN_DE1,1),size_fct)= ...
                                    DN_DE1(:,1);
                                clear DN_DE1 e_moy couche1
                                size_fct = size_fct+1;
                            end
            else
            %}

            
            calcul_spectre_test % calcul le spectre en comparant un fit theorique et les points exp.
            
            if strcmp(type_spectre,'spectre_par_str')
                energie=energie_str';
                %nbr_tir=nomData(16:end-4);
                graph_name = strcat('E0=',E0,', E1=',E1,', Emax=',energie_max,'MeV, erreur='...
                    ,erreur);
                nbr_tir=strcat(nbr_tir,'par_str');
                ylabel('dN/dEd\Omega (part./MeV/sr)');
                xlabel('Energy (MeV)');
                axis([0 60  1e6 1e11])
                set(gca, 'XScale', 'linear')
                set(gca, 'YScale', 'log')
                cd(dirfichier)
                save DN_DE_str.txt DN_DE1 -ASCII -DOUBLE -TABS;
                cd(dircur)
            
            elseif strcmp(type_spectre,'spectre_moy_pxl')
                energie=energie_moy_pxl;
                %nbr_tir=nomData(16:end-4);
                graph_name = strcat('E0=',E0,', E1=',E1,', Emax=',energie_max,'MeV, erreur='...
                    ,erreur);
                nbr_tir=strcat(nbr_tir,'moy_par_str');
                ylabel('nbr de proton moyen par str');
                cd(dirfichier)
                save DN_DE_moy_str.txt DN_DE1 -ASCII -DOUBLE -TABS;
                cd(dircur)
            else
                %nbr_tir=nomData(12:end-4);
                Edepos_J=num2str(Edepos_J,'%10.4f');
                graph_name = strcat('E0=',E0,', E1=',E1,', Emax=',energie_max,'MeV, erreur='...
                    ,erreur,', Energie-dep=',Edepos_J,'J' );
                ylabel('dN/dE (part./MeV)');
                %label=strcat('energie (Mev)                       Tir',nbr_tir);
                xlabel('Energy (MeV)');
                axis([0 60  1e6 1e11])
                set(gca, 'XScale', 'linear')
                set(gca, 'YScale', 'log')
                cd(dirfichier)
                save DN_DE.txt DN_DE1 -ASCII -DOUBLE -TABS;
                cd(dircur)
                
            end
            
            title({strcat('Proton spectrum - Shot-',nbr_tir),graph_name})
            cd(dirfichier)
            saveas (h1,strcat('spectre-tir',nbr_tir),'fig');
            save_string=strcat('spectre-tir',nbr_tir,'.tif');
            print('-dtiff', save_string);
            cd(dircur)
            
        end
        etude = questdlg('Start again ?','start again','Yes','No','Yes');
    end