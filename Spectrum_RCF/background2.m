% Just to set the unexposed Background
% used in main line 128

%---------------------------------------------------------------------------------------------
%outsoucring of this part, in our special case
% h=figure;
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% imagesc(a);
% axis equal
%
% global H
% global j1
% global i1
%


n_filmb=num2str(n_film);
H = warndlg('Click on five points on the RCF, where the film is not exposed (reference of the transmission)',strcat('RCF #',n_filmb));
waitfor(H) %attend la fermeture de la boite
[j1,i1]=(ginput(5));%io=row=y jo=column=x !!


clear aaa;
aaa= zeros(length(j1),1);
for i=1:max(size(j1))
    aaa(i)=a(floor(i1(i)),floor(j1(i)));
end

T0=mean(aaa);     %floor cuts the numbers after the comma
a=double(a);
clear [j1,i1];
a=a./T0;

close;

u=find(a>1);
a(u)=1;
pic=double(max(max(a)));

current_background_type=RCF_only{n_film};
