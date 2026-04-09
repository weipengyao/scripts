% script to 
% 1. load the TS lineout data
% 2. convert the wavelength axis into nm centered at lmd0
% 3. normalize the signal axis


opts = detectImportOptions('/Users/yao/Nextcloud/PROJECTS/LULI2000/cross_talk_Nov_2023/TS_data_lineouts/shot80.txt');
M = readmatrix('/Users/yao/Nextcloud/PROJECTS/LULI2000/cross_talk_Nov_2023/TS_data_lineouts/shot80.txt',opts);

res = 0.0037; % nm/px
lmd0 = 526.5; % probe laser wavelength

noise = (M(1,2)+M(end,2))/2.;

secondary_peak = 0.336;

plot(M(:,1)*res-max(M(:,1)*res)/2.+lmd0, ...
    (M(:,2)-noise)/max(M(:,2)-noise)/secondary_peak, ...
    'sk',LineWidth=1);
xlabel('\lambda0 (nm)');
ylabel('A. U.');
xlim([lmd0-2.5,lmd0+2.5]);
ylim([0,1]);
set(gca,'FontSize',20);

