% test script for Thomson code
% Makes a plot of a Thomson scattering spectrum

lmd0 = 526.5;
ne= 1e26; % electron density m^{-3}
Te=150; % electron temperature eV
Ti=100; % ion temperature eV
Z=1; % average ionization 
A=1; % atomic mass 
vI=6e4; % ion velocity [m s^-1]
J=-6e7;% plasma current [A cm^-3];
fract=1; % fractional contribution to plasma (multi spieces plasma)
lambda0 = lmd0; % Probe wavelength
lambdaRange = (lmd0-2.5:0.01:lmd0+2.5); % wavelengths to calculate the cross section at
scatAngle = 90; % Scattering angle
order =1; % Order of calcualtion (1st or 2nd order relativistic correction)
cgsUnits = 0; % CGS unit flag (0 defaults to SI)
nmLambda =1; % nm wavelength units flag

%% Calculate the underlying spectrum
 [scatteringCrossSection,formFactor,totalScatteringRatio]  ...
    = thomsonCrossSection(lambdaRange,lambda0,ne,Te,fract,Ti,A,Z,vI,J,scatAngle,order,cgsUnits,nmLambda);

% need to compare alpha with (Z*Te/3./Ti - 1.)^(-0.5)
% to check if collective or non-collective
% alpha = 1 / (k*lmd_D) [J. Ross, et al., RSI 81, 10D523 (2010)]
% sprintf('%.2f',alpha)

figure('Color','white');
nmLambda =1;
plot(lambdaRange,scatteringCrossSection/max(scatteringCrossSection),'-r');
xlabel('Wavelength [nm]');
ylabel('Cross Section [m^-1 nm^-1 sr^-1]')
set(gca,'YScale','linear');



%% Apply Arbitary Spectrometer Broadening
instrumentWidth = 0.04; % gaussian stanrdard deviation of instrument fn in nm.
scatteringCrossSection =gaussianBroadening(lambdaRange,scatteringCrossSection,instrumentWidth);

% figure('Color','white');

plot(lambdaRange,scatteringCrossSection/max(scatteringCrossSection), ...
    '-b',LineStyle='--',LineWidth=2);
xlabel('Wavelength [nm]');
% ylabel('Cross Section [m^-1 nm^-1 sr^-1]');
ylabel('A.U.');
set(gca,'YScale','linear');
