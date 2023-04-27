%close all
%clear all
%clc

H_Glacier = 980;

T = linspace(-10,0,100)+273.15;
Q = 150e3;
R = 8.31446261815324;
Tref = 273.15;
A0 = 5e-24;
V = -1.74e-5;

T_to_A = @(T,p) A0*exp(-(Q+p*V)./R * (1/T-1/Tref));

for i=1:length(T)
	A_A0(i) = T_to_A(T(i),0);
end

figure
plot(T-273.15, A_A0)

TProfileData = importdata('foxxTemperatureProfile.csv');
TProfile.X = TProfileData.data(:,2);
TProfile.T = TProfileData.data(:,1);

%transforming to ground-up coordinates
HData = -min(TProfile.X);
TProfile.X = (TProfile.X+HData)*H_Glacier/HData;
TProfile.X = [-500;TProfile.X]; %adding rock T for interpolating purposes
TProfile.T = [TProfile.T(1);TProfile.T];

T_Ice = griddedInterpolant(TProfile.X, TProfile.T,'linear');


YPlot = linspace(0,H_Glacier,1000);
TPlot = T_Ice(YPlot);
for i=1:length(TPlot)
	A_Plot(i) = T_to_A(TPlot(i)+273.15, 910*9.81*(980-YPlot(i)));
end

figure
subplot(1,2,1)
plot(TPlot, YPlot,'k-')
hold on
plot(TProfile.T(2:end), TProfile.X(2:end),'r*')
xlabel('T [degC]')
ylabel('y [m]')

subplot(1,2,2)
plot(A_Plot, YPlot,'k-')
xlabel('A [Pa^{-3} s^{-1}]')
ylabel('y [m]')

save("TProfile.mat","T_Ice")

