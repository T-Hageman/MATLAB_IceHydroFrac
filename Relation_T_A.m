close all
clear all
clc

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

% figure
% plot(T-273.15, A_A0)

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
	ft_Plot(i) = 2.0-0.0068*(TPlot(i)+273.15);   %https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2012JE004101
end

f = figure
tiledlayout(1,3,"Padding","tight","TileSpacing","tight")
nexttile
plot(TPlot, YPlot,'k-','LineWidth',1)
%hold on
%plot(TProfile.T(2:end), TProfile.X(2:end),'r*')
xlabel('$T\;[^\circ \mathrm{C}]$','Interpreter','latex')
ylabel('$y\;[\mathrm{m}]$','Interpreter','latex')
ylim([0 980])
ax = gca;
ax.FontSize = 8;

nexttile
plot(A_Plot*1e24, YPlot,'k-','LineWidth',1)
xlabel('$A\;[10^{-24}\;\mathrm{Pa}^{-3} \mathrm{s}^{-1}]$','Interpreter','latex')
ylim([0 980])
ax = gca;
ax.FontSize = 8;
ax.YTick = []

nexttile
plot(ft_Plot, YPlot,'k-','LineWidth',1)
xlabel('$f_t\;[\mathrm{MPa}]$','Interpreter','latex')
ylim([0 980])
ax = gca;
ax.FontSize = 8;
ax.YTick = []

%saveFigNow(f,"Figures/TempProfile",6,false,false)

save("TProfile.mat","T_Ice")

function saveFigNow(fg, sname, HFig, SaveDouble, hasColorbar, cb)
	fprintf(sname+"  ")

	fg.Units = 'centimeters';
	if (SaveDouble)
		fg.Position = [2 2 17.8 HFig];
	else
		fg.Position = [2 2 8.7 HFig];
	end


	drawnow();
	print(fg, sname+".png",'-dpng','-r1200'); fprintf(".png  ")
	print(fg, sname+".jpg",'-djpeg','-r1200'); fprintf(".jpg  ")
	print(fg, sname+".eps",'-depsc','-r1200'); fprintf(".eps  ")
	print(fg, sname+".emf",'-dmeta','-r1200'); fprintf(".emf\n")
end
