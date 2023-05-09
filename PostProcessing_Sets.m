close all
clear all
clc

addpath(genpath('./Models'))
addpath(genpath('./Shapes'))

global maxplotStep defscale saveIndividualP saveIndividualS;
maxplotStep = 3600;
defscale = 1000;
saveIndividualP = false;
saveIndividualS = false;

for set=[4]
	switch set
		case 1 %Viscous
			SetName = "Das_Visc"

			ind = 0;
			cstyles = { "k", "r", "c", "b", "m", ....
						     "r-.", "c-.", "b-.", "m-.", ....
						     "r:", "c:", "b:", "m:"};
			for Ti=1:8
				T = 1-Ti;
				ind = ind+1;
				fnames{ind} = "Results_Visc/"+string(T);
				snames{ind} = "T"+string(T);
				lnames{ind} = "T="+string(T);
				lstyles{ind} = cstyles{Ti};
			end
		case 2 %Lin-Elastic
			SetName = "Das_LinEl"
			cstyles = { "k", "r", "c", "b", "m", ....
			     "r-.", "c-.", "b-.", "m-.", ....
			     "r:", "c:", "b:", "m:"};

			ind = 0;

			for Ti=1:11
				T = 1-Ti;
				ind = ind+1;
				fnames{ind} = "Results_NoVisc/"+string(T);
				snames{ind} = "T"+string(T);
				lnames{ind} = "T="+string(T);
				lstyles{ind} = cstyles{Ti};
			end
		case 3 %Composite
			SetName = "Das_Composite"
			cstyles = { "k", "r", "c", "b", "m", ....
						     "r-.", "c-.", "b-.", "m-.", ....
						     "r:", "c:", "b:", "m:"};
		
			fnames{2} = "Results_NoVisc/0";
			snames{2} = "NoVisc";
			lnames{2} = "Linear-Elastic";
			lstyles{2} = cstyles{2};

			fnames{1} = "Results_Visc/0";
			snames{1} = "Visc";
			lnames{1} = "Visco-Elastic";
			lstyles{1} = cstyles{1};

		case 4 %Single
			SetName = "Das_Single"
			cstyles = { "k", "r", "c", "b", "m", ....
						     "r-.", "c-.", "b-.", "m-.", ....
						     "r:", "c:", "b:", "m:"};

			fnames{1} = "Results/Visc_FrictionFactor";
			snames{1} = "Visc";
			lnames{1} = "Visco-Elastic";
			lstyles{1} = cstyles{1};
	end

	ProcessingForSet(fnames, snames, lnames, lstyles, SetName)
end



function ProcessingForSet(fnames, snames, lnames, lstyles, setName)
	T0 = 20;

	Das_QData = [0.028143727632721438, -0.1951255992949945
				0.332136882451348, -0.2014419462029151
				0.671285708235394, -0.2267073338345984
				0.998026650149292, -0.27723810909796454
				1.3268355727082146, -0.6246371890336087
				1.6659843984922604, -1.951070039696977
				1.8272868888041849, -3.9217702749682672]; 
	Das_QData(:,1) = Das_QData(:,1)*60-T0; Das_QData(:,2) = Das_QData(:,2)-Das_QData(1,2);
	for i=1:length(Das_QData)-1
		Das_dQData(2*i-1,1)=Das_QData(i,1);
		Das_dQData(2*i  ,1)=Das_QData(i+1,1);
		Das_dQData(2*i-1,2)=(Das_QData(i,2)-Das_QData(i+1,2))/(Das_QData(i+1,1)-Das_QData(i,1))*60;
		Das_dQData(2*i  ,2)=(Das_QData(i,2)-Das_QData(i+1,2))/(Das_QData(i+1,1)-Das_QData(i,1))*60;
	end

	Das_HData = [	0.0010022909507449995, -0.19111111111111123
					0.17998281786941628, -0.20351190476190517
					0.33748568155784686, -0.19359126984127029
					0.5093069873997715, -0.17623015873015913
					0.7026059564719362, -0.1539087301587303
					0.881586483390608, -0.12910714285714286
					1.0713058419243993, -0.10182539682539682
					1.2502863688430703, -0.07454365079365122
					1.4364261168384886, -0.024940476190476346
					1.5796105383734256, 0.044503968253968296
					1.6869988545246284, 0.1313095238095232
					1.7621706758304705, 0.20819444444444413
					1.8266036655211915, 0.2801190476190474
					1.8981958762886606, 0.36692460317460296
					1.9697880870561288, 0.4834920634920632
					2.0413802978235975, 0.5802182539682539
					2.098654066437572, 0.6769444444444443
					2.173825887743414, 0.7786309523809523
					2.248997709049256, 0.8703968253968253
					2.3993413516609396, 0.9497619047619046
					2.5747422680412377, 0.974563492063492
					2.721506300114548, 0.9596825396825397
					2.864690721649486, 0.9249603174603174
					2.9971363115693017, 0.8852777777777776
					3.1295819014891184, 0.8480753968253967
					3.233390607101948, 0.8034325396825395
					3.3228808705612836, 0.7563095238095237
					3.4517468499427273, 0.7042261904761904
					3.5698739977090503, 0.6719841269841269
					3.698739977090493, 0.642222222222222
					3.856242840778924, 0.6000595238095237
					4.013745704467354, 0.5554166666666664
					4.146191294387172, 0.5107738095238092
					4.296534936998857, 0.4611706349206348
					4.4790950744559, 0.40908730158730133
					4.640177548682704, 0.38676587301587295
					4.797680412371134, 0.37684523809523807
					4.958762886597938, 0.36196428571428574];
	Das_HData(:,1) = Das_HData(:,1)*60-T0; Das_HData(:,2) = Das_HData(:,2)-Das_HData(1,2);


	global maxplotStep defscale saveIndividualP saveIndividualS;

		LFrac = 3.2e3;
		ALake = 5.6e6; 
		m3_msec_to_m_hour = LFrac/ALake*3600;
		m3_to_mWaterHeigt = LFrac/ALake;

	f1 = figure();
	f2 = figure();
	f3 = figure();
	f4 = figure();

	f5 = figure();
	f6 = figure();
	f7 = figure();
	f8 = figure();



	mkdir("./Figures/"+setName)
	for f=1:length(fnames)
		folder = fnames{f}+"/";
		Files=dir(folder);
		nmax = 10*(length(Files)-2)-10;
		plotHalf = false;
		if nmax>maxplotStep
			nmax = maxplotStep;
			%plotHalf = true;
		end
		fname = folder+string(nmax)+".mat";
		load(fname, "mesh","physics","solver","dt","n_max", "TimeSeries");

		figure(f1)
			plot(TimeSeries.tvec/60, TimeSeries.Lfrac,lstyles{f},'LineWidth',1,'DisplayName',lnames{f});
			hold on

		figure(f2)
			semilogy(TimeSeries.tvec/60, physics.models{8}.QTotal,lstyles{f},'LineWidth',1,'DisplayName',lnames{f});
			hold on

		figure(f3)
			plot(TimeSeries.tvec/60, TimeSeries.Qvec(:,1)/335000/910,lstyles{f},'LineWidth',1,'DisplayName',lnames{f});
			hold on

		figure(f4)
			plot(TimeSeries.tvec/60, TimeSeries.Qvec(:,2)/335000/910,lstyles{f},'LineWidth',1,'DisplayName',lnames{f});
			hold on

			if (setName=="Das_Single")
				figure(f5)
					yyaxis left
					plot(TimeSeries.tvec/60, TimeSeries.Qinflow*m3_to_mWaterHeigt,'k-','LineWidth',1,'DisplayName',"$Q$");
					hold on
		
					plot(TimeSeries.tvec/60, TimeSeries.qCurrent*m3_msec_to_m_hour,'k-.','LineWidth',1,'DisplayName',"$\dot{Q}$");
					hold on

					plot(Das_QData(:,1),-Das_QData(:,2),'r-','LineWidth',1,'DisplayName',"$Q_{Das}$");
					plot(Das_dQData(:,1), Das_dQData(:,2),'r-.','LineWidth',1,'DisplayName',"$\dot{Q}_{Das}$")

					yyaxis right
					hold on
				figure(f7)
					plot(TimeSeries.tvec/60, -2*TimeSeries.upLift(:,1),'k-','LineWidth',1,'DisplayName',"$h$");
					hold on
					plot(TimeSeries.tvec/60, TimeSeries.upLift(:,2),'k-.','LineWidth',1,'DisplayName',"$u_y$");
					plot(Das_HData(:,1), Das_HData(:,2),'r-.','LineWidth',1,'DisplayName',"$u_{y,das}$");
					
			else
				figure(f5)
					yyaxis left
					plot(TimeSeries.tvec/60, TimeSeries.Qinflow*m3_to_mWaterHeigt,'LineWidth',1,'DisplayName',lnames{f});
					hold on
		
					yyaxis right
					plot(TimeSeries.tvec/60, TimeSeries.qCurrent*m3_msec_to_m_hour,'LineWidth',1,'DisplayName',lnames{f});
					hold on

				figure(f7)
					yyaxis left
					plot(TimeSeries.tvec/60, -2*TimeSeries.upLift(:,1),'LineWidth',1,'DisplayName',lnames{f});
					hold on

					yyaxis right
					plot(TimeSeries.tvec/60, TimeSeries.upLift(:,2),'LineWidth',1,'DisplayName',lnames{f});
					hold on
			end


		kmax = 1;
		if plotHalf
			kmax = 4;
		end
		for k=1:kmax
			suffix = "";
			if (k>1.5)
				fname = folder+string(450*(k-1))+".mat";
				load(fname, "mesh","physics","solver","dt","tvec","n_max","Lfrac","Qvec","MeltqVec");
				suffix = "__"+string(15*(k-1))+"min";
			end
	
			if (saveIndividualP)
				fg = figure;
					physics.PlotNodal("pd", defscale, "Fracture", 1e-6);
    				physics.models{4}.plotTotHeight(physics, defscale);
					zoomLim = plot_boundaries(physics, defscale);
					xlim([-3000 3000]);
					cb = colorbar;
					cb.Title.String = {"$p$", "[$\mathrm{MPa}$]"};
					cb.Title.Interpreter='latex';
					%axis off
					HFig = (max(ylim)-min(ylim))/300;
					saveFigNow(fg, "./Figures/"+setName+"/L_"+snames{f}+suffix,HFig, false, true, cb)
		
					%xlim([-10 10])
					%saveFigNow(fg, "./Figures/"+setName+"/LZ_"+snames{f},HFig, true, true, cb)
		
					close(fg);
			end
			if (saveIndividualS)
				for i=1:6
					fg = figure;
						sstring = PlotStresses(physics, defscale, i);
						view(2);
						cb = colorbar;
						cb.Title.String = {"$\sigma_{"+sstring+"}$", "[$\mathrm{MPa}$]"};
						cb.Title.Interpreter='latex';	
						xlim([-3000 3000]);
						if i==1
							caxis([-0.1 0.1]);
						elseif i==2
							caxis([0 0.5]);
							cb.Ticks = 0:0.1:0.5;
						else
	
						end
						saveFigNow(fg, "./Figures/"+setName+"/s"+sstring+"_"+snames{f}+suffix,HFig, false, true, cb)
		
					close(fg)
				end
			end
		end
	end


	figure(f1)
	xlim([0 120])
	ylim([0 5000])
	legend('Location','southoutside','NumColumns',3,'Interpreter','latex')
	xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
	ylabel('$L_{frac}\;[\mathrm{m}]$','Interpreter','latex')
	saveFigNow(f1, "./Figures/"+setName+"LFrac",6, false, false)

	figure(f2)
	legend('Location','southoutside','NumColumns',3,'Interpreter','latex')
	xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
	ylabel('$Q_{total}\;[\mathrm{m^3}/\mathrm{m}]$','Interpreter','latex')
	xlim([0 120])
	ylim([1e-1 inf])
	ax = gca;
	ax.YTick = logspace(-2,6,9);
	saveFigNow(f2, "./Figures/"+setName+"QInflow",6, false, false)

	figure(f3)
	xlim([0 120])
	legend('Location','southoutside','NumColumns',3,'Interpreter','latex')
	xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
	ylabel('$Q_{ice}/\rho_i\mathcal{L}\;[\mathrm{m}^3]$','Interpreter','latex')
	saveFigNow(f3, "./Figures/"+setName+"QICE",6, false, false)

	figure(f4)
	xlim([0 120])
	legend('Location','southoutside','NumColumns',3,'Interpreter','latex')
	xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
	ylabel('$Q_{flow}/\rho_i\mathcal{L}\;[\mathrm{m}^3]$','Interpreter','latex')
	saveFigNow(f4, "./Figures/"+setName+"QFLOW",6, false, false)

	if (setName=="Das_Single")
		figure(f5)
		yyaxis left
		xlim([0 120])
		ax = gca;
		ax.YAxis(1).Color = [0 0 0];
		ax.YAxis(2).Color = [0 0 0];
		ylims = [0 1.5];
		ylim(ylims);
		legend('Location','southoutside','Interpreter','latex','NumColumns',3)
		xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
		ylabel(["$Q_{flow} [\mathrm{m}_{lakeHeight}]$","$\dot{Q}_{flow} [\mathrm{m}_{lakeHeight}/\mathrm{hour}]$"],'Interpreter','latex')
	
		yyaxis right
		xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
		ylabel(["$Q_{flow} [\mathrm{m}^3/\mathrm{m}]$","$\dot{Q}_{flow} [\mathrm{m}^3/\mathrm{m}\mathrm{hour}]$"],'Interpreter','latex')
		xlim([0 120])
		ylim(ylims/m3_to_mWaterHeigt)

		saveFigNow(f5, "./Figures/"+setName+"QLake",6, false, false)
	
		figure(f7)
		xlim([0 120])
		legend('Location','southoutside','Interpreter','latex','NumColumns',3)
		xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
		ylabel('$h\;[\mathrm{m}]$  $u_y\;[\mathrm{m}]$','Interpreter','latex')
		saveFigNow(f7, "./Figures/"+setName+"uy",6, false, false)
	else
		figure(f5)
		yyaxis left
		xlim([0 120])
		legend('Location','southoutside','NumColumns',2)
		xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
		ylabel('$Q_{flow} [\mathrm{m}_{lakeHeight}]$','Interpreter','latex')
	
		yyaxis right
		xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
		ylabel('$q [\mathrm{m}_{lakeHeight}/\mathrm{hour}]$','Interpreter','latex')
		xlim([0 120])
		saveFigNow(f5, "./Figures/"+setName+"QLake",6, false, false)
	
		figure(f7)
		yyaxis left
		xlim([0 120])
		legend('Location','southoutside','NumColumns',2,'Interpreter','latex')
		xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
		ylabel('$u_x\;[\mathrm{m}]$','Interpreter','latex')
	
		yyaxis right
		xlim([0 120])
		legend('Location','southoutside','NumColumns',2,'Interpreter','latex')
		xlabel('$\mathrm{time}\;[\mathrm{minutes}]$','Interpreter','latex')
		ylabel('$u_y\;[\mathrm{m}]$','Interpreter','latex')
		saveFigNow(f7, "./Figures/"+setName+"uy",6, false, false)
	end
end


function saveFigNow(fg, sname, HFig, SaveDouble, hasColorbar, cb)
	fprintf(sname+"  ")
	ax = gca;

	fg.Units = 'centimeters';
	if (SaveDouble)
		fg.Position = [2 2 16 HFig];
	else
		fg.Position = [2 2 8.6 HFig];
	end
	ax.FontSize = 9;
	if (hasColorbar)
		axP = ax.Position;
		cb.Position(1) = 0.90;
		cb.Position(2) = 0.05;
		cb.Position(3) = 0.03;
		cb.Position(4) = 0.55;
		cb.Ticks = 0:1:10;
		ax.FontSize = 7;
		cb.FontSize = 9;
		ax.Position(1) = 0.07;
		ax.Position(2) = axP(2);
		ax.Position(3) = 0.78;
		cb.Units = "centimeters";
		cb.Position(4) = HFig-cb.Position(3)-0.80;
	end



	drawnow();
	print(fg, sname+".png",'-dpng','-r1200'); fprintf(".png  ")
	print(fg, sname+".jpg",'-djpeg','-r1200'); fprintf(".jpg  ")
	print(fg, sname+".eps",'-depsc','-r1200'); fprintf(".eps  ")
	print(fg, sname+".emf",'-dmeta','-r1200'); fprintf(".emf\n")
end

function zoomLim = plot_boundaries(physics, sc)
	Egroups = [2 3 4 5 7 8];

    ipc = physics.mesh.ipcount1D;
    Svec = physics.StateVec;

	zoomLim = -2000;
	for e=1:length(Egroups)
        for n_el=1:size(physics.mesh.Elementgroups{Egroups(e)}.Elems, 1)
            Elem_Nodes   = physics.mesh.getNodes(Egroups(e), n_el);
            [N, ~, ~]    = physics.mesh.getVals(Egroups(e), n_el);
            
            dofsX = physics.dofSpace.getDofIndices(physics.models{4}.dofTypeIndices(1), Elem_Nodes);
            dofsY = physics.dofSpace.getDofIndices(physics.models{4}.dofTypeIndices(2), Elem_Nodes);

            X = Svec(dofsX);
            Y = Svec(dofsY);

            xy = physics.mesh.getIPCoords(Egroups(e), n_el);
            for ip=1:ipc
                coords(ip,1) = xy(1, ip)+N(ip,:)*X(1:3)*sc;
                coords(ip,2) = xy(2, ip)+N(ip,:)*Y(1:3)*sc;

				if (e==5)
					zoomLim = max(zoomLim, coords(ip,1));
				end
            end
            lines(n_el,:,:) = coords;
        end
        hold on
        plot(squeeze(lines(:,:,1))', squeeze(lines(:,:,2))', 'k');
	end

	zoomLim = abs(zoomLim);
end


function outName = PlotStresses(physics, defscale, itype)
	switch itype
		case 1
			varName = 'sxy';
			outName = "xy";
		case 2
			varName = 'sdev';
			outName = "dev";
		case 3
			varName = 'svol';
			outName = "vol";
		case 4
			varName = 'sxx';
			outName = "xx";
		case 5
			varName = 'syy';
			outName = "yy";
		case 6
			varName = 'szz';
			outName = "zz";
	end

	for e=1:size(physics.mesh.Elementgroups{1}.Elems, 1)
		Elem_Nodes = physics.mesh.getNodes(1, e);
        IPvar = physics.Request_Info(varName, e, "Interior");
        xy = physics.mesh.getIPCoords(1, e);
		[N, ~, ~] = physics.mesh.getVals(1, e);
		dofsX = physics.dofSpace.getDofIndices(1, Elem_Nodes);
        dofsY = physics.dofSpace.getDofIndices(2, Elem_Nodes);
        
		for ip=1:length(IPvar)
			col = mod(ip,5)+1;
			row = floor((ip-1)/5)+1;
        	x(row,col) = xy(1,ip) + defscale*N(ip,:)*physics.StateVec(dofsX);
        	y(row,col) = xy(2,ip) + defscale*N(ip,:)*physics.StateVec(dofsY);
        	z(row,col) = IPvar(ip);
		end
		surface(x,y,z*1e-6,'EdgeColor','none','FaceColor','interp')
		hold on
	end
end

