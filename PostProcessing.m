close all
clear all
clc

addpath(genpath('./Models'))
addpath(genpath('./Shapes'))


%Define file name for inputs, and file name for refernce displacements
fileZero = "Reference150.mat";
fileName = "Results1950.mat";

%scale for deformations
defScale = 1000;
iceHeight = 300;

%loading from file
load(fileZero, "mesh","physics","solver", "TimeSeries");
physicsOffset = physics;
load(fileName, "mesh","physics","solver", "TimeSeries");


%% plot pressures and deformations
f = figure(1);
	tiledlayout('flow')
	nexttile
	plotDeformationsPressure(physics, physicsOffset, defScale, false);
	axis equal tight

	%clim([-0.5 0.5])
	ax = gca; 
	axis off

    cb = colorbar;
    cb.Layout.Tile = 'east'; 
    cb.Title.String = {"$p-\rho_{\mathrm{i}}gy$", "[$\mathrm{MPa}$]"};
    cb.Title.Interpreter='latex';
    cb.FontSize = 8;

    saveFigNow(f, "Figures/PSurf", 3, false)

 %% plotting displacement over time
 f = figure(2);

	[~,iStart]=find(TimeSeries.tvec<0); iStart=max(iStart); if(isempty(iStart)) iStart=1; end
	[~,iBasal]=find(TimeSeries.Lfrac>iceHeight); iBasal=min(iBasal);

	lcs = [360, 257, 205, 154, 102];
	xLoc = [0, 500, 1000, 1500, 2000];

	cstyles = { "k", "c", "m", "b","g","y"};  %

	cc = turbo(15);

	m = ["x","o","d","s","*"]; m = [m m];
	for l=1:length(lcs)
		hp = plot(TimeSeries.tvec(iStart:end)/60, TimeSeries.SurfaceDisp(iStart:end,lcs(l),2)-TimeSeries.SurfaceDisp(iStart,1,2),'Color',cc(l,:),'LineWidth',1,'DisplayName',"$x="+string(xLoc(l))+"\;\mathrm{m}$");
		hold on
		plotsparsemarkers(hp, [], m(l), 6, false);
	end

	[tVec,yVec] = GetRef("uy",0);
	plot(tVec, yVec-yVec(1), "r-.",'LineWidth',1,'DisplayName',"Observations (Das et al.)")

	xlim([0 max(TimeSeries.tvec)/60]);
	ylim([0 max(TimeSeries.SurfaceDisp(iStart:end,lcs(1),2)-TimeSeries.SurfaceDisp(iStart,1,2))])
	plot([TimeSeries.tvec(iBasal) TimeSeries.tvec(iBasal)]/60,[-0.1 1],'k:','HandleVisibility','off')

	legend('Location','southoutside','Interpreter','latex','NumColumns',2)
	xlabel('time $[\mathrm{min}]$','Interpreter','latex')
	ylabel('Surface uplift $[\mathrm{m}]$','Interpreter','latex')

	saveFigNow(f, "Figures/Uplift", 6, false)

%% thermal fluxes
f = figure(3);

	QScale = 1/(910*335000);
	hp = plot(TimeSeries.tvec(iStart:end)/60, QScale*TimeSeries.Qvec(iStart:end,1),'k','LineWidth',1,'DisplayName',"Conduction");
	hold on
	plotsparsemarkers(hp, [], "x", 6, false);

	hp = plot(TimeSeries.tvec(iStart:end)/60, QScale*TimeSeries.Qvec(iStart:end,2),'k','LineWidth',1,'DisplayName',"Friction");
	hold on
	plotsparsemarkers(hp, [], "d", 6, false);

	hp = plot(TimeSeries.tvec(iStart:end)/60, QScale*TimeSeries.Qvec(iStart:end,3),'k','LineWidth',1,'DisplayName',"Freezing");
	hold on
	plotsparsemarkers(hp, [], "*", 6, false);

	xlim([0 max(TimeSeries.tvec)/60]);
	plotylim = max(max(QScale*abs(TimeSeries.Qvec(iStart:end,:))));
	ylim([-plotylim plotylim])

	plot([TimeSeries.tvec(iBasal) TimeSeries.tvec(iBasal)]/60,[-6 6],'k:','HandleVisibility','off')

	legend('Location','southoutside','Interpreter','latex','NumColumns',2)
	xlabel('time $[\mathrm{min}]$','Interpreter','latex')
	ylabel('Thermal energy/$\rho_i \mathcal{L} \;[\mathrm{m}^3/\mathrm{m}]$','Interpreter','latex')

	saveFigNow(f, "Figures/Thermals", 8, false)

%% fluid fluxes
f = figure(4);
	tiledlayout('flow')
	LFrac = 3.2e3;
	ALake = 5.6e6; 
	m3_msec_to_m_hour = LFrac/ALake*3600;
	m3_to_mWaterHeigt = LFrac/ALake;

	nexttile
	yyaxis left
	plot(TimeSeries.tvec(iStart:end)/60, TimeSeries.Qinflow(iStart:end)*m3_to_mWaterHeigt,"-k",'LineWidth',1,'DisplayName',"Simulation");
	hold on
	[tVec,yVec] = GetRef("Q",0);
	plot(tVec, yVec,"-ro",'LineWidth',1,'DisplayName',"Reference");
	hold on

	xlim([0 max(TimeSeries.tvec)/60]);
	xlabel({'time $[\mathrm{min}]$','(a)'},'Interpreter','latex')
	ylabel({'Water level change [$\mathrm{m}$]'},'Interpreter','latex')	
	%legend('Location','southoutside','Interpreter','latex','NumColumns',2)

	plot([TimeSeries.tvec(iBasal) TimeSeries.tvec(iBasal)]/60,[0 2],'k:','HandleVisibility','off')

	ax = gca;
	ax.YColor = 'k';
	ylim([0 max(TimeSeries.Qinflow(iStart:end))*m3_to_mWaterHeigt])

	yyaxis right
	ylabel({'Total inflow [$\mathrm{m}^3/\mathrm{m}$]'},'Interpreter','latex')	
	ax = gca;
	ax.YColor = 'k';
	ylim([0 max(TimeSeries.Qinflow(iStart:end))*m3_to_mWaterHeigt]/m3_to_mWaterHeigt)

	ax = gca;
	ax.FontSize = 8;

	nexttile
	yyaxis left
	plot(TimeSeries.tvec(iStart:end)/60, TimeSeries.qCurrent(iStart:end)*m3_msec_to_m_hour,"k-",'LineWidth',1,'DisplayName',"Simulation");
	hold on
	[tVec,yVec] = GetRef("dQ",0);

    plot(tVec(1:2), yVec(1:2),"r-",'LineWidth',1,'DisplayName',"Observations");
    for iii=3:2:length(tVec)-2
	    plot(tVec(iii:iii+1), yVec(iii:iii+1),"r-",'LineWidth',1,'DisplayName',"Reference",'HandleVisibility','off');
        plot(tVec(iii+1:iii+2), yVec(iii+1:iii+2),"r:",'LineWidth',1,'DisplayName',"Reference",'HandleVisibility','off');
    end
	plot(tVec, yVec,"ro",'LineWidth',1,'DisplayName',"Reference",'HandleVisibility','off');

	xlim([0 max(TimeSeries.tvec)/60]);
	xlabel({'time $[\mathrm{min}]$','(b)'},'Interpreter','latex')
	ylabel({'Change rate [$\mathrm{m}/\mathrm{hour}$]'},'Interpreter','latex')	
	l = legend('Interpreter','latex','NumColumns',3);

	plot([TimeSeries.tvec(iBasal) TimeSeries.tvec(iBasal)]/60,[0 2],'k:','HandleVisibility','off')

	ax = gca;
	ax.YColor = 'k';
	ylim([0 max(TimeSeries.qCurrent(iStart:end)*m3_msec_to_m_hour)])

	yyaxis right
	ylabel({'Inflow rate [$\mathrm{m}^3/\mathrm{m}/\mathrm{hour}$]'},'Interpreter','latex')	
	ax = gca;
	ax.YColor = 'k';
	ylim([0 max(TimeSeries.qCurrent(iStart:end)*m3_msec_to_m_hour)]/m3_to_mWaterHeigt)
	l.Layout.Tile = "south";

	saveFigNow(f, "Figures/FluidFluxes", 9, false)

%% opening vs melting
f = figure(5);
	semilogy([-1],[1],'w','DisplayName','\textbf{Deformation:}');
	hold on
	hp = semilogy(TimeSeries.tvec(iStart:end)/60, -2*TimeSeries.upLift(iStart:end,1),'k-','LineWidth',1,'DisplayName',"Results");
	plotsparsemarkers(hp, [], "d", 6, false);

	semilogy([-1],[1],'w','DisplayName','\textbf{Melting:}');
	hp = semilogy(TimeSeries.tvec(iStart:end)/60, physics.models{6}.hMeltHist(iStart:end),'k-.','LineWidth',1,'DisplayName',"Results");
	plotsparsemarkers(hp, [], "d", 6, false);

	plot([TimeSeries.tvec(iBasal) TimeSeries.tvec(iBasal)]/60,[1e-3 1],'k:','HandleVisibility','off')

	xlim([0 max(TimeSeries.tvec)/60]);
	ylim([1e-4 1])
	xlabel('time $[\mathrm{min}]$','Interpreter','latex')
	ylabel('Opening width [$\mathrm{m}$]','Interpreter','latex')	
	legend('Location','southoutside','Interpreter','latex','NumColumns',2)

	saveFigNow(f, "Figures/OpeningHeights", 6, false)

%% uplift region
f = figure(6);
	LBase = TimeSeries.Lfrac-iceHeight; LBase(LBase<=0)=nan;
	LCrev = TimeSeries.Lfrac-iceHeight; LCrev(LCrev>0)=nan; LCrev=LCrev+iceHeight;
	hp = plot(TimeSeries.tvec(iStart:end)/60, LBase(iStart:end)*1e-3,'k-','LineWidth',1,'DisplayName',"Viscoplastic");
	hold on
	plotsparsemarkers(hp, [], "d", 6, false);
	hp = plot(TimeSeries.tvec(iStart:end)/60, LCrev(iStart:end)*1e-3,'k-.','LineWidth',1,'DisplayName',"Viscoplastic",'HandleVisibility','off');
	plotsparsemarkers(hp, [], "d", 6, false);

	plot([TimeSeries.tvec(iBasal) TimeSeries.tvec(iBasal)]/60,[0 4],'k:','HandleVisibility','off')

	xlim([0 max(TimeSeries.tvec)/60]);
	ylim([0 max(max(LBase), max(LCrev))]*1e-3);
	xlabel('time $[\mathrm{min}]$','Interpreter','latex')
	ylabel({'Crevasse depth [$\mathrm{km}$]';'Hor. crack length [$\mathrm{km}$]'},'Interpreter','latex')	
	legend('Location','southoutside','Interpreter','latex','NumColumns',2)

	saveFigNow(f, "Figures/FractureLengths", 6, false)


function saveFigNow(fg, sname, HFig, SaveDouble)
	fprintf(sname+"  ")
	ax = gca;

	fg.Units = 'centimeters';
	if (SaveDouble)
		fg.Position = [2 2 16 HFig];
	else
		fg.Position = [2 2 8 HFig];
	end
	ax.FontSize = 8;

	drawnow();
	print(fg, sname+".png",'-dpng','-r1200'); fprintf(".png  ")
	print(fg, sname+".jpg",'-djpeg','-r1200'); fprintf(".jpg  ")
	%print(fg, sname+".tif",'-dtiffn','-r1200'); fprintf(".tif  ")
	print(fg, sname+".svg",'-dsvg','-r1200'); fprintf(".svg  ")
	print(fg, sname+".eps",'-depsc','-r1200'); fprintf(".eps  ")
	print(fg, sname+".emf",'-dmeta','-r1200'); fprintf(".emf\n")
end

function [tVec,yVec] = GetRef(whichToPlot,T0)
	if (whichToPlot=="uy")
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

		tVec =  Das_HData(:,1)*60-T0;
		yVec =	Das_HData(:,2);
	end

	if (whichToPlot=="Q" || whichToPlot=="dQ")
		Das_QData = [0.028143727632721438, -0.1951255992949945
				0.332136882451348, -0.2014419462029151
				0.671285708235394, -0.2267073338345984
				0.998026650149292, -0.27723810909796454
				1.3268355727082146, -0.6246371890336087
				1.6659843984922604, -1.951070039696977
				1.8272868888041849, -3.9217702749682672]; 
		Das_QData(:,2) = Das_QData(:,2)-Das_QData(1,2);

		tVec =  Das_QData(:,1)*60-T0;
		yVec =	-(Das_QData(:,2)-Das_QData(1,2));

		if (whichToPlot=="dQ")
			for i=1:length(Das_QData)-1
				Das_dQData(2*i-1,1)=Das_QData(i,1);
				Das_dQData(2*i  ,1)=Das_QData(i+1,1);
				Das_dQData(2*i-1,2)=(Das_QData(i,2)-Das_QData(i+1,2))/(Das_QData(i+1,1)-Das_QData(i,1));
				Das_dQData(2*i  ,2)=(Das_QData(i,2)-Das_QData(i+1,2))/(Das_QData(i+1,1)-Das_QData(i,1));
			end

			tVec =  Das_dQData(:,1)*60-T0;
			yVec =	Das_dQData(:,2);
		end
	end


end

function plotDeformationsPressure(physics, physicsOffset, sc, fillICE)
	HIce = max(physics.mesh.Nodes(:,2));
	% Boundaries
	Egroups = [2 3 4 5 6 7 8];

    ipc = physics.mesh.ipcount1D;
    Svec = physics.StateVec;
	OffsetVec = physicsOffset.StateVec;

	zoomLim = -2000;
	for e=1:length(Egroups)
        for n_el=1:size(physics.mesh.Elementgroups{Egroups(e)}.Elems, 1)
            Elem_Nodes   = physics.mesh.getNodes(Egroups(e), n_el);
            [N, ~, ~]    = physics.mesh.getVals(Egroups(e), n_el);
            
            dofsX = physics.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(1), Elem_Nodes);
            dofsY = physics.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(2), Elem_Nodes);

            X = Svec(dofsX);
            Y = Svec(dofsY);

	        dofsX = physicsOffset.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(1), Elem_Nodes);
            dofsY = physicsOffset.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(2), Elem_Nodes);

            XOffset = OffsetVec(dofsX);
            YOffset = OffsetVec(dofsY);		

            xy = physics.mesh.getIPCoords(Egroups(e), n_el);
            for ip=1:ipc
                coords(ip,1) = xy(1, ip)+N(ip,:)*(X(1:3)-XOffset(1:3))*sc;
                coords(ip,2) = xy(2, ip)+N(ip,:)*(Y(1:3)-YOffset(1:3))*sc;

				if (e==5)
					zoomLim = max(zoomLim, coords(ip,1));
				end
            end
            lines(n_el,:,:) = coords;
        end
        hold on
        plot(squeeze(lines(:,:,1))', squeeze(lines(:,:,2))', 'k');
	end

	% Interior Colouring
	if (fillICE)
		Egroup = 1;
		XAreasIce = []; XAreasRock = [];
		YAreasIce = []; YAreasRock = []; 
		C = [];
		for n_el=1:size(physics.mesh.Elementgroups{Egroup}.Elems, 1)
        	Elem_Nodes   = physics.mesh.getNodes(Egroup, n_el);
        	Elem_NodesOffet   = physicsOffset.mesh.getNodes(Egroup, n_el);
	
        	dofsX = physics.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(1), Elem_Nodes);
        	dofsY = physics.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(2), Elem_Nodes);
	
			dofsXOffset = physicsOffset.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(1), Elem_NodesOffet);
        	dofsYOffset = physicsOffset.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(2), Elem_NodesOffet);
	
			XN = physics.mesh.Nodes(Elem_Nodes,1);
			YN = physics.mesh.Nodes(Elem_Nodes,2);
	
        	X = sc*(Svec(dofsX)-OffsetVec(dofsXOffset)) + XN;
        	Y = sc*(Svec(dofsY)-OffsetVec(dofsYOffset)) + YN;
	
			if (mean(YN)<0)
				XAreasRock(end+1,:) = X([1 3 9 7 1]);
				YAreasRock(end+1,:) = Y([1 3 9 7 1]);
			else
				XAreasIce(end+1,:) = X([1 3 9 7 1]);
				YAreasIce(end+1,:) = Y([1 3 9 7 1]);
			end
		end
		fill(XAreasIce', YAreasIce', [1 1 1]*0.9,'EdgeColor','none');
		fill(XAreasRock', YAreasRock', [1 1 1]*0.4,'EdgeColor','none');
	end

	% pressures & crack
	Egroup = 9;

	xlineLeft = [];
	ylineLeft = [];
	xlineRight = [];
	ylineRight = [];
	xlineCentre = [];
	ylineCentre = [];
	plineCentre = [];

	for n_el=1:size(physics.mesh.Elementgroups{Egroup}.Elems, 1)
	    Elem_Nodes   = physics.mesh.getNodes(Egroup, n_el);
		Elem_NodesOffset = [min(Elem_Nodes([1 4])) min(Elem_Nodes([2 5])) min(Elem_Nodes([3 6]))];

    	dofsX = physics.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(1), Elem_Nodes);
    	dofsY = physics.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(2), Elem_Nodes);
		dofsP = physics.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(3), Elem_Nodes(1:3));

		try
		dofsXOffset = physicsOffset.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(1), Elem_NodesOffset);
    	dofsYOffset = physicsOffset.dofSpace.getDofIndices(physics.models{5}.dofTypeIndices(2), Elem_NodesOffset);

		dofsXOffset = [dofsXOffset;dofsXOffset];
		dofsYOffset = [dofsYOffset;dofsYOffset];

		XN = physics.mesh.Nodes(Elem_Nodes,1);
		YN = physics.mesh.Nodes(Elem_Nodes,2);

    	X = sc*(Svec(dofsX)-OffsetVec(dofsXOffset)) + XN;
    	Y = sc*(Svec(dofsY)-OffsetVec(dofsYOffset)) + YN;
		P = Svec(dofsP);

		xlineLeft(end+1,:) = X(1:3);
		ylineLeft(end+1,:) = Y(1:3);
		xlineRight(end+1,:) = X(4:6);
		ylineRight(end+1,:) = Y(4:6);
		xlineCentre(end+1,:) = [0.5*(X(1:3)+X(4:6));NaN];
		ylineCentre(end+1,:) = [0.5*(Y(1:3)+Y(4:6));NaN];
		plineCentre(end+1,:) = 1e-6*[P - 910*9.81*(HIce-YN(1:3));NaN];
		catch ME

		end

	end
    plot(xlineLeft', ylineLeft', 'k');
	plot(xlineRight', ylineRight', 'k');

    patch(xlineCentre',ylineCentre',plineCentre','EdgeColor','interp','FaceColor','None','Marker','None','MarkerFaceColor','flat','LineWidth',2);
end



