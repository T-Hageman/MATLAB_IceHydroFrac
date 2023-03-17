close all
clear all
clc

delete(gcp('nocreate'))
parpool('threads')

addpath(genpath('./Models'))
addpath(genpath('./Shapes'))

%% pressure/opening plots
plotType = {"Sxy","Sdev","SVol","PFrac"};
folder = "Results_Visc";

for irun=1:11
		T = 1-irun;
		file = folder+"/"+string(T)+"/"
		leg_title = "T="+string(T);

		defscale = 1000;
		for i=4:4
			vidfile = VideoWriter("Animations/"+folder+"/"+plotType{i}+"_"+leg_title+".mp4",'MPEG-4');
				vidfile.FrameRate = 3;
				open(vidfile);
			try
				Files=dir(file);
				nmax = min(3600, 10*(length(Files)-2-1));
				fg = figure;
				fg.Units = 'centimeters';
				fg.Position = [-30 2 28 16];
				drawnow();
				for j=30:30:nmax
					clf
		
					load(file+string(j)+".mat", "physics","TimeSeries");
					if i==4
							physics.PlotNodal("pd", defscale, "Fracture", 1e-6);
        					physics.models{4}.plotTotHeight(physics, defscale);
							plot_boundaries(physics, defscale)
							xlabel('x [m]')
							ylabel('y [m]')
							title('t='+string(round(TimeSeries.tvec(end)))+' s');
							ax = gca;
							cb=ax.Colorbar;
							cb.Title.String="p [MPa]";
					
						fg.Units = 'centimeters';
						fg.Position = [-30 2 28 16];
						ax.Position = [0.09 0.17 0.8 0.75];
						cb.Position(1) = 0.90;
						xlim([-3000 3000]);
						ylim([-200 1000]);
						drawnow()
						disp("irun="+string(irun)+", i="+string(i)+", Time= "+string(TimeSeries.tvec(end)))
					else
						sstring = PlotStresses(physics, defscale, i);
						xlabel('x [m]')
						ylabel('y [m]')
						title('t='+string(round(tvec(end)))+' s');
						view(2);
						cb = colorbar;
						cb.Title.String = {"$\sigma_{"+sstring+"}$", "[$MPa$]"};
						cb.Title.Interpreter='latex';	
						xlim([-2000 2000]);
						if i==1
							clim([-0.1 0.1]);
						elseif i==2
							clim([0 0.5]);
						else
		
						end
						ax = gca;
						fg.Units = 'centimeters';
						fg.Position = [-30 2 28 16];
						ax.Position = [0.09 0.17 0.8 0.75];
						cb.Position(1) = 0.90;
						xlim([-3000 3000]);
						ylim([-200 1000]);
						drawnow();
						disp("irun="+string(irun)+", i="+string(i)+", Time= "+string(TimeSeries.tvec(end)))
					end
					frame = getframe(gcf);
					writeVideo(vidfile, frame);
				end
			catch ME
				rethrow(ME)
			end
			close(vidfile);
			close all
		end
end


function plot_boundaries(physics, sc)
	Egroups = [2 3 4 5 7 8];

    ipc = physics.mesh.ipcount1D;
    Svec = physics.StateVec;

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
            end
            lines(n_el,:,:) = coords;
        end
        hold on
        plot(squeeze(lines(:,:,1))', squeeze(lines(:,:,2))', 'k');
	end
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
	end

	x = zeros(size(physics.mesh.Elementgroups{1}.Elems, 1),5,5);
	y = zeros(size(physics.mesh.Elementgroups{1}.Elems, 1),5,5);
	z = zeros(size(physics.mesh.Elementgroups{1}.Elems, 1),5,5);
	parfor e=1:size(physics.mesh.Elementgroups{1}.Elems, 1)
		Elem_Nodes = physics.mesh.getNodes(1, e);
        IPvar = physics.Request_Info(varName, e, "Interior");
        xy = physics.mesh.getIPCoords(1, e);
		[N, ~, ~] = physics.mesh.getVals(1, e);
		dofsX = physics.dofSpace.getDofIndices(1, Elem_Nodes);
        dofsY = physics.dofSpace.getDofIndices(2, Elem_Nodes);
        
		xloc = zeros(5,5);
		yloc = zeros(5,5);
		zloc = zeros(5,5);
		for ip=1:length(IPvar)
			col = mod(ip,5)+1;
			row = floor((ip-1)/5)+1;
        	xloc(row,col) = xy(1,ip) + defscale*N(ip,:)*physics.StateVec(dofsX);
        	yloc(row,col) = xy(2,ip) + defscale*N(ip,:)*physics.StateVec(dofsY);
        	zloc(row,col) = IPvar(ip);
		end

		x(e,:,:) = xloc;
		y(e,:,:) = yloc;
		z(e,:,:) = zloc;
	end

	for e=1:size(physics.mesh.Elementgroups{1}.Elems, 1)
		surface(squeeze(x(e,:,:)),squeeze(y(e,:,:)),squeeze(z(e,:,:))*1e-6,'EdgeColor','none','FaceColor','interp');
		hold on
	end


end