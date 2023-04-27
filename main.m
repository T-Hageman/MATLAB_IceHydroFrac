%function main(irun, hasVisc)

	addpath(genpath('./Models'))
	addpath(genpath('./Shapes'))

	hasVisc=1;
	HIce = 980;
	dx_elem = 2.5;
	T_Ice = "profile"; 
	meshName = "mesh_Das_25.mphtxt";

	flowtype = "FrictionFactor";   %"CubicLaw";"FrictionFactor"
	if (hasVisc == 1)
		savefolder = "./Results/Visc_";
	else
		savefolder = "./Results/NoVisc_";
	end
	savefolder=savefolder+flowtype;
	mkdir(savefolder);
	savefolder = savefolder + "/";
	
	fprintf('Starting job\n')
	maxNumCompThreads(16);
	
	delete(gcp('nocreate'))
	parpool('threads')
	
	tmr = tic;
	

	%% input properties
	Files=dir(savefolder);
	nmax = 10*(length(Files)-2);
	restartfrom = nmax;
	if restartfrom>0
		restart = true;
		restart_num = restartfrom;
	else
		restart = false;
		restart_num = 0;
	end
	%restart = false;
	if restart == false
    	% mesh properties
		mesh_in.type = "File";
		mesh_in.FileName = meshName;
		mesh_in.nfrac = 30/dx_elem;
		mesh_in.ipcount1D = 3;
		mesh_in.zeroWeight = true;
	
    	%physics models
    	dt0 = 2;
	
    	physics_in{1}.type = "ViscoElastic";
    	physics_in{1}.Egroup = "Internal";
    	physics_in{1}.young = [9e9; 20e9]; %[9e9; 20e9];
    	physics_in{1}.poisson = [0.33; 0.25]; %[0.33; 0.25];
		physics_in{1}.A = hasVisc*5e-24;
		physics_in{1}.n = 3;
    	physics_in{1}.Hmatswitch = 0;
	
    	physics_in{2}.type = "SelfWeight";
    	physics_in{2}.Egroup = "Internal";
    	physics_in{2}.density = [910; 2500];
	
    	physics_in{3}.type = "FractureCZM";
    	physics_in{3}.Egroup = "Fracture";
    	physics_in{3}.tensile = 0.3e6;
    	physics_in{3}.energy = (5.0e-3)*physics_in{3}.tensile;
    	physics_in{3}.dummy = 1e10; %1e10
    	physics_in{3}.Hmatswitch = 0;
	
    	physics_in{4}.type = "FractureFluid";
    	physics_in{4}.Egroup = "Fracture";
    	physics_in{4}.QTip = 0;
    	physics_in{4}.visc = 1.0e-3;
    	physics_in{4}.Kf   = 1.0e9;
    	physics_in{4}.FlowModel = flowtype;%"CubicLaw";"FrictionFactor"
    	physics_in{4}.melt = true;
		physics_in{4}.freeze = true;
		physics_in{4}.rho_ice = 910;
		physics_in{4}.rho_water = 1000;
		physics_in{4}.cp_ice = 2115;
		physics_in{4}.melt_heat = 335000;
		physics_in{4}.T_ice = T_Ice;
		physics_in{4}.k_ice = 2;
	
    	physics_in{5}.type = "Constrainer";
    	physics_in{5}.Egroup = "Bottom";
    	physics_in{5}.dofs = {"dy"};
    	physics_in{5}.conVal = [0];
	
    	physics_in{6}.type = "Constrainer";
    	physics_in{6}.Egroup = "Left";
    	physics_in{6}.dofs = {"dx"};
    	physics_in{6}.conVal = 0*[-1e-5]*4000*0.5;
	
    	physics_in{7}.type = "Constrainer";
    	physics_in{7}.Egroup = "Right";
    	physics_in{7}.dofs = {"dx"};
    	physics_in{7}.conVal = 0*[1e-5]*4000*0.5;
	
    	physics_in{8}.type = "LakeBoundary";
    	physics_in{8}.Egroup = "Fracture";
    	physics_in{8}.p0 = 1.0e5;
    	physics_in{8}.dummy = 1e-6;
		physics_in{8}.Reference = 1;
	
		physics_in{9}.type = "Inertia";
		physics_in{9}.Egroup = "Internal";
    	physics_in{9}.density = 950;
		physics_in{9}.beta    = 0.4;
		physics_in{9}.gamma   = 0.75;
	
    	%% solver inputs
    	solver_in.maxIt = 100;
    	solver_in.Conv = 1e-6;
    	solver_in.tiny = 0.1;%1e-6;
    	solver_in.linesearch = true;
    	solver_in.linesearchLims = [0.1 1.0];
	
    	%% initialization
    	mesh = Mesh(mesh_in);
    	%mesh.plot(true, true, false, false);
    	mesh.check();
	
    	physics = Physics(mesh, physics_in, dt0);
    	physics.time = 0.0;
	
    	n_max = 3600;
    	solver = Solver(physics, solver_in);
    	TimeSeries.tvec = 0;
    	TimeSeries.Lfrac = [mesh.Area(9)];
		TimeSeries.Qvec = [0,0,0];
		TimeSeries.MeltqVec = 0;
		TimeSeries.Qinflow = 0;
		TimeSeries.qCurrent = 0;
		TimeSeries.upLift = [0,0];
    	startstep = 1;
	else
    	filename = savefolder+string(restart_num);
    	load(filename, "mesh","physics","solver","dt","n_max", "TimeSeries");
    	startstep = restart_num+1;
	end
	
	fprintf('Starting timesteps\n')
	
	for tstep = startstep:n_max
    	disp("Step: "+string(tstep));
        dt = physics.dt;

    	solver.Solve();
    	
    	TimeSeries.tvec(end+1) = physics.time;
    	TimeSeries.Lfrac(end+1) = mesh.Area(9);
		TimeSeries.Qvec(end+1,:) = physics.models{4}.QMeltTot;
		TimeSeries.Qinflow(end+1) = physics.models{8}.QTotal(end);
		TimeSeries.qCurrent(end+1) = physics.models{8}.qCurrent(end);
		TimeSeries.upLift(end+1,1) = physics.models{8}.dxCurrent(end);
		TimeSeries.upLift(end,2) = physics.models{8}.dyCurrent(end);

    	disp("time: "+physics.time);
	
    	%close all
		if mod(tstep, 10) == 0
			if (false)
				plotres(physics, TimeSeries);
			end
        	
        	filename = savefolder+string(tstep);
        	save(filename, "mesh","physics","solver","dt","n_max","TimeSeries");
		end
	
		frozen = physics.models{4}.checkFrozen(physics);
		if frozen
			fprintf("Fracture frozen, exiting")
			break
		end
	end
	
	if (false)
		plotres(physics, TimeSeries);
	end
	filename = savefolder+"End";
	save(filename, "mesh","physics","solver","dt","n_max","TimeSeries");
	
	toc(tmr)
%end


	function plotres(physics, TimeSeries)
		LFrac = 3.2e3;
		ALake = 5.6e6; 
		m3_msec_to_m_hour = LFrac/ALake*3600;
		m3_to_mWaterHeigt = LFrac/ALake;

		figure(71)
		clf
		subplot(2,1,1)
    		%physics.PlotNodal("dx", 1000, "Internal",1);
			title("t="+string(physics.time));
		subplot(2,1,2)
    		physics.PlotNodal("pd", 1000, "Fracture",1);
    		physics.models{4}.plotTotHeight(physics, 1000);
			plot_boundaries(physics, 1000)
	
		figure(72)
		clf
		subplot(3,1,1)
			plot(TimeSeries.tvec, TimeSeries.Qvec(:,1));
			hold on
			plot(TimeSeries.tvec, TimeSeries.Qvec(:,2));
			plot(TimeSeries.tvec, TimeSeries.Qvec(:,3));
			legend('Ice desorbtion', 'Flow produced', 'Melting process');
			xlabel('t [sec]')
			ylabel('Q_{thermal} [J]')
	
		subplot(3,1,2)
    		yyaxis left
    		plot(TimeSeries.tvec, TimeSeries.Qinflow*m3_to_mWaterHeigt);
			ylabel('Q_{flow} [m_{lakeHeight}]')
    		yyaxis right
    		plot(TimeSeries.tvec, TimeSeries.Lfrac);
			xlabel('t [sec]')
			ylabel('L_{frac} [m]')

		subplot(3,1,3)
			yyaxis left
			plot(TimeSeries.tvec, TimeSeries.qCurrent*m3_msec_to_m_hour);
			xlabel('t [sec]')
			ylabel('q [m_{lakeHeight}/hour]')
			yyaxis right
			plot(TimeSeries.tvec, TimeSeries.upLift(:,1));
			hold on
			plot(TimeSeries.tvec, TimeSeries.upLift(:,2));
			ylabel('dx dy  [m]')

		drawnow();
	end

function plot_boundaries(physics, sc)
	Egroups = [2 3 4 5];

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
