function main()
	addpath(genpath('./Models'))
	addpath(genpath('./Shapes'))

	%flag indicating whether running locally or on a high-performance
	%computing cluster (surpresses visual outputs)
	hpc = false;

	%folder to which the results will be saved
	savefolder = "./Results";

	%Temperature profile used
	if true
		T_Ice = @(y) min(0.0,1.0e-4*y*(y-250)); 
	else
		load("TProfile.mat","T_Ice");
	end

	%file containing the mesh, and size of interface elements (must match
	%the size used within the mesh file)
	dx_elem = 2.5;
	meshName = "Meshes/mesh_3000x300m.mphtxt"; %mesh_Das_25.mphtxt

	%start parrallel pool
	mkdir(savefolder);
	savefolder = savefolder + "/";
	
	fprintf('Starting job\n')
	if hpc
		maxNumCompThreads(32);
	else
		maxNumCompThreads(8);
	end
	
	delete(gcp('nocreate'))
	parpool('threads')
	
	tmr = tic;

	% check if simulation needs to be restarted from previous files
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
	
	%% input properties
	if restart == false
    	% mesh properties
		mesh_in.type = "File";		
		mesh_in.FileName = meshName;	
		mesh_in.nfrac = 30/dx_elem;		%number of initial interface elements
		mesh_in.ipcount1D = 3;			%number of integration points per direction
		mesh_in.zeroWeight = true;		%must be set to true, add zero-weight integration points to boundary of elements 
	
    	%physics models
    	dt0 = 2;
	
		%Interior of ice and rock: Momentum balance
    	physics_in{1}.type = "ViscoElastic";
    	physics_in{1}.Egroup = "Internal";
    	physics_in{1}.young = [9e9; 20e9];		%Youngs modulus of Ice and Rock [Pa]
    	physics_in{1}.poisson = [0.33; 0.25];	%Poisson ratio of ice and rock [-]
		physics_in{1}.A0 = 5e-24;				%Glenns law creep coefficient [Pa^-3 s^-1]
		physics_in{1}.Q = 150e3;				%Energy for scaling creep with temperature [J]
		physics_in{1}.TRef = 273.15;			%Reference temperature at which A0 is determined [K]
		physics_in{1}.n = 3;					%Glenns law creep exponent [-]
		physics_in{1}.T_Ice = T_Ice;			%Temperature profile for the ice
    	physics_in{1}.Hmatswitch = ["FLeft", "FRight"];			%Element groups indicating ice-rock boundary

		%Interior of ice and rock: Contribution due to inertia
		physics_in{2}.type = "Inertia";
		physics_in{2}.Egroup = "Internal";
    	physics_in{2}.density = 950;		%Density [kg/m^3]
		physics_in{2}.beta    = 0.4;		%Time discretisation constant (Newmark scheme)
		physics_in{2}.gamma   = 0.75;		%Time discretisation constant (Newmark scheme)
	
		%Interior of ice and rock: Contribution due to gravity
    	physics_in{3}.type = "SelfWeight";
    	physics_in{3}.Egroup = "Internal";
    	physics_in{3}.density = [910; 2500];	%Density [kg/m^3]
	
		%Fracture interface: Cohesive zone model and propagation
    	physics_in{4}.type = "FractureCZM";
    	physics_in{4}.Egroup = "Fracture";
    	physics_in{4}.energy = 10;			%Fracture release energy [J/m^2]
    	physics_in{4}.dummy = 0*1e10;		%Dummy stiffness to prevent walls from penetrating
    	physics_in{4}.Hmatswitch = ["FLeft", "FRight"];		%Depth of ice-rock interface
		physics_in{4}.T_ice = T_Ice;		%Temperature profile of ice
	
		%Fracture interface: Fluid flow and melting
    	physics_in{5}.type = "FractureFluid";
    	physics_in{5}.Egroup = "Fracture";
    	physics_in{5}.visc = 1.0e-3;		%water viscosity [Pa s]
    	physics_in{5}.Kf   = 1.0e9;			%Water compressibility [Pa]
    	physics_in{5}.FlowModel = "FrictionFactor";   %"CubicLaw";"FrictionFactor", Model used for fluid flow within crevasse
    	physics_in{5}.melt = true;			%Include wall melting
		physics_in{5}.freeze = true;		%Include wall freezing
		physics_in{5}.rho_ice = 910;		%Density of ice [kg/m^3]
		physics_in{5}.rho_water = 1000;		%Density of water [kg/m^3]
		physics_in{5}.cp_ice = 2115;		%heat capacity of ice [J/kg]
		physics_in{5}.melt_heat = 335000;	%latent heat of ice-water tansition [J/kg]
		physics_in{5}.T_ice = T_Ice;		%Temperature profile
		physics_in{5}.k_ice = 2;			%Thermal conductivity of ice [J/K/m^2]

		%Fracture inlet: Imposed inlet pressure
		physics_in{6}.type = "LakeBoundary";
    	physics_in{6}.Egroup = "Fracture";
		physics_in{6}.Surface = "Top";
    	physics_in{6}.p0 = 1.0e5;		%Reference pressure imposed at the top [Pa]
    	physics_in{6}.dummy = 1e-6;		%dummy constant used to enforce this pressure

		%mechanical boundary conditions:
    	physics_in{7}.type = "Constrainer";
    	physics_in{7}.Egroup = "Bottom";
    	physics_in{7}.dofs = {"dy"};
    	physics_in{7}.conVal = [0];
	
    	physics_in{8}.type = "Constrainer";
    	physics_in{8}.Egroup = "Left";
    	physics_in{8}.dofs = {"dx"};
    	physics_in{8}.conVal = 0;
	
    	physics_in{9}.type = "Constrainer";
    	physics_in{9}.Egroup = "Right";
    	physics_in{9}.dofs = {"dx"};
    	physics_in{9}.conVal = 0;
	
    	%% solver inputs
    	solver_in.maxIt = 100;		%max number of iterations for nonlinear solver
    	solver_in.Conv = 1e-6;		%relative energy error below which system is considered converged
    	solver_in.tiny = 1e-3;		%absolute energy error below which system is considered converged
    	solver_in.linesearch = true;%use a linear line-search scheme
    	solver_in.linesearchLims = [0.1 1.0];%limits for line-search to be performed within
	
    	%% initialization
		%load mesh from file
    	mesh = Mesh(mesh_in);
    	mesh.plot(true, true, false, false);
    	mesh.check();
	
		%initialize physics
    	physics = Physics(mesh, physics_in, dt0);

		%time solver settings
    	physics.time = -24*3600; %initial time (negative for initialization period) [s]
    	t_max = 3600*2;	%max simulation time [s]

		%initialize solvers
    	solver = Solver(physics, solver_in);

		%structure for saving output time series
    	TimeSeries.tvec = physics.time;		%times of outputs[s]
    	TimeSeries.Lfrac = [mesh.Area(9)];	%crevasse depth/length
		TimeSeries.Qvec = [0,0,0];			%thermal fluxes ('Ice desorbtion', 'Flow produced', 'Melting process')
		TimeSeries.Qinflow = 0;				%total volume of fluid that has entered the crevasse
		TimeSeries.qCurrent = 0;			%current inflow rate
		TimeSeries.upLift = [0,0];			%surface uplift at the centre of the surface
		physics.models{6}.updateSurfaceElevation(physics);	
		TimeSeries.SurfaceCoords(:,1) = physics.models{6}.surface_X;	%full uplift profile at top surface (x coordinates of points)
		TimeSeries.SurfaceCoords(:,2) = physics.models{6}.surface_Y;	%full uplift profile at top surface (y coordinates of points)
		TimeSeries.SurfaceDisp(1,:,1) = physics.models{6}.surface_dX;	%full uplift profile at top surface	(horizontal displacement)
		TimeSeries.SurfaceDisp(1,:,2) = physics.models{6}.surface_dY;	%full uplift profile at top surface (vertical displacement)

		startstep = 1;
	else
    	filename = savefolder+string(restart_num);
    	load(filename, "mesh","physics","solver","t_max", "TimeSeries");
    	startstep = restart_num+1;
	end
	
	%start of actual simultation
	fprintf('Starting timesteps\n')
	n_max = 1e6;
	for tstep = startstep:n_max
    	disp("Step: "+string(tstep)+", Time: "+string(physics.time));

		%Set time increment based on either being in initial or propagation
		%period
		if (tstep<10)
			physics.dt = 1;
		elseif (physics.time<-600)
			physics.dt = 600;
		elseif (physics.time<-0.001)
			physics.dt = -physics.time;
		else
			physics.dt = 2;
		end

		%solve current time increment
    	solver.Solve();
    	
		%save timedata 
		physics.models{6}.updateSurfaceElevation(physics);
    	TimeSeries.tvec(end+1) = physics.time;%times of outputs[s]
    	TimeSeries.Lfrac(end+1) = mesh.Area(9);%crevasse depth/length
		TimeSeries.Qvec(end+1,:) = physics.models{5}.QMeltTot;%thermal fluxes ('Ice desorbtion', 'Flow produced', 'Melting process')
		TimeSeries.Qinflow(end+1) = physics.models{6}.QTotal(end);%total volume of fluid that has entered the crevasse
		TimeSeries.qCurrent(end+1) = physics.models{6}.qCurrent(end);%current inflow rate
		TimeSeries.upLift(end+1,1) = physics.models{6}.dxCurrent(end);%surface uplift at the centre of the surface
		TimeSeries.upLift(end,2) = physics.models{6}.dyCurrent(end);%surface uplift at the centre of the surface
		TimeSeries.SurfaceDisp(end+1,:,1) = physics.models{6}.surface_dX;%full uplift profile at top surface	(horizontal displacement)
		TimeSeries.SurfaceDisp(end,:,2) = physics.models{6}.surface_dY;%full uplift profile at top surface (vertical displacement)
	
		%plot current state
		if (hpc == false)
			plotres(physics, TimeSeries);
		end

		%save outputs for post-processing and restarting
		if mod(tstep, 10) == 0
        	filename = savefolder+string(tstep);
        	save(filename, "mesh","physics","solver","t_max","TimeSeries");
		end
	
		%stopping criteria
		frozen = physics.models{5}.checkFrozen(physics);
		if frozen
			fprintf("Fracture frozen, exiting")
			break
		end
		if (physics.time>t_max)
			break
		end
	end
	
	if (hpc == false)
		plotres(physics, TimeSeries);
	end
	filename = savefolder+"End";
	save(filename, "mesh","physics","solver","t_max","TimeSeries");
	
	toc(tmr)
end


function plotres(physics, TimeSeries)
	% function for plotting outputs during simulations
	t0 = find(TimeSeries.tvec==0);
	if (isempty(t0))
		t0=1;
	end

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
		physics.models{5}.plotTotHeight(physics, 1000);
		plot_boundaries(physics, 1000)

	figure(72)
	clf
	subplot(3,1,1)
		plot(TimeSeries.tvec(t0:end), TimeSeries.Qvec(t0:end,1));
		hold on
		plot(TimeSeries.tvec(t0:end), TimeSeries.Qvec(t0:end,2));
		plot(TimeSeries.tvec(t0:end), TimeSeries.Qvec(t0:end,3));
		legend('Ice desorbtion', 'Flow produced', 'Melting process');
		xlabel('t [sec]')
		ylabel('Q_{thermal} [J]')

	subplot(3,1,2)
		yyaxis left
		plot(TimeSeries.tvec(t0:end), TimeSeries.Qinflow(t0:end)*m3_to_mWaterHeigt);
		ylabel('Q_{flow} [m_{lakeHeight}]')
		yyaxis right
		plot(TimeSeries.tvec(t0:end), TimeSeries.Lfrac(t0:end));
		xlabel('t [sec]')
		ylabel('L_{frac} [m]')

	subplot(3,1,3)
		yyaxis left
		plot(TimeSeries.tvec(t0:end), TimeSeries.qCurrent(t0:end)*m3_msec_to_m_hour);
		xlabel('t [sec]')
		ylabel('q [m_{lakeHeight}/hour]')
		yyaxis right
		plot(TimeSeries.tvec(t0:end), TimeSeries.upLift(t0:end,1));
		hold on
		plot(TimeSeries.tvec(t0:end), TimeSeries.upLift(t0:end,2));
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
