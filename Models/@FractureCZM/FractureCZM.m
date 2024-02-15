classdef FractureCZM < BaseModel
    %FractureCZM Implements cohesive fracture and its propagation crierium,
	%input parameters required:
    % physics_in{4}.type = "FractureCZM";
	% physics_in{4}.Egroup = "Fracture";
	% physics_in{4}.energy = 10;			%Fracture release energy [J/m^2]
	% physics_in{4}.dummy = 0*1e10;		%Dummy stiffness to prevent walls from penetrating
	% physics_in{4}.Hmatswitch = 0;		%Depth of ice-rock interface
	% physics_in{4}.T_ice = T_Ice;		%Temperature profile of ice
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        energy
        dummy
        
        hist %history parameter (element, ip)
        histOld
		T_ice
		ft
    end
    
    methods
        function obj = FractureCZM(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "FractureCZM";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType({"dx","dy"});
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            obj.energy = inputs.energy;
            obj.dummy = inputs.dummy;
			obj.T_ice = inputs.T_ice;
			obj.ft = @(T) 2.0-0.0068*(T+273.15);
            
            obj.hist = 1e10+zeros(size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount);
            obj.histOld = obj.hist;
        end
        
        function Commit(obj, physics, commit_type)
            if (commit_type == "Pathdep")
                obj.histOld = obj.hist;
            end
        end
        
		function [fc, dofsXY] = get_fc(obj, physics)
			[direction, elem, ips] = obj.mesh.getNextFracture();
			stresses = physics.Request_Info("stresses",elem,"Interior");
			dstresses = physics.Request_Info("dstressesdxy",elem,"Interior");
			nds = obj.mesh.getNodes(1, elem);

	        dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), nds);
            dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), nds);
            dofsXY = [dofsX; dofsY];

			if (direction == "Horizontal")
				y = 0;
			else
				y = mean(obj.mesh.Nodes(nds,2));
			end
        
			fc = -1e12;
			for i=1:length(ips)
				stress = squeeze(stresses(1, ips(i), :));
	
				Ice_temp = obj.T_ice(y);
				ts = obj.ft(Ice_temp);
	
				if (direction == "Horizontal") %compare yy component
					sel = [0;1;0;0];
				else %compare xx component
					sel = [1;0;0;0];
				end
	
				fc  = max(fc,sel'*stress-ts);
			end
		end

        function Irr = Irreversibles(obj, physics)
            Irr = false;
            
			[fc, ~] = obj.get_fc(physics);
			if (fc>=0 && physics.time>=0)
				obj.mesh.Propagate_Disc_New(physics);
				Irr = true;
			end
            fprintf("Fracture criterium: " + string(fc) + "\n" );

        end
        
        function getKf(obj, physics)
            while (size(obj.hist, 1)<size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1))
                obj.histOld(end+1,:) = 1e-4;
                obj.hist(end+1,:) = 1e-4;
            end
            
            fprintf("        FractureCZM get Matrix:")
            t = tic;
            
            Svec = physics.StateVec;

            dofmatX = [];
            dofmatY = [];
            kmat = [];
            fvec = [];
            dofvec = [];

            newHist = zeros(size(obj.hist));
            ipc = obj.mesh.ipcount1D;

            parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                Elem_Nodes   = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, ~, w]    = obj.mesh.getVals(obj.myGroupIndex, n_el);
                [nvec, tvec] = obj.mesh.getNormals(obj.myGroupIndex, n_el);
                xy = obj.mesh.getIPCoords(obj.myGroupIndex, n_el);
                
                dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                dofsXY = [dofsX; dofsY];
                
                X = Svec(dofsX);
                Y = Svec(dofsY);
                XY = [X;Y];
                
                f_el = zeros(length(dofsXY), 1);
                K_el = zeros(length(dofsXY));
                
                C_Lumped = zeros(size(N, 2),1);
                for ip=1:ipc
                    Nd = obj.getNd(N(ip,:));
                    
                    h = nvec(ip,:)*Nd*XY;

					Ice_temp = obj.T_ice(xy(2,ip));
					f_t = obj.ft(Ice_temp);
                    
                    %Lumped integration of dummy stiffness
                    C_Lumped = C_Lumped + N(ip,:)'*w(ip);

                    %CZM STUFF
                    hstOld = obj.histOld(n_el, ip);
                    tau = zeros(2,1);
                    dtaudh = zeros(2,2);
                    if (h>=hstOld)
                         hloc = h;
                         tau(1) = f_t * exp(-f_t*h/obj.energy);
                         dtaudh(1,1) = -f_t.^2/obj.energy * exp(-f_t*h/obj.energy);
					elseif (h>0)
                         hloc = hstOld;
                         tau(1) = f_t * exp(-f_t*hstOld/obj.energy)*h/hstOld;
                         dtaudh(1,1) = f_t * exp(-f_t*hstOld/obj.energy)/hstOld;
					else
						hloc = hstOld;
						%tractions via no-pen condition
					end
                    newHist(n_el, ip) = hloc;
                    R = zeros(2,2);
                    R(1,:) = nvec(ip,:); R(2,:) = tvec(ip,:);
                    f_el = f_el + w(ip)*Nd'*R'*tau;
                    K_el = K_el + w(ip)*Nd'*R'*dtaudh*R*Nd;
                end
                
                for cp=1:length(C_Lumped)
                    NL = zeros(size(N, 2), 1);
                    NL(cp) = 1.0;
                    Nd = obj.getNd(NL);
                    
                    n_est = nvec(1,:);
                    if (n_est*Nd*XY < 0)
                    	f_el = f_el + C_Lumped(cp)*obj.dummy*(Nd'*(n_est'*n_est)*Nd)*XY;
                    	K_el = K_el + C_Lumped(cp)*obj.dummy*(Nd'*(n_est'*n_est)*Nd);
					end
				end

                fvec = [fvec; f_el];
                dofvec = [dofvec; dofsXY];

                [dofmatxloc,dofmatyloc] = ndgrid(dofsXY,dofsXY);
                dofmatX = [dofmatX; dofmatxloc(:)];
                dofmatY = [dofmatY; dofmatyloc(:)];
                kmat = [kmat; K_el(:)];
            end 
            
            obj.hist = newHist;

            physics.fint = physics.fint + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint), 1);
            physics.K = physics.K + sparse(dofmatX, dofmatY, kmat, length(physics.fint),length(physics.fint));		

            tElapsed = toc(t);
            fprintf("            (Assemble time:"+string(tElapsed)+")\n");
        end
        
        function Nd = getNd(~, vals)
            cp_count = length(vals);
            Nd = zeros(2, cp_count*2);
            for ii = 1:cp_count %using plane strain e_zz = 0
				%dx
				Nd(1, ii) = vals(ii);
				Nd(1, ii + cp_count) = -vals(ii);

				%dy
				Nd(2, ii + 2*cp_count) = vals(ii);
				Nd(2, ii + 3*cp_count) = -vals(ii);
            end
        end

        
    end
end

