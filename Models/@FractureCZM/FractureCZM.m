classdef FractureCZM < BaseModel
    %FractureCZM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        tensile
        energy
        dummy
        Hmatswitch
        
        hist %history parameter (element, ip)
        histOld
        
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
            
            obj.tensile = inputs.tensile;
            obj.energy = inputs.energy;
            obj.dummy = inputs.dummy;
            obj.Hmatswitch = inputs.Hmatswitch;
            
            obj.hist = 1e10+zeros(size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount);
            obj.histOld = obj.hist;
        end
        
        function Commit(obj, physics, commit_type)
            if (commit_type == "Pathdep")
                obj.histOld = obj.hist;
            end
        end
        
		function [fc, f_t, K_tu, K_tt, dofsXY] = get_fc(obj, physics)
			if physics.ArcTime.Enable
            	TDof = physics.dofSpace.getDofIndices(physics.tType, 1);
            	time    = physics.StateVec(TDof);
            	timeOld = physics.StateVec_Old(TDof);
            	dt = time - timeOld;
			end

			[direction, elem, ips] = obj.mesh.getNextFracture();
			stresses = physics.Request_Info("stresses",elem,"Interior");
			dstresses = physics.Request_Info("dstressesdxy",elem,"Interior");
			nds = obj.mesh.getNodes(1, elem);

	        dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), nds);
            dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), nds);
            dofsXY = [dofsX; dofsY];
        
			stress = squeeze(stresses(1, ips, :));
            dstress = squeeze(dstresses(1, ips, :, :));
			if (direction == "Horizontal") %compare yy component
				sel = [0;1;0;0];
				ts = obj.tensile;
			else %compare xx component
				sel = [1;0;0;0];
				ts = obj.tensile;
			end

			f_t = 0;
			K_tu = zeros(1, length(dofsXY));
			K_tt = 0;

			fc = sel'*stress-ts;
			if physics.ArcTime.Enable
            	if (dt < physics.ArcTime.TMin && fc>0)
                	f_t = f_t + physics.ArcTime.KDummy * (time - timeOld - physics.ArcTime.TMin+1e-9);
                	K_tt = K_tt + physics.ArcTime.KDummy;
				elseif (dt < physics.ArcTime.TMin)
                	f_t = f_t + physics.ArcTime.KDummy * (time - timeOld - physics.ArcTime.TMin-1e-9);
                	K_tt = K_tt + physics.ArcTime.KDummy;
            	elseif (dt >= physics.ArcTime.TMax && fc<0)
                	f_t = f_t + physics.ArcTime.KDummy * (time - timeOld - physics.ArcTime.TMax-1e-9);
                	K_tt = K_tt + physics.ArcTime.KDummy;
				elseif (dt >= physics.ArcTime.TMax && fc>0)
                	f_t = f_t + physics.ArcTime.KDummy * (time - timeOld - physics.ArcTime.TMax+1e-3);
                	K_tt = K_tt + physics.ArcTime.KDummy;
				else
                	f_t = f_t + sel'*stress-ts;
                	K_tu = K_tu + sel'*dstress;
                	K_tt = K_tt + physics.ArcTime.KStab;
	% 				if (fc<0)
	% 					f_t = f_t + physics.ArcTime.KStab * fc/1e6*((time-timeOld) - 1.1*(time-timeOld));
	% 					K_tt = K_tt + physics.ArcTime.KStab * fc/1e6*-0.1;
	% 				else
	% 					f_t = f_t + physics.ArcTime.KStab * fc/1e6*((time-timeOld) - 0.9*(time-timeOld));
	% 					K_tt = K_tt + physics.ArcTime.KStab * fc/1e6*0.1;
	% 				end
				end
			else
				f_t = f_t + sel'*stress-ts;
			end

		end

        function Irr = Irreversibles(obj, physics)
            Irr = false;
            
            sf = 0;
            if (physics.ArcTime.Enable)
                sf = physics.ArcTime.Tol;
            end

			[fc, ~, ~, ~, ~] = obj.get_fc(physics);
			if (fc+sf>=0)
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
                    
                    %if (h < 0 )
                        %Lumped integration of dummy stiffness
                        C_Lumped = C_Lumped + N(ip,:)'*w(ip);
                    %end

                    %CZM STUFF
                    hstOld = obj.histOld(n_el, ip);
                    tau = zeros(2,1);
                    dtaudh = zeros(2,2);
                     if (xy(2,end-1)>obj.Hmatswitch+1e-2 || true)  %vertical prop
                         if (h>=hstOld)
                             hloc = h;
                             tau(1) = obj.tensile * exp(-obj.tensile*h/obj.energy);
                             dtaudh(1,1) = -obj.tensile.^2/obj.energy * exp(-obj.tensile*h/obj.energy);
						 elseif (h>0)
                             hloc = hstOld;
                             tau(1) = obj.tensile * exp(-obj.tensile*hstOld/obj.energy)*h/hstOld;
                             dtaudh(1,1) = obj.tensile * exp(-obj.tensile*hstOld/obj.energy)/hstOld;
						else
							hloc = hstOld;
						end
                     else
                         hloc = h;
                         tau(1) = 0.0;
                         dtaudh(1,1) = 0.0;
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
                
%                 physics.fint(dofsXY) = physics.fint(dofsXY) + f_el;
%                 physics.K(dofsXY, dofsXY) = physics.K(dofsXY, dofsXY) + K_el;

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


            %% arctime stuff
			if (physics.ArcTime.Enable)
                tscale = 0.00001;
                TDof = physics.dofSpace.getDofIndices(physics.tType, 1);
                
				[fc, f_t, K_tu, K_tt, dofsXY] = obj.get_fc(physics);

                physics.fint(TDof) = physics.fint(TDof) + tscale*f_t;
                physics.K(TDof, dofsXY) = physics.K(TDof, dofsXY) + tscale*K_tu;
                physics.K(TDof, TDof) = physics.K(TDof, TDof) + tscale*K_tt;

				fprintf("Fracture criterium: " + string(fc) + "\n" );
			end


			

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

