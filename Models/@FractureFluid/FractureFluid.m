classdef FractureFluid < BaseModel
    %FractureFluid Fluid flows and thermals within the crevasse. Input
	%parameters:
    % physics_in{5}.type = "FractureFluid";
	% physics_in{5}.Egroup = "Fracture";
	% physics_in{5}.visc = 1.0e-3;		%water viscosity [Pa s]
	% physics_in{5}.Kf   = 1.0e9;			%Water compressibility [Pa]
	% physics_in{5}.FlowModel = "FrictionFactor";   %"CubicLaw";"FrictionFactor", Model used for fluid flow within crevasse
	% physics_in{5}.melt = true;			%Include wall melting
	% physics_in{5}.freeze = true;		%Include wall freezing
	% physics_in{5}.rho_ice = 910;		%Density of ice [kg/m^3]
	% physics_in{5}.rho_water = 1000;		%Density of water [kg/m^3]
	% physics_in{5}.cp_ice = 2115;		%heat capacity of ice [J/kg]
	% physics_in{5}.melt_heat = 335000;	%latent heat of ice-water tansition [J/kg]
	% physics_in{5}.T_ice = T_Ice;		%Temperature profile
	% physics_in{5}.k_ice = 2;			%Thermal conductivity of ice [J/K/m^2]
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        FlowModel
        visc
        Kf

        hmelt
        hmeltOld
		qxSaved
		qxSavedOld

		melt
		freeze

		rho_ice
		rho_water
		cp_ice
		melt_heat
		T_ice
		k_ice

		QMeltTot
		QMeltTotOld

		qMeltTotOld
		qMeltTot

		tfrac
    end
    
    methods
        function obj = FractureFluid(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "FractureFluid";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType({"dx","dy","pd"});
            
            halfnodes = [];
            for n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                Elem_Nodes   = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                halfnodes(end+1:end+length(Elem_Nodes)/2) = Elem_Nodes(1:length(Elem_Nodes)/2);
            end
            halfnodes = unique(halfnodes);
            
            obj.dofSpace.addDofs(obj.dofTypeIndices, halfnodes);
            
            obj.FlowModel = inputs.FlowModel;
            obj.visc = inputs.visc;
            obj.Kf = inputs.Kf;

            obj.melt = inputs.melt;
			obj.freeze = inputs.freeze;

			obj.rho_ice = inputs.rho_ice;
			obj.rho_water = inputs.rho_water;
			obj.cp_ice = inputs.cp_ice;
			obj.melt_heat = inputs.melt_heat;
			obj.T_ice = inputs.T_ice;
			obj.k_ice = inputs.k_ice;

            obj.hmelt = 1e-3+zeros(size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount);
            obj.hmeltOld = obj.hmelt;
			obj.qxSaved = zeros(size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount);
			obj.qxSavedOld = obj.qxSaved;

			obj.tfrac = -1e10+zeros(size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount);

			obj.QMeltTotOld = [0, 0, 0];
			obj.QMeltTot = [0, 0, 0];

			obj.qMeltTot = 0;
			obj.qMeltTotOld = 0;
        end
        
        function Commit(obj, physics, commit_type)
            if (commit_type == "Timedep")
                obj.hmeltOld = obj.hmelt;
				obj.qxSavedOld = obj.qxSaved;

				obj.QMeltTotOld = obj.QMeltTot;
				obj.qMeltTotOld = obj.qMeltTot;
            end
        end

		function frozen = checkFrozen(obj, physics)
            while (size(obj.hmelt, 1)<size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1))
                obj.hmeltOld(end+1,:) = 0.0;
                obj.hmelt(end+1,:) = 0.0;
            end

            ipc = obj.mesh.ipcount1D;
            Svec = physics.StateVec;

			LengthFrozen = 0;
			LengthNotFrozen = 0;
			for n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                Elem_Nodes   = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, G, w]    = obj.mesh.getVals(obj.myGroupIndex, n_el);
                [nvec, tvec] = obj.mesh.getNormals(obj.myGroupIndex, n_el);
                
                dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);

                X = Svec(dofsX);
                Y = Svec(dofsY);
                XY = [X;Y];

				for ip=1:ipc
                    Nd = obj.getNd(N(ip,:));
                    ujump = max(0,nvec(ip,:)*Nd*XY);
					mlt = obj.hmelt(n_el,ip);

					if (ujump+mlt>=-1e-9)
						LengthNotFrozen = LengthNotFrozen + w(ip);
					else
						LengthFrozen = LengthFrozen + w(ip);
					end
				end
			end
			
			FrozenPercentage = LengthFrozen/(LengthFrozen+LengthNotFrozen);
			fprintf("Fracture "+string(FrozenPercentage*100)+"%% frozen\n")
			if FrozenPercentage>0.50
				frozen = true;
			else
				frozen = false;
			end

		end

        function getKf(obj, physics)
            fprintf("        FractureFluid get Matrix:")
            t = tic;

            dt = physics.dt;
			tOld = physics.time;


			while (size(obj.qxSaved, 1)<size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1))
				obj.qxSaved(end+1,:) = 0.0;
				obj.qxSavedOld(end+1,:) = 0.0;
				obj.tfrac(end+1,:) = tOld;
			end
			while (size(obj.hmelt, 1)<size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1))
                obj.hmeltOld(end+1,:) = 0.0;
                obj.hmelt(end+1,:) = 0.0;
			end
            
            if (length(obj.mesh.Tjunction)>0) %T-junction node constrain
                TNodes = obj.mesh.Tjunction;
                
                Tdofs = physics.dofSpace.getDofIndices(obj.dofTypeIndices(3), TNodes );
                Ptips = physics.StateVec(Tdofs);
                ktip = 1e5; %penalty approach constraints
                
                q_el = zeros(2, 1);
                K_el = zeros(2, 2);
                q_el(1) = ktip*(Ptips(2) - Ptips(1));
                q_el(2) = ktip*(Ptips(1) - Ptips(2));
                K_el = [-ktip ktip; ktip -ktip];
                
                physics.fint(Tdofs) = physics.fint(Tdofs) + q_el;
                physics.K(Tdofs, Tdofs) = physics.K(Tdofs, Tdofs) + K_el;
            end
            
            dofmatX = [];
            dofmatY = [];
            kmat = [];
            fvec = [];
            dofvec = [];

            hmeltNew = zeros(size(obj.hmelt,1), size(obj.hmelt,2));
			qx_ips = zeros(size(obj.hmelt,1), size(obj.hmelt,2));

            ipc = obj.mesh.ipcount1D;

            Svec = physics.StateVec;
            SvecOld = physics.StateVec_Old;

			Q_Update = [0, 0, 0];
			MeltProd = 0;

            parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)  %parfor
                Elem_Nodes   = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, G, w]    = obj.mesh.getVals(obj.myGroupIndex, n_el);
                [nvec, tvec] = obj.mesh.getNormals(obj.myGroupIndex, n_el);
				xy = obj.mesh.getIPCoords(obj.myGroupIndex, n_el);
                
                dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                dofsPD = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes(1:length(Elem_Nodes)/2));
                dofsXY = [dofsX; dofsY];
                
                X = Svec(dofsX);
                Y = Svec(dofsY);
                PD= Svec(dofsPD);
                XY = [X;Y];
                
                XOld = SvecOld(dofsX);
                YOld = SvecOld(dofsY);
                PDOld= SvecOld(dofsPD);
                XYOld = [XOld;YOld];
                
                f_u = zeros(length(dofsXY), 1);
                K_up = zeros(length(dofsXY), length(dofsPD));
                
                f_p = zeros(length(dofsPD), 1);
                K_pu = zeros(length(dofsPD), length(dofsXY));
                K_pp = zeros(length(dofsPD), length(dofsPD));
                
                g = [0;-9.81];
				C_Lumped = zeros(size(N, 2),1);
                for ip=1:ipc
                    Nd = obj.getNd(N(ip,:));
                    ujump = nvec(ip,:)*Nd*XY;
                    ujumpOld = nvec(ip,:)*Nd*XYOld;
                    dp_dx = G(ip,:,1)*PD-tvec(ip,:)*g*obj.rho_water;

                    %% Forces
                    %f_u = f_u - w(ip)*Nd'*nvec(ip,:)'*N(ip,:)*PD;
                    %K_up = K_up - w(ip)*Nd'*nvec(ip,:)'*N(ip,:);
					C_Lumped = C_Lumped + N(ip,:)'*w(ip);
                    
					%% displacement opening fluid flux
                    %f_p  = f_p  - w(ip)*N(ip,:)'*(ujump-ujumpOld)/dt;
                    %K_pu = K_pu - 0.5*w(ip)*N(ip,:)'*nvec(ip,:)*Nd/dt;
					

					%% Return mapping scheme
					Ice_temp = obj.T_ice(xy(2,ip));
					hMelt_hist = obj.hmeltOld(n_el, ip);
					tfr = obj.tfrac(n_el, ip);
					%initGuess = [obj.qxSavedOld(n_el, ip); obj.hmeltOld(n_el, ip); obj.hmeltOld(n_el, ip)+ujump];
					initGuess = [obj.qxSaved(n_el, ip); obj.hmelt(n_el, ip); obj.hmelt(n_el, ip)+ujump];
					[qx, hmelt_ip, h, derivs, Qprod] = obj.get_qx_hmelt_h(dp_dx, ujump, tfr, tOld, dt, hMelt_hist, initGuess, Ice_temp);
				

					Q_Update = Q_Update + Qprod*w(ip)*dt;
					hmeltNew(n_el, ip) = hmelt_ip;
					qx_ips(n_el, ip) = qx;

					%% wall melting
					f_p = f_p   - w(ip)*(1-obj.rho_ice/obj.rho_water)/dt*N(ip,:)'*(hmelt_ip-hMelt_hist);

					K_pu = K_pu - w(ip)*(1-obj.rho_ice/obj.rho_water)/dt*derivs.hmelt_du*N(ip,:)'*nvec(ip,:)*Nd;
                    K_pp = K_pp - w(ip)*(1-obj.rho_ice/obj.rho_water)/dt*derivs.hmelt_dpdx*N(ip,:)'*G(ip,:,1);

					MeltProd = MeltProd + obj.rho_ice/obj.rho_water*w(ip)*(hmelt_ip-hMelt_hist);

					%% fluid flux
                    f_p  = f_p + w(ip)*G(ip,:,1)'*qx;

                    K_pu = K_pu + w(ip)*G(ip,:,1)'*derivs.qx_du*nvec(ip,:)*Nd;
                    K_pp = K_pp + w(ip)*G(ip,:,1)'*derivs.qx_dpdx*G(ip,:,1);

                    %% compressibility
                    f_p  = f_p - w(ip) * h/obj.Kf * N(ip,:)'*N(ip,:)*(PD-PDOld)/dt;

                    K_pu = K_pu - w(ip)*N(ip,:)'*derivs.h_du*nvec(ip,:)*Nd * 1/obj.Kf * (N(ip,:)*(PD-PDOld))/dt ;
                    K_pp = K_pp - w(ip)*N(ip,:)'*(1/obj.Kf*derivs.h_dpdx*G(ip,:,1)*(N(ip,:)*(PD-PDOld)/dt) + h/obj.Kf * N(ip,:)/dt);			
                end
                
                for cp=1:length(C_Lumped)
                    NL = zeros(size(N, 2), 1)'; 
                    NL(cp) = 1.0;
                    Nd = obj.getNd(NL);
                    
                    n_est = nvec(1,:);
                    f_u  = f_u  - C_Lumped(cp)*(Nd'*n_est'*NL)*PD;
                    K_up = K_up - C_Lumped(cp)*(Nd'*n_est'*NL);

					f_p  = f_p  - C_Lumped(cp)*NL'*n_est*Nd*(XY-XYOld)/dt;
                    K_pu = K_pu - 0.5*C_Lumped(cp)*NL'*n_est*Nd/dt;

					K_pp(cp,cp) = K_pp(cp,cp) - C_Lumped(cp)*1e-9; %some extra damping/stabilisation
                end

				%forces
				fvec = [fvec; f_u];
				dofvec = [dofvec; dofsXY];
				
				[dofmatxloc,dofmatyloc] = ndgrid(dofsXY,dofsPD);
				dofmatX = [dofmatX; dofmatxloc(:)];
				dofmatY = [dofmatY; dofmatyloc(:)];
				kmat = [kmat; K_up(:)];
				
				%fluxes
				fvec = [fvec; f_p];
				dofvec = [dofvec; dofsPD];
				
				[dofmatxloc,dofmatyloc] = ndgrid(dofsPD,dofsXY);
				dofmatX = [dofmatX; dofmatxloc(:)];
				dofmatY = [dofmatY; dofmatyloc(:)];
				kmat = [kmat; K_pu(:)];
				
				[dofmatxloc,dofmatyloc] = ndgrid(dofsPD,dofsPD);
				dofmatX = [dofmatX; dofmatxloc(:)];
				dofmatY = [dofmatY; dofmatyloc(:)];
				kmat = [kmat; K_pp(:)];
            end 

            physics.fint = physics.fint + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint), 1);
            physics.K = physics.K + sparse(dofmatX, dofmatY, kmat, length(physics.fint),length(physics.fint));

            obj.hmelt = hmeltNew;
			obj.QMeltTot = obj.QMeltTotOld + Q_Update;
			obj.qxSaved = qx_ips;
			obj.qMeltTot = obj.qMeltTotOld + MeltProd;
            
            tElapsed = toc(t);
            fprintf("        (Assemble time:"+string(tElapsed)+")\n");
		end

		function [qx, hmelt, h, derivs, Qprod] = get_qx_hmelt_h(obj, dp_dx, u, tfr, tOld, dt, hMelt_hist, initGuess, Ice_temp)
			% if (u<max(0,hMelt_hist))
			% 	u=max(0,hMelt_hist);
			% end

			stop = false;
			sol = initGuess;
			%ext = [dpdx; u; dt];

			it=0;

			[f, C, D] = obj.getFracK(sol, dp_dx, u, tfr, tOld, dt, hMelt_hist, Ice_temp);
			while (stop == false)
				if (sum(isnan(sol)+isinf(sol))>0)
					something wrong here
				end

				[E,F,G] = equilibrate(C);

				dsol = -G*((F*E*C*G)\(F*E*f));
				%dsol = -C\f;
				sol = sol+dsol;

				fOld = f;
				[f, C, D] = obj.getFracK(sol, dp_dx, u, tfr, tOld, dt, hMelt_hist, Ice_temp);

				%linesearch
				LineSearch = -fOld'*dsol/((f'-fOld')*dsol);
				LineSearch = min(max(LineSearch, 1e-1), 1);

				sol = sol-(1-LineSearch)*dsol;

				[f, C, D] = obj.getFracK(sol, dp_dx, u, tfr, tOld, dt, hMelt_hist, Ice_temp);

				% error
				if (it==0)
					err0 = sum(abs(dsol.*f))+eps;
					err = 1;
				else
					err = sum(abs(dsol.*f))/err0;
				end

				% convergence check
				if (err<1e-12 || err*err0<1e-18)
					stop = true;
				end
				it=it+1;
				if (it>1000)
					stop = true;
					fprintf("        (ip not converged)"+string(err)+"/"+string(err*err0)+"\n");
				end
			end

			%results
			qx = sol(1); 
			hmelt = sol(2);
			h = max(1e-4, sol(3));

			%derivatives
			[f, C, D, Qprod] = obj.getFracK(sol, dp_dx, u, tfr, tOld, dt, hMelt_hist, Ice_temp);

			%[E,F,G] = equilibrate(C);
			%drvs = -G*((F*E*C*G)\(F*E*D));
			drvs = -C\D;

			derivs.qx_dpdx = drvs(1,1);
			derivs.qx_du = drvs(1,2);

			derivs.hmelt_dpdx = drvs(2,1);
			derivs.hmelt_du = drvs(2,2);

			derivs.h_dpdx = drvs(3,1);
			derivs.h_du = drvs(3,2);
		end
        
		function [f, C, D, Qres] = getFracK(obj, sol, dp_dx, u, tfr, tOld, dt, hMelt_hist, Ice_temp)
			f = [0;0;0];
			C = zeros(3,3);
			D = zeros(3,3);

			smallH = 1e-4;

			dH = 1.0;
			hFlow = sol(3);
			if (hFlow<smallH)
				hFlow = smallH;
				dH = 0.0;
			end
			[qxFlow, dqx_dh, dqx_dpdx] = obj.getFlow(dp_dx, hFlow);

			ThermCap = obj.rho_ice*obj.melt_heat; 
			Qfreeze = 2*obj.k_ice^0.5*Ice_temp*(obj.rho_ice*obj.cp_ice)^0.5*pi^(-0.5)*(tOld+dt-tfr)^(-0.5);

			if (obj.melt)
				Qmelt = -sol(1)*dp_dx;
				dQmelt_ds = -dp_dx;
				dQmelt_dp = -sol(1);
			else
				Qmelt = 0;
				dQmelt = 0;
			end

			Q = Qfreeze + Qmelt;
			Qres(1) = Qfreeze;
			Qres(2) = Qmelt;
			Qres(3) = -(sol(2)-hMelt_hist)*ThermCap/dt;

			f(1) = sol(1)-qxFlow; 
			f(2) = (sol(2)-hMelt_hist)*ThermCap/dt-Q;
			f(3) = sol(3) - (sol(2)+u);

			% qx, hmelt, h
			C(1,1) = 1; 
			C(1,2) = 0;
			C(1,3) = -dqx_dh*dH;
			C(2,1) = -dQmelt_ds;
			C(2,2) = ThermCap/dt;
			C(2,3) = 0;
			C(3,1) = 0; 
			C(3,2) = -1;
			C(3,3) = 1;

% 			if (sol(3)<smallH)
% 				f(3) = f(3) + 1e5*(sol(3)-smallH)^2;
% 				C(3,3) = C(3,3) + 2*1e5*(sol(3)-smallH);
% 			end

			% dpdx, ujump
			D(1,1) = -dqx_dpdx;
			D(1,2) = 0;
			D(2,1) = -dQmelt_dp;
			D(2,2) = 0;
			D(3,1) = 0;
			D(3,2) = -1;
		end

		function [qxFlow, dqx_dh, dqx_dpdx] = getFlow(obj, dp_dx, h)
			if (h<0) 
				h=0; 
			end
        	if (obj.FlowModel == "CubicLaw")
            	qxFlow = -h^3/(12*obj.visc)*dp_dx;
            	dqx_dh = -3*(h)^2/(12*obj.visc)*dp_dx;
            	dqx_dpdx = -(h)^3/(12*obj.visc);
        	end
        	if (obj.FlowModel == "FrictionFactor")
            	k = 1e-2; f0 = 0.143;
				preFac = obj.rho_water^(-0.5)*k^(-1/6)*sqrt(4)*f0^(-0.5);
            	qxFlow   = -preFac    *(max(1,abs(dp_dx)))^(-0.5)*dp_dx*h^(5/3);
            	dqx_dh   = -preFac*5/3*(max(1,abs(dp_dx)))^(-0.5)*dp_dx*h^(2/3);
            	dqx_dpdx = -preFac*0.5*(max(1,abs(dp_dx)))^(-0.5)*h^(5/3);
        	end
		end

        function Nd = getNd(~, vals)
            cp_count = length(vals);
            Nd = zeros(2, cp_count*2);
            for ii = 1:cp_count 
				%dx
				Nd(1, ii) = vals(ii);
				Nd(1, ii + cp_count) = -vals(ii);

				%dy
				Nd(2, ii + 2*cp_count) = vals(ii);
				Nd(2, ii + 3*cp_count) = -vals(ii);
            end
        end

        function plotTotHeight(obj, physics, defscale)
            while (size(obj.hmelt, 1)<size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1))
                obj.hmeltOld(end+1,:) = 0.0;
                obj.hmelt(end+1,:) = 0.0;
            end

            ipc = obj.mesh.ipcount1D;
            Svec = physics.StateVec;

            lines = []; lines_melt = [];
            lines2 = []; lines2_melt = [];

            for n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                Elem_Nodes   = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, G, w]    = obj.mesh.getVals(obj.myGroupIndex, n_el);
                [nvec, tvec] = obj.mesh.getNormals(obj.myGroupIndex, n_el);
                
                dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                dofsXY = [dofsX; dofsY];

                X = Svec(dofsX);
                Y = Svec(dofsY);
                XY = [X;Y];

                xy = obj.mesh.getIPCoords(obj.myGroupIndex, n_el);
                for ip=1:ipc
                    if (abs(nvec(ip,1))>0.1)
                        hm(1) = 0.5*obj.hmelt(n_el,ip);
                        hm(2) = hm(1);
                        sc = [defscale, defscale];
                    else
                        hm(1) = 0.0;
                        hm(2) = obj.hmelt(n_el,ip);
                        sc = [defscale, defscale];
                    end
                    coords(ip,1) = xy(1, ip)+N(ip,:)*X(1:3)*sc(1);
                    coords(ip,2) = xy(2, ip)+N(ip,:)*Y(1:3)*sc(2);
                    coords2(ip,1) = xy(1, ip)+N(ip,:)*X(4:6)*sc(1);
                    coords2(ip,2) = xy(2, ip)+N(ip,:)*Y(4:6)*sc(2);

                    coords_melt(ip,:) = coords(ip,:) + nvec(ip,:)*hm(1).*sc; 
                    coords_melt2(ip,:) = coords2(ip,:) - nvec(ip,:)*hm(2).*sc; 
                end
                lines(n_el,:,:) = coords;
                lines2(n_el,:,:)= coords2;
                lines_melt(n_el,:,:) = coords_melt;
                lines2_melt(n_el,:,:) = coords_melt2;
            end
            hold on
            plot(squeeze(lines(:,:,1))', squeeze(lines(:,:,2))', 'k');
            plot(squeeze(lines2(:,:,1))', squeeze(lines2(:,:,2))', 'k');
            plot(squeeze(lines_melt(:,:,1))', squeeze(lines_melt(:,:,2))', 'r');
            plot(squeeze(lines2_melt(:,:,1))', squeeze(lines2_melt(:,:,2))', 'r');
		end

    end
end

