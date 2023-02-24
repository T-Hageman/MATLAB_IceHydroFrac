classdef ViscoElastic < BaseModel
    %ViscoElastic Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        poisson
        young
        D_el
        D_el2
        Hmatswitch
        myK
		A
		n
		P
		FVisc;
		strain_visc;
		strain_viscOld;
    end
    
    methods
        function obj = ViscoElastic(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "ViscoElastic";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType({"dx","dy"});
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            %% constuct plane-strain elastic stiffness matrix
            obj.poisson = inputs.poisson;
            obj.young = inputs.young;
            obj.Hmatswitch = inputs.Hmatswitch;
			obj.A = inputs.A;
			obj.n = inputs.n;
            
            D_el = zeros(4,4);
            a = obj.young(1) / ((1.0 + obj.poisson(1)) * (1.0 - 2.0*obj.poisson(1)));

            D_el(1, 1) = a * (1.0 - obj.poisson(1));
            D_el(1, 2) = a * obj.poisson(1);
            D_el(1, 3) = a * obj.poisson(1);
            D_el(2, 1) = a * obj.poisson(1);
            D_el(2, 2) = a * (1.0 - obj.poisson(1));
            D_el(2, 3) = a * obj.poisson(1);
            D_el(3, 1) = a * obj.poisson(1);
            D_el(3, 2) = a * obj.poisson(1);
            D_el(3, 3) = a * (1.0 - obj.poisson(1));
            D_el(4, 4) = a * 0.5 * (1.0 - 2.0*obj.poisson(1));
            
            obj.D_el = D_el;
            
			D_el2 = zeros(4,4);
            a = obj.young(2) / ((1.0 + obj.poisson(2)) * (1.0 - 2.0*obj.poisson(2)));

            D_el2(1, 1) = a * (1.0 - obj.poisson(2));
            D_el2(1, 2) = a * obj.poisson(2);
            D_el2(1, 3) = a * obj.poisson(2);
            D_el2(2, 1) = a * obj.poisson(2);
            D_el2(2, 2) = a * (1.0 - obj.poisson(2));
            D_el2(2, 3) = a * obj.poisson(2);
            D_el2(3, 1) = a * obj.poisson(2);
            D_el2(3, 2) = a * obj.poisson(2);
            D_el2(3, 3) = a * (1.0 - obj.poisson(2));
            D_el2(4, 4) = a * 0.5 * (1.0 - 2.0*obj.poisson(2));
            obj.D_el2 = D_el2;

			obj.P = zeros(4,4); 
			obj.P(1,1)=2/3; obj.P(2,2)=2/3; obj.P(3,3) = 2/3; obj.P(4,4)=1;
			obj.P(2,1)=-1/3; obj.P(3,1)=-1/3; obj.P(1,2)=-1/3; obj.P(3,2)=-1/3; obj.P(1,3)=-1/3; obj.P(2,3)=-1/3;
            
            obj.myK = sparse(0,0);
			obj.strain_visc = zeros(size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount, 4);
			obj.strain_viscOld = obj.strain_visc;
        end
        
        function Commit(obj, physics, commit_type)
            if (commit_type == "Timedep")
                obj.strain_viscOld = obj.strain_visc;
			end
			if (commit_type == "Irrevirsibles")
                obj.UpdateVisc(physics);
			end
			if (commit_type == "StartIt")
				if (obj.A>0)
					obj.UpdateVisc(physics);
				end
            end
        end

		function [stress, strain_ve_new] = return_mapping(obj, strain_total, strain_ve_old, A, D, dt)
			sol(1:4,1) = D*(strain_total - strain_ve_old);
			sol(5:8,1) = strain_ve_old;

			[K,f] = obj.return_mapping_getKF(sol, strain_total, strain_ve_old, A, D, dt);
			conv = false;
			it = 0;
			while (conv == false)
				[E,F,G] = equilibrate(K);
				dsol = -G*((F*E*K*G)\(F*E*f));
				sol = sol +dsol;

				fOld = f;
				[K,f] = obj.return_mapping_getKF(sol, strain_total, strain_ve_old, A, D, dt);

				LineSearch = -fOld'*dsol/((f'-fOld')*dsol);
				LineSearch = min(max(LineSearch, 1e-1), 1);

				sol = sol-(1-LineSearch)*dsol;

				[K,f] = obj.return_mapping_getKF(sol, strain_total, strain_ve_old, A, D, dt);
				err = sum(abs(f(5:8)));
				%disp(err);

				it=it+1;
				if (it>100)
					conv = true;
					disp(err);
				end
				if (err<1e-6)
					conv = true;
				end
			end
			stress = sol(1:4);
			strain_ve_new = sol(5:8);

		end

		function [K,f] = return_mapping_getKF(obj, sol, strain_total, strain_ve_old, A, D, dt)
			f(1:4,1) = sol(1:4) - D*(strain_total-sol(5:8));
			f(5:8,1) = sol(5:8) - strain_ve_old - dt*A*(sol(1:4)'*obj.P'*obj.P*sol(1:4))^((obj.n-1)/2)*obj.P*sol(1:4);

			K(1:4,1:4) = eye(4,4); K(1:4,5:8) = D;
			K(5:8,5:8) = eye(4,4);
			K(5:8,1:4) = - dt*A*obj.n*(sol(1:4)'*obj.P'*obj.P*sol(1:4))^((obj.n-1)/2)*obj.P;
		end

		function UpdateVisc(obj, physics)
			fprintf("        ViscoElastic UpdateVisc:")
			t = tic;

			fvec = [];
            dofvec = [];
			SVec = physics.StateVec;

			if (physics.ArcTime.Enable)
                TDof = physics.dofSpace.getDofIndices(physics.tType, 1);
                dt = physics.StateVec(TDof) - physics.StateVec_Old(TDof);
			else
                dt = physics.dt;
            end
			strain_viscTosave = [];
			parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [~, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
				xy = obj.mesh.getIPCoords(obj.myGroupIndex, n_el);

               	dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                dofsXY = [dofsX; dofsY];

				X = SVec(dofsX);
				Y = SVec(dofsY);
                XY = [X;Y];

				f_el = zeros(length(dofsXY), 1);
				Strain_VP = zeros(length(w), 4);
				for ip=1:length(w)
                    B = obj.getB(G(ip,:,:));
                    strain = B*XY;
                        if (xy(2,ip) < obj.Hmatswitch)
                            Aloc = 0;
                            D = obj.D_el2;                       
                        else
                            Aloc = obj.A;
                            D = obj.D_el;
                        end
                    %stress = obj.D_el*(strain-squeeze(obj.strain_viscOld(n_el,ip,1:4)));

					%dVPStrain_dt = Aloc*(max(1e-6,0.5*stress'*obj.P'*obj.P*stress))^((obj.n-1)/2)*obj.P*stress;
					%Strain_VP(ip,1:4) = squeeze(obj.strain_viscOld(n_el,ip,1:4)) + dVPStrain_dt*dt;
					[~, Strain_VP(ip,1:4)] = obj.return_mapping(strain, squeeze(obj.strain_viscOld(n_el,ip,1:4)), Aloc, D, dt);

                    f_el = f_el - w(ip)*B'*D*Strain_VP(ip,1:4)';

					if (sum(isnan(Strain_VP(ip,1:4)))>0 || sum(abs(Strain_VP(ip,1:4)))>100)
						disp(Strain_VP);
						disp(stress);
						disp(XY);
					end
				end
				strain_viscTosave(n_el,:,:) = Strain_VP;
				
				fvec = [fvec; f_el];
                dofvec = [dofvec; dofsXY];
			end
			obj.strain_visc = strain_viscTosave;
			obj.FVisc = sparse(dofvec, 0*dofvec+1, fvec, physics.dofSpace.NDofs, 1);
			
            tElapsed = toc(t);
            fprintf("            (Assemble time:"+string(tElapsed)+")\n");
		end

        function getKf(obj, physics)
            fprintf("        ViscoElastic get Matrix:")
            t = tic;
            
            if (length(obj.myK) == length(physics.fint))
                recalc = false;
            else
                recalc = true;
            end
            
            
            if (recalc == true)
                dofmatX = [];
                dofmatY = [];
                kmat = [];
                parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                    [~, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
                    xy = obj.mesh.getIPCoords(obj.myGroupIndex, n_el);

                    dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                    dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                    dofsXY = [dofsX; dofsY];

                    %X = physics.StateVec(dofsX);
                    %Y = physics.StateVec(dofsY);
                    %XY = [X;Y];

                    %f_el = zeros(length(dofsXY), 1);
                    K_el = zeros(length(dofsXY));
                    for ip=1:length(w)
                        B = obj.getB(G(ip,:,:));
                        %strain = B*XY;

                        if (xy(2,ip) < obj.Hmatswitch)
                            %stress = obj.D_el2*strain;
                            %f_el = f_el + B'*stress*w(ip);
                            K_el = K_el + B'*obj.D_el2*B*w(ip);                       
                        else
                            %stress = obj.D_el*strain;
                            %f_el = f_el + B'*stress*w(ip);
                            K_el = K_el + B'*obj.D_el*B*w(ip);
                        end
                    end

                    [dofmatxloc,dofmatyloc] = ndgrid(dofsXY,dofsXY);
                    dofmatX = [dofmatX; dofmatxloc(:)];
                    dofmatY = [dofmatY; dofmatyloc(:)];
                    kmat = [kmat; K_el(:)];
                end 
                obj.myK = sparse(dofmatX, dofmatY, kmat, length(physics.fint),length(physics.fint));
            end
            physics.fint = physics.fint + obj.myK*physics.StateVec ;
			if (obj.A>0)
				physics.fint = physics.fint + obj.FVisc;
			end
            physics.K = physics.K + obj.myK;
            
            
            tElapsed = toc(t);
            fprintf("            (Assemble time:"+string(tElapsed)+")\n");
        end
        
        function B = getB(~, grads)
            cp_count = size(grads, 2);
            B = zeros(4, cp_count*2);
            for ii = 1:cp_count %using plane strain e_zz = 0
				%dx
				B(1, ii) = grads(1,ii, 1);
				B(4, ii) = grads(1,ii, 2);

				%dy
				B(2, ii + cp_count) = grads(1,ii, 2);
				B(4, ii + cp_count) = grads(1,ii, 1);
            end
        end

        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
           hasInfo = false;
           provided = [];
           
           if (var == "stresses" || var=="sxx" || var=="syy" || var=="szz" || var=="sxy" || var == "dstressesdxy" || var == "sdev" || var == "svol") %elem, ip, componenets
               hasInfo = true;
               if (var == "stresses")
                   provided = zeros(length(elems), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount, 4);
               elseif (var == "dstressesdxy")  %elem, ip, comp, dxy
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, elems(1));
                    provided = zeros(length(elems), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount, 4, length(Elem_Nodes)*2 );
               else
                   provided = zeros(length(elems), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount);
               end
                for el=1:length(elems)
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, elems(el));
                    [~, G, w] = obj.mesh.getVals(obj.myGroupIndex, elems(el));
                    xy = obj.mesh.getIPCoords(obj.myGroupIndex, elems(el));

                    dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                    dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);

                    X = physics.StateVec(dofsX);
                    Y = physics.StateVec(dofsY);
                    XY = [X;Y];

                    for ip=1:length(w)
                        B = obj.getB(G(ip,:,:));
                        strain = B*XY-squeeze(obj.strain_viscOld(elems(el),ip,1:4));
                        
                        if (xy(2,ip) < obj.Hmatswitch)
                            stress = obj.D_el2*strain; 
                            dStress = obj.D_el2*B;
                        else
                            stress = obj.D_el*strain;
                            dStress = obj.D_el*B;
                        end

                        if (loc == "Interior")
                            switch var 
                                case "stresses"
                                    provided(el, ip, :) = stress;
                                case "sxx"
                                    provided(el, ip) = stress(1);
                                case "syy"
                                    provided(el, ip) = stress(2);
                                case "szz"
                                    provided(el, ip) = stress(3);
                                case "sxy"
                                    provided(el, ip) = stress(4);
                                case "dstressesdxy"
                                    provided(el, ip, :, :) = dStress;
								case "sdev"
                                    provided(el, ip) = sqrt(stress'*obj.P'*obj.P*stress);
								case "svol"
                                    provided(el, ip) = sqrt(stress'*(eye(4)-obj.P)'*(eye(4)-obj.P)*stress);
                            end
                        end
                        
                    end
                end
           end
        end
        
    end
end

