classdef Inertia < BaseModel
    %LINEARELASTIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        density
		beta
		gamma
        myK
    end
    
    methods
        function obj = Inertia(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "Inertia";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType({"dx","dy"});
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            %% constuct plane-strain elastic stiffness matrix
            obj.density = inputs.density;
			obj.beta = inputs.beta;
			obj.gamma = inputs.gamma;
            
            obj.myK = sparse(0,0);
        end
        
        function getKf(obj, physics)
            fprintf("        Inertia get Matrix:")
            t = tic;
            
            dt = physics.dt;

            if (length(obj.myK) == length(physics.fint))
                recalc = false;
            else
                recalc = true;
            end
            
			physics.VAvec(:,1) = obj.gamma/(obj.beta*dt)*(physics.StateVec-physics.StateVec_Old)-(obj.gamma/obj.beta-1)*physics.VAvecOld(:,1)-dt*(obj.gamma/obj.beta/2-1)*physics.VAvecOld(:,2);
            physics.VAvec(:,2) = 1/(obj.beta*dt^2)*(physics.StateVec-physics.StateVec_Old)-1/(obj.beta*dt)*physics.VAvecOld(:,1)-(1/obj.beta/2-1)*physics.VAvecOld(:,2);

            if (recalc == true)
                dofmatX = [];
                dofmatY = [];
                kmat = [];
                parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                    [N, ~, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
                    xy = obj.mesh.getIPCoords(obj.myGroupIndex, n_el);

                    dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                    dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                    dofsXY = [dofsX; dofsY];

                    %X = physics.StateVec(dofsX);
                    %Y = physics.StateVec(dofsY);
                    %XY = [X;Y];

                    %f_el = zeros(length(dofsXY), 1);
                    K_el = zeros(length(dofsX));
                    for ip=1:length(w)
						K_el = K_el + w(ip)*obj.density*N(ip,:)'*N(ip,:);
                    end

                    [dofmatxloc,dofmatyloc] = ndgrid(dofsX,dofsX);
                    dofmatX = [dofmatX; dofmatxloc(:)];
                    dofmatY = [dofmatY; dofmatyloc(:)];
                    kmat = [kmat; K_el(:)];

					[dofmatxloc,dofmatyloc] = ndgrid(dofsY,dofsY);
                    dofmatX = [dofmatX; dofmatxloc(:)];
                    dofmatY = [dofmatY; dofmatyloc(:)];
                    kmat = [kmat; K_el(:)];
                end 
                obj.myK = sparse(dofmatX, dofmatY, kmat, length(physics.fint),length(physics.fint));
            end
            physics.fint = physics.fint + obj.myK*physics.VAvec(:,2);
            physics.K = physics.K + obj.myK*(1/obj.beta/dt^2);
            
            tElapsed = toc(t);
            fprintf("            (Assemble time:"+string(tElapsed)+")\n");
		end
    end
end

