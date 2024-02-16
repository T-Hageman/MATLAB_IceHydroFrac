classdef LakeBoundary < BaseModel
    %FractureFluid Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
		dispDofIndices
        
        p0
		kdummy

        QTotal
		qCurrent
		dxCurrent
		dyCurrent
        hMeltHist

		refY

		SurfaceGroup
		surface_dX
		surface_dY
		surface_X
		surface_Y
    end
    
    methods
        function obj = LakeBoundary(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "LakeBoundary";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
			obj.SurfaceGroup = obj.mesh.getGroupIndex(inputs.Surface);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType({"pd"});
			obj.dispDofIndices = obj.dofSpace.addDofType({"dx","dy"});
            
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.Elementgroups{obj.myGroupIndex}.Elems(1,1));
            
            obj.p0 = inputs.p0;
            obj.kdummy = inputs.dummy;

			obj.QTotal = 0.0;
			obj.qCurrent = 0.0;
			obj.dxCurrent = 0.0;
			obj.dyCurrent = 0.0;
			obj.hMeltHist = 0.0;
			obj.refY = 0;
        end
        
        function getKf(obj, physics)
            fprintf("        LakeBoundary get Matrix:")
            t = tic;
            
            Elem_Node   = obj.mesh.Elementgroups{obj.myGroupIndex}.Elems(1,1);
            dofsPD = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Node);
            PD = physics.StateVec(dofsPD);

            time = physics.time+physics.dt;
            ptarget = obj.p0;
            
            q = obj.kdummy * (ptarget-PD); 
            k = -obj.kdummy;
            
            physics.fint(dofsPD) = physics.fint(dofsPD) + q;
            physics.K(dofsPD, dofsPD) = physics.K(dofsPD, dofsPD) + k;

            tElapsed = toc(t);
            fprintf("        (Assemble time:"+string(tElapsed)+")\n");
		end

		function updateSurfaceElevation(obj, physics)
			obj.surface_dX = [];
			obj.surface_dY = [];
			obj.surface_X  = [];
			obj.surface_Y  = [];

			for n_el=1:size(obj.mesh.Elementgroups{obj.SurfaceGroup}.Elems, 1)
                Elem_Nodes = obj.mesh.getNodes(obj.SurfaceGroup, n_el);
                [N, ~, w] = obj.mesh.getVals(obj.SurfaceGroup, n_el);
                xy = obj.mesh.getIPCoords(obj.SurfaceGroup, n_el);
	
                dofsX = obj.dofSpace.getDofIndices(obj.dispDofIndices(1), Elem_Nodes);
                dofsY = obj.dofSpace.getDofIndices(obj.dispDofIndices(2), Elem_Nodes);
                X = physics.StateVec(dofsX);
                Y = physics.StateVec(dofsY);
				for ip=1:length(w)
					obj.surface_X = [obj.surface_X xy(1,ip)];
					obj.surface_Y = [obj.surface_Y xy(2,ip)];

					obj.surface_dX = [obj.surface_dX N(ip,:)*X];
					obj.surface_dY = [obj.surface_dY N(ip,:)*Y];
				end

			end

		end
       
        function Commit(obj, physics, commit_type)
            if (commit_type == "Timedep")
                Elem_Node   = obj.mesh.Elementgroups{obj.myGroupIndex}.Elems(1,1);
                dofsPD = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Node);
				dofsX = obj.dofSpace.getDofIndices(obj.dispDofIndices(1), Elem_Node);
                dofsY = obj.dofSpace.getDofIndices(obj.dispDofIndices(2), Elem_Node);

                X = physics.StateVec(dofsX);
                Y = physics.StateVec(dofsY);
                PD = physics.StateVec(dofsPD);

                q = obj.kdummy * (obj.p0-PD);
                obj.QTotal(end+1) = obj.QTotal(end)+q*physics.dt;
				obj.qCurrent(end+1) = q;
				obj.dxCurrent(end+1)= X;
				obj.dyCurrent(end+1)= Y-obj.refY;
				obj.hMeltHist(end+1) = physics.models{5}.hmelt(1,1);

                fprintf("Total fluid flow:"+string(obj.QTotal(end))+"\n");

				if (physics.time < 0)
					obj.refY = Y(1);
				end
            end
        end
        
    end
end

