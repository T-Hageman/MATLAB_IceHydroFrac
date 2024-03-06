classdef Physics < handle
    %PHYSICS Class that handles the definitions and assembly of stiffness
	%matrices, constraints, and statevectors
    
    properties
        mesh
        models
        dofSpace
        
        K
        fint
        StateVec
        StateVec_Old
		VAvec
		VAvecOld
        time
        dt
        
        condofs
        convals
        convals_corr
        conMat
        unconMat
		nonz
    end
    
    methods
        PlotNodal(obj, dofName, dispscale, plotloc, zscale)
        PlotIP(obj, varName, plotloc)
        
        function obj = Physics(mesh, inputs, dt0)
            obj.mesh = mesh;
            obj.dofSpace = DofSpace(obj.mesh);
			obj.time = 0;
            
            for i=1:length(inputs)
                f = str2func(inputs{i}.type);
                obj.models{i} = f(mesh, obj, inputs{i});
            end
            
            obj.dt = dt0;
            
            dofcount = obj.dofSpace.NDofs;
            obj.StateVec = zeros(dofcount, 1);
            obj.StateVec_Old = obj.StateVec;
			obj.nonz = 0;

			obj.VAvec = zeros(dofcount,2);
			obj.VAvecOld = zeros(dofcount,2);
            
            obj.K = sparse(dofcount, dofcount);
            obj.fint = zeros(dofcount,1);
        end
        
        function Assemble(obj, cdt)
            dofcount = obj.dofSpace.NDofs;

			obj.condofs = [];
			obj.convals = [];

			if isempty(obj.K)

			else
				obj.nonz = round(nnz(obj.K));
			end
            obj.K = spalloc(dofcount, dofcount, obj.nonz);
            obj.fint = zeros(dofcount, 1);

            disp("    Assembling:")
            for m=1:length(obj.models)
                obj.models{m}.getKf(obj);
			end
        end
       
        function Commit(obj, commit_type)
            for m=1:length(obj.models)
                obj.models{m}.Commit(obj, commit_type);
            end
            
            if (commit_type == "Timedep")
				obj.time = obj.time + obj.dt;

                obj.StateVec_Old = obj.StateVec;
				obj.VAvecOld = obj.VAvec;
				obj.K = [];
            end
        end
        
        function anyIrr = Irreversibles(obj)
            anyIrr = false;
            for m=1:length(obj.models)
                anyIrr = anyIrr + obj.models{m}.Irreversibles(obj);
            end
        end
            
        function Constrain(obj)
            obj.convals_corr = obj.convals - obj.StateVec(obj.condofs);
            basemat = speye(size(obj.K));
            obj.unconMat = basemat;
            obj.unconMat(:, obj.condofs) = [];
            obj.conMat = basemat(:, obj.condofs);

            obj.fint = obj.unconMat'*obj.fint + obj.unconMat'*obj.K*obj.conMat*obj.convals_corr;
            obj.K    = obj.unconMat'*obj.K*obj.unconMat;
        end
        
        function dx = Update(obj, dx)
			dx_uncon = obj.unconMat*dx;

            obj.StateVec = obj.StateVec + dx_uncon + obj.conMat*obj.convals_corr;
            dx = obj.unconMat'*dx_uncon;
        end
        
        function info = Request_Info(obj, var, elems, loc)
            info = false;
            for m=1:length(obj.models)
                [hasInfo, provided] = obj.models{m}.Provide_Info(obj, var, elems, loc);
                if (hasInfo)
                    info = provided;
                end
            end   
        end

    end
end

