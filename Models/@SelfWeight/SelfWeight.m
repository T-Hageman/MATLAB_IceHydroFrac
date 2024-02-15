classdef SelfWeight < BaseModel
    %SelfWeight Adds the gravity term to the momentum balance. Inputs:
	% physics_in{3}.type = "SelfWeight";
	% physics_in{3}.Egroup = "Internal";
	% physics_in{3}.density = [910; 2500];	%Density [kg/m^3]
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        density
        gravity
        fvec
    end
    
    methods
        function obj = SelfWeight(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "SelfWeight";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType({"dy"});
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            %% other properties
            obj.density = inputs.density;
            obj.gravity = -9.81;
            obj.fvec = [];
        end
        
        function getKf(obj, physics)
            fprintf("        SelfWeight get Matrix:")
            t = tic;
            
            if (length(obj.fvec) == length(physics.fint))
                recalc = false;
            else
                recalc = true;
            end
            
            if (recalc == true)
                fvecl = [];
                dofvec = [];

                parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                    [N, ~, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
					xy = obj.mesh.getIPCoords(obj.myGroupIndex, n_el);

                    dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices, Elem_Nodes);

                    f_el = zeros(length(dofsY), 1);
                    for ip=1:length(w)
						if (xy(2,ip) < 0)
							rho = obj.density(2); %rock
						else
							rho = obj.density(1); %ice
						end
                        f_el = f_el - rho*obj.gravity*N(ip,:)'*w(ip);
                    end             
                    fvecl = [fvecl; f_el];
                    dofvec = [dofvec; dofsY];
                end 

                obj.fvec = sparse(dofvec, 0*dofvec+1, fvecl, length(physics.fint), 1);
            end
            physics.fint = physics.fint + obj.fvec;
            
            tElapsed = toc(t);
            fprintf("            (Assemble time:"+string(tElapsed)+")\n");
        end
        
    end
end

