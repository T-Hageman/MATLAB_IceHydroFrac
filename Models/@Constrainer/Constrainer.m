classdef Constrainer < BaseModel
    %CONSTRAINER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        myName
        mesh
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        conVal
    end
    
    methods
        function obj = Constrainer(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "Constrainer";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType(inputs.dofs);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            obj.conVal = inputs.conVal;
        end
        
        function getKf(obj, physics)
            allNodes = obj.mesh.GetAllNodesForGroup(obj.myGroupIndex);
            
            for i=1:length(obj.dofTypeIndices)
                newcons = obj.dofSpace.getDofIndices(obj.dofTypeIndices(i), allNodes);
                
                physics.condofs = [physics.condofs; newcons];
                physics.convals = [physics.convals; newcons*0+obj.conVal(i)];
            end
        end
    end
end

