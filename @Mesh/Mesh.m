classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Area
        Elementgroups
        Nodegroups
        Nodes
        ipcount1D
        zeroWeight
        Tjunction
    end
    
    methods
        [Nodes, Egroups, Ngroups, Area, rect] = TFrac_Generator(obj, props);
        f1 = plot(obj, plotnodes, plotelems, plotnames, plotnodenames);
        xy = getIPCoords(obj, group, elems);
        check(obj);
        [N, G, w] = getVals(obj, group, elem);
        [n, t] = getNormals(obj, group, elem);
        elems = getConnected(obj, groupnum, node);
        Propagate_Disc(obj, physics, tipnode, direction);
		Propagate_Disc_New(obj, physics);

		[direction, elem, ips] = getNextFracture(obj);
        
        function obj = Mesh(inProps)
            %MESH Construct an instance of this class
            %   Detailed explanation goes here
			if (inProps.type == "File")
				[obj.Nodes, obj.Elementgroups, obj.Nodegroups, obj.Area, rectangular] = obj.MeshFromFile(inProps);
                
                obj.Tjunction = [];
			end

            obj.zeroWeight = inProps.zeroWeight;
            obj.ipcount1D = inProps.ipcount1D;
            for g=1:length(obj.Elementgroups)
               f = str2func(obj.Elementgroups{g}.type);
               obj.Elementgroups{g}.ShapeFunc = f(obj.ipcount1D, rectangular, obj.zeroWeight);
            end
            if (obj.zeroWeight)
                obj.ipcount1D = obj.ipcount1D+2;
            end
            
        end
        
        function groupIndex = getGroupIndex(obj, groupname)
            groupIndex = -1;
            for i=1:length(obj.Elementgroups)
               if (groupname == obj.Elementgroups{i}.name)
                   groupIndex = i;
               end
            end
            if (groupIndex == -1)
                msg = "Elementgroup "+groupname+" is not defined.";
                error(msg);
            end
        end
        
        function myNode = getNodes(obj, group, elem)
            myNode = obj.Elementgroups{group}.Elems(elem,:);
        end
        
        function nds = GetAllNodesForGroup(obj, groupindex)
            myElems = obj.Elementgroups{groupindex}.Elems;
            nds = unique(reshape(myElems, [], 1));
        end
        
        
    end
end

