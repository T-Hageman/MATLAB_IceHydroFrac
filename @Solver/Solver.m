classdef Solver
    %SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        physics
        
        maxIt
        Conv
        tiny
        linesearch
        linesearchLims
    end
    
    methods
        Solve(obj);
        
        function obj = Solver(physics, inputs)
            obj.physics = physics;
            
            obj.maxIt = inputs.maxIt;
            obj.Conv = inputs.Conv;
            obj.tiny = inputs.tiny;
            obj.linesearch = inputs.linesearch;
            if (obj.linesearch)
                obj.linesearchLims = inputs.linesearchLims;
            end
        end
        

    end
end

