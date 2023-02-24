function check(obj)

    disp("Checking mesh");

    newArea = 0.0*obj.Area;
    for g=1:length(obj.Elementgroups)
		ThisArea = 0;
        parfor ielem = 1:size(obj.Elementgroups{g}.Elems, 1)
            [N, G, w] = obj.getVals(g, ielem);
            ThisArea = ThisArea + sum(w);
        end
        newArea(g) = ThisArea;
        disp("    "+obj.Elementgroups{g}.name+ ": Old Area: "+string(obj.Area(g))+"    new Area: "+string(newArea(g)))
    end


    obj.Area = newArea;
end

