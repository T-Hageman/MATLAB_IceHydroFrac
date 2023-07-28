function [direction, elem, ips] = getNextFracture(obj)

	downGroup = 6;
	leftGroup = 7;
	rightGroup= 8;

	if (length(obj.Elementgroups{downGroup}.Elems)>0)	%propagate downwards
		direction = "Vertical";
		Nds = obj.Elementgroups{downGroup}.Elems(end,:);
		middleNode = Nds(2);
		[elem_ind,~] = find(obj.Elementgroups{1}.Elems==middleNode);
		xc = obj.Nodes(obj.Elementgroups{1}.Elems(elem_ind,5),1);
		if (xc(1)>xc(2)) %evaluating based on right element
			elem = elem_ind(1);
		else
			elem = elem_ind(2);
		end
	else %propagate sideways
		direction = "Horizontal";
		Nds = obj.Elementgroups{rightGroup}.Elems(1,:);
		middleNode = Nds(2);
		[elem_ind,~] = find(obj.Elementgroups{1}.Elems==middleNode);
		yc = obj.Nodes(obj.Elementgroups{1}.Elems(elem_ind,5),2);
		if (yc(1)>yc(2)) %evaluating based on top element
			elem = elem_ind(1);
		else
			elem = elem_ind(2);
		end

	end

	Elem_Nodes = obj.Elementgroups{1}.Elems(elem, :);
	if (middleNode == Elem_Nodes(2)) %bottom of element
		ips = [2:obj.ipcount1D-1];
	elseif (middleNode == Elem_Nodes(4)) %left of element
		ips = [obj.ipcount1D+1:obj.ipcount1D:obj.ipcount1D*(obj.ipcount1D-2)+1];
	elseif (middleNode == Elem_Nodes(6)) %right of element
		ips = [2*obj.ipcount1D:obj.ipcount1D:obj.ipcount1D*(obj.ipcount1D-1)];
	else %top of element
		ips = [obj.ipcount1D*(obj.ipcount1D-1)+2:obj.ipcount1D*obj.ipcount1D-1];
	end

	if (direction == "Horizontal")
		ipc = obj.getIPCoords(1, elem);
		ipc = ipc(1,ips);
		[~,ind] = min(ipc);
		%[~,ind] = sort(ipc);
		ind=ind(1);
		ip = ips(ind);
	else
		ipc = obj.getIPCoords(1, elem);
		ipc = ipc(2,ips);
		[~,ind] = max(ipc);
		%[~,ind] = sort(ipc);
		ind = ind(1);
		ip = ips(ind);
	end

end

