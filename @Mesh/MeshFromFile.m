function [Nodes, Elementgroups, Nodegroups, Area, rect] = MeshFromFile(obj, props)

	input_fileName = props.FileName;
	textfile = fopen(input_fileName);
	
	%% preamble
	stop = false;
	while stop==false
		line = fgetl(textfile);
		if (line == "# Mesh vertex coordinates")
			stop = true;
		end
	end

	%% nodes
	stop = false;
	c=0;
	while stop==false
		line = fgetl(textfile);
		if (line == "")
			stop = true;
		else
			c=c+1;
			N = str2num(line);
			Nodes(c,:) = N;
		end
	end	
	
	%% boundaries
	stop = false;
	while stop==false
		line = fgetl(textfile);
		if (line == "3 # number of vertices per element")
			stop = true;
			line = fgetl(textfile);
			line = fgetl(textfile);
		end
	end

	stop = false;
	c=0;
	while stop==false
		line = fgetl(textfile);
		if (line == "")
			stop = true;
			line = fgetl(textfile);
			line = fgetl(textfile);
		else
			c=c+1;
			B = str2num(line);
			B_Elems(c,:) = B+1;
		end
	end		

	stop = false;
	c=0;
	while stop==false
		line = fgetl(textfile);
		if (line == "")
			stop = true;
			line = fgetl(textfile);
			line = fgetl(textfile);
		else
			c=c+1;
			G = str2double(line);
			G_Elems(c) = G+1;
		end
	end			

	num_groups = max(G_Elems);
	c = zeros(num_groups,1);
	for i=1:length(G_Elems)
		grp = G_Elems(i);
		c(grp,:) = c(grp,:)+1;
		BE{grp}(c(grp),:)=B_Elems(i,:);
	end

	if false
    	Elementgroups{2}.name = "Bottom";
    	Elementgroups{2}.type = "L3B";
		Elementgroups{2}.Elems = [BE{2}; BE{7}];
	
    	Elementgroups{3}.name = "Top";
    	Elementgroups{3}.type = "L3B";
		Elementgroups{3}.Elems = [BE{5}; BE{10}];
	
		Elementgroups{4}.name = "Left";
    	Elementgroups{4}.type = "L3B";
		Elementgroups{4}.Elems = [BE{1}; BE{3}];
	
		Elementgroups{5}.name = "Right";
    	Elementgroups{5}.type = "L3B";
		Elementgroups{5}.Elems = [BE{11}; BE{12}];
		
		Elementgroups{6}.name = "FDown";
    	Elementgroups{6}.type = "L3B";
		Elementgroups{6}.Elems = BE{8};
	
		Elementgroups{7}.name = "FLeft";
    	Elementgroups{7}.type = "L3B";
		Elementgroups{7}.Elems = BE{4};	
	
		Elementgroups{8}.name = "FRight";
    	Elementgroups{8}.type = "L3B";
		Elementgroups{8}.Elems = BE{9};	
	else
		bg = [];
		for i=1:length(BE)
			ycoords = Nodes(BE{i},2);
			if (sum(ycoords<-150)==length(ycoords))
				bg(end+1) = i;
			end
		end
    	Elementgroups{2}.name = "Bottom";
    	Elementgroups{2}.type = "L3B";
		Elementgroups{2}.Elems = [BE{bg(1)}; BE{bg(2)}];
	
		tg = [];
		for i=1:length(BE)
			ycoords = Nodes(BE{i},2);
			if (sum(ycoords>200)==length(ycoords))
				tg(end+1) = i;
			end
		end

     	Elementgroups{3}.name = "Top";
     	Elementgroups{3}.type = "L3B";
 		Elementgroups{3}.Elems = [BE{tg(1)}; BE{tg(2)}];

		lg = [];
		for i=1:length(BE)
			xcoords = Nodes(BE{i},1);
			if (sum(xcoords<-500)==length(xcoords))
				lg(end+1) = i;
			end
		end	
	
		Elementgroups{4}.name = "Left";
    	Elementgroups{4}.type = "L3B";
		Elementgroups{4}.Elems = [BE{lg(1)}; BE{lg(2)}];

		rg = [];
		for i=1:length(BE)
			xcoords = Nodes(BE{i},1);
			if (sum(xcoords>500)==length(xcoords))
				rg(end+1) = i;
			end
		end	
	
		Elementgroups{5}.name = "Right";
    	Elementgroups{5}.type = "L3B";
		Elementgroups{5}.Elems = [BE{rg(1)}; BE{rg(2)}];

		UsedGroups = [lg rg tg bg];
		AvailableGroups = 1:length(BE);
		AvailableGroups(UsedGroups) = [];

		for i=1:length(AvailableGroups)
			XMean(i) = mean(Nodes(BE{AvailableGroups(i)},1));
			YMean(i) = mean(Nodes(BE{AvailableGroups(i)},2));
		end

		[~,downGroup] = max(YMean);
		downGroup = AvailableGroups(downGroup);

		[~,leftGroup] = min(XMean);
		leftGroup = AvailableGroups(leftGroup);

		[~,rightGroup] = max(XMean);
		rightGroup = AvailableGroups(rightGroup);

		
		Elementgroups{6}.name = "FDown";
    	Elementgroups{6}.type = "L3B";
		Elementgroups{6}.Elems = BE{downGroup};
	
		Elementgroups{7}.name = "FLeft";
    	Elementgroups{7}.type = "L3B";
		Elementgroups{7}.Elems = BE{leftGroup};	
	
		Elementgroups{8}.name = "FRight";
    	Elementgroups{8}.type = "L3B";
		Elementgroups{8}.Elems = BE{rightGroup};	
	end

	for EG=2:length(Elementgroups)
		clear firstCoords Elems_Temp

		Elems_Old = Elementgroups{EG}.Elems;
		for i=1:size(Elems_Old,1)
			crds = Nodes(Elems_Old(i,:),:);
			if (abs(crds(1,1)-crds(2,1))<1e-6) %vertical
				[~, ord] = sort(crds(:,2));
			else
				[~, ord] = sort(crds(:,1));
			end
			Elems_Temp(i,:) = Elems_Old(i,ord);
			firstCoords(i,:) = Nodes(Elems_Temp(i,1),:);
		end

		if (abs(firstCoords(1,1)-firstCoords(2,1))<1e-6)
			[~, ord] = sort(firstCoords(:,2));
		else
			[~, ord] = sort(firstCoords(:,1));
		end
		Elementgroups{EG}.Elems = Elems_Temp(ord,:);
	end

	%% interior
    Elementgroups{1}.name = "Internal";
    Elementgroups{1}.type = "Q9B";

	stop = false;
	while stop==false
		line = fgetl(textfile);
		if (line == "9 # number of vertices per element")
			stop = true;
			line = fgetl(textfile);
			line = fgetl(textfile);
		end
	end

	stop = false;
	c=0;
	while stop==false
		line = fgetl(textfile);
		if (line == "")
			stop = true;
			line = fgetl(textfile);
			line = fgetl(textfile);
		else
			c=c+1;
			G = str2num(line);
			Elementgroups{1}.Elems(c,:) = G([1;5;2;6;7;8;3;9;4])+1;
		end
	end		

	%%Insert Interfaces
	if true
		Elementgroups{9}.name = "Fracture";
    	Elementgroups{9}.type = "LI6B";
		for i=1:props.nfrac
			Nds = Elementgroups{6}.Elems(end,:);
			nnodes = length(Nodes);
	
			ToSplit = Nds(2:3);
			NewNodeNums = [nnodes+1; nnodes+2];
			Nodes(NewNodeNums(1),:) = Nodes(ToSplit(1),:);
			Nodes(NewNodeNums(2),:) = Nodes(ToSplit(2),:);
	
			Elementgroups{6}.Elems(end,:) = [];
			Elementgroups{9}.Elems(i,:) = [Nds(3:-1:1), NewNodeNums(2:-1:1)', Nds(1)];
			if (i>1)
				Elementgroups{9}.Elems(i-1,6) = NewNodeNums(2);
			end
			
			for j=1:length(ToSplit)
				for EG = 1:length(Elementgroups)-1
					for e=1:length(Elementgroups{EG}.Elems)
						nodeloc = find(Elementgroups{EG}.Elems(e,:)==ToSplit(j));
						if (length(nodeloc) >0)
							if EG==1
								centre = 5;
							else
								centre = 2;
							end
							if (Nodes(Elementgroups{EG}.Elems(e,centre),1)>Nodes(ToSplit(1),1))
								for k=1:length(nodeloc)
                                	Elementgroups{EG}.Elems(e,nodeloc(k)) = NewNodeNums(j);
								end
							end
						end
					end
				end
			end
		end
	end


	for g=1:length(Elementgroups)
        Nodegroups{g}.name = Elementgroups{g}.name;
        Nodegroups{g}.Nodes = unique(reshape(Elementgroups{g}.Elems,[],1));
	end
	Area = zeros(g,1);
	rect = false;
end

