function Propagate_Disc_New(obj, physics)

	downGroup = 6;
	leftGroup = 7;
	rightGroup= 8;

	if (size(obj.Elementgroups{downGroup}.Elems,1)>0)	%propagate downwards
		ifrac = length(obj.Elementgroups{9}.Elems)+1;
		nnodes = length(obj.Nodes);

		Nds = obj.Elementgroups{downGroup}.Elems(end,:);

		ToSplit = Nds(2:3);

		newnodes = [nnodes+1; nnodes+2];
		obj.Nodes(newnodes(1),:) = obj.Nodes(ToSplit(1),:);
		obj.Nodes(newnodes(2),:) = obj.Nodes(ToSplit(2),:);

        dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  
        physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodes);
        for i=1:2
            physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), ToSplit(1:2)));
            physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), ToSplit(1:2)));  

			physics.VAvec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes),:) = physics.VAvec(physics.dofSpace.getDofIndices(dofTypeIndices(i), ToSplit(1:2)),:);
            physics.VAvecOld(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes),:) = physics.VAvecOld(physics.dofSpace.getDofIndices(dofTypeIndices(i), ToSplit(1:2)),:);  
        end

        physics.dofSpace.addDofs(dofTypeIndices(3), Nds(1:2) );
        physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), Nds(1:2) )) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), Nds(3) ));
        physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), Nds(1:2) )) = 0;

		physics.VAvec(physics.dofSpace.getDofIndices(dofTypeIndices(3), Nds(1:2) ),:) = 0;
		physics.VAvecOld(physics.dofSpace.getDofIndices(dofTypeIndices(3), Nds(1:2) ),:) = 0;

		obj.Elementgroups{downGroup}.Elems(end,:) = [];
		obj.Elementgroups{9}.Elems(ifrac,:) = [Nds(3:-1:1), newnodes(2:-1:1)', Nds(1)];
		obj.Elementgroups{9}.Elems(ifrac-1,6) = newnodes(2);
		
		for j=1:length(ToSplit)
			for EG = 1:length(obj.Elementgroups)-1
				for el=1:size(obj.Elementgroups{EG}.Elems,1)
					nodeloc = find(obj.Elementgroups{EG}.Elems(el,:)==ToSplit(j));
					if (length(nodeloc) >0)
						if EG==1
							centre = 5;
						else
							centre = 2;
						end
						if (obj.Nodes(obj.Elementgroups{EG}.Elems(el,centre),1)>obj.Nodes(ToSplit(1),1))
							for k=1:length(nodeloc)
                                obj.Elementgroups{EG}.Elems(el,nodeloc(k)) = newnodes(j);
							end
						end
					end
				end
			end
		end
	elseif (isempty(obj.Tjunction))  % T Junction
		ifrac = length(obj.Elementgroups{9}.Elems)+1;
		nnodes = length(obj.Nodes);

		NdsLeft = obj.Elementgroups{leftGroup}.Elems(end,:);
		NdsRight = obj.Elementgroups{rightGroup}.Elems(1,:);

		ToSplit = [NdsLeft(2), NdsLeft(3), NdsRight(2)];
		newnodes = [nnodes+1; nnodes+2; nnodes+3; nnodes+4];
		oldnodes = [ToSplit(1); ToSplit(2); ToSplit(2); ToSplit(3)];

		obj.Nodes(newnodes,:) = obj.Nodes(oldnodes,:);

		obj.Elementgroups{leftGroup}.Elems(end,:) = [];
		obj.Elementgroups{rightGroup}.Elems(1,:) = [];

		obj.Elementgroups{9}.Elems(ifrac,:)   = [NdsLeft(1:3), NdsLeft(1), newnodes(1:2)'];   %NewLeft
		obj.Elementgroups{9}.Elems(ifrac+1,:) = [NdsRight(1:3), newnodes(3:4)', NdsRight(3)]; %NewRight

		obj.Elementgroups{9}.Elems(ifrac-1,3) = newnodes(2);
		obj.Elementgroups{9}.Elems(ifrac-1,6) = newnodes(3);

		obj.Tjunction = [newnodes(2), NdsLeft(3)];

        dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  
        physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodes);
        for i=1:2
            physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), oldnodes));
            physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), oldnodes));  

			physics.VAvec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes),:) = physics.VAvec(physics.dofSpace.getDofIndices(dofTypeIndices(i), oldnodes),:);
            physics.VAvecOld(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes),:) = physics.VAvecOld(physics.dofSpace.getDofIndices(dofTypeIndices(i), oldnodes),:);  
        end

        physics.dofSpace.addDofs(dofTypeIndices(3), [NdsLeft(1:2) NdsRight(2:3)] );
        physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), [NdsLeft(1:2) NdsRight(2:3)] )) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), NdsLeft(3) ));
        physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), [NdsLeft(1:2) NdsRight(2:3)] )) = 0;

		physics.VAvec(physics.dofSpace.getDofIndices(dofTypeIndices(3), [NdsLeft(1:2) NdsRight(2:3)] ),:) = 0.0;
		physics.VAvecOld(physics.dofSpace.getDofIndices(dofTypeIndices(3), [NdsLeft(1:2) NdsRight(2:3)] ),:) = 0.0;
        
        physics.dofSpace.addDofs(dofTypeIndices(3), newnodes(2) );
        physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), newnodes(2) )) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), NdsLeft(3) ));
        physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), newnodes(2) )) = 0;
		
		physics.VAvec(physics.dofSpace.getDofIndices(dofTypeIndices(3), newnodes(2) ),:) = 0.0;
		physics.VAvecOld(physics.dofSpace.getDofIndices(dofTypeIndices(3), newnodes(2) ),:) = 0.0;

		%renumbering
		for j=1:length(oldnodes)
			for EG = 1:length(obj.Elementgroups)-1
				if (EG ~= leftGroup && EG ~= rightGroup)
					for el=1:size(obj.Elementgroups{EG}.Elems,1)
						nodeloc = find(obj.Elementgroups{EG}.Elems(el,:)==oldnodes(j));
						if (length(nodeloc) >0)
							if EG==1
								centre = 5;
							else
								centre = 2;
							end
							if (obj.Nodes(obj.Elementgroups{EG}.Elems(el,centre),2)>obj.Nodes(oldnodes(j),2))
								for k=1:length(nodeloc)
									if (round(nodeloc/2)==nodeloc/2) %centre node, valid for 1 and 4
										obj.Elementgroups{EG}.Elems(el,nodeloc(k)) = newnodes(j);
									else
										if (obj.Nodes(obj.Elementgroups{EG}.Elems(el,centre),1)>obj.Nodes(oldnodes(2),1)) %right element 3
											obj.Elementgroups{EG}.Elems(el,nodeloc(k)) = newnodes(3);
										else %left element 2
											obj.Elementgroups{EG}.Elems(el,nodeloc(k)) = newnodes(2);
										end
									end
								end
							end
						end
					end
				end
			end
		end
		obj.Propagate_Disc_New(physics);
	else  %SideWays
		ifrac = length(obj.Elementgroups{9}.Elems)+1;
		nnodes = length(obj.Nodes);

		NdsLeft = obj.Elementgroups{leftGroup}.Elems(end,:);
		NdsRight = obj.Elementgroups{rightGroup}.Elems(1,:);

		ToSplit = [NdsLeft(2), NdsLeft(3), NdsRight(1), NdsRight(2)];
		newnodes = [nnodes+1; nnodes+2; nnodes+3; nnodes+4];

		obj.Nodes(newnodes,:) = obj.Nodes(ToSplit,:);

		obj.Elementgroups{leftGroup}.Elems(end,:) = [];
		obj.Elementgroups{rightGroup}.Elems(1,:) = [];

		obj.Elementgroups{9}.Elems(ifrac,:)   = [NdsLeft(1:3), NdsLeft(1), newnodes(1:2)'];   %NewLeft
		obj.Elementgroups{9}.Elems(ifrac+1,:) = [NdsRight(1:3), newnodes(3:4)', NdsRight(3)]; %NewRight

		obj.Elementgroups{9}.Elems(ifrac-2,4) = newnodes(2);
		obj.Elementgroups{9}.Elems(ifrac-1,6) = newnodes(3);

        dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  
        physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodes);
        for i=1:2
            physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), ToSplit));
            physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), ToSplit));  

			physics.VAvec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes),:) = physics.VAvec(physics.dofSpace.getDofIndices(dofTypeIndices(i), ToSplit),:);
            physics.VAvecOld(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes),:) = physics.VAvecOld(physics.dofSpace.getDofIndices(dofTypeIndices(i), ToSplit),:);  
        end

        physics.dofSpace.addDofs(dofTypeIndices(3), [NdsLeft(1:2) NdsRight(2:3)] );
        physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), [NdsLeft(1:2) NdsRight(2:3)] )) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), NdsLeft(3) ));
        physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), [NdsLeft(1:2) NdsRight(2:3)] )) = 0;	

		physics.VAvec(physics.dofSpace.getDofIndices(dofTypeIndices(3), [NdsLeft(1:2) NdsRight(2:3)] ),:) = 0.0;
		physics.VAvecOld(physics.dofSpace.getDofIndices(dofTypeIndices(3), [NdsLeft(1:2) NdsRight(2:3)] ),:) = 0.0;

		for j=1:length(ToSplit)
			for EG = 1:length(obj.Elementgroups)-1
				if (EG ~= leftGroup && EG ~= rightGroup)
					for el=1:size(obj.Elementgroups{EG}.Elems,1)
						nodeloc = find(obj.Elementgroups{EG}.Elems(el,:)==ToSplit(j));
						if (length(nodeloc) >0)
							if EG==1
								centre = 5;
							else
								centre = 2;
							end
							if (obj.Nodes(obj.Elementgroups{EG}.Elems(el,centre),2)>obj.Nodes(ToSplit(j),2))
								for k=1:length(nodeloc)
									obj.Elementgroups{EG}.Elems(el,nodeloc(k)) = newnodes(j);
								end
							end
						end
					end
				end
			end
		end




	end


% 
% 
% 
%     neighbours = obj.getConnected("Internal", tipnode);
% 
%     if (abs(direction(1))<1e-6 && abs(direction(2)--1)<1e-6) % propagate downwards
%         elem_nodes = obj.Elementgroups{1}.Elems(neighbours(1),:);
%         SideNodes = elem_nodes([9 6 3]);
% 
%         newnodes = substitude(obj, SideNodes(1:2), "Vertical");
% 
%         new_nodes_int = [SideNodes(1:3),newnodes,SideNodes(3)];
%         obj.Elementgroups{6}.Elems(end+1,:) = new_nodes_int';
% 
%         dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  
% 
%         physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodes);
%         for i=1:2
%             physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(1:2)));
%             physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(1:2)));  
%         end
% 
%         physics.dofSpace.addDofs(dofTypeIndices(3), SideNodes(2:3) );
%         physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(2:3) )) = 0.0;
%         physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(2:3) )) = 0.0;
%     elseif (abs(direction(1)--1)<1e-6 && abs(direction(2))<1e-6) % propagate left
%         elem_nodes = obj.Elementgroups{1}.Elems(neighbours(1),:);
%         if (size(neighbours,1) == 2)
%             SideNodes = elem_nodes([7 8 9]);
%         else
%             SideNodes = elem_nodes([7 8 9]);
%         end
%         
%         
%         newnodes = substitude(obj, SideNodes(2:3), "Horizontal");
%         
%         new_nodes_int = [SideNodes(1:3),SideNodes(1),newnodes];
%         obj.Elementgroups{6}.Elems(end+1,:) = new_nodes_int';
%         
%         dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  
%         physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodes);
%         for i=1:2
%             physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(2:3)));
%             physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(2:3)));  
%         end
%         
%         physics.dofSpace.addDofs(dofTypeIndices(3), SideNodes(1:2) );
%         physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(1:2) )) = 0.0;
%         physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(1:2) )) = 0.0;
%         
%     elseif (abs(direction(1)-1)<1e-6 && abs(direction(2))<1e-6) % propagate right
%         elem_nodes = obj.Elementgroups{1}.Elems(neighbours(3),:);
%         SideNodes = elem_nodes([7 8 9]);
%         
%         newnodes = substitude(obj, SideNodes(1:2), "Horizontal");
%         new_nodes_int = [SideNodes(1:3),newnodes,SideNodes(3)];
%         obj.Elementgroups{6}.Elems(end+1,:) = new_nodes_int';
%         
%         dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  
%         physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodes);
%         for i=1:2
%             physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(1:2)));
%             physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodes)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodes(1:2)));  
%         end
%         
%         physics.dofSpace.addDofs(dofTypeIndices(3), SideNodes(2:3) );
%         physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(2:3) )) = 0.0;
%         physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodes(2:3) )) = 0.0;
%      elseif (abs(direction(1)-42)<1e-6 && abs(direction(2)--42)<1e-6) % T-junction
%         %split tip 
%         newnode = substitude(obj, tipnode, "Tip");
%         dofTypeIndices = physics.dofSpace.getDofType({"dx","dy","pd"})  ;  
%         physics.dofSpace.addDofs(dofTypeIndices(1:2), newnode);
%         for i=1:2
%             physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnode)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), tipnode));
%             physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnode)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), tipnode));  
%         end        
%         
%         
%         %new elements
%         elem_nodes  = obj.Elementgroups{1}.Elems(neighbours(1),:);
%         elem_nodes2 = obj.Elementgroups{1}.Elems(neighbours(3),:);
%         
%         SideNodesLeft  = elem_nodes([7 8 9]);
%         SideNodesRight = elem_nodes2([7 8 9]);
%         
%         newnodesLeft = substitude(obj, SideNodesLeft(2:3), "Horizontal");
%         new_nodes_intLeft = [SideNodesLeft(1:3),SideNodesLeft(1), newnodesLeft];
%         obj.Elementgroups{6}.Elems(end+1,:) = new_nodes_intLeft';
%         
%         newnodesRight = substitude(obj, SideNodesRight(2), "Horizontal");
%         new_nodes_intRight = [SideNodesRight(1:3), newnode, newnodesRight, SideNodesRight(3)];
%         obj.Elementgroups{6}.Elems(end+1,:) = new_nodes_intRight';
%         
%         physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodesLeft);
%         for i=1:2
%             physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodesLeft)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodesLeft(2:3)));
%             physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodesLeft)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodesLeft(2:3)));  
%         end  
%         
%         physics.dofSpace.addDofs(dofTypeIndices(1:2), newnodesRight);
%         for i=1:2
%             physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodesRight)) = physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodesRight(2:2)));
%             physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), newnodesRight)) = physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(i), SideNodesRight(2:2)));  
%         end
%         
%         
%         physics.dofSpace.addDofs(dofTypeIndices(3), SideNodesLeft(1:3) );
%         physics.dofSpace.addDofs(dofTypeIndices(3), SideNodesRight(2:3) );
%         
%         physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodesLeft(1:3) )) = 0.0;
%         physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodesLeft(1:3) )) = 0.0;
%         
%         physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodesRight(2:3) )) = 0.0;
%         physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), SideNodesRight(2:3) )) = 0.0;
%         
%         % last corrections
%         obj.Tjunction = [new_nodes_intLeft(3) new_nodes_intLeft(6)];
%         obj.Elementgroups{6}.Elems(end-2,3) = new_nodes_intLeft(6);
%         
%         physics.dofSpace.addDofs(dofTypeIndices(3), new_nodes_intLeft(6) );
%         physics.StateVec(physics.dofSpace.getDofIndices(dofTypeIndices(3), new_nodes_intLeft(6) )) = 0.0;
%         physics.StateVec_Old(physics.dofSpace.getDofIndices(dofTypeIndices(3), new_nodes_intLeft(6) )) = 0.0;
%     else
%         
%        error("still to implement direction for propagation"); 
%         
%     end


    obj.check();
end
% 
% function newnodes = substitude(obj, oldnodes, direction)
%     ToSplit = oldnodes;
%     for j=1:length(ToSplit)
%         nnodes = size(obj.Nodes, 1) +1;
%         obj.Nodes(nnodes,:) = obj.Nodes(ToSplit(j),:);
%         for g=1:6
%             for el=1:size(obj.Elementgroups{g}.Elems, 1)
%                 nodeloc = find(obj.Elementgroups{g}.Elems(el,:)==ToSplit(j));
%                 if (length(nodeloc) >0)
%                     for k=1:length(nodeloc)
%                         if (direction == "Vertical")
%                             if (obj.Elementgroups{g}.type == "Q9" && (nodeloc(k) == 1 || nodeloc(k) == 4 || nodeloc(k) == 7))
%                                 obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
%                             end
%                             if (obj.Elementgroups{g}.type == "L3" && nodeloc(k) == 3 )
%                                 obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
%                             end
%                             if (obj.Elementgroups{g}.type == "LI6" && nodeloc(k) == 6)
%                                 obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
%                             end
%                         elseif (direction == "Horizontal")
%                             if (obj.Elementgroups{g}.type == "Q9" && (nodeloc(k) == 1 || nodeloc(k) == 2 || nodeloc(k) == 3))
%                                 obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
%                             end
%                             if (obj.Elementgroups{g}.type == "L3" && (nodeloc(k) == 3 || nodeloc(k) == 1))
%                                 obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
%                             end
%                             if (obj.Elementgroups{g}.type == "LI6" && (nodeloc(k) == 6 || nodeloc(k) == 4))
%                                 obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
%                             end
%                         elseif (direction == "Tip")
%                             if (obj.Elementgroups{g}.type == "Q9" && nodeloc(k) == 1)
%                                 obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
%                             end
%                             if (obj.Elementgroups{g}.type == "LI6" && (nodeloc(k) == 6 ))
%                                 obj.Elementgroups{g}.Elems(el,nodeloc(k)) = nnodes;
%                             end
%                         else
%                            error("Fracture propagation error") 
%                         end
%                     end
%                 end
%             end
%         end
%         newnodes(j) = nnodes;
%     end
% 
% end