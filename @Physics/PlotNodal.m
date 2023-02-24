function PlotNodal(obj, dofName, dispscale, plotloc, zscale)
    dxTypes = obj.dofSpace.getDofType({"dx";"dy";dofName});


    for g=1:length(obj.mesh.Elementgroups)
        if (obj.mesh.Elementgroups{g}.name == plotloc)

            if (obj.mesh.Elementgroups{g}.type == "Q9")
                for el=1:size(obj.mesh.Elementgroups{g}.Elems, 1)
                    elnodes = obj.mesh.Elementgroups{g}.Elems(el,:);

                    xdofs = obj.dofSpace.getDofIndices(dxTypes(1), elnodes);
                    ydofs = obj.dofSpace.getDofIndices(dxTypes(2), elnodes);
                    zdofs = obj.dofSpace.getDofIndices(dxTypes(3), elnodes);

                    order = [1 2 3 6 9 8 7 4];
                    X(el,:) = obj.mesh.Nodes(elnodes(order),1)+dispscale*obj.StateVec(xdofs(order));
                    Y(el,:) = obj.mesh.Nodes(elnodes(order),2)+dispscale*obj.StateVec(ydofs(order));
                    Z(el,:) = obj.StateVec(zdofs(order));
                end
                patch(X',Y',Z'*zscale,'EdgeColor','interp','FaceColor','interp');
                hold on
                colorbar
			end
			if (obj.mesh.Elementgroups{g}.type == "Q9B")
                for el=1:size(obj.mesh.Elementgroups{g}.Elems, 1)
                    elnodes = obj.mesh.Elementgroups{g}.Elems(el,:);

                    xdofs = obj.dofSpace.getDofIndices(dxTypes(1), elnodes);
                    ydofs = obj.dofSpace.getDofIndices(dxTypes(2), elnodes);
                    zdofs = obj.dofSpace.getDofIndices(dxTypes(3), elnodes);

                    order = [1 3 9 7];
                    X(el,:) = obj.mesh.Nodes(elnodes(order),1)+dispscale*obj.StateVec(xdofs(order));
                    Y(el,:) = obj.mesh.Nodes(elnodes(order),2)+dispscale*obj.StateVec(ydofs(order));
                    Z(el,:) = obj.StateVec(zdofs(order));
                end
                patch(X',Y',Z'*zscale,'EdgeColor','interp','FaceColor','interp');
                hold on
                colorbar
            end
            if (obj.mesh.Elementgroups{g}.type == "L3" || obj.mesh.Elementgroups{g}.type == "LI6")
                for el=1:size(obj.mesh.Elementgroups{g}.Elems, 1)
                    elnodes = obj.mesh.Elementgroups{g}.Elems(el,:);

                    zdofs = obj.dofSpace.getDofIndices(dxTypes(3), elnodes(1:3));
					xdofs = obj.dofSpace.getDofIndices(dxTypes(1), elnodes);
                    ydofs = obj.dofSpace.getDofIndices(dxTypes(2), elnodes);

                    order = [1 2 3];
					if (obj.mesh.Elementgroups{g}.type == "L3")
                    	X(el,:) = [obj.mesh.Nodes(elnodes(order),1)+dispscale*obj.StateVec(xdofs(order));NaN];
                    	Y(el,:) = [obj.mesh.Nodes(elnodes(order),2)+dispscale*obj.StateVec(ydofs(order));NaN];
                    	Z(el,:) = [obj.StateVec(zdofs(order));NaN];
					else
                    	X(el,:) = [obj.mesh.Nodes(elnodes(order),1)+0*0.5*dispscale*obj.StateVec(xdofs(order))+2*0.5*dispscale*obj.StateVec(xdofs(order+3));NaN];
                    	Y(el,:) = [obj.mesh.Nodes(elnodes(order),2)+0*0.5*dispscale*obj.StateVec(ydofs(order))+2*0.5*dispscale*obj.StateVec(ydofs(order+3));NaN];
                    	Z(el,:) = [obj.StateVec(zdofs(order));NaN];
					end

                end
                patch(X',Y',Z'*zscale,'EdgeColor','interp','FaceColor','None','Marker','None','MarkerFaceColor','flat','LineWidth',4);
                hold on
                colorbar
			end
			 if (obj.mesh.Elementgroups{g}.type == "L3B" || obj.mesh.Elementgroups{g}.type == "LI6B")
                for el=1:size(obj.mesh.Elementgroups{g}.Elems, 1)
                    elnodes = obj.mesh.Elementgroups{g}.Elems(el,:);

                    zdofs = obj.dofSpace.getDofIndices(dxTypes(3), elnodes(1:3));
					xdofs = obj.dofSpace.getDofIndices(dxTypes(1), elnodes);
                    ydofs = obj.dofSpace.getDofIndices(dxTypes(2), elnodes);

                    order = [1 3];
					if (obj.mesh.Elementgroups{g}.type == "L3B")
                    	X(el,:) = [obj.mesh.Nodes(elnodes(order),1)+dispscale*obj.StateVec(xdofs(order));NaN];
                    	Y(el,:) = [obj.mesh.Nodes(elnodes(order),2)+dispscale*obj.StateVec(ydofs(order));NaN];
                    	Z(el,:) = [obj.StateVec(zdofs(order));NaN];
					else
                    	X(el,:) = [obj.mesh.Nodes(elnodes(order),1)+0.5*dispscale*obj.StateVec(xdofs(order))+0.5*dispscale*obj.StateVec(xdofs(order+3));NaN];
                    	Y(el,:) = [obj.mesh.Nodes(elnodes(order),2)+0.5*dispscale*obj.StateVec(ydofs(order))+0.5*dispscale*obj.StateVec(ydofs(order+3));NaN];
                    	Z(el,:) = [obj.StateVec(zdofs(order));NaN];
					end

                end
                patch(X',Y',Z'*zscale,'EdgeColor','interp','FaceColor','None','Marker','None','MarkerFaceColor','flat','LineWidth',4);
                hold on
                colorbar
            end
        end
    end
end


