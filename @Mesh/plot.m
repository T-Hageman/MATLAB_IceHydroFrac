function f1 = plot(obj, plotnodes, plotelems, plotnames, plotnodenames)
    f1 = figure();
    if (plotnodes)
        plot(obj.Nodes(:,1),obj.Nodes(:,2), 'k*');
		if plotnodenames
			for i=1:length(obj.Nodes)
				textloc(i,:) = [obj.Nodes(i,1)+0.05*rand() obj.Nodes(i,2)+0.05*rand()];
                textlabel(i) = "N"+string(i);
			end
			text(textloc(:,1), textloc(:,2), textlabel);
		end
    end 
    hold on
    if (plotelems)
        for g=1:length(obj.Elementgroups)
			clear X Y order textloc textlabel
            if (obj.Elementgroups{g}.type == "Q9" || obj.Elementgroups{g}.type == "Q9B")
                for el=1:size(obj.Elementgroups{g}.Elems, 1)
                    elnodes = obj.Elementgroups{g}.Elems(el,:);
                    order = [1 2 3 6 9 8 7 4];
                    X(el,:) = obj.Nodes(elnodes(order),1);
                    Y(el,:) = obj.Nodes(elnodes(order),2);

                    textloc(el,:) = [obj.Nodes(elnodes(5),1) obj.Nodes(elnodes(5),2)];
                    textlabel(el) = "E"+string(el);
                end
                patch('XData',X','YData',Y','EdgeColor','black','FaceColor','none');
                hold on
                if (plotnames)
                   text(textloc(:,1), textloc(:,2), textlabel);
                end
			end
            if (obj.Elementgroups{g}.type == "L3" || obj.Elementgroups{g}.type == "L3B")
                for el=1:size(obj.Elementgroups{g}.Elems, 1)
                    elnodes = obj.Elementgroups{g}.Elems(el,:);
                    order = [1 2 3];
                    X(el,:) = obj.Nodes(elnodes(order),1);
                    Y(el,:) = obj.Nodes(elnodes(order),2);

					textloc(el,:) = [obj.Nodes(elnodes(2),1) obj.Nodes(elnodes(2),2)];
                    textlabel(el) = "E"+string(el);
				end
				if (el>0)
                	plot(X,Y,'r')
                	hold on
					if (plotnames)
						text(textloc(:,1), textloc(:,2), textlabel,'Color','red');
					end
				end
            end
            if (obj.Elementgroups{g}.type == "LI6" || obj.Elementgroups{g}.type == "LI6B")
                for el=1:size(obj.Elementgroups{g}.Elems, 1)
                    elnodes = obj.Elementgroups{g}.Elems(el,:);
                    order = [1 2 3 6 5 4];
                    X(el,:) = obj.Nodes(elnodes(order),1);
                    Y(el,:) = obj.Nodes(elnodes(order),2);

					textloc(el,:) = [obj.Nodes(elnodes(2),1) obj.Nodes(elnodes(2),2)];
                    textlabel(el) = "E"+string(el);
                end
                plot(X',Y','b')
                hold on
				if (plotnames)
                   text(textloc(:,1), textloc(:,2), textlabel,'Color','blue');
                end
            end
        end
    end
end

