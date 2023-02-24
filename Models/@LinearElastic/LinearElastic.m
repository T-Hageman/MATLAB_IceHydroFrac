classdef LinearElastic < BaseModel
    %LINEARELASTIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        poisson
        young
        D_el
        D_el2
        Hmatswitch
        myK
    end
    
    methods
        function obj = LinearElastic(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "LinearElastic";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType({"dx","dy"});
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            %% constuct plane-strain elastic stiffness matrix
            obj.poisson = inputs.poisson;
            obj.young = inputs.young;
            obj.Hmatswitch = inputs.Hmatswitch;
            
            D_el = zeros(4,4);
            a = obj.young(1) / ((1.0 + obj.poisson(1)) * (1.0 - 2.0*obj.poisson(1)));

            D_el(1, 1) = a * (1.0 - obj.poisson(1));
            D_el(1, 2) = a * obj.poisson(1);
            D_el(1, 3) = a * obj.poisson(1);
            D_el(2, 1) = a * obj.poisson(1);
            D_el(2, 2) = a * (1.0 - obj.poisson(1));
            D_el(2, 3) = a * obj.poisson(1);
            D_el(3, 1) = a * obj.poisson(1);
            D_el(3, 2) = a * obj.poisson(1);
            D_el(3, 3) = a * (1.0 - obj.poisson(1));
            D_el(4, 4) = a * 0.5 * (1.0 - 2.0*obj.poisson(1));
            
            obj.D_el = D_el;
            

			D_el2 = zeros(4,4);
            a = obj.young(2) / ((1.0 + obj.poisson(2)) * (1.0 - 2.0*obj.poisson(2)));

            D_el2(1, 1) = a * (1.0 - obj.poisson(2));
            D_el2(1, 2) = a * obj.poisson(2);
            D_el2(1, 3) = a * obj.poisson(2);
            D_el2(2, 1) = a * obj.poisson(2);
            D_el2(2, 2) = a * (1.0 - obj.poisson(2));
            D_el2(2, 3) = a * obj.poisson(2);
            D_el2(3, 1) = a * obj.poisson(2);
            D_el2(3, 2) = a * obj.poisson(2);
            D_el2(3, 3) = a * (1.0 - obj.poisson(2));
            D_el2(4, 4) = a * 0.5 * (1.0 - 2.0*obj.poisson(2));
            obj.D_el2 = D_el2;
            
            obj.myK = sparse(0,0);
        end
        
        function getKf(obj, physics)
            fprintf("        LinearElastic get Matrix:")
            t = tic;
            
            if (length(obj.myK) == length(physics.fint))
                recalc = false;
            else
                recalc = true;
            end
            
            
            if (recalc == true)
                dofmatX = [];
                dofmatY = [];
                kmat = [];
                parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                    [~, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
                    xy = obj.mesh.getIPCoords(obj.myGroupIndex, n_el);

                    dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                    dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                    dofsXY = [dofsX; dofsY];

                    %X = physics.StateVec(dofsX);
                    %Y = physics.StateVec(dofsY);
                    %XY = [X;Y];

                    %f_el = zeros(length(dofsXY), 1);
                    K_el = zeros(length(dofsXY));
                    for ip=1:length(w)
                        B = obj.getB(G(ip,:,:));
                        %strain = B*XY;

                        if (xy(2,ip) < obj.Hmatswitch)
                            %stress = obj.D_el2*strain;
                            %f_el = f_el + B'*stress*w(ip);
                            K_el = K_el + B'*obj.D_el2*B*w(ip);                       
                        else
                            %stress = obj.D_el*strain;
                            %f_el = f_el + B'*stress*w(ip);
                            K_el = K_el + B'*obj.D_el*B*w(ip);
                        end
                    end

                    [dofmatxloc,dofmatyloc] = ndgrid(dofsXY,dofsXY);
                    dofmatX = [dofmatX; dofmatxloc(:)];
                    dofmatY = [dofmatY; dofmatyloc(:)];
                    kmat = [kmat; K_el(:)];
                end 
                obj.myK = sparse(dofmatX, dofmatY, kmat, length(physics.fint),length(physics.fint));
            end
            physics.fint = physics.fint + obj.myK*physics.StateVec;
            physics.K = physics.K + obj.myK;
            
            
            tElapsed = toc(t);
            fprintf("            (Assemble time:"+string(tElapsed)+")\n");
        end
        
        function B = getB(~, grads)
            cp_count = size(grads, 2);
            B = zeros(4, cp_count*2);
            for ii = 1:cp_count %using plane strain e_zz = 0
				%dx
				B(1, ii) = grads(1,ii, 1);
				B(4, ii) = grads(1,ii, 2);

				%dy
				B(2, ii + cp_count) = grads(1,ii, 2);
				B(4, ii + cp_count) = grads(1,ii, 1);
            end
        end

        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
           hasInfo = false;
           provided = [];
           
           if (var == "stresses" || var=="sxx" || var=="syy" || var=="szz" || var=="sxy" || var == "dstressesdxy") %elem, ip, componenets
               hasInfo = true;
               if (var == "stresses")
                   provided = zeros(length(elems), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount, 4);
               elseif (var == "dstressesdxy")  %elem, ip, comp, dxy
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, elems(1));
                    provided = zeros(length(elems), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount, 4, length(Elem_Nodes)*2 );
               else
                   provided = zeros(length(elems), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount);
               end
                for el=1:length(elems)
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, elems(el));
                    [~, G, w] = obj.mesh.getVals(obj.myGroupIndex, elems(el));
                    xy = obj.mesh.getIPCoords(obj.myGroupIndex, elems(el));

                    dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                    dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);

                    X = physics.StateVec(dofsX);
                    Y = physics.StateVec(dofsY);
                    XY = [X;Y];

                    for ip=1:length(w)
                        B = obj.getB(G(ip,:,:));
                        strain = B*XY;
                        
                        if (xy(2,ip) < obj.Hmatswitch)
                            stress = obj.D_el2*strain; 
                            dStress = obj.D_el2*B;
                        else
                            stress = obj.D_el*strain;
                            dStress = obj.D_el*B;
                        end

                        if (loc == "Interior")
                            switch var 
                                case "stresses"
                                    provided(el, ip, :) = stress;
                                case "sxx"
                                    provided(el, ip) = stress(1);
                                case "syy"
                                    provided(el, ip) = stress(2);
                                case "szz"
                                    provided(el, ip) = stress(3);
                                case "sxy"
                                    provided(el, ip) = stress(4);
                                case "dstressesdxy"
                                    provided(el, ip, :, :) = dStress;
                            end
                        end
                        
                    end
                end
           end
        end
        
    end
end

