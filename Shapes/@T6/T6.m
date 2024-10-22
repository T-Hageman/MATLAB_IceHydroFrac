classdef T6
    %Q9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ipcount1D
        ipcount
        
        rectangular
        
        Nbase   % ip, shapefunc
        Gbase   % ip, shapefunc, dx/dy
        
        wbase   % ip
        xbase   % ip
        ybase   % ip
    end
    
    methods
        function obj = T6(ipcount1D, rect, zeroWeight)

            [x, w] = obj.getIpscheme(ipcount1D, zeroWeight);
			obj.rectangular = false;
            
            if (zeroWeight)
                ipcount1D = ipcount1D+2;
            end
            obj.ipcount1D = ipcount1D;
            obj.ipcount = obj.ipcount1D*obj.ipcount1D;

			xbase = x(:,1);
			ybase = x(:,2);
			wbase = w;

            syms s t 
            for i=1:length(xbase)
                x = xbase(i); y = ybase(i);
                x1 = 0; x2 = 1/2; x3 = 1;

				%Ntriangle6nodes=[2*(1-s-t)*((1/2)-s-t) 2*s*(s-(1/2)) 2*t*(t-(1/2)) 4*s*(1-s-t) 4*s*t 4*t*(1-s-t)];
				Ntriangle6nodes=0.5*[2*(1-s-t)*((1)-s-t) 2*s*(s-(0)) 2*t*(t-(0)) 4*s*(1-s-t) 4*s*t 4*t*(1-s-t)];

				Nbase(i,:) = double(subs(subs(Ntriangle6nodes, x), y));
				Gbase(i,:,1) = double(subs(subs(diff(Ntriangle6nodes,s), x), y));
                Gbase(i,:,2) = double(subs(subs(diff(Ntriangle6nodes,t), x), y));
            end
            
            obj.xbase = xbase;
            obj.ybase = ybase;
            obj.wbase = wbase;
            obj.Nbase = Nbase;
            obj.Gbase = Gbase;
            
        end
        
        function [N, G, w] = getVals(obj, X, Y)
            N = obj.Nbase;
            
            dXdXi = obj.Gbase(:,:,1)*X; dXdEta = obj.Gbase(:,:,2)*X;
            dYdXi = obj.Gbase(:,:,1)*Y; dYdEta = obj.Gbase(:,:,2)*Y;

            J(:,1,1) = dXdXi; J(:,2,1) = dXdEta;
            J(:,1,2) = dYdXi; J(:,2,2) = dYdEta;

            G = 0.0*obj.Gbase;
            for i=1:size(obj.Nbase, 1)
                Jinv = inv(squeeze(J(i,:,:)));
                for j=1:size(obj.Nbase, 2)
                    G(i,j,:) = Jinv*squeeze(obj.Gbase(i,j,:));
                end

                w(i) = obj.wbase(i)*abs(det(squeeze(J(i,:,:))));
            end
            

        end
        
        function xy = getIPGlobal(obj, X,Y)
            
            if obj.rectangular 
                xy(1,:) = X(1) + (0.5*(obj.xbase+1))*(X(3)-X(1));
                xy(2,:) = Y(1) + (0.5*(obj.ybase+1))*(Y(7)-Y(1));
            else
                %% not implemented
            end

        end
 
    end
    
    
    methods (Access = private)
        function [x, w] = getIpscheme(obj, ipcount1D, zeroWeight)
            if (ipcount1D == 1)
                x = [1/3, 1/3];
                w = 0.5;
            elseif (ipcount1D == 2)
                x = [2/3, 1/6; 1/6, 1/6; 1/6, 2/3];
				w = [1/3; 1/3; 1/3]*0.5;
            elseif (ipcount1D == 3)
                x = [1/3, 1/3; 3/5, 1/5; 1/5, 1/5; 1/5, 3/5];
                w = [-27/48; 25/48; 25/48; 25/48]*0.5;                     
            else
                error("Higer order ip schemes not implemented in Shapes.Q9");
            end
            
            if (zeroWeight)
               x = [x; 0, 0; 0, 1; 1, 0];
               w = [w; 0; 0; 0];
            end
            
        end
    end
end

