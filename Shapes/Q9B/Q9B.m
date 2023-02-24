classdef Q9B
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
        function obj = Q9B(ipcount1D, rect, zeroWeight)

            [x1D, w1D] = obj.getIpscheme(ipcount1D, zeroWeight);
            
            if (zeroWeight)
                ipcount1D = ipcount1D+2;
            end
            obj.ipcount1D = ipcount1D;
            obj.ipcount = obj.ipcount1D*obj.ipcount1D;
            obj.rectangular = rect;
            
            k = 0;
            for j=1:obj.ipcount1D
                for i=1:obj.ipcount1D
                    k=k+1;
                    xbase(k) = x1D(i);
                    ybase(k) = x1D(j);
                    wbase(k) = w1D(i)*w1D(j);
                end
            end
            
            for i=1:length(xbase)
                x = xbase(i); y = ybase(i);
                Nx = [(1-x).^2; 
                      2*(1-x)*x; 
                      x.^2];
                Dx = [-2*(1-x); 
                      2*(1-x)-2*x; 
                      2*x];
                Ny = [(1-y).^2; 
                      2*(1-y)*y; 
                      y.^2];
                Dy = [-2*(1-y); 
                      2*(1-y)-2*y; 
                      2*y];
                
                Nbase(i,:) = kron(Ny, Nx);
                Gbase(i,:,1) = kron(Ny, Dx);
                Gbase(i,:,2) = kron(Dy, Nx);
            end
            
            obj.xbase = xbase;
            obj.ybase = ybase;
            obj.wbase = wbase;
            obj.Nbase = Nbase;
            obj.Gbase = Gbase;
            
        end
        
        function [N, G, w] = getVals(obj, X, Y)
            N = obj.Nbase;
            
            if obj.rectangular 
            	G(:,:,1) = obj.Gbase(:,:,1)/(X(3)-X(1));
            	G(:,:,2) = obj.Gbase(:,:,2)/(Y(7)-Y(1));
                w = obj.wbase * (X(3)-X(1)) * (Y(7)-Y(1));
            else %% probably requires correcting
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

        end
        
        function xy = getIPGlobal(obj, X,Y)
            
            if obj.rectangular 
                xy(1,:) = X(1) + ((obj.xbase))*(X(3)-X(1));
                xy(2,:) = Y(1) + ((obj.ybase))*(Y(7)-Y(1));
            else
                xy(1,:) = obj.Nbase*X;
				xy(2,:) = obj.Nbase*Y;
            end

        end
 
    end
    
    
    methods (Access = private)
        function [x1D, w1D] = getIpscheme(obj, ipcount1D, zeroWeight)
            if (ipcount1D == 1)
                x1D = 0;
                w1D = 2;
            elseif (ipcount1D == 2)
                x1D = [-1/sqrt(3); 1/sqrt(3)];
                w1D = [1; 1];                
            elseif (ipcount1D == 3)
                x1D = [-sqrt(3/5); 0; sqrt(3/5)];
                w1D = [5/9; 8/9; 5/9];                     
            elseif (ipcount1D ==4)
                x1D = [-sqrt(3/7+2/7*sqrt(6/5)); -sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7+2/7*sqrt(6/5))];
                w1D = [(18-sqrt(30))/36; (18+sqrt(30))/36; (18+sqrt(30))/36; (18-sqrt(30))/36];                     
            elseif (ipcount1D == 5)
                x1D = [-1/3*sqrt(5+2*sqrt(10/7)); -1/3*sqrt(5-2*sqrt(10/7)); 0; 1/3*sqrt(5-2*sqrt(10/7)); 1/3*sqrt(5+2*sqrt(10/7))];
                w1D = [(322-13*sqrt(70))/900; (322+13*sqrt(70))/900; 128/225; (322+13*sqrt(70))/900; (322-13*sqrt(70))/900];                     
            else
                error("Higer order ip schemes not implemented in Shapes.Q9");
            end
            
            if (zeroWeight)
               x1D = [-1; x1D; 1];
               w1D = [0; w1D; 0];
            end
            x1D = (x1D+1)/2;
			w1D = w1D/2;
        end
    end
end

