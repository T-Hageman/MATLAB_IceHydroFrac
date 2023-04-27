function Solve(obj)

    stop = false;
    it = 0;
 
	obj.physics.Commit("StartIt");
    obj.physics.Assemble(true);
    obj.physics.Constrain();
    recalc_pre=true;
    En_err0 = -1;
    curr_max_it = obj.maxIt;
	con_dt = true;
	first_noncon = true;

    while(stop == false)
        fprintf("    Solving it:" + string(it) + "      ");
        tsolve = tic;

		if (con_dt == false && first_noncon == true)
			obj.physics.Assemble(con_dt);
			obj.physics.Constrain();
			first_noncon = false;
		end
        
        recalc_pre = true;
        if (recalc_pre)
            [P,R,C] = equilibrate(obj.physics.K);
            recalc_pre = false;
        end

        if true
            d = -R*P*obj.physics.fint;
            B = R*P*obj.physics.K*C;
            
            if true
                dy = B\d;
            else
                [L,U] = ilu(B,struct('type','ilutp','droptol',1e-6));
                dy = gmres(B,d,[],1e-6,200,L,U);
            end
            
            dx = C*dy;
        else
            dx = -obj.physics.K\obj.physics.fint;
		end
		if (sum(isnan(dx)+isinf(dx))>0)
			something wrong in solver
		end
        tsolve = toc(tsolve);
        fprintf("        (Solver time:"+string(tsolve)+")\n");

        if (obj.linesearch && it>0)
            e0 = obj.physics.fint'*dx;
            dx = obj.physics.Update(dx);
            
            obj.physics.Assemble(con_dt);
            obj.physics.Constrain();
            
			e1 = obj.physics.fint'*dx;

            factor = -e0/(e1-e0);
            factor = max(obj.linesearchLims(1), min(obj.linesearchLims(2), factor));
            obj.physics.Update(-(1-factor)*dx);
            fprintf("    Linesearch: " + string(e0) + " -> " + string(e1) + ":  eta=" + string(factor) +"\n");
        else
            obj.physics.Update(dx);
        end
        
        % convergence
        obj.physics.Assemble(con_dt);
        obj.physics.Constrain();
        if (En_err0 < 0)
            En_err0 = sum(abs(obj.physics.fint.*dx));
            En_err = En_err0;
		else
			En_err = sum(abs(obj.physics.fint.*dx));
        end
        En_err_n = En_err/En_err0;
		         
        fprintf("    Residual:" + string(En_err_n) + "   ("+string(En_err)  +") \n");
        if (obj.physics.ArcTime.Enable)
			[~,fc_error,~,~] = obj.physics.models{3}.get_fc(obj.physics);
			fc_error = abs(fc_error);
            fprintf("    dt:" + obj.physics.dt + "\n");
			fprintf("    fc error:" + string(fc_error)  +"\n");

			if (En_err_n<1e-4)
				con_dt = false;
			end
		else
			fc_error = 0;
		end
        
        it=it+1;
        if (it>curr_max_it || ((En_err_n<obj.Conv || En_err<obj.tiny) && fc_error<obj.physics.ArcTime.Tol))
            obj.physics.Commit("Pathdep");
            irr = obj.physics.Irreversibles();
            if (irr == false || obj.physics.ArcTime.Enable == true)
                stop = true;
			else
				obj.physics.Commit("Irrevirsibles");
                obj.physics.Assemble(con_dt);
                obj.physics.Constrain();
                recalc_pre = true;
                En_err0 = -1;
                curr_max_it = it + obj.maxIt;
            end
        end
    end
    
    obj.physics.Commit("Timedep");
    
end

