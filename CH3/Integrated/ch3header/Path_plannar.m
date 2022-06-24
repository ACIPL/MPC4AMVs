classdef Path_plannar < handle
    properties
        N
        model %Model = AUV(ones(11,1),3,zeros(6,1),zeros(3,1))
        stroke
        mlti
    end
    
    methods
        function obj = Path_plannar(N,model,stroke,mlti)
            obj.N = N;
            obj.model = model;
            obj.stroke = stroke;
            obj.mlti = mlti;
        end
        
        function obj = calc_upperbound(obj,upperbound)
            if ~isempty(upperbound)
                m = length(upperbound);
                Thrust_max = zeros(m*obj.N,1);
                for i=1:1:obj.N
                    Thrust_max(1+m*(i-1):m*i,1) = upperbound;
                end
                obj.upperbound = Thrust_max;
            else
                obj.upperbound = [];
            end
        end
        
        
        function cost = planning_cost(obj,a,bbrks,dt )
            mlti = obj.mlti;
            a=a';
            knots=augknt(bbrks,mlti);
            sp=spmak(knots,a);
            
            fprime=fnder(sp);
            f2prime=fnder(sp,2);
            f3prime=fnder(sp,3);
            
            nbrk=length(bbrks);
            
            cost=0;
            
            for x=bbrks(1):dt:bbrks(nbrk)
                
                kprimex=(fnval(f3prime,x)*(1+fnval(fprime,x)^2)^(3/2)-fnval(f2prime,x)*(3/2)*(1+fnval(fprime,x)^2)^(1/2)*2*fnval(fprime,x)*fnval(f2prime,x))/(1+fnval(fprime,x)^2)^3;
                kprimex=kprimex^2;
                f0x=kprimex/(1+fnval(fprime,x)^2)^(1/2);
                cost=cost+dt*f0x;
                
            end
            
        end
        
        
        function [ y,ceq ] = planning_init_contraints( obj,a,bbrks,dt,X,Y)
            strok = obj.stroke;
            mlti = obj.mlti;
            a=a';
            knots=augknt(bbrks,mlti);
            sp=spmak(knots,a);
            
            nbrk=length(bbrks);
            
            x=bbrks(1):dt:bbrks(nbrk);
            x=x';
            nx=length(x);
            y=zeros(2*nx,1);
            ceq=zeros(2*nx,1);
            
            for i=1:1:nx
                
                y(i)=fnval(sp,x(i))-plygfun2(x(i),X,Y,strok);
                
            end
            
            for j=(nx+1):1:(2*nx)
                
                y(j)=plygfun1(x(j-nx),X,Y)-fnval(sp,x(j-nx));
            end
            
        end
        function [y, ceq] = planning_constraints(obj, a,bbrks,T,X,Y,y0,y0d,y02d,eps )
            mlti1 = obj.mlti;
            strok = obj.stroke;
            a=a';
            knots=augknt(bbrks,mlti1);
            sp=spmak(knots,a);
            
            nbrk=length(bbrks);
            
            x=bbrks(1):T:bbrks(nbrk);
            x=x';
            nx=length(x);
            y=zeros(2*nx+6,1);
            ceq=zeros(2*nx,1);
            
            for i=1:1:nx
                
                y(i)=fnval(sp,x(i))-plygfun2(x(i),X,Y,strok);
                
            end
            
            for j=(nx+1):1:(2*nx)
                
                y(j)=plygfun1(x(j-nx),X,Y)-fnval(sp,x(j-nx));
            end
            
            
            y(2*nx+1)=fnval(sp,X(1))-y0-eps;
            y(2*nx+2)=y0-fnval(sp,X(1))-eps;
            
            y(2*nx+3)=fnval(fnder(sp,1),X(1))-y0d-eps;
            y(2*nx+4)=y0d-fnval(fnder(sp,1),X(1))-eps;
            
            y(2*nx+5)=fnval(fnder(sp,2),X(1))-y02d-eps;
            y(2*nx+6)=y02d-fnval(fnder(sp,2),X(1))-eps;
            
        end
        
        
        function [x, fval] = calc_path1(obj,x0,bbrks,X,Y,dt)
%             options = optimset('Algorithm','interior-point');
            [x, fval] = fmincon(@(x) obj.planning_cost(x,bbrks,dt),x0,[],[],...
                [],[],[],[],...
                @(x)obj.planning_init_contraints(x,bbrks,dt,X,Y));
        end
        
        function [x, fval] = calc_path2(obj,x0,bbrks,X,Y,y0,y0d,y02d,eps,dt)
%             if everdone < 30
%                 options = optimset('Algorithm','interior-point','MaxFunEvals',1200);
%             else
                options = optimset('Algorithm','sqp','MaxFunEvals',1200);
%             end
%             options = optimset('Algorithm','sqp'); 
            [x, fval] = fmincon(@(x) obj.planning_cost(x,bbrks,dt),x0,[],[],...
                [],[],[],[],...
                @(x)obj.planning_constraints(x,bbrks,dt,X,Y,y0,y0d,y02d,eps),options);
        end
        %%%%%%%%%%%%%%%%%%
        function [eta_Ref,P] = reference(obj,t,Spline_map2,SpAry,i,M,T_plt,y_scale,dt)
            N = obj.N;
            
            xR=obj.xRef(t);
            xRdot=obj.xRefdot(t);
            
            i_spl=floor((i-1)/M)+1;
            Spline_map2(i,1)=(i-1)*T_plt;
            Spline_map2(i,2)=y_scale*fnval(SpAry{i_spl,1},Spline_map2(i,1));
            Spline_map2(i,3)=y_scale*fnval(fnder(SpAry{i_spl,1},1),Spline_map2(i,1));
            Spline_map2(i,4)=y_scale*fnval(fnder(SpAry{i_spl,1},2),Spline_map2(i,1));
            Spline_map2(i,5)=y_scale*fnval(fnder(SpAry{i_spl,1},3),Spline_map2(i,1));
            yR=Spline_map2(i,2);
            yRdot=Spline_map2(i,3);
            
            psiR=obj.psiRef(xRdot,yRdot);
            eta_Ref=[xR;yR;psiR];
            t1=t;
            for j=1:1:N % Generate reference upto N steps ahead
                xR=obj.xRef(t1);
                xRdot=obj.xRefdot(t1);
                xRddot=obj.xRefddot(t1);
                xRdddot=obj.xRefdddot(t1);
                
                i_spl=floor((i+j-1)/M)+1;
                Spline_map2(i+j,1)=(i+j-1)*T_plt;
                Spline_map2(i+j,2)=y_scale*fnval(SpAry{i_spl,1},Spline_map2(i+j,1));
                Spline_map2(i+j,3)=y_scale*fnval(fnder(SpAry{i_spl,1},1),Spline_map2(i+j,1));
                Spline_map2(i+j,4)=y_scale*fnval(fnder(SpAry{i_spl,1},2),Spline_map2(i+j,1));
                Spline_map2(i+j,5)=y_scale*fnval(fnder(SpAry{i_spl,1},3),Spline_map2(i+j,1));
                
                yR=Spline_map2(i+j-1,2);
                yRdot=Spline_map2(i+j-1,3);
                yRddot=Spline_map2(i+j-1,4);
                yRdddot=Spline_map2(i+j-1,5);
                
                psiR=obj.psiRef(xRdot,yRdot);
                psiRdot=obj.psiRefdot(xRdot,xRddot,yRdot,yRddot);
                psiRddot=obj.psiRefddot(xRdot,xRddot,xRdddot,yRdot,yRddot,yRdddot);
                
                uR=sqrt(xRdot^2+yRdot^2);
                vR=0;
                rR=psiRdot;
                
                uRdot=(xRdot^2+yRdot^2)^(-1/2)*(xRdot*xRddot+yRdot*yRddot);
                vRdot=0;
                rRdot=psiRddot;
                
                P(:,j)=[xR;yR;psiR;uR;vR;rR;uRdot;vRdot;rRdot];
                
                t1=t1+dt;
            end
        end
        function [ xR ] = xRef( obj,t )
            a=1;
            xR=a*t;
        end
        
        function [ xRdot ] = xRefdot(obj, t )
            a=1;
            xRdot=a;
        end
        
        function [ xddot ] = xRefddot(obj, t )
            xddot=0;
        end
        
        function [ xRdddot ] = xRefdddot( obj,t )
            xRdddot=0;   
        end
        
        function [ psiR ] = psiRef( obj,xRdot,yRdot )
            psiR=atan2(yRdot,xRdot);
        end
        
        function [ psiRdot ] = psiRefdot(obj, xRdot,xRddot,yRdot,yRddot )
            psiRdot=(xRdot*yRddot-yRdot*xRddot)/(xRdot^2+yRdot^2);
        end
        
  
        function [ psiRddot ] = psiRefddot(obj, xRdot,xRddot,xRdddot,yRdot,yRddot,yRdddot )
            psiRddot=(1/(xRdot^2+yRdot^2)^2)*((xRddot*yRddot+xRdot*yRdddot-yRddot*xRddot-yRdot*xRdddot)*(xRdot^2+yRdot^2)-(2*xRdot*xRddot+2*yRdot*yRddot)*(xRdot*yRddot-yRdot*xRddot));
        end
        
        function get_params(obj)
            fprintf('weights: \n')
            disp(obj.weights)
        end
    end
end


