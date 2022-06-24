classdef (Abstract) Trajectory2D < handle
    properties
        t
        xR
        xRd
        xRdd
        xRddd
        yR
        yRd
        yRdd
        yRddd
        psiR
        psiRd
        psiRdd       
        kappa
    end
    methods(Abstract)
        xR = Xr(obj,t)
        xRd = Xrdot(obj,t)
        xRdd = Xrddot(obj,t)
        xRddd = Xrdddot(obj,t)
        yR = Yr(obj,t)
        yRd = Yrdot(obj,t)
        yRdd = Yrddot(obj,t)
        yRddd = Yrdddot(obj,t)
%         psiR = PSIr(obj,t)
%         psiRd = PSIrdot(obj,t)
%         psiRdd = PSIrddot(obj,t)
    end
    methods
        % curvature
        function k = Kappa(obj,t)
           xd = obj.Xrdot(t);
           xdd = obj.Xrddot(t);
           yd = obj.Yrdot(t);
           ydd = obj.Yrddot(t);
           if (xd~= 0 || yd ~= 0)
               k = (xd*ydd - yd*xdd)/(xd^2 + yd^2)^(3/2);   
           else
               k = 0;
           end
        end
            
        function psiR = PSIr(obj,t)
           xd = obj.Xrdot(t);
           yd = obj.Yrdot(t);
           psiR = atan2(yd,xd);
        end
        
        function psiRd = PSIrdot(obj,t)
           xd = obj.Xrdot(t);
           xdd = obj.Xrddot(t);
           yd = obj.Yrdot(t);
           ydd = obj.Yrddot(t);
           
           if (xd~= 0 || yd ~= 0)
               psiRd = (xd*ydd - yd*xdd)/(xd^2 + yd^2);   
           else
               psiRd = 0;
           end
        end
        
        function psiRdd = PSIrddot(obj,t)
           xd = obj.Xrdot(t);
           xdd = obj.Xrddot(t);
           xddd = obj.Xrdddot(t);
           yd = obj.Yrdot(t);
           ydd = obj.Yrddot(t);
           yddd = obj.Yrdddot(t);
           if (xd~= 0 || yd ~= 0)
               psiRdd =(1/(xd^2+yd^2)^2)*((xdd*ydd+xd*yddd-ydd*xdd-yd*xddd)*(xd^2+yd^2)-(2*xd*xdd+2*yd*ydd)*(xd*ydd-yd*xdd));
           else
               psiRdd = 0;
           end
        end
        
        function obj = update(obj,dt)
            obj.t = obj.t + dt;
            obj.xR = obj.Xr(obj.t);
            obj.xRd = obj.Xrdot(obj.t);
            obj.xRdd = obj.Xrddot(obj.t);
            obj.xRddd = obj.Xrdddot(obj.t);
            obj.yR = obj.Yr(obj.t);
            obj.yRd = obj.Yrdot(obj.t);
            obj.yRdd = obj.Yrddot(obj.t);
            obj.yRddd = obj.Yrdddot(obj.t);
            obj.psiR = obj.PSIr(obj.t);
            obj.psiRd = obj.PSIrdot(obj.t);
            obj.psiRdd = obj.PSIrddot(obj.t);
            obj.kappa = obj.Kappa(obj.t);
        end
        
        function obj = setTime(obj,t)
            obj.t = t;
        end
        
        function obj = init(obj,t0)
            obj = setTime(obj,t0);
            obj = update(obj,0);
        end
        
        function plot(obj,timeInterval,sampleNum)
            assert(length(timeInterval) == 2,'timeInterval is 2x1 vector [t0; tf]');
            t0 = timeInterval(1);
            tf = timeInterval(2);
            assert(t0 <= tf,'timeInterval is 2x1 vector [t0; tf] and t0 <= tf')
            if sampleNum <= 0
                sampleNum = 0;
            else
                sampleNum = floor(sampleNum);
            end
            dt = (tf - t0)/sampleNum;
            Data = zeros(sampleNum + 1, 2);
            Kappa = zeros(sampleNum + 1, 1);
            s = t0;
            for i = 1:1:sampleNum + 1
                Data(i, 1) = obj.Xr(s);
                Data(i, 2) = obj.Yr(s);
                Kappa(i) = obj.Kappa(s);
                s = s + dt;
            end
            fig1=figure(1);
            set(fig1,'name','Reference Trajectory');
            plot(Data(:, 1),Data(:, 2));
            fprintf('Kappa: \n')
            disp(Kappa);
        end
    end
end