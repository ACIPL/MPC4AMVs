classdef sine2D < Trajectory2D
    properties
        mag
        freq
    end
    methods
        function obj = sine2D(m, f, t0)
            if nargin == 3
                obj.mag = m;
                obj.freq = f;
                obj.t = t0;
                obj = init(obj,t0);
            else
                errID = 'myComponent:inputError';
                msg = 'Three input arguments (freq, mag, t0) should be provided.';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end
        function xR = Xr(obj,t)
            b = obj.freq;
            xR = b*t;
        end
        function xRd = Xrdot(obj,t)
            b = obj.freq;
            xRd = b;
        end
        function xRdd = Xrddot(obj,t)
            xRdd = 0;
        end
        function xRddd = Xrdddot(obj,t)
            xRddd = 0;
        end
        function yR = Yr(obj,t)
            a = obj.mag;
            b = obj.freq;
            yR = a*sin(b*t);
        end
        function yRd = Yrdot(obj,t)
            a = obj.mag;
            b = obj.freq;
            yRd = a*b*cos(b*t);
        end
        function yRdd = Yrddot(obj,t)
            a = obj.mag;
            b = obj.freq;
            yRdd = -a*b^2*sin(b*t);
        end
        function yRddd = Yrdddot(obj,t)
            a = obj.mag;
            b = obj.freq;
            yRddd = -a*b^3*cos(b*t);
        end
       
    end
      
end