classdef eight_shape2D < Trajectory2D
    properties
        mag
        freq1
        freq2
    end
    methods
        function obj = eight_shape2D(m, f1, f2, t0)
            if nargin == 4
                obj.mag = m;
                obj.freq1 = f1;
                obj.freq2 = f2;
                obj.t = t0;
                obj = init(obj,t0);
            else
                errID = 'myComponent:inputError';
                msg = 'Four input arguments (freq, mag, t0) should be provided.';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end
        function xR = Xr(obj,t)
            a = -obj.mag;
            b = obj.freq1;
            xR = a*sin(b*t);
        end
        function xRd = Xrdot(obj,t)
            a = -obj.mag;
            b = obj.freq1;
            xRd = a*b*cos(b*t);
        end
        function xRdd = Xrddot(obj,t)
            a = -obj.mag;
            b = obj.freq1;
            xRdd = -a*b^2*sin(b*t);
        end
        function xRddd = Xrdddot(obj,t)
            a = -obj.mag;
            b = obj.freq1;
            xRddd = -a*b^3*cos(b*t);
        end
        function yR = Yr(obj,t)
            a = obj.mag;
            b = obj.freq2;
            yR = a*sin(b*t);
        end
        function yRd = Yrdot(obj,t)
            a = obj.mag;
            b = obj.freq2;
            yRd = a*b*cos(b*t);
        end
        function yRdd = Yrddot(obj,t)
            a = obj.mag;
            b = obj.freq2;
            yRdd = -a*b^2*sin(b*t);
        end
        function yRddd = Yrdddot(obj,t)
            a = obj.mag;
            b = obj.freq2;
            yRddd = -a*b^3*cos(b*t);
        end
       
    end
      
end