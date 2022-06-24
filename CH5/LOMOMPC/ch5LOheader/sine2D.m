classdef sine2D < Path2D
    properties
        mag
        freq
    end
    methods
        function obj = sine2D(m, f)
            if nargin == 3
                obj.mag = m;
                obj.freq = f;
            else
                errID = 'myComponent:inputError';
                msg = 'Three input arguments (freq, mag, s) should be provided.';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end
        function xR = px(obj,s)
            b = obj.freq;
            xR = b*s;
        end
        function xRd = px_dot(obj,s)
            b = obj.freq;
            xRd = b;
        end
        function xRdd = px_ddot(obj,s)
            xRdd = 0;
        end
        function xRddd = px_dddot(obj,s)
            xRddd = 0;
        end
        function yR = py(obj,s)
            a = obj.mag;
            b = obj.freq;
            yR = a*sin(b*s);
        end
        function yRd = py_dot(obj,s)
            a = obj.mag;
            b = obj.freq;
            yRd = a*b*cos(b*s);
        end
        function yRdd = py_ddot(obj,s)
            a = obj.mag;
            b = obj.freq;
            yRdd = -a*b^2*sin(b*s);
        end
        function yRddd = py_dddot(obj,s)
            a = obj.mag;
            b = obj.freq;
            yRddd = -a*b^3*cos(b*s);
        end
    end
      
end