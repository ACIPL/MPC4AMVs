classdef point2D < Trajectory2D 
    methods    
        function obj = point2D(x,y,psi)
            if nargin == 3
                obj.xR = x;
                obj.yR = y;
                obj.psiR = psi;
            else
                errID = 'myComponent:inputError';
                msg = 'Three input arguments (x, y, psi) should be provided.';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end
        
        function xR = Xr(obj,t)
            xR = obj.xR;
        end
        function xRd = Xrdot(obj,t)
            xRd = 0;
        end
        function xRdd = Xrddot(obj,t)
            xRdd = 0;
        end
        function xRddd = Xrdddot(obj,t)
            xRddd = 0;
        end
        function yR = Yr(obj,t)
            yR = obj.yR;
        end
        function yRd = Yrdot(obj,t)
            yRd = 0;
        end
        function yRdd = Yrddot(obj,t)
            yRdd = 0;
        end
        function yRddd = Yrdddot(obj,t)
            yRddd = 0;
        end
    end
      
end

