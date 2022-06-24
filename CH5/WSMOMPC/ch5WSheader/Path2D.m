classdef (Abstract) Path2D < handle
    properties
        s
        dt
        xR
        xRd
        xRdd
        yR
        yRd
        yRdd
        psiR
    end
    methods(Abstract)
        xR = px(obj,s)
        xRd = px_dot(obj,s)
        xRdd = px_ddot(obj,s)
        yR = py(obj,s)
        yRd = py_dot(obj,s)
        yRdd = py_ddot(obj,s)
    end
    methods
     
        function psiR = PSIr(s)
            psiR = atan(cos(s));
        end
        
        function s_plus = path_evol(s,v_s,dt)
            lambda=0;
            s_dot=-lambda*s + v_s;
            s_plus= s + s_dot*dt;
        end
        
    end
end