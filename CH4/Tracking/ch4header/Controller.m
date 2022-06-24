classdef (Abstract) Controller < handle
    methods(Abstract)
        control = calc_control(obj,ref,state)
        get_params(obj)
    end 
end

