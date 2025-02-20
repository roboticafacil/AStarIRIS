classdef (Abstract) duplicableObject < handle
    methods (Abstract)   
        newObj=duplicate(obj);
    end
end