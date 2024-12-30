classdef Model < handle
    %SHORTRATEMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %simulation profile
        timeGrid
    end
    
    methods (Abstract)
        generateScenario(obj)
    end
    
end

