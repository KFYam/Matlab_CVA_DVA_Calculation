classdef PDTermStructure < TermStructure
    %PDTERMSTRUCTURE Summary of this class goes here
    %   Detailed explanation goes here
    properties
        name
        %attributes
        dates %dates array
        times %times array
        defaultProbs
        survivalProbs
        recoveryRate
        %currency only support backward piecewise constant interpolation
    end
    
    %Public methods
    methods 
        % Constructor
        function obj = PDTermStructure(dayCount, refDate)
            if nargin < 2 
                error('PDTermStructure:Too few parameters');
            end
            
            % call super class constructor
            obj = obj@TermStructure(dayCount, refDate);
        end % end constructor

        % survival probability
        % default probability
        % default density
        % hazard rate
        % (Conditional) PD curve
        function c = CumPD(obj, d)
            c = CumPDByTime(obj, timeFromReference(obj,d));
        end
        
        function c = CumPDByTime(obj, T)
            c = PD(obj, T);
        end
        
        function c = SurvByTime(obj, T)
            c = Survival(obj,T);
        end
    end

    % Abstract methods
    methods (Abstract, Access = protected)
        Survival(obj, T)
        PD(obj, T)
    end  
end

