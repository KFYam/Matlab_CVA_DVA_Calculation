classdef ConstantRateModel < Model
%   Model with constant instant rate and 'partially' diffused FX using USD
%   simulated factor
%   r(t) = r(0)
%   st = s0*capit_USD/irTS

    properties
        name
        irTS   %Yield temr structure
                       
        %fx property
        HW2Model_dom
        s0 %fx spot              
        st %simulated FX
        xt
        yt
    end
    
    methods
        function obj = ConstantRateModel(name,irTS,s0, HW2Model_dom)
            if ~isa(irTS, 'YieldTermStructure')
                error(message('HW2ModelConstructor:need a yield curve term structure as input'));
            end
            obj.name = name;
            obj.irTS = irTS;
            obj.s0 = s0;
            obj.HW2Model_dom = HW2Model_dom;            
        end
        
        function c = P(obj, t, T, xt, yt)      
            c = discountFactorByTime(obj.irTS, T)';
            c = repmat(c,1,size(xt,2));
        end
        
        function obj = generateScenario(obj, timeGrid, dw)
            NoTimeGrid = length(timeGrid);
            NoPath = size(obj.HW2Model_dom.capit, 2);
            obj.st = zeros(NoTimeGrid, NoPath);
            obj.xt = zeros(NoTimeGrid, NoPath);
            obj.yt = zeros(NoTimeGrid, NoPath);
            obj.timeGrid = timeGrid;            
            capitPrev = 1;            
            for iTimeGrid = 1 : NoTimeGrid                                
                capit = discountFactorByTime(obj.irTS,timeGrid(iTimeGrid))/capitPrev;
                obj.st(iTimeGrid,:) = (obj.HW2Model_dom.capit(iTimeGrid,:)/capit);
                if iTimeGrid == 1
                    obj.st(iTimeGrid,:) = obj.s0*obj.st(iTimeGrid,:);
                else                    
                    obj.st(iTimeGrid,:) = obj.st(iTimeGrid-1,:).*obj.st(iTimeGrid,:);                    
                end
                capitPrev = discountFactorByTime(obj.irTS,timeGrid(iTimeGrid));
            end % for each time grid point                    
        end
    end
end