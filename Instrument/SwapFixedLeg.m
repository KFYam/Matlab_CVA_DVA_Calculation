classdef SwapFixedLeg < Instrument
    %SWAPFIXEDLEG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected, GetAccess = public)
        startDate
        maturityDate
        dayCountConvention
        notional
        fixedRate
        paymentFrequency
        currency
        businessDayConvention
        calendar
        endOfMonthRule
        notionalAtStart
        notionalAtMaturity
        swapType
        payAtMaturity  %if coupon payment delayed to maturity date
        exAmtNotional %whether exchange armotized notional
        
        % leg schedule
        CFAmounts
        CFDates
        CFTimes
        
        %pricing info
        %discount curve handle
%         discountHandle
    end
        
    methods

        function obj = SwapFixedLeg(swapType, startDate, maturityDate, dayCountConvention, notional, ...
                                    fixedRate, paymentFrequency, currency, varargin)
            
            obj = obj@Instrument();
                                
            if nargin < 8
                error(message('SwapFixedLeg:too few input parameters'));
            end
            
            if maturityDate < startDate
                error(message('SwapFixedLeg:maturity date cant be larger than start date'));
            end
            
            if ~ismember(paymentFrequency, [0, 1, 2, 3, 4, 6, 12])
                error(message('SwapFixedLeg:wrong payment frequency'));
            end
            
            if ~ismember(swapType, ['IRS', 'CCS'])
                error(message('SwapFixedLeg:wrong swap type input, only support IRS or CCS'));
            end
            
            obj.swapType = swapType;
            obj.startDate = startDate;
            obj.maturityDate = maturityDate;
            obj.dayCountConvention = dayCountConvention;
            obj.notional =  notional;
            obj.fixedRate = fixedRate;
            obj.paymentFrequency = paymentFrequency;
            obj.currency = currency;            
            
            if ~isempty(varargin) && ischar(varargin{1})
                p = inputParser;
   
                p.addParamValue('Calendar',[]);
                p.addParamValue('BusinessDayConvention',{'actual'},...
                    @(x) all(ismember(x,{'actual','follow','previous','modifiedfollow',...
                    'modifiedprevious'})))
                p.addParamValue('EndOfMonthRule', false);
                p.addParamValue('notionalAtStart', false);
                p.addParamValue('notionalAtMaturity', false);
                p.addParamValue('payAtMaturity', false);
                p.addParamValue('exAmtNotional', false);
                
                try
                    p.parse(varargin{:});
                catch ME
                    newMsg = message('SwapFixedleg:Can not parse');
                    newME = MException(newMsg.Identifier,getString(newMsg));
                    newME = addCause(newME,ME);
                    throw(newME)
                    throw(newME)
                end      
                obj.calendar = p.Results.Calendar;
                obj.businessDayConvention = p.Results.BusinessDayConvention;
                obj.endOfMonthRule = p.Results.EndOfMonthRule;
                obj.notionalAtStart = p.Results.notionalAtStart;
                obj.notionalAtMaturity = p.Results.notionalAtMaturity;
                obj.payAtMaturity = p.Results.payAtMaturity;
                obj.exAmtNotional = p.Results.exAmtNotional;
            else
                obj.calendar = [];
                obj.businessDayConvention = 'actual';
                obj.endOfMonthRule = false;
                obj.notionalAtStart = false;
                obj.notionalAtMaturity = false;
                obj.payAtMaturity = false;
                obj.exAmtNotional = false;
            end %end if
            
            % generate leg schedule
            if length(notional) > 1
                %for armotizing swap
                if ~isvector(startDate) || ~isvector(startDate) || size(startDate, 2) ~=1 || size(maturityDate, 2) ~=1
                    error(message('FloatingLeg:Need vertical vector of dates for armotizing swap'));
                end
                
                if size(notional, 2)~=1 || size(notional, 1) ~= size(startDate, 1) || size(maturityDate, 1) ~= size(startDate, 1)
                    error(message('FloatingLeg:start date, end date and notional profile size need to be the same'));
                end
                
                obj.CFDates = startDate;
                obj.CFDates(length(startDate)+1, 1) = maturityDate(length(maturityDate));
                obj.CFTimes = yearfrac(startDate(1),obj.CFDates,dayCountConvention);
                yfs = obj.CFTimes;
                yfs = diff(yfs);
                obj.CFAmounts(1) = 0;
                obj.CFAmounts(2:length(yfs)+1,1) = fixedRate .* obj.notional .* yfs;
                obj.CFAmounts(end) = obj.CFAmounts(end) + obj.notional(end);
                
                %ser the notional profile
                obj.notional(1) = notional(1);
                obj.notional(2:length(notional)+1,1) = notional;
            else
                [obj.CFAmounts, obj.CFDates] = cfamounts(fixedRate, startDate, maturityDate, ...
                 'Period',paymentFrequency, 'Basis', dayCountConvention, 'EndMonthRule', obj.endOfMonthRule, ...
                 'Face', notional, 'StartDate', startDate, 'BusinessDayConvention', obj.businessDayConvention, 'holiday', holidays);
             
                obj.CFDates = obj.CFDates';
                obj.CFAmounts = obj.CFAmounts';
                obj.CFTimes = yearfrac(obj.CFDates(1),obj.CFDates, dayCountConvention);
                obj.CFAmounts(1)= 0; %no taken into account of accrualed interest
                
                 obj.notional(length(obj.CFTimes),1) = notional;
            end
             
             % adjust the notional amount since cfamount default it to bond
             % convention (principle only paid at the end
             %adjust notional at beginning
             if obj.notionalAtStart == true
                 obj.CFAmounts(1) = obj.CFAmounts(1) - obj.notional(1);
             end
             
             %adjust notiona at end  
             if obj.notionalAtMaturity == false 
                obj.CFAmounts(length(obj.CFAmounts)) = obj.CFAmounts(length(obj.CFAmounts)) - obj.notional(end);
             end
             
             %handle exchange of notional 
             if obj.exAmtNotional == true
                notionalAmt = diff(obj.notional(2:end));
                if size(obj.CFAmounts, 1) == length(obj.CFDates)
                    %all cash flow not expired, apply the notional
                    %exchange from second cashflow
                    obj.CFAmounts(2:end - 1) = obj.CFAmounts(2:end - 1, :) - notionalAmt;
                else
                    obj.CFAmounts(1:end - 1, :) = obj.CFAmounts(1:end - 1, :) - ...
                        notionalAmt(end - length(obj.CFAmounts, 1) + 1:end);
                end
            end
                        
        end %end function
        
        %calculate MTM swap leg price
        function c = calculate(obj, valueDate, env)
            
            if ~isa(env, 'Environment')
                error(message('SwapFixedLeg:Need environment obj'));
            end
                
            modelAssignment = GlobalSetting.getModelAssignment('IR_IRS_MTM');
            
            switch modelAssignment
                case 'ANALYTICAL'
                    c = swapAnalyticalEngine(obj, valueDate, env);                
                otherwise
                    error(message('Fixedleg:pricing model not supported'));
            end %End switch
          
            obj.npv = c;            
        end
        
        %Analytic swap fixed rate pricer
        function c = swapAnalyticalEngine(obj, valueDate, env)
                       
           disCurveName = GlobalSetting.getCurveAssignment(obj.swapType, obj.currency, 'DISCOUNT', '');
           discountHandle = getIRCurveHandle(env, disCurveName);
           
           expired = obj.CFDates < valueDate;
           
           unexpiredCFAmount = obj.CFAmounts;
           unexpiredCFDates = obj.CFDates;
           
           %set expired cash flow to 0
           unexpiredCFAmount(expired) = 0;
           unexpiredCFDates(expired) = discountHandle.referenceDate;
           
           dfs = forwardDiscountFactor(discountHandle, valueDate, unexpiredCFDates');
           c = dfs * unexpiredCFAmount;
        end
        
        %Path pricing
        function [mtmTimeGrid, results] = pathPricer(obj, valueDate, env, pricingGrid)
%             modelAssignment = GlobalSetting.getModelAssignment('IR_IRS_CVA');
%             
%             switch modelAssignment
%                 case 'G2'
%                     modelHandle = getIRModelHandle(env, strcat(obj.currency, '_', modelAssignment));
%                 otherwise
%                     error(message('SwapFixedLeg:model not supported'));
%             end

            if isKey(env.IRModelCollection,  strcat(obj.currency, '_G2'))
                modelHandle = getIRModelHandle(env, strcat(obj.currency, '_G2'));
            else 
                modelHandle = getIRModelHandle(env, strcat(obj.currency, '_CRM'));
            end
            
            timeGrid = modelHandle.timeGrid;
            
            if isempty(timeGrid)
                error(message('SwapFixedLeg:Model need a simulated path'));
            end
            
            nTimeGrid = length(timeGrid);
            nPath = size(modelHandle.xt, 2);
            results = zeros(nTimeGrid, nPath);
            
            %for iPath = 1:nPath
            for iTimeGrid = 1:nTimeGrid
                %if the dela already expired
                if timeGrid(iTimeGrid) >= timeFromReference(modelHandle.irTS, obj.maturityDate)
                    results(iTimeGrid:end, :) = 0;
                    break
                end
               
                %if not
                xt = modelHandle.xt(iTimeGrid, :);
                yt = modelHandle.yt(iTimeGrid, :);  
                
                if length(modelHandle)~= 0
                    st = modelHandle.st(iTimeGrid, :);  
                else
                    st = ones(1, nPath);
                end
                
                unexpired =  timeFromReference(modelHandle.irTS, obj.CFDates) > timeGrid(iTimeGrid); 
                
                dfs = P(modelHandle, timeGrid(iTimeGrid),  obj.CFTimes(unexpired)', xt, yt);
                results(iTimeGrid, :) = (dfs' * obj.CFAmounts(unexpired))'.* st;           
            end
            
            %add exposure at time 0
            currentMTM  = swapAnalyticalEngine(obj, valueDate, env);
            obj.npv = currentMTM; 
            rowOne = zeros(1, nPath);
            rowOne(:) = currentMTM * getFXSpot(env, obj.currency);
            results = [rowOne; results];
            mtmTimeGrid = [0; timeGrid];
            % all mtm profile in USD
            
        end %end of path pricer
    
        function c = yyPrice(obj, vYF, discMatrix)
            %assume a disc matrix is ready to be used here to discount all
            %cash flow and do the libor projection
            [iNbRow,iNbIteration] = size(discMatrix);

            iNbCashFlow = length(obj.CFAmounts);
            c = sum(obj.CFAmounts(iNbCashFlow - iNbRow + 1:end) * ones(1,iNbIteration) .*discMatrix,2)/ iNbIteration;

        end

        function mC = yyPriceForProfile(obj, vYF, discMatrix)
            %assume a disc matrix is ready to be used here to discount all
            %cash flow and do the libor projection
            [iNbRow,iNbIteration] = size(discMatrix);

            iNbCashFlow = length(obj.CFAmounts);
            mC = obj.CFAmounts(iNbCashFlow - iNbRow + 1:end) * ones(1,iNbIteration) .*discMatrix;

        end

    end %end methods

end

