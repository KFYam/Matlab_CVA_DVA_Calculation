classdef HW2Model < Model
%   G2MODEL Summary of this class goes here
%   r(t) = x(t) + y(t) + phi(t)
%   dx(t) = -a*x(t)dt + sigma * dW_1(t), x(0) = 0
%   dy(t) = -b*y(t)dt + eta * dW_2(t), y(0) = 0
%   dW_1(t)*dW_2(t) = rho
%   dst/st = (r_dom(t)-r(t)+(correlation(1)*sigma + correlation(2)*eta)*volTS)dt+volTS*dW_FX(t)
%
    properties
        name
        irTS   %Yield temr structure
        a
        b
        sigma  
        eta    
        volIRYF
        rho
        coeff
                
        %fx property
        isDomesticCCY
        HW2Model_dom
        volYF
        volTS   % determinstic vol term structure
        s0 %fx spot       
        
        %change of measure 
        correlation % correlation matrix with the format [IR1 IR2 IR1_dom, IR2_dom FX]
        %for change of measure only the last column will be used
                
        %calibration instrument
        %%%%include object calibration instrument
        vSwaption
        vCaplet
        vFXCall
        calibType
        %simulated path
        xt
        yt
        st %simulated FX
        phi
        ft  
        capit
        
    end
    
    methods
        function obj = HW2Model(name,irTS, a, b, sigma, eta,volIRYF, rho, coeff, HW2Model_dom,volYF,volTS,s0,correlation,isDomesticCCY,varargin)
            if ~isa(irTS, 'YieldTermStructure')
                error(message('HW2ModelConstructor:need a yield curve term structure as input'));
            end
            
            if ~isscalar(a),error(message('HW2Model:invalidA')),end
            if ~isscalar(b),error(message('HW2Model:invalidB')),end
            if ~isscalar(rho),error(message('HW2Model:invalidRho')),end
            if ~isempty(varargin) && ischar(varargin{1})
                p = inputParser;
                p.addParamValue('calibtype',[]);                
                try
                    p.parse(varargin{:});
                    catch ME
                        newMsg = message('HW2Model:Can not parse');
                        newME = MException(newMsg.Identifier,getString(newMsg));
                        newME = addCause(newME,ME);
                        throw(newME)
                        throw(newME)
                    end      
                obj.calibType = p.Results.calibtype;                    
            else
                obj.calibType = 'none';                    
            end
            obj.name = name;
            obj.a = a;
            obj.b = b;
            obj.irTS = irTS;
            obj.sigma = sigma;
            obj.eta = eta;
            obj.coeff = coeff;
            obj.volIRYF = volIRYF;
            obj.rho = rho;
            obj.isDomesticCCY = isDomesticCCY;
            obj.HW2Model_dom = HW2Model_dom;
            obj.volYF = volYF;            
            obj.volTS = volTS;
            obj.s0 = s0;
            obj.correlation = correlation;
            obj.vSwaption = [];
            obj.vCaplet = [];
            obj.vFXCall = [];
        end
        
        function obj = addSwaption(obj,swaption)
            obj.vSwaption = [obj.vSwaption;swaption];
        end
        function obj = addCaplet(obj, caplet)
            obj.vCaplet = [obj.vCaplet;caplet];
        end
        function obj = addFXCall(obj, fxCall)
            obj.vFXCall = [obj.vFXCall;fxCall];
        end
        
        %calculate zero coupon bond price at time t with simulated xt yt
        %t is scalar, T can be a vertical array, xt, yt are both vertical
        %array
        function c = P(obj, t, T, xt, yt)        
            vDetermP = yyCalcDetermP(obj, t, T)';
            B = @(z,t,T) (1 - exp(-bsxfun(@times,z,(T-t))))/z;                      
            vDetermP = repmat(vDetermP, 1,size(xt, 2));
            c =  vDetermP .* exp(-B(obj.a, t,T)' * xt - B(obj.b, t, T)' * yt);           
        end
%         
        function c = PtT(obj, t, T, xt, yt)
            vDetermP = yyCalcDetermP(obj, t, t+T)';
            B = @(z,T) (1 - exp(-z*T))/z;                      
            c =  vDetermP .* exp(-B(obj.a,T)' * mean(xt,2) - B(obj.b, T)' * mean(yt,2));           
        end
        
        function c = yyCalcDetermP(obj, t, T)
        %compute the deterministc part of the zero coupon price formula
            iNb_t = length(t);
            iNb_T = length(T);
            PM = @(t)discountFactorByTime(obj.irTS, t);
            
            V = @(t,T) yyVtT(obj,t,T);
        
            c = zeros(iNb_t, iNb_T);
            %the change of measure is alreading handled in the stochastic part, no
            %need to do it here.
            for i = 1:iNb_t
                for j = 1:iNb_T
                    %drift = sum(yyCalcDriftIR(obj,t(i),T(j)));
                    c(i,j) = PM(T(j))/PM(t(i))*exp(0.5*V(0,t(i))-0.5*V(0,T(j))+0.5*V(t(i),T(j)));                                
                end
            end
        end
        
        function obj = generateScenario(obj, timeGrid, dw)
            %dw1 dw2 are matrix of single row represent a single gimt e
            %grid point and row represent a path            
            dw1 = dw(:, :, 1);
            dw2 = dw(:, :, 2);
            if ~obj.isDomesticCCY
                dw3 = dw(:, :, 3);
            end
            NoTimeGrid = length(timeGrid);
            NoPath = size(dw1, 2);
            
            if size(dw1, 1) ~= size(dw2, 1) || size(dw1, 2) ~= size(dw2, 2) || NoTimeGrid ~= size(dw2, 1)
               error(message('G2Model:Size of the simulated matrix must be the same'));
            end
            
            obj.xt = zeros(NoTimeGrid, NoPath);
            obj.yt = zeros(NoTimeGrid, NoPath);
            obj.capit = zeros(NoTimeGrid, NoPath);
            obj.st = zeros(NoTimeGrid, NoPath);
            obj.phi = zeros(NoTimeGrid, 1);
            obj.ft = zeros(NoTimeGrid, 1);
            obj.timeGrid = timeGrid;
            tPrev = 0;
            vDiffTimeGrid = diff(timeGrid);
            for iTimeGrid = 1 : NoTimeGrid                
                obj.ft(iTimeGrid) = -(log(discountFactorByTime(obj.irTS, timeGrid(iTimeGrid))) - log(discountFactorByTime(obj.irTS, timeGrid(iTimeGrid)-0.0027)))/0.0027;
                obj.phi(iTimeGrid) = obj.ft(iTimeGrid) + calcPhi(obj,0,timeGrid(iTimeGrid));
                vVarXY = yyVarianceX(obj,tPrev, timeGrid(iTimeGrid));
                vVarCapit = yyVarianceForSimCapit(obj,tPrev, timeGrid(iTimeGrid));                
                [dDetermCapit,dDriftX,dDriftY] = yyCalcDetermForSimIR(obj,tPrev, timeGrid(iTimeGrid));
                VarFX = yyVarStoFX(obj,tPrev, timeGrid(iTimeGrid));
                tPrev = timeGrid(iTimeGrid);
                
                obj.xt(iTimeGrid,:) = sqrt(vVarXY(1))* dw1(iTimeGrid,:) + dDriftX;
                obj.yt(iTimeGrid,:) = sqrt(vVarXY(2))* dw2(iTimeGrid,:) + dDriftY;
                    
                if iTimeGrid == 1
                    obj.capit(iTimeGrid,:) = dDetermCapit*exp(sqrt(vVarCapit(1))*dw1(iTimeGrid,:) + sqrt(vVarCapit(2))*dw2(iTimeGrid,:) ...
                                                        - obj.xt(iTimeGrid,:)/obj.a - obj.yt(iTimeGrid,:)/obj.b);                    
                    if ~obj.isDomesticCCY
                        obj.st(iTimeGrid,:) = (obj.HW2Model_dom.capit(iTimeGrid,:)./obj.capit(iTimeGrid,:)) .* exp(-0.5*VarFX + sqrt(VarFX)*dw3(iTimeGrid,:));
                        obj.st(iTimeGrid,:) = obj.s0*obj.st(iTimeGrid,:);
                    
                    else
                        obj.st = ones(NoTimeGrid,NoPath);
                    end
                else
                    obj.xt(iTimeGrid,:) = obj.xt(iTimeGrid,:) + obj.xt(iTimeGrid-1,:)*exp(-obj.a*vDiffTimeGrid(iTimeGrid-1));             
                    obj.yt(iTimeGrid,:) = obj.yt(iTimeGrid,:) + obj.yt(iTimeGrid-1,:)*exp(-obj.b*vDiffTimeGrid(iTimeGrid-1));             
                    
                    obj.capit(iTimeGrid,:) = dDetermCapit*exp(sqrt(vVarCapit(1))*dw1(iTimeGrid,:) + sqrt(vVarCapit(2))*dw2(iTimeGrid,:) ...
                                                                +(obj.xt(iTimeGrid-1,:) - obj.xt(iTimeGrid,:))/obj.a ...
                                                                +(obj.yt(iTimeGrid-1,:) - obj.yt(iTimeGrid,:))/obj.b);

                    if ~obj.isDomesticCCY
                        obj.st(iTimeGrid,:) = (obj.HW2Model_dom.capit(iTimeGrid,:)./obj.capit(iTimeGrid,:)) .* exp(-0.5*VarFX + sqrt(VarFX)*dw3(iTimeGrid,:));
                        obj.st(iTimeGrid,:) = obj.st(iTimeGrid-1,:).*obj.st(iTimeGrid,:);
                    end
                end
            end % for each time grid point        
            
        end% end simulation       
        function phi = calcPhi(obj,t,T)
            iNbVolDates =length(obj.volIRYF);    
            if iNbVolDates > 1 
                iNumVol =1;    
                phi = 0;
                tPrev = t;            
                while obj.volIRYF(iNumVol)<=t && iNumVol < iNbVolDates
                    iNumVol = iNumVol + 1;
                end
                B = @(x,t, T) (exp(-x*t) - exp(-x*T))/x;

                while obj.volIRYF(iNumVol) <=T && iNumVol <= iNbVolDates                
                    phi = phi + obj.sigma(iNumVol)^2 * B(obj.a,tPrev, obj.volIRYF(iNumVol))^2 /2+ ...
                                        obj.eta(iNumVol)^2 * B(obj.b,tPrev, obj.volIRYF(iNumVol))^2 /2 + ...
                                        obj.rho * obj.sigma(iNumVol) * obj.eta(iNumVol) * B(obj.a,tPrev, obj.volIRYF(iNumVol))* B(obj.b,tPrev, obj.volIRYF(iNumVol));
                    tPrev = obj.volIRYF(iNumVol);
                    if iNumVol == iNbVolDates
                        break;
                    else
                        iNumVol= iNumVol + 1;
                    end        
                end
                if tPrev < T        
                    phi = phi + obj.sigma(iNumVol)^2 * B(obj.a,tPrev, T)^2 /2+ ...
                        obj.eta(iNumVol)^2 * B(obj.b,tPrev, T)^2 /2 + ...
                        obj.rho * obj.sigma(iNumVol) * obj.eta(iNumVol) * B(obj.a,tPrev, T)* B(obj.b,tPrev, T);
                end
            else
               phi = obj.sigma^2 * (1 - exp(-obj.a * timeGrid(iTimeGrid)))^2 /(2 * obj.a^2)+ ...
                        obj.eta^2 * (1 - exp(-obj.b * timeGrid(iTimeGrid)))^2 /(2 * obj.b^2) + ...
                        obj.rho * obj.sigma * obj.eta * (1 - exp(-obj.a * timeGrid(iTimeGrid)))* (1 - exp(-obj.b * timeGrid(iTimeGrid)))/(obj.a*obj.b);
            end
            

        end
        function V = yyVtT(obj,t,T)
            iNumVol =1;    
            V = 0;
            iNbVolDates =length(obj.volIRYF);    
            dDiff2 = T - t;            
            while obj.volIRYF(iNumVol)<=t && iNumVol < iNbVolDates
                iNumVol = iNumVol + 1;
            end
            B = @(x,diff1, diff2) (exp(-x*diff1) - exp(-x*diff2))/x;
            
            while obj.volIRYF(iNumVol) <= T && iNumVol <= iNbVolDates
                dDiff1 = T - obj.volIRYF(iNumVol);                
                V = V + obj.sigma(iNumVol)^2*(dDiff2 - dDiff1 - 2*B(obj.a,dDiff1,dDiff2) + B(2*obj.a,dDiff1,dDiff2))/obj.a/obj.a + ...
                        obj.eta(iNumVol)^2*(dDiff2 - dDiff1  - 2*B(obj.b,dDiff1,dDiff2) + B(2*obj.b,dDiff1,dDiff2))/obj.b/obj.b +...
                        2*obj.rho*obj.sigma(iNumVol)*obj.eta(iNumVol)*(dDiff2 - dDiff1  - B(obj.a,dDiff1,dDiff2) - B(obj.b,dDiff1,dDiff2) +...
                                                                                            B(obj.a+obj.b,dDiff1,dDiff2))/obj.a/obj.b;
                dDiff2 = dDiff1;
                if iNumVol== iNbVolDates
                    break;
                else
                    iNumVol= iNumVol + 1;
                end        
            end
            if dDiff2 >0        
                V = V + obj.sigma(iNumVol)^2*(dDiff2 - 2*B(obj.a,0,dDiff2) + B(2*obj.a,0,dDiff2))/obj.a/obj.a + ...
                        obj.eta(iNumVol)^2*(dDiff2 - 2*B(obj.b,0,dDiff2) + B(2*obj.b,0,dDiff2))/obj.b/obj.b +...
                        2*obj.rho*obj.sigma(iNumVol)*obj.eta(iNumVol)*(dDiff2 - B(obj.a,0,dDiff2) - B(obj.b,0,dDiff2) +...
                                                                                            B(obj.a+obj.b,0,dDiff2))/obj.a/obj.b;
            end
        end
                   
        function [c,driftX,driftY] = yyCalcDetermForSimIR(obj, t, T)
        %compute the deterministc part of the zero coupon price formula                  	
            PM = @(t)discountFactorByTime(obj.irTS, t);
            V = @(t,T) yyVtT(obj,t,T);
            
            dDrift = 0;
            driftX = 0;
            driftY = 0;
            if ~obj.isDomesticCCY 
                %dDrift = yyCalcDriftCapit2(obj,t,T);
                if length(obj.volIRYF) > 1
                    vYF = union(obj.volYF, obj.volIRYF);                    
                    vFX = yyGenerateVolVector(vYF,obj.volTS, obj.volYF);                
                    vSigma = yyGenerateVolVector(vYF,obj.sigma, obj.volIRYF);
                    vEta = yyGenerateVolVector(vYF,obj.eta, obj.volIRYF);     
                else
                    vYF = obj.volYF;                
                    vFX = obj.volTS;                
                    vSigma = obj.sigma*ones(length(obj.volYF));
                    vEta = obj.eta*ones(length(obj.volYF));
                end
                dDrift = yyCalcDriftCapit(obj,vSigma, vEta,vFX,vYF,t,T);
                [driftX,driftY] = yyCalcDriftIR(obj,vSigma, vEta,vFX,vYF, t,T);
            end 
            c = PM(t)/PM(T)*exp(-0.5*V(0,t) + 0.5*V(0,T) + dDrift);
        end
        function drift = yyCalcDriftCapit(obj,sigma,eta,vfx,yf,t,T)         
            iNumVol =1;    
            drift = 0;
            iNbVolDates =length(yf);    
            dDiff2 = T - t;
            vCorrelIRFX = obj.correlation(1:2,5);
            while yf(iNumVol)<=t && iNumVol < iNbVolDates
                iNumVol = iNumVol + 1;
            end
            while yf(iNumVol) <=T && iNumVol <= iNbVolDates
                dDiff1 = T - yf(iNumVol);                
                drift = drift + vCorrelIRFX(1)*sigma(iNumVol) * vfx(iNumVol)*(dDiff2 - dDiff1)/obj.a;
                drift = drift + vCorrelIRFX(2)*eta(iNumVol) * vfx(iNumVol)*(dDiff2 - dDiff1)/obj.b;
                dDiff2 = dDiff1;
                if iNumVol== iNbVolDates
                    break;
                else
                    iNumVol= iNumVol + 1;
                end        
            end
            if dDiff2 >0        
                drift = drift + vCorrelIRFX(1)*sigma(iNumVol) * vfx(iNumVol)*dDiff2/obj.a;        
                drift = drift + vCorrelIRFX(2)*eta(iNumVol) * vfx(iNumVol)*dDiff2/obj.b;                        
            end
        end

        
        function [driftX,driftY] = yyCalcDriftIR(obj,sigma, eta,vfx,yf, t,T)
            %calculation
            if obj.isDomesticCCY
                driftX = 0;
                driftY = 0;
            else
                iNumVol =1;    
                driftX = 0;
                driftY = 0;
                iNbVolDates =length(yf);    
                dDiff2 = T - t;
                vCorrelIRFX = obj.correlation(1:2,5);
                while yf(iNumVol)<=t && iNumVol < iNbVolDates
                    iNumVol = iNumVol + 1;
                end
                while yf(iNumVol) <=T && iNumVol <= iNbVolDates
                    dDiff1 = T - yf(iNumVol);
                    driftX = driftX +  vCorrelIRFX(1)*sigma(iNumVol) * vfx(iNumVol)*(exp(-obj.a*dDiff1) - exp(-obj.a*dDiff2))/obj.a;        
                    driftY = driftY +  vCorrelIRFX(2)*eta(iNumVol)* vfx(iNumVol)*(exp(-obj.b*dDiff1) - exp(-obj.b*dDiff2))/obj.b;        
                    dDiff2 = dDiff1;
                    if iNumVol== iNbVolDates
                        break;
                    else
                        iNumVol= iNumVol + 1;
                    end        
                end
                if dDiff2 >0 
                    driftX = driftX +  vCorrelIRFX(1)*sigma(iNumVol) * vfx(iNumVol)*(1 - exp(-obj.a*dDiff2))/obj.a;        
                    driftY = driftY +  vCorrelIRFX(2)*eta(iNumVol)* vfx(iNumVol)*(1 - exp(-obj.b*dDiff2))/obj.b;                        
                end
            end
        end
        function var = yyVarStoFX(obj,t, T)
            if obj.isDomesticCCY
                var = 0;
            else
                iNumVol =1;    
                var = 0;
                iNbVolDates = length(obj.volYF);                    
                dTprev = t;
                while obj.volYF(iNumVol)<=t && iNumVol < iNbVolDates
                    iNumVol = iNumVol + 1;
                end
                while obj.volYF(iNumVol) <=T && iNumVol <= iNbVolDates
                    var = var + obj.volTS(iNumVol)^2 * (obj.volYF(iNumVol)-dTprev);
                    dTprev = obj.volYF(iNumVol);
                    if iNumVol== iNbVolDates
                        break;
                    else
                        iNumVol= iNumVol + 1;
                    end        
                end
                if dTprev < T
                    var = var + obj.volTS(iNumVol)^2 * (T - dTprev);
                end                
            end
        end 
        function vVarX = yyVarianceX(obj,t, T)
            %compute the variance and covariance of the IR factors X1 and X2 
            iNumVol =1;    
            vVarX = [0;0;0];
            iNbVolDates = length(obj.volIRYF);
            if iNbVolDates > 1 
                dDiff2 = T - t;
                while obj.volIRYF(iNumVol)<=t && iNumVol < iNbVolDates
                    iNumVol = iNumVol + 1;
                end
                while obj.volIRYF(iNumVol) <=T && iNumVol <= iNbVolDates
                    dDiff1= T - obj.volIRYF(iNumVol);
                    vVarX(1) = vVarX(1) + obj.sigma(iNumVol)*obj.sigma(iNumVol)*(exp(-2*obj.a*dDiff1) - exp(-2*obj.a*dDiff2))/(2*obj.a);
                    vVarX(2) = vVarX(2) + obj.eta(iNumVol)*obj.eta(iNumVol)*(exp(-2*obj.b*dDiff1) - exp(-2*obj.b*dDiff2))/(2*obj.b);
                    vVarX(3) = vVarX(3) + 2 * obj.rho * obj.sigma(iNumVol)*obj.eta(iNumVol)*(exp(-(obj.a+obj.b)*dDiff1) - exp(-(obj.a+obj.b)*dDiff2))/(obj.a+obj.b);            
                    dDiff2 = dDiff1;
                    if iNumVol== iNbVolDates
                        break;
                    else
                        iNumVol= iNumVol + 1;
                    end        
                end
                if dDiff2 > 0
                    vVarX(1) = vVarX(1) + obj.sigma(iNumVol)*obj.sigma(iNumVol)*(1 - exp(-2*obj.a*dDiff2))/(2*obj.a);
                    vVarX(2) = vVarX(2) + obj.eta(iNumVol)*obj.eta(iNumVol)*(1 - exp(-2*obj.b*dDiff2))/(2*obj.b);
                    vVarX(3) = vVarX(3) + 2 * obj.rho * obj.sigma(iNumVol)*obj.eta(iNumVol)*(1 - exp(-(obj.a+obj.b)*dDiff2))/(obj.a+obj.b);                    
                end                
            else
                diff = T - t;
                vVarX(1) = obj.sigma*obj.sigma*(1 - exp(-2*obj.a*diff))/(2*obj.a);
                vVarX(2) = obj.eta*obj.eta*(1 - exp(-2*obj.b*diff))/(2*obj.b);
                vVarX(3) = 2 * obj.rho * obj.sigma*obj.eta*(1 - exp(-(obj.a+obj.b)*diff))/(obj.a+obj.b);            
            end
        end
        function varForCapit = yyVarianceForSimCapit(obj,t, T)
            %compute the variance of the stochastic term for the capit factor
            iNumVol =1;    
            varForCapit = [0;0;0];
            iNbVolDates = length(obj.volIRYF);
            if iNbVolDates > 1 
                tPrev = t;            
                while obj.volIRYF(iNumVol)<=t && iNumVol < iNbVolDates
                    iNumVol = iNumVol + 1;
                end
                while obj.volIRYF(iNumVol) <=T && iNumVol <= iNbVolDates
                    diff = obj.volIRYF(iNumVol) - tPrev;
                    varForCapit(1) = varForCapit(1) + obj.sigma(iNumVol)*obj.sigma(iNumVol)/obj.a/obj.a * diff;
                    varForCapit(2) = varForCapit(2) + obj.eta(iNumVol)*obj.eta(iNumVol)/obj.b/obj.b * diff;
                    varForCapit(3) = varForCapit(3) + 2*obj.rho*obj.sigma(iNumVol)*obj.eta(iNumVol)/obj.a/obj.b * diff;             
                    tPrev = obj.volIRYF(iNumVol);
                    if iNumVol== iNbVolDates
                        break;
                    else
                        iNumVol= iNumVol + 1;
                    end        
                end
                if tPrev < T
                    diff = T - tPrev;
                    varForCapit(1) = varForCapit(1) + obj.sigma(iNumVol)*obj.sigma(iNumVol)/obj.a/obj.a * diff;
                    varForCapit(2) = varForCapit(2) + obj.eta(iNumVol)*obj.eta(iNumVol)/obj.b/obj.b * diff;
                    varForCapit(3) = varForCapit(3) + 2*obj.rho*obj.sigma(iNumVol)*obj.eta(iNumVol)/obj.a/obj.b * diff;            
                end                
            else
                diff = T - t;
                varForCapit(1) = obj.sigma*obj.sigma/obj.a/obj.a * diff;
                varForCapit(2) = obj.eta*obj.eta/obj.b/obj.b * diff;
                varForCapit(3) = 2*obj.rho*obj.sigma*obj.eta/obj.a/obj.b * diff;                            
            end
            
        end
        function obj = calibrate(obj)
            if length(obj.vSwaption) + length(obj.vCaplet) < 2
                error(message('G2Model:calibration needs at least 3 instrument'));
            end
            options = optimset('disp','iter-detailed','TolFun',1e-6,'UseParallel','always','Diagnostics', 'on','MaxFunEvals',10000);    
            vWeigth = [obj.vSwaption.weight]'/length([obj.vSwaption.weight]);
            isConstantVol = false;
            if length(obj.volIRYF) == 0
                obj.volIRYF = unique([obj.vSwaption.maturityYF]'.*(vWeigth>0));
                obj.volIRYF = obj.volIRYF(obj.volIRYF>0);
                isConstantVol = true;
            end
            vMarketPrice = [obj.vSwaption.marketPrice]';                   
            isCapletIncluded = false;
            if length(obj.vCaplet)>1
                isCapletIncluded = true;
                if ~isConstantVol 
                    vtmp = unique([obj.vCaplet.maturityYF]'.*([obj.vCaplet.weight]'>0));
                    vtmp = vtmp(vtmp>0);                    
                    obj.volIRYF = unique(union(obj.volIRYF, vtmp));
                    clear vtmp;
                end                
                vWeigth = [vWeigth ;[obj.vCaplet.weight]'/length([obj.vCaplet.weight])]; 
                vMarketPrice = [vMarketPrice;[obj.vCaplet.marketPrice]'];                       
            end
           
            %%Optimisation  - depending on the calibration type
            % Initial point and constraint          
            iNbVol = length(obj.volIRYF);            
            if strcmp(obj.calibType, 'MeanReversionInput')
                x0 = [0.001*ones(iNbVol,1);0.001*ones(iNbVol,1)];
                ub = [ones(iNbVol,1);ones(iNbVol,1)];
                lb = [0.0001*ones(iNbVol,1);0.0001*ones(iNbVol,1)];         
                % Objective function
                if isCapletIncluded 
                    HW2Fobjfun = @(x) vWeigth .*(1 - [ModelPrice2(obj.vSwaption,obj.a,obj.b,x(1:iNbVol),x(iNbVol+1:end),obj.rho,obj.volIRYF);ModelPrice(obj.vCaplet,obj,x)]./vMarketPrice);
                else                
                    HW2Fobjfun = @(x) vWeigth .*(1 - ModelPrice2(obj.vSwaption,obj.a,obj.b,x(1:iNbVol),x(iNbVol+1:end),obj.rho,obj.volIRYF)./vMarketPrice);                   
                end
            elseif strcmp(obj.calibType, 'CoeffInput')
                x0 = [0.08;0.005;0.001*ones(iNbVol,1)];
                ub = [1;0.1;ones(iNbVol,1)];
                lb = [0.05;0.00001;0.0001*ones(iNbVol,1)];
                % Objective function
                if isCapletIncluded 
                    HW2Fobjfun = @(x) vWeigth .*(1 - [ModelPrice2(obj.vSwaption,x(1),x(2),x(3:end),obj.coeff*x(3:end),obj.rho,obj.volIRYF);ModelPrice(obj.vCaplet,obj,x)]./vMarketPrice);       
                else
                    HW2Fobjfun = @(x) vWeigth .*(1 - ModelPrice2(obj.vSwaption,x(1),x(2),x(3:end),obj.coeff*x(3:end),obj.rho,obj.volIRYF)./vMarketPrice);
                end
            elseif strcmp(obj.calibType, 'all')
                x0 = 0.001*ones(iNbVol,1);
                ub = ones(iNbVol,1);ones(iNbVol,1);
                lb = 0.0001*ones(iNbVol,1);
                % Objective function
                if isCapletIncluded 
                    HW2Fobjfun = @(x) vWeigth .*(1 - [ModelPrice2(obj.vSwaption,x(1),obj.a,obj.b,x,obj.coeff*x,obj.rho,obj.volIRYF);ModelPrice(obj.vCaplet,obj,x)]./vMarketPrice);       
                else                
                    HW2Fobjfun = @(x) vWeigth .*(1 - ModelPrice2(obj.vSwaption,obj.a,obj.b,x,obj.coeff*x,obj.rho,obj.volIRYF)./vMarketPrice);                  
                end
            else
                x0 = [0.08;0.005;0.001*ones(iNbVol,1);0.001*ones(iNbVol,1)];
                ub = [1;0.1;ones(iNbVol,1);ones(iNbVol,1)];
                lb = [0.05;0.00001;0.0001*ones(iNbVol,1);0.0001*ones(iNbVol,1)];         
                % Objective function
                if isCapletIncluded 
                    HW2Fobjfun = @(x) vWeigth .*(1 - [ModelPrice2(obj.vSwaption,x(1),x(2),x(3:3+iNbVol-1),x(3+iNbVol:end),obj.rho,obj.volIRYF);ModelPrice(obj.vCaplet,obj,x)]./vMarketPrice);       
                else                
                    HW2Fobjfun = @(x) vWeigth .*(1 - ModelPrice2(obj.vSwaption,x(1),x(2),x(3:3+iNbVol-1),x(3+iNbVol:end),obj.rho,obj.volIRYF)./vMarketPrice);
                end
            end
            
            fprintf('IR Calibration %s \n', char(obj.name));
            tic;      
            % Min finding
            HW2F_Param = lsqnonlin(HW2Fobjfun,x0,lb,ub,options);
            time = toc;
            time = time/60;
            fprintf('End IR calibration %s, time spent = %fmin \n', char(obj.name),time);
            
            if strcmp(obj.calibType, 'MeanReversionInput')
                obj.sigma = HW2F_Param(1:iNbVol);
                obj.eta = HW2F_Param(iNbVol+1:end);  
            elseif strcmp(obj.calibType, 'CoeffInput')
                obj.a = HW2F_Param(1);
                obj.b = HW2F_Param(2);
                obj.sigma = HW2F_Param(3:end);
                obj.eta = obj.coeff*HW2F_Param(3:end);
            elseif strcmp(obj.calibType, 'all')
                obj.sigma = HW2F_Param;
                obj.eta = obj.coeff*HW2F_Param;  
            else
                obj.a = HW2F_Param(1);
                obj.b = HW2F_Param(2);
                obj.sigma = HW2F_Param(3:3+iNbVol-1);
                obj.eta = HW2F_Param(3+iNbVol:end);  
            end
            
            if ~obj.isDomesticCCY
                %%FX Calibration
                iNbFXVol = length(obj.volYF);                
                if iNbFXVol  == 0
                    obj.volYF = unique([obj.vFXCall.maturityYF]'.*([obj.vFXCall.weight]'>0));
                    obj.volYF = obj.volYF(obj.volYF>0);
                    iNbFXVol = length(obj.volYF);
                end
                
                lb = 0.0001*ones(iNbFXVol,1);
                ub = ones(iNbFXVol,1);
                x0 = 0.05*ones(iNbFXVol,1);

                FXobjfun = @(x) [obj.vFXCall.weight]'.*(1 - (ModelPrice(obj.vFXCall,obj,x)./[obj.vFXCall.marketPrice]'));

                fprintf('FX Calibration %s \n', char(obj.name));
                tic;                                        
                FX_Param = lsqnonlin(FXobjfun,x0,lb,ub,options);
                time = toc;
                time = time/60;
                fprintf('End FX calibration %s, time spent = %fmin \n', strcat(char(obj.name),char(obj.HW2Model_dom.name)),time);

                obj.volTS = FX_Param;
            end
        end % end of calibration               
        
        function var = yyVarFX(obj,t, T)
            if obj.isDomesticCCY
                var = 0;
            else
                iNumVol =1;    
                var = 0;
                vCorrelIRFX = obj.correlation(1:4,5);
                dTprev = t;
                cov = [0;0];
                vYF = obj.volYF;
                if length(obj.volIRYF) > 1
                    vYF = union(vYF, obj.volIRYF);
                end
                if length(obj.HW2Model_dom.volIRYF) > 1
                	vYF = union(vYF, obj.HW2Model_dom.volIRYF);
                end
                iNbVolDates = length(vYF);    
                vFX = yyGenerateVolVector(vYF,obj.volTS, obj.volYF);                
                if length(obj.volIRYF) > 1
                    vSigma = yyGenerateVolVector(vYF,obj.sigma, obj.volIRYF);
                    vEta = yyGenerateVolVector(vYF,obj.eta, obj.volIRYF);                
                else
                    vSigma = obj.sigma*ones(iNbVolDates ,1);
                    vEta = obj.eta*ones(iNbVolDates ,1);
                end
                if length(obj.HW2Model_dom.volIRYF) > 1
                	vSigma_dom = yyGenerateVolVector(vYF,obj.HW2Model_dom.sigma, obj.HW2Model_dom.volIRYF);
                    vEta_dom = yyGenerateVolVector(vYF,obj.HW2Model_dom.eta, obj.HW2Model_dom.volIRYF);                
                else
                    vSigma_dom = obj.HW2Model_dom.sigma*ones(iNbVolDates ,1);
                    vEta_dom = obj.HW2Model_dom.eta*ones(iNbVolDates ,1);                
                end
                
                dDiff2 = T - t;                
                B = @(x,Diff1,Diff2) (exp(-x*Diff1) - exp(-x*Diff2))/x;
                
                while vYF(iNumVol)<=t && iNumVol < iNbVolDates
                    iNumVol = iNumVol + 1;
                end
                while vYF(iNumVol) <= T && iNumVol <= iNbVolDates
                    dDiff1 = T - vYF(iNumVol);
                    cov(1) = cov(1) + vCorrelIRFX(1)*vSigma(iNumVol)*vFX(iNumVol)*(dDiff2 - dDiff1 - B(obj.a,dDiff1,dDiff2))/obj.a;
                    cov(1) = cov(1) + vCorrelIRFX(2)*vEta(iNumVol)*vFX(iNumVol)*(dDiff2 - dDiff1 - B(obj.b,dDiff1,dDiff2))/obj.b;            
                    
                    cov(2) = cov(2) + vCorrelIRFX(3)*vSigma_dom(iNumVol) * vFX(iNumVol)*(dDiff2 - dDiff1 - B(obj.HW2Model_dom.a,dDiff1,dDiff2))/obj.HW2Model_dom.a;
                    cov(2) = cov(2) + vCorrelIRFX(4)*vEta_dom(iNumVol) *vFX(iNumVol)*(dDiff2 - dDiff1 - B(obj.HW2Model_dom.b,dDiff1,dDiff2))/obj.HW2Model_dom.b;            
                    
                    dDiff2 = dDiff1;
                    
                    var = var + vFX(iNumVol)^2 * (vYF(iNumVol) - dTprev);
                    dTprev = vYF(iNumVol);
                    if iNumVol == iNbVolDates
                        break;
                    else
                        iNumVol= iNumVol + 1;
                    end        
                end
                if dTprev < T
                    var = var + vFX(iNumVol)^2 * (T - dTprev);
                    cov(1) = cov(1) + vCorrelIRFX(1)*vSigma(iNumVol)*vFX(iNumVol)*(dDiff2 - B(obj.a,0,dDiff2))/obj.a;
                    cov(1) = cov(1) + vCorrelIRFX(2)*vEta(iNumVol)*vFX(iNumVol)*(dDiff2 - B(obj.b,0,dDiff2))/obj.b;            
                    
                    cov(2) = cov(2) + vCorrelIRFX(3)*vSigma_dom(iNumVol) * vFX(iNumVol)*(dDiff2 - B(obj.HW2Model_dom.a,0,dDiff2))/obj.HW2Model_dom.a;
                    cov(2) = cov(2) + vCorrelIRFX(4)*vEta_dom(iNumVol) * vFX(iNumVol)*(dDiff2 -  B(obj.HW2Model_dom.b,0,dDiff2))/obj.HW2Model_dom.b;            

                end
                %%the variance    
                covIRIR = CovarianceIRIR(obj,t, T, obj.HW2Model_dom);                
                varIRDom = CovarianceIRIR(obj.HW2Model_dom,t, T,obj.HW2Model_dom);
                varIRFor = CovarianceIRIR(obj,t, T, obj);                            

                var = var + varIRDom + varIRFor - 2*covIRIR + 2*cov(2)- 2*cov(1);
            end
        end 
        

        %%covariance between IR volatilities
        function cov = CovarianceIRIR(obj,t, T, H2Fobj_i)
            if strcmp(obj.name,H2Fobj_i.name)
                %variance of the economy
                correl = [1 obj.rho; obj.rho 1];
                %correl = obj.correlation(1:2,1:2);
            else
                %covariance foreign-dom
                correl = obj.correlation(1:2,3:4);
            end
            B = @(x,tPrev,tNext, T) (exp(-x*(T-tNext)) - exp(-x*(T-tPrev)))/x;
            vYF = [];
            if length(obj.volIRYF) > 1
                vYF = union(vYF, obj.volIRYF);
            end
            if length(H2Fobj_i.volIRYF) > 1
                vYF = union(vYF, H2Fobj_i.volIRYF);
            end
            iNbVolDates = max(1,length(vYF));            
            if length(obj.volIRYF) > 1
                vSigma = yyGenerateVolVector(vYF,obj.sigma, obj.volIRYF);
                vEta = yyGenerateVolVector(vYF,obj.eta, obj.volIRYF);                
            else
                vSigma = obj.sigma*ones(iNbVolDates ,1);
                vEta = obj.eta*ones(iNbVolDates ,1);
            end
            if length(H2Fobj_i.volIRYF) > 1
                vSigma_i = yyGenerateVolVector(vYF,H2Fobj_i.sigma, H2Fobj_i.volIRYF);
                vEta_i = yyGenerateVolVector(vYF,H2Fobj_i.eta, H2Fobj_i.volIRYF);                
            else
                vSigma_i = H2Fobj_i.sigma*ones(iNbVolDates ,1);
                vEta_i = H2Fobj_i.eta*ones(iNbVolDates ,1);                
            end
           
            iNumVol = 1;
            iNbVolDates = length(vYF);
            tPrev = t;
            while vYF(iNumVol)<=t && iNumVol < iNbVolDates
                iNumVol = iNumVol + 1;
            end
            while vYF(iNumVol) <=T && iNumVol <= iNbVolDates
                cov = (correl(1,1)/obj.a/H2Fobj_i.a)*vSigma(iNumVol)*vSigma_i(iNumVol)*(vYF(iNumVol) - tPrev - B(obj.a,tPrev,vYF(iNumVol),T)- B(H2Fobj_i.a,tPrev,vYF(iNumVol),T) ...
                                                                                        + B(obj.a + H2Fobj_i.a,tPrev,vYF(iNumVol),T));
                cov = cov + (correl(1,2)/obj.a/H2Fobj_i.b)*vSigma(iNumVol)*vEta_i(iNumVol)*(vYF(iNumVol) - tPrev - B(obj.a,tPrev,vYF(iNumVol),T) - B(H2Fobj_i.b,tPrev,vYF(iNumVol),T) ...
                                                                                        + B(obj.a + H2Fobj_i.b,tPrev,vYF(iNumVol),T));
                cov = cov + (correl(2,1)/obj.b/H2Fobj_i.a)*vEta(iNumVol) * vSigma_i(iNumVol) *(vYF(iNumVol) - tPrev - B(obj.b,tPrev,vYF(iNumVol),T) - B(H2Fobj_i.a,tPrev,vYF(iNumVol),T) ...
                                                                                        + B(obj.b + H2Fobj_i.a,tPrev,vYF(iNumVol),T));
                cov = cov + (correl(2,2)/obj.b/H2Fobj_i.b)*vEta(iNumVol) *vEta_i(iNumVol) *(vYF(iNumVol) - tPrev - B(obj.b,tPrev,vYF(iNumVol),T) - B(H2Fobj_i.b,tPrev,vYF(iNumVol),T) ...
                                                                                        + B(obj.b + H2Fobj_i.b,tPrev,vYF(iNumVol),T));            
                tPrev = vYF(iNumVol);
                if iNumVol== iNbVolDates
                    break;
                else
                    iNumVol= iNumVol + 1;
                end        
            end
            if tPrev < T
                cov = (correl(1,1)/obj.a/H2Fobj_i.a)*vSigma(iNumVol)*vSigma_i(iNumVol)*(T - tPrev - B(obj.a,tPrev,T,T)- B(H2Fobj_i.a,tPrev,T,T) ...
                                                                                        + B(obj.a + H2Fobj_i.a,tPrev,T,T));
                cov = cov + (correl(1,2)/obj.a/H2Fobj_i.b)*vSigma(iNumVol)*vEta_i(iNumVol)*(T - tPrev - B(obj.a,tPrev,T,T) - B(H2Fobj_i.b,tPrev,T,T) ...
                                                                                        + B(obj.a + H2Fobj_i.b,tPrev,T,T));
                cov = cov + (correl(2,1)/obj.b/H2Fobj_i.a)*vEta(iNumVol) * vSigma_i(iNumVol) *(T - tPrev - B(obj.b,tPrev,T,T) - B(H2Fobj_i.a,tPrev,T,T) ...
                                                                                        + B(obj.b + H2Fobj_i.a,tPrev,T,T));
                cov = cov + (correl(2,2)/obj.b/H2Fobj_i.b)*vEta(iNumVol) *vEta_i(iNumVol) *(T - tPrev - B(obj.b,tPrev,T,T) - B(H2Fobj_i.b,tPrev,T,T) ...
                                                                                        + B(obj.b + H2Fobj_i.b,tPrev,T,T));                            
            end
        end
        
        function cov = yyCovarianceIRFX(obj,t, T)            
            if obj.isDomesticCCY 
                cov = 0;
            else
                vCorrelIRFX = obj.correlation(1:4,5);
                iNumVol =1;    
                cov = [0;0];
                
                vYF = obj.volYF;
                if length(obj.volIRYF) > 1
                    vYF = union(vYF, obj.volIRYF);
                end
                if length(obj.HW2Model_dom.volIRYF) > 1
                	vYF = union(vYF, obj.HW2Model_dom.volIRYF);
                end
                iNbVolDates = length(vYF);    
                vFX = yyGenerateVolVector(vYF,obj.volTS, obj.volYF);                
                if length(obj.volIRYF) > 1
                    vSigma = yyGenerateVolVector(vYF,obj.sigma, obj.volIRYF);
                    vEta = yyGenerateVolVector(vYF,obj.eta, obj.volIRYF);                
                else
                    vSigma = obj.sigma*ones(iNbVolDates ,1);
                    vEta = obj.eta*ones(iNbVolDates ,1);
                end
                if length(obj.HW2Model_dom.volIRYF) > 1
                	vSigma_dom = yyGenerateVolVector(vYF,obj.HW2Model_dom.sigma, obj.HW2Model_dom.volIRYF);
                    vEta_dom = yyGenerateVolVector(vYF,obj.HW2Model_dom.eta, obj.HW2Model_dom.volIRYF);                
                else
                    vSigma_dom = obj.HW2Model_dom.sigma*ones(iNbVolDates ,1);
                    vEta_dom = obj.HW2Model_dom.eta*ones(iNbVolDates ,1);                
                end
                
                iNbVolDates = length(vYF);    
                dDiff2 = T - t;                
                B = @(x,Diff1,Diff2) (exp(-x*Diff1) - exp(-x*Diff2))/x;
                
                while vYF(iNumVol)<=t && iNumVol < iNbVolDates
                    iNumVol = iNumVol + 1;
                end
                while vYF(iNumVol) <=T && iNumVol <= iNbVolDates
                    dDiff1 = T - obj.volYF(iNumVol);
                    cov(1) = cov(1) + vCorrelIRFX(1)*vSigma(iNumVol)*vFX(iNumVol)*(dDiff2 - dDiff1 - B(obj.a,dDiff1,dDiff2))/obj.a;
                    cov(1) = cov(1) + vCorrelIRFX(2)*vEta(iNumVol)*vFX(iNumVol)*(dDiff2 - dDiff1 - B(obj.b,dDiff1,dDiff2))/obj.b;            
                    
                    cov(2) = cov(2) + vCorrelIRFX(3)*vSigma_dom(iNumVol) * vFX(iNumVol)*(dDiff2 - dDiff1 - B(obj.HW2Model_dom.a,dDiff1,dDiff2))/obj.HW2Model_dom.a;
                    cov(2) = cov(2) + vCorrelIRFX(4)*vEta_dom(iNumVol) *vFX(iNumVol)*(dDiff2 - dDiff1 - B(obj.HW2Model_dom.b,dDiff1,dDiff2))/obj.HW2Model_dom.b;            
                    
                    dDiff2 = dDiff1;
                    if iNumVol== iNbVolDates
                        break;
                    else
                        iNumVol= iNumVol + 1;
                    end        
                end
                if dDiff2 >0
                    cov(1) = cov(1) + vCorrelIRFX(1)*vSigma*vFX(iNumVol)*(dDiff2 - B(obj.a,0,dDiff2))/obj.a;
                    cov(1) = cov(1) + vCorrelIRFX(2)*vEta*vFX(iNumVol)*(dDiff2 - B(obj.b,0,dDiff2))/obj.b;            
                    
                    cov(2) = cov(2) + vCorrelIRFX(3)*vSigma_dom * vFX(iNumVol)*(dDiff2 - B(obj.HW2Model_dom.a,0,dDiff2))/obj.HW2Model_dom.a;
                    cov(2) = cov(2) + vCorrelIRFX(4)*vEta_dom *vFX(iNumVol)*(dDiff2 -  B(obj.HW2Model_dom.b,0,dDiff2))/obj.HW2Model_dom.b;            
                end
            end
        end       
    end %end all methods
end

% 
% function drift = yyCalcDriftCapit2(obj,t,T)         
%             iNumVol =1;    
%             drift = 0;
%             iNbVolDates =length(obj.volYF);    
%             dDiff2 = T - t;
%             vCorrelIRFX = obj.correlation(1:2,5);
%             while obj.volYF(iNumVol)<=t && iNumVol < iNbVolDates
%                 iNumVol = iNumVol + 1;
%             end
%             while obj.volYF(iNumVol) <=T && iNumVol <= iNbVolDates
%                 dDiff1 = T - obj.volYF(iNumVol);                
%                 drift = drift + vCorrelIRFX(1)*obj.sigma * obj.volTS(iNumVol)*(dDiff2 - dDiff1)/obj.a;        
%                 drift = drift + vCorrelIRFX(2)*obj.eta * obj.volTS(iNumVol)*(dDiff2 - dDiff1)/obj.b;                        
%                 dDiff2 = dDiff1;
%                 if iNumVol== iNbVolDates
%                     break;
%                 else
%                     iNumVol= iNumVol + 1;
%                 end        
%             end
%             if dDiff2 >0        
%                 drift = drift + vCorrelIRFX(1)*obj.sigma * obj.volTS(iNumVol)*dDiff2/obj.a;        
%                 drift = drift + vCorrelIRFX(2)*obj.eta * obj.volTS(iNumVol)*dDiff2/obj.b;                        
%             end
%         end
%         
%         function var = yyVarFX(obj,t, T)
%             if obj.isDomesticCCY
%                 var = 0;
%             else
%                 iNumVol =1;    
%                 var = 0;
%                 iNbVolDates = length(obj.volYF);    
%                 vCorrelIRFX = obj.correlation(1:4,5);
%                 dTprev = t;
%                 cov = [0;0];
%                 dDiff2 = T - t;
%                 B = @(x,Diff1,Diff2) (exp(-x*Diff1) - exp(-x*Diff2))/x;
%                 while obj.volYF(iNumVol)<=t && iNumVol < iNbVolDates
%                     iNumVol = iNumVol + 1;
%                 end
%                 while obj.volYF(iNumVol) <=T && iNumVol <= iNbVolDates
%                     dDiff1 = T - obj.volYF(iNumVol);
%                     cov(1) = cov(1) + vCorrelIRFX(1)*obj.sigma*obj.volTS(iNumVol)*(dDiff2 - dDiff1 - B(obj.a,dDiff1,dDiff2))/obj.a;
%                     cov(1) = cov(1) + vCorrelIRFX(2)*obj.eta*obj.volTS(iNumVol)*(dDiff2 - dDiff1 - B(obj.b,dDiff1,dDiff2))/obj.b;            
%                     
%                     cov(2) = cov(2) + vCorrelIRFX(3)*obj.HW2Model_dom.sigma * obj.volTS(iNumVol)*(dDiff2 - dDiff1 - B(obj.HW2Model_dom.a,dDiff1,dDiff2))/obj.HW2Model_dom.a;
%                     cov(2) = cov(2) + vCorrelIRFX(4)*obj.HW2Model_dom.eta *obj.volTS(iNumVol)*(dDiff2 - dDiff1 - B(obj.HW2Model_dom.b,dDiff1,dDiff2))/obj.HW2Model_dom.b;            
%                     
%                     dDiff2 = dDiff1;
%                     
%                     var = var + obj.volTS(iNumVol)^2 * (obj.volYF(iNumVol)-dTprev);
%                     dTprev = obj.volYF(iNumVol);
%                     if iNumVol== iNbVolDates
%                         break;
%                     else
%                         iNumVol= iNumVol + 1;
%                     end        
%                 end
%                 if dTprev < T
%                     var = var + obj.volTS(iNumVol)^2 * (T - dTprev);
%                     cov(1) = cov(1) + vCorrelIRFX(1)*obj.sigma*obj.volTS(iNumVol)*(dDiff2 - B(obj.a,0,dDiff2))/obj.a;
%                     cov(1) = cov(1) + vCorrelIRFX(2)*obj.eta*obj.volTS(iNumVol)*(dDiff2 - B(obj.b,0,dDiff2))/obj.b;            
%                     
%                     cov(2) = cov(2) + vCorrelIRFX(3)*obj.HW2Model_dom.sigma * obj.volTS(iNumVol)*(dDiff2 - B(obj.HW2Model_dom.a,0,dDiff2))/obj.HW2Model_dom.a;
%                     cov(2) = cov(2) + vCorrelIRFX(4)*obj.HW2Model_dom.eta *obj.volTS(iNumVol)*(dDiff2 -  B(obj.HW2Model_dom.b,0,dDiff2))/obj.HW2Model_dom.b;            
% 
%                 end
%                 %%the variance    
%                 covIRIR = CovarianceIRIR(obj,t, T, obj.HW2Model_dom);                
%                 varIRDom = CovarianceIRIR(obj.HW2Model_dom,t, T,obj.HW2Model_dom);
%                 varIRFor = CovarianceIRIR(obj,t, T, obj);                            
% 
%                 var = var + varIRDom + varIRFor - 2*covIRIR + 2*cov(2)- 2*cov(1);
%             end
%         end 