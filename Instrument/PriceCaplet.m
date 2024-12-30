function CapletPrice = PriceCaplet(vInCurve,X,vYFFixing,vYFMaturity,a,b,sigma,eta,rho)

    if ~isafin(vInCurve,'RateSpec') && ~isa(vInCurve,'IRDataCurve')
        error(message('fininst:swaptionbylg2f:invalidCurve'));
    end

    if ~isscalar(a),error(message('fininst:swaptionbylg2f:invalidA')),end 
            
    nCaplets = length(vYFMaturity);
    CapletPrice = zeros(nCaplets,1);
    
    for cplidx = 1:nCaplets                       
        curve = vInCurve(cplidx);        
        PM = @(t)discountFactorByTime(curve, t);
        vTZC = [vYFFixing(cplidx);vYFMaturity(cplidx)];
        P0 = PM(vTZC);       
        dDiff = diff(vTZC);
        dK = 1/(1 + dDiff*X(cplidx));
        
        %dV = VarZC(a, vSigma,vYFVolDates,vTZC);
        vV = yyVarianceX(a,b,sigma,eta,rho,vTZC(1), vTZC(2));
        dV = vV(1)*((1-exp(-a*diff(vTZC)))/a)^2 + vV(2)*((1-exp(-b*diff(vTZC)))/b)^2+ ...
                                            vV(3)*(1-exp(-a*diff(vTZC)))*(1-exp(-b*diff(vTZC)))/a/b;
        dD1 = (log(P0(2)/P0(1)/dK)+0.5*dV)/sqrt(dV);
        dD2 = dD1 - sqrt(dV);
        
        CapletPrice(cplidx) = P0(1)*normcdf(dD1)-P0(2)*dK*normcdf(dD2);
    end
end

function vV = yyVarianceX(a,b,sigma,eta,rho,t, T)
    diff = T - t;
    vV(1) = sigma*sigma*(1 - exp(-2*a*diff))/(2*a);
    vV(2) = eta*eta*(1 - exp(-2*b*diff))/(2*b);
    vV(3) = 2 * rho * sigma * eta*(1 - exp(-(a+b)*diff))/(a+b);
end
