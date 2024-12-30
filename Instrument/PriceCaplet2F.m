function CapletPrice = PriceCaplet2F(vInCurve,HW2F_obj,X,vYFFixing,vYFMaturity)

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
        vV = yyVarianceX(HW2F_obj,vTZC(1), vTZC(2));
        dV = vV(1)*((1-exp(-HW2F_obj.a*diff(vTZC)))/HW2F_obj.a)^2 + vV(2)*((1-exp(-HW2F_obj.b*diff(vTZC)))/HW2F_obj.b)^2+ ...
                                            vV(3)*(1-exp(-HW2F_obj.a*diff(vTZC)))*(1-exp(-HW2F_obj.b*diff(vTZC)))/HW2F_obj.a/HW2F_obj.b;
        dD1 = (log(P0(2)/P0(1)/dK)+0.5*dV)/sqrt(dV);
        dD2 = dD1 - sqrt(dV);
        
        CapletPrice(cplidx) = P0(1)*normcdf(dD1)-P0(2)*dK*normcdf(dD2);
    end
end