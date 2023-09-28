% Bisect Function

lamb1 = 10^-10;
lamb2 = 10^10;

while (lamb2-lamb1)/(lamb1+lamb2) > eps
    lambmid = (lamb1+lamb2)/2;
    for e=1:ne
        rho_newe = rho(e)*(-fp(e)/(lambmid*gp(e)))^eta;
        if rho_newe <= rho_min
            rho_newe = rho_min;
        elseif rho_newe >= 1
            rho_newe = 1;
        else
            rho_newe = rho_newe;
        end
        rho_new(e) = rho_newe;
    end
    g = rho_new'*v-V_con;
    if g > 0
        lamb1 = lambmid;
    else 
        lamb2 = lambmid;
    end
end
