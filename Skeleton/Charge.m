%% Function that calculates the charge of superparticles and background density
function [Q_1,Q_2,rho_back] = Charge(QM1,QM2,IB,N1,N2,L,WP)
    Q_1 = (WP^2*L)/(N1*QM1);
    if QM2>0
        Q_2 = -Q_1*(N1/N2);
    end
    if QM2<0
        Q_2 = Q_1*(N1/N2);
    end
    if QM2 == 0
        Q_2 = 0;
    end
    rho_back = (-Q_1/L)*(IB);   %Background density
end
