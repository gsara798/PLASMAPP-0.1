% Function for field equations
function [Phi,Eg] = Field(Field_method,charge_density,Ng,dx,L,Ax,kap)
    switch Field_method
        case 'Finite Difference Method'
            Phi = Ax\(-charge_density(1:Ng-1)*dx^2);
            Phi = [Phi;0];
            Eg = ([Phi(Ng); Phi(1:Ng-1)]-[Phi(2:Ng);Phi(1)])/(2*dx);
        case 'Fast Fourier Transform (FFT)'
            rho_k = fft(charge_density(1:Ng));
            phi_k = rho_k./kap.^2;
            Phi = real(ifft(phi_k,'symmetric'));
            Phi(Ng) = 0;
            Eg_k = -1i.*kap.*phi_k;
            Eg =  real(ifft(Eg_k,'symmetric'));
        case 'Direct Integration'
            G0 = (charge_density(1)+charge_density(Ng))*dx/2;
            G = G0;
            G1(1) = G0;
            for k = 2:(Ng)
                G = G + (charge_density(k)+charge_density(k-1))*dx/2;
                G1(k) = G;
            end

            F0 = (G1(1)+ G1(end))*dx/2;
            F = F0;
            F1(1) = F0;
            for k = 2:(Ng)
                F = F + (G1(k)+G1(k-1))*dx/2;
                F1(k) = F;
            end

            Phi_prime = (1/L)*F1(end);
            Eg = -Phi_prime +  G1;
            Eg = Eg';
            
            Phi0 = (Eg(1)+ Eg(end))*dx/2;
            Phi1 = Phi0;
            Phi(1) = Phi0;
            for k = 2:(Ng)
                Phi1 = Phi1 + (Eg(k)+Eg(k-1))*dx/2;
                Phi(k) = -Phi1;
            end
            Phi = Phi';
    end
end
