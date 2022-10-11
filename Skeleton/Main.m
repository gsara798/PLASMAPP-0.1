clear
% Parameters
L = 64; dt = 0.1; Nt = 1000; Ng = 256; Num_beams = 2;
% Beam 1
N1 = 10000; V01 = 5; Vth1 = 1; QM1 = -1; XP1 = 0; Mode1 = 0; WP = 1;
% Beam 2
N2 = 10000; V02 = -5; Vth2 = 1; QM2 = -1; XP2 = 0; Mode2 = 0;
% Ions in the background
IB = 20000;

% Methods
Motion_method = 'Leapfrog';  % 'Leapfrog' 'Runge Kutta (RK4)' 'Euler method'
Field_method = 'Finite Difference Method';  % 'Finite Difference Method' 'Fast Fourier Transform (FFT)' 'Direct Integration'
Interpolation_method = 'Cloud in Cell (CIC)'; % 'Nearest Grid Point (NGP)'

% Size of the cell
dx = L/(Ng);
time = 0:dt:(Nt-1)*dt;

% Charge
[Q1,Q2,rho_back] = Charge(QM1,QM2,IB,N1,N2,L,WP);

% Initial loading
[xp1,vp1] = InitialLoading(N1,V01,Vth1,XP1,Mode1,L);
[xp2,vp2] = InitialLoading(N2,V02,Vth2,XP2,Mode2,L);


% Auxiliarity vectors
p1 = AuxVector(N1); p2 = AuxVector(N2);

% Poisson's equation preparation
un = ones(Ng-1, 1); % Ng-1 * 1
Ax = spdiags([un -2*un un], [-1 0 1], Ng-1, Ng-1); % Matrix that gives the indices of the Poisson equation Ng-1 * Ng-1.
kap = (2*pi/L)*(-Ng/2:Ng/2-1);
kap = fftshift(kap');
kap(1) = 1;

% Others
mat1 = 0; mat2 = 0; Eg = 0; Phi =0;

% Computational cycle
for it = 1:Nt

    mom(it) = (Q1/QM1)*(sum(vp1))+(Q2/QM2)*(sum(vp2));
    E_kin(it) = 0.5*abs(Q1/QM1)*(sum(vp1.^2)) + 0.5*abs(Q2/QM2)*(sum(vp2.^2));
    E_pot(it) = 0.5*sum(Eg.^2)*dx;

    switch Motion_method
        case 'Leapfrog'
            vp1 = MotionV(vp1,QM1,mat1,Eg,N1,Motion_method,dt);
            vp2 = MotionV(vp2,QM2,mat2,Eg,N2,Motion_method,dt);
    end

    %Updating positions
    xp1 = MotionX(xp1,vp1,Motion_method,dt);
    xp2 = MotionX(xp2,vp2,Motion_method,dt);

    %Periodic boundary conditions:
    xp1 = PeriodicBC(xp1,0,L);
    xp2 = PeriodicBC(xp2,0,L);

    %Interpolation functions
    mat1 = InterpolationF(Interpolation_method,dx,Ng,xp1,N1,p1);
    mat2 = InterpolationF(Interpolation_method,dx,Ng,xp2,N2,p2);

    %Charge density:
    rho1 = Charge_density(Q1,mat1,dx);
    rho2 = Charge_density(Q2,mat2,dx);
    rhot = rho1 + rho2 + rho_back;

    % Field equations:
    [Phi,Eg] = Field(Field_method,rhot,Ng,dx,L,Ax,kap);

    %Updating velocity
    vp1 = MotionV(vp1,QM1,mat1,Eg,N1,Motion_method,dt);
    vp2 = MotionV(vp2,QM2,mat2,Eg,N2,Motion_method,dt);

    scatter(xp1,vp1,'.')
    hold on
    scatter(xp2,vp2,'.')
    pause(0.00001)
    hold off
end

figure;
plot(time,mom)
xlabel('Time')
ylabel('Momentum')

figure;
hold on
plot(time,E_kin,'k')
plot(time,E_pot,'r')
plot(time,E_kin+E_pot,'b')
xlabel('Time')
ylabel('Energy')

