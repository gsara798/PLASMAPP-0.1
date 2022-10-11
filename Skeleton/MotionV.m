%% Funtion that moves the particles
function velocity = MotionV(velocity,QM,interp,Efieldg,N_particles,Motion_method,dt)
    if N_particles == 0
        velocity = 0;
    else
        switch Motion_method
            case 'Euler method'
                velocity = velocity + QM*interp*Efieldg*dt;
            case 'Leapfrog'
                    velocity = velocity + 0.5*QM*interp*Efieldg*dt;
            case 'Runge Kutta (RK4)'
                l1 = (QM*interp*Efieldg)*dt;
                l2 = (QM*interp*Efieldg+(1/2)*l1)*dt;
                l3 = (QM*interp*Efieldg+(1/2)*l2)*dt;
                l4 = (QM*interp*Efieldg+l3)*dt;
                velocity = velocity + (1/6)*(l1+2*l2+2*l3+l4);
        end
    end
end
