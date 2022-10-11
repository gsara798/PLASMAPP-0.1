%% Function for calculate position of particles
function position = MotionX(position,velocity,Motion_method,dt)
    switch Motion_method
        case 'Euler method'
            position = position + velocity*dt;
        case 'Leapfrog'
            position = position + velocity*dt;
        case 'Runge Kutta (RK4)'
            k1 = velocity*dt;
            k2 = (velocity+(1/2)*k1)*dt;
            k3 = (velocity+(1/2)*k2)*dt;
            k4 = (velocity+k3)*dt;
            position = position + (1/6)*(k1+2*k2+2*k3+k4);
    end
end
