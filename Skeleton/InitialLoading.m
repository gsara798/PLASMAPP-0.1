function [position,velocity] = InitialLoading(N_particles,V0,Vth,Amplitude,Mode,L)
    position = transpose(linspace(0,L-L/N_particles,N_particles));
    % Maxwellian distribution of velocities of the Beam.
    velocity = Vth*randn(N_particles,1);
    % Adding the mean velocity.
    velocity = velocity+V0;
    % Perturbation
    if Amplitude == 0
    else
        position = position+Amplitude*cos(2*pi*Mode*position/L);
    end
end
