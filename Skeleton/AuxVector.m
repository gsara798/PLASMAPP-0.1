function p = AuxVector(N_particles)
    if N_particles ~= 0
        p = 1:N_particles; p = [p p];
    end
    if N_particles == 0
        p = 0;
    end
end
