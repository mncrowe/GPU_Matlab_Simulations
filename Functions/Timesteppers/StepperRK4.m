function [u, t] = StepperRK4(u, RHS, t, dt, filter, Nt)
% Fourth-order Runge-Kutta solver
%
% Inputs:
% - u: current solution
% - RHS: right-hand side function
% - t: current time
% - dt: timestep
% - filter: array, solution is multiplied by filter each timestep
% - Nt: number of timesteps

    arguments
        u (:,:,:)
        RHS function_handle
        t (1,1) double
        dt (1,1) double
        filter (:,:)
        Nt (1,1) double
    end

    for i = 1:Nt

        k1 = RHS(u, t);
        k2 = RHS(u + k1 * dt/2, t + dt/2);
        k3 = RHS(u + k2 * dt/2, t + dt/2);
        k4 = RHS(u + k3 * dt, t + dt);

        u = filter .* (u + (k1 + 2*k2 + 2*k3 + k4) * dt/6);
        t = t + dt;
 
    end

end