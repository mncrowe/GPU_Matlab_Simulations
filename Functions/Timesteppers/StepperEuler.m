function [u, t] = StepperEuler(u, RHS, t, dt, filter, Nt)
% First-order Euler solver
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

        u = filter .* (u + RHS(u, t) * dt);
        t = t + dt;

    end

end