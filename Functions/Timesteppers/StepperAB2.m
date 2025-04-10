function [u, t, RHS1] = StepperAB2(u, RHS, t, dt, filter, Nt, RHS1)
% Two-step Adams-Bashforth solver
%
% Inputs:
% - u: current solution
% - RHS: right-hand side function
% - t: current time
% - dt: timestep
% - filter: array, solution is multiplied by filter each timestep
% - Nt: number of timesteps
% - RSH1: previous value of RHS(u, t), Euler timestep used if not available

    arguments
        u (:,:,:)
        RHS function_handle
        t (1,1) double
        dt (1,1) double
        filter (:,:)
        Nt (1,1) double
        RHS1 (:,:,:)        = NaN
    end

    for i = 1:Nt

        % use the Euler method for initial timestep, then use multistep:

        if isnan(RHS1)

            RHS1 = RHS(u, t);
            u = filter .* (u + RHS1 * dt);

        else

            RHS2 = RHS1;
            RHS1 = RHS(u, t);

            u = filter .* (u + (3/2 * RHS1  - 1/2 * RHS2) * dt);

        end

        t = t + dt;

    end

end