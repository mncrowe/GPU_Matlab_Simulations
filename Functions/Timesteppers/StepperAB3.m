function [u, t, RHS1, RHS2] = StepperAB3(u, RHS, t, dt, filter, Nt, RHS1, RHS2)
% Three-step Adams-Bashforth solver
%
% Inputs:
% - u: current solution
% - RHS: right-hand side function
% - t: current time
% - dt: timestep
% - filter: array, solution is multiplied by filter each timestep
% - Nt: number of timesteps
% - (RSH1, RHS2): previous values of RHS(u, t), Euler timestep used if not available

    arguments
        u (:,:,:)
        RHS function_handle
        t (1,1) double
        dt (1,1) double
        filter (:,:)
        Nt (1,1) double
        RHS1 (:,:,:)      = NaN
        RHS2 (:,:,:)      = NaN
    end

    for i = 1:Nt

        % use the Euler method for first two timesteps, then use multistep:

        if isnan(RHS1)

            RHS1 = RHS(u, t);
            u = filter .* (u + RHS1 * dt);

        else

            if isnan(RHS2)

                RHS2 = RHS(u, t);
                u = filter .* (u + RHS2 * dt);

            else

                RHS3 = RHS2;
                RHS2 = RHS1;
                RHS1 = RHS(u, t);
    
                u = filter .* (u + (23 * RHS1  - 16 * RHS2 + 5 * RHS3) * dt/12);

            end

        end

        t = t + dt;

    end

end