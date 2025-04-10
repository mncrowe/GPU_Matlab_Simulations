function [prob, U, T, runtime] = StepForward(prob, Nt, options)
% Evolves the ODE or PDE given by 'prob' forward in time.
%
% Inputs:
% - prob: problem structure created using 'CreateProblem'
% - options: structure containing;
%       - Nt: number of timesteps (default: 1)
%       - timestepper, string describing timestepper (default: "RK4")
%       - Ns: number of saves, excluding initial save (default: 1)
%       - GPU_out: if true, output U is moved to a GPU (default: false)
%
% Outputs:
% - prob: updated problem
% - U: array containing saved solution values
% - T: array containing times corresponding to U saves
% - runtime: time taken to run timestepping (ms)

    arguments
        prob struct
        Nt (1,1) double            = 1
        options.timestepper string = "RK4"
        options.Ns                 = 1
        options.GPU_out            = false
    end

    % set Ns = Nt if Ns > Nt:

    if options.Ns > Nt
        Ns = Nt;
    else
        Ns = options.Ns;
    end

    % get number of variables in u (i.e. length of u = (u1, u2, ...)):

    Su = size(prob.u);

    if length(Su) > 2
        Nu = Su(3);
    else
        Nu = 1;
    end

    % get problem size from grid:

    Nx = prob.grid.Nx;
    Ny = prob.grid.Ny;

    % get problem fields:

    u = prob.u;
    t = prob.t;
    dt = prob.dt;
    filter = prob.filter;

    % create save arrays:

    U = zeros(Nx, Ny, Nu, Ns + 1);
    T = zeros(1, Ns + 1);

    if options.GPU_out
        U = gpuArray(U);
    end

    % set initial value of save arrays to current solution:

    U(:, :, :, 1) = u;
    T(1) = t;

    % define iteration intervals between each save:

    Nb = ceil(Nt/Ns);
    I = [1:Nb:Nt Nt+1];

    % create right-hand side function:

    if Nu > 1
        RHS_shape = @(F) squeeze(permute(reshape(F, [Nx Nu Ny]), [1 3 2]));
        RHS = @(u, t) RHS_shape(prob.N(u)) + squeeze(sum(prob.L .* u, 3));
    else
        RHS = @(u, t) prob.N(u) + prob.L .* u;
    end

    % start timer:   

    tic

    % perform timestepping for Nt iterations:

    for i = 1:Ns

        if isequal(options.timestepper, "RK2")
            [u, t] = StepperRK2(u, RHS, t, dt, filter, I(i+1) - I(i));
            U(:, :, :, i+1) = u;
            T(i+1) = t;
        end

        if isequal(options.timestepper, "RK4")
            [u, t] = StepperRK4(u, RHS, t, dt, filter, I(i+1) - I(i));
            U(:, :, :, i+1) = u;
            T(i+1) = t;
        end

        if isequal(options.timestepper, "Euler")
            [u, t] = StepperEuler(u, RHS, t, dt, filter, I(i+1) - I(i));
            U(:, :, :, i+1) = u;
            T(i+1) = t;
        end

        if isequal(options.timestepper, "AB2")
            if ~exist('RHS1', 'var'); RHS1 = NaN; end
            [u, t, RHS1] = StepperAB2(u, RHS, t, dt, filter, I(i+1) - I(i), RHS1);
            U(:, :, :, i+1) = u;
            T(i+1) = t;
        end

        if isequal(options.timestepper, "AB3")
            if ~exist('RHS1', 'var'); RHS1 = NaN; end
            if ~exist('RHS2', 'var'); RHS2 = NaN; end
            [u, t, RHS1, RHS2] = StepperAB3(u, RHS, t, dt, filter, I(i+1) - I(i), RHS1, RHS2);
            U(:, :, :, i+1) = u;
            T(i+1) = t;
        end

        % print progress to screen:

        disp(['Iteration: ' num2str(I(i + 1) - 1)])

        % identify NaN in solution and raise error:

        if max(isnan(u),[],"all")
            error("NaN detected in u.")
        end

    end

    % finish timing:

    runtime = toc;

    % remove dimensions of length 1 from U:

    U = squeeze(U);

    % update solution in prob:

    prob.u = u;
    prob.t = t;
    
end