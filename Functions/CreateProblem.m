function prob = CreateProblem(grid, options)
% Creates a structure containing an ODE or PDE.
%
% Inputs:
% - grid: grid structure created using 'CreateGrid'
% - options: structure containing;
%       - L: linear part of equation, array (default: zeros(Nx, Ny))
%       - N: nonlinear part of equation, function (default: @(u) 0*u)
%       - u: value of function at current time (default: zeros(Nx, Ny)
%       - filter: mask used during timestepping (default: 1)
%       - dt: timestep (default: 0.1)
%       - t: current time (default: 0)

    arguments
       grid struct
       options.L {DoubleOrFunction} = zeros(grid.Nx, grid.Ny)
       options.N function_handle    = @(u) 0*u
       options.u (:,:,:)            = zeros(grid.Nx, grid.Ny)
       options.filter (:,:)         = 1
       options.dt (1,1) double      = 0.1
       options.t (1,1) double       = 0
    end

    % get (Nx, Ny) from grid:

    Nx = grid.Nx;
    Ny = grid.Ny;

    % get linear and nonlinear parts of equation:

    L = options.L;
    N = options.N;

    % raise error if size of L is inconsistent with grid:

    if ~isa(L, "function_handle") && ~isequal(size(L(:,:,1,1)), [Nx, Ny])
        if ~isequal(size(L(:,:,1,1)), [1, 1])
            error(['L must be of size ' num2str(Nx) ' x ' num2str(Ny)])
        end
    end

    % calculate L if entered as a function in Fourier space:

    if isa(L, "function_handle")
        L = L(grid.k, grid.l, grid.K2inv);
    end

    % define problem fields:

    prob.L = L;
    prob.grid = grid;
    prob.N = N;
    prob.u = options.u;
    prob.filter = options.filter;
    prob.t = options.t;
    prob.dt = options.dt;

    % move problem arrays to GPU if grid.GPU = true:

    if grid.GPU
        prob.L = gpuArray(prob.L);
        prob.u = gpuArray(prob.u);
        prob.filter = gpuArray(prob.filter);
    end

end

% validation function for L input:

function DoubleOrFunction(a)

    assert(isa(a, "function_handle") || isa(a, "double") || isa(a, "gpuArray"), ...
        "Input must be a double array or function handle.")

end