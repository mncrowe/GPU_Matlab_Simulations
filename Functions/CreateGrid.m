function grid = CreateGrid(Nx, Ny, Lx, Ly, options)
% Creates a periodic 2D grid in both real space and Fourier space.
%
% Inputs:
% (Nx, Ny): number of gridpoints in x and y directions (default: 128)
% (Lx, Ly): domain limits in x and y, scalars or vectors (default: [-1 1])
% - options: structure containing;
%       - GPU: creates grid on a GPU if set to true

    arguments
        Nx (1,1) double = 128
        Ny (1,1) double = 128
        Lx (1,:) double = [-1 1]
        Ly (1,:) double = [-1 1]
        options.GPU (1,1) logical = false
    end

    % set Lx and Ly to vectors for scalar inputs:

    if length(Lx) == 1
        Lx = [0 Lx];
    end

    if length(Ly) == 1
        Ly = [0 Ly];
    end

    % calculate grid spacing:

    dx = (Lx(2) - Lx(1)) / Nx;
    dy = (Ly(2) - Ly(1)) / Ny;

    % define grid and create x, y, k, l:

    grid.GPU = options.GPU;

    grid.Nx = Nx;
    grid.Ny = Ny;

    grid.Lx = Lx;
    grid.Ly = Ly;
    
    grid.x = (Lx(1):dx:(Lx(2)-dx))';
    grid.y = Ly(1):dy:(Ly(2)-dy);

    grid.k = 2 * pi / (Lx(2) - Lx(1)) * [0:Nx/2-1 -Nx/2:-1]';
    grid.l = 2 * pi / (Ly(2) - Ly(1)) * [0:Ny/2-1 -Ny/2:-1];

    % calculate K^2 = |(k, l)|^2 and 1/K^2:

    grid.K2 = grid.k .^2 + grid.l.^2;

    grid.K2inv = 1 ./ grid.K2;
    grid.K2inv(1, 1) = 0;

    % move grid arrays to a GPU if GPU = true:

    if options.GPU
        grid.K2 = gpuArray(grid.K2);
        grid.K2inv = gpuArray(grid.K2inv);
        grid.k = gpuArray(grid.k);
        grid.l = gpuArray(grid.l);
    end

end