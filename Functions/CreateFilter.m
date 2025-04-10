function filter = CreateFilter(grid, options)
% Creates an array containing an exponential filter to damp high wavenumber
% oscillations and perform dealiasing.
%
% Inputs:
% - grid: grid structure created using 'CreateGrid'
% - options: structure containing;
%       - order: order of filter (default: 4)
%       - Klims: filtered region of wavenumber space (default: [2/3 1])
%       - tol: value of filter at K = Klims(2)

    arguments
        grid struct
        options.order (1,1) double = 4
        options.Klims (2,1) double = [2/3 1]
        options.tol (1,1) double = 1e-15
    end

    % create array of nondimensional wavenumbers K:

    dx = (grid.Lx(2) - grid.Lx(1)) / grid.Nx;
    dy = (grid.Ly(2) - grid.Ly(1)) / grid.Ny;

    K = sqrt((grid.k / pi * dx).^2 + (grid.l / pi * dy).^2);

    % set decay scale so filter = tol at K = Klims(2):

    decay = -log(options.tol) / (options.Klims(2) - options.Klims(1)) ^ options.order;

    % calculate filter:

    filter = exp(- decay * (K - options.Klims(1)) .^ options.order);
    filter(K < options.Klims(1)) = 1;

end