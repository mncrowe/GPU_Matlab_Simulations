% Solve the wave equation in Fourier space
%
% This can be written as:
%
% ḏ  [u] _ [    0        1][u]
% dt [v] ‾ [c^2*d^2/dx^2 0][v]
%
% where we've used a first-order reduction by writing v = du/dt.
%
% This can easily be solved in Fourier space using c^2*d^2/dx^2 = -c^2*K^2,
% where K^2 = k^2 + l^2 is the magnitude of the wave-vector in Fourier
% space.

% Define grid and set parameters:

c = 1;
GPU = true;

Nx = 512;
Ny = 512;

dt = 0.01;
Nt = 1000;

Lx = [-5, 5];
Ly = [-5, 5];

grid = CreateGrid(Nx, Ny, Lx, Ly, GPU = GPU);

% Define linear operator for L in Fourier space:

L = zeros(Nx, Ny, 2, 2);

L(:, :, 1, 1) = zeros(Nx, Ny);
L(:, :, 1, 2) = ones(Nx, Ny);
L(:, :, 2, 1) = -c^2*grid.K2;
L(:, :, 2, 2) = zeros(Nx, Ny);

% Set initial condition in real space:

u = zeros(Nx, Ny, 2);
u(:, :, 1) = exp(-4*(grid.x.^2 + grid.y.^2));
u(:, :, 2) = zeros(Nx, Ny);

% Move to Fourier space:

u = fft2(u);

% Create problem:

prob = CreateProblem(grid, L = L, u = u, dt = dt);

% Timestep problem in Fourier space:

[prob, U, T] = StepForward(prob, Nt, Ns = 100, timestepper = "RK4");

% Move solution to real space:

U = real(ifft2(U));