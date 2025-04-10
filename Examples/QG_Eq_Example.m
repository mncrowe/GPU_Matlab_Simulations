% Solve the QG equation in Fourier space:
%
% dQ/dt + J[psi, Q] = nu*Lap[Q]
%
% for Q = Lap[psi].

% Define grid and set parameters:

nu = 0;
GPU = false; %true;

Nx = 512;
Ny = 512;

t_end = 1;
dt = 0.001;
Nt = t_end/dt;
Ns = 10;

Lx = [-pi, pi];
Ly = [-pi, pi];

grid = CreateGrid(Nx, Ny, Lx, Ly, GPU = GPU);

% Define wavenumber parameters for RHS function and create filter:

k = grid.k;
l = grid.l;
K2inv = grid.K2inv;
K2 = grid.K2;
K = sqrt(grid.K2);

filter = CreateFilter(grid);

% Define linear operator for L in Fourier space:

L = - nu*grid.K2;

% Set initial condition:

psih = (K2 .* (1 + K.^4)).^(-1/2) .* (randn(Nx, Ny) + 1i*randn(Nx, Ny));
psih(1, 1) = 0;

psih = fft2(real(ifft2(psih)));
E = sum(0.5 * K2(:) .* abs(psih(:)).^2 / (Nx * Ny)^2);
psih = psih * sqrt(0.5/E);

qh = -K2 .* psih .* filter;

% Define problem:

prob = CreateProblem(grid, L = L, N = @(u) RHS(u, K2inv, k, l), ...
    u = qh, filter = filter, dt = dt);

% Timestep problem in Fourier space:

[prob, Q, T, runtime] = StepForward(prob, Nt, Ns = Ns, timestepper = "AB3");

disp(['Time per timestep: ', num2str(runtime/Nt*1000, '%1.3f'), ' ms']);

% Move solution to real space:

Q = real(ifft2(Q));

% Create nonlinear RHS function:

function qh = RHS(qh, K2inv, k, l)

    psih = -K2inv .* qh;
    uh = - 1i * l .* psih;
    vh =   1i * k .* psih;
    
    u = real(ifft2(uh));
    v = real(ifft2(vh));
    q = real(ifft2(qh));
    
    qh = -1i * k .* (fft2(u .* q)) - 1i * l .* (fft2(v .* q));

end