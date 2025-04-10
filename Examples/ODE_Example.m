% Example for the ODE: y' = -y, y(0) = 1.
%
% This script solves the ODE using various timesteppers then compares the
% accuracy of the resulting numerical solutions.

N = @(u) -u;

dt = 0.01;
t_end = 1;

Nt = t_end/dt;

grid = CreateGrid(1, 1);
prob = CreateProblem(grid, N = N, u = 1, dt = dt);

[~, U_Euler, ~] = StepForward(prob, Nt, timestepper="Euler", Ns = 10);
[~, U_RK4, ~] = StepForward(prob, Nt, timestepper="RK4", Ns = 10);
[~, U_RK2, ~] = StepForward(prob, Nt, timestepper="RK2", Ns = 10);
[~, U_AB2, ~] = StepForward(prob, Nt, timestepper="AB2", Ns = 10);
[~, U_AB3, T] = StepForward(prob, Nt, timestepper="AB3", Ns = 10);

U = exp(-T)';

semilogy(T, abs(U - U_Euler), T, abs(U - U_RK4), T, ...
    abs(U - U_RK2), T, abs(U - U_AB2), T, abs(U - U_AB3))

legend('Euler', 'RK4', 'RK2', 'AB2', 'AB3')