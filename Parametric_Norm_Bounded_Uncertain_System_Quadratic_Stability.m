%The system xdot = Ax + Mp, p = delta*q, q = Nx + Qp, is quadratically
%stable under parametric norm-bounded uncertainty if the following LMI in
%A,M,N,Q is feasible.

clear all;

A = [1 4 2; 0 1 1; 3 8 5];
M = [6; 0; 2];
N = [0 0 4; 9 9 1; 3 5 4];
Q = [0; 3; 7];

eps = 0.0001;

P = sdpvar(size(A,1),size(A,2)); 
mu = sdpvar(1,1);

Constraints = P >= eps*eye(size(P));
Constraints = [Constraints, mu >= 0];
Mat = [A*P+P*A' P*N'; N*P zeros(size(P*N'))] + mu*[M*M' M*Q'; Q*M' Q*Q'-eye(size(Q*Q'))];
Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

sol = optimize(Constraints);

if (sol.problem == 1)
    'Problem Infeasible: Not Quadratically Stable'
else
    'Problem Feasible! Quadratically Stable!'
end