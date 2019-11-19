%The system xdot = Ax + Bu + Mp, p = delta*q, q = Nx + Qp + D_12 u, is quadratically
%stable under parametric norm-bounded uncertainty if the following LMI in
%A,B,M,N,Q,D_12 is feasible. The controller u = Kx is formed with K = ZP^-1

clear all;

A = [1 4 2; 0 1 1; 3 8 5];
B = [3; 5.5; 9];
M = [6; 0; 2];
N = [0 0 4; 9 9 1; 3 5 4];
Q = [0; 3; 7];
D12 = [0; 1; 1];

eps = 0.0001;

P = sdpvar(size(A,1),size(A,2)); mu = sdpvar(1,1);
Z = sdpvar(size(B,2),size(B,1));

Constraints = P >= eps*eye(size(P));
Constraints = [Constraints, mu >= 0];
Mat = [A*P+B*Z+P*A'+Z'*B' P*N'+Z'*D12'; N*P+D12*Z zeros(size(P*N'))] + mu*[M*M' M*Q'; Q*M' Q*Q'-eye(size(Q*Q'))];
Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

sol = optimize(Constraints);

if (sol.problem == 1)
    'Problem Infeasible: Not Quadratically Stable'
else
    'Problem Feasible! Quadratically Stable!'
end

Z = value(Z);
P = value(P);
K = Z*inv(P)