%There exists a theta such that ||theta*M*theta^(-1)]]^2 < gamma
%if the following LMI is feasible

clear all;

A = [1 4 2; 0 1 1; 3 8 5];
B = [3; 5.5; 9];
C = [0 1 1];
D = 0;
eps = 0.0001;

X = sdpvar(size(A,1),size(A,2)); gamma = sdpvar(1,1);
theta = sdpvar(1,1);

Constraints = X >= eps*eye(size(X));
Mat = [A'*X+X*A X*B; B'*X -theta] + (1/gamma^2)*[C'; D']*theta*[C D];
Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

sol = bisection(Constraints,gamma);