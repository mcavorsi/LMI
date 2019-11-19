%Minimize the maximum eigenvalue of a matrix that depends affinely on a
%variable

clear all;

A = [1 4 2; 0 1 1; 3 8 5];
B = [3; 5.5; 9];
C = [0 1 1];
eps = 0.0001;

P = sdpvar(size(A,1),size(A,2)); gamma = sdpvar(1,1);

Constraints = P >= eps*eye(size(P));
Mat = [-A'*P-P*A-C'*C P*B; B'*P gamma];
Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

sol = optimize(Constraints,gamma);