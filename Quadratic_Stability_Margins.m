%The quadratic stability margin of the system xdot = Ax + B_p*p is defined
%as the largest alpha >= 0 for which the system is quadratically stable.
%Find the largest alpha using the following code

clear all;

A = [4 3 6; 8.2 7 1.3; 0 9 2.4];
Bp = [6; 4.3; 1];
Cq = [1 1 5];

eps = 0.0001;

P = sdpvar(size(A,1),size(A,1));
lambda = sdpvar(1,1);
alpha = sdpvar(1,1);

beta = alpha*alpha;

Constraints = alpha >= eps;
Mat = [A'*P+P*A+beta*lambda*Cq'*Cq P*Bp; Bp'*P -lambda];
Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

sol = optimize(Constraints,-alpha);

Alpha = value(alpha)