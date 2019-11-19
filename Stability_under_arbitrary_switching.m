%The systems A1*x(t) and A2*x(t) are stable under arbitrary switching if
%the following LMI is feasible.

clear all;

A1 = [4 3 6; 8.2 7 1.3; 0 9 2.4];
A2 = [2 9 0; 3.1 0 8; 4 1.6 7.5];

eps = 0.0001;

P = sdpvar(size(A1,1),size(A1,1));

Constraints = P >= eps*eye(size(P));
Constraints = [Constraints, A1'*P+P*A1 <= eps*eye(size(P))];
Constraints = [Constraints, A2'*P+P*A2 <= eps*eye(size(P))];

sol = optimize(Constraints);

if (sol.problem ~= 0)
    'Switched system not stable under arbitrary switching'
else
    'Switched system is stable under arbitrary switching'
end