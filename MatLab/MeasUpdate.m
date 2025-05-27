function [K, x, P] = MeasUpdate(x, z, g, s, G, P, n)

m = length(z);
disp(m);
Inv_W = zeros(m,m);

for i=1:m
    Inv_W(i,i) = s(i)*s(i);    % Inverse weight (measurement covariance)
end

% Kalman gain
K = P*G'*inv(Inv_W+G*P*G');
disp(K);
disp(z-g);
% State update
x = x + K*(z-g);
disp(x)
% Covariance update
P = (eye(n)-K*G)*P;

