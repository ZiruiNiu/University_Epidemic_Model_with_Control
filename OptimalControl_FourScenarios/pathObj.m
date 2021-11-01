function J = pathObj(x, u, x_tar, Q, R)
% Define the Lagrangian term in the cost function
% It's a quadratic form in our design
    J = sum(Q*((x - x_tar).^2)) + sum(R*(u.^2));
end