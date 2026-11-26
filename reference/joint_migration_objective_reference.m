function [Q, gradQ, hessQ, loss] = joint_migration_objective_reference(Y, power, A, dA, ddA)
% Fixed-power objective matching the existing maximizing Newton convention.
% All inputs must already have identical extraction, compensation, local crop.
% Q = 2*power*<Y,A> - power^2*<A,A> = ||Y||_F^2 - loss.
% This Hessian holds power fixed; profiling power requires its Schur correction.
    validateattributes(power, {'numeric'}, {'real', 'scalar', 'finite', 'nonnegative'});
    assert(isreal(Y) && isreal(A) && isequal(size(Y), size(A)), 'Y and A must be equally sized real surfaces.');
    e = Y - power*A;
    loss = sum(e.^2, 'all');
    Q = 2*power*sum(Y.*A, 'all') - power^2*sum(A.^2, 'all');
    gradQ = zeros(3, 1);
    hessQ = zeros(3, 3);
    for i = 1:3
        Ai = dA(:, :, i);
        gradQ(i) = 2*power*sum(e.*Ai, 'all');
        for j = i:3
            Aj = dA(:, :, j);
            Aij = ddA(:, :, i, j);
            hessQ(i, j) = 2*power*sum(e.*Aij, 'all') - 2*power^2*sum(Ai.*Aj, 'all');
            hessQ(j, i) = hessQ(i, j);
        end
    end
end
