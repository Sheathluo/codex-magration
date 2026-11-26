function [A, dA, ddA] = joint_migration_template_reference(X, theta, gamma, grid_step, cfg, delay_sign)
% Three-parameter power-template reference; independent of the production pipeline.
% theta = [tau0; nu0; eta], eta = normalized Doppler increment per symbol.
% tau(k) = tau0 + delay_sign*gamma*(k*nu0 + 0.5*k^2*eta).
% nu(k)  = nu0 + k*eta, k = 0,...,K-1.
% delay_sign = -1: approaching-positive Doppler, received phase exp(+j*phi).
% delay_sign = +1: reproduces the existing positive delay-drift algebra;
%                 carrier phase / Doppler convention must be made consistent.
% X: real baseline POWER surface on the complete uniform template grid.
% grid_step = [h_tau, h_nu] in normalized delay and normalized Doppler units.
% cfg.symbol_power_weight: K fixed power weights (production convention: mean 1).
% cfg.migration_compensate_bin: 2-by-K fixed forward drifts in SEARCH-grid bins.
% cfg.echo_rows, cfg.echo_cols: extraction indices into the full template.
% Processing order: full-grid Fourier shift -> extract -> search-grid circshift.
% A: search-grid template; dA(:,:,i); ddA(:,:,i,j), i,j = 1:3.
% Fourier interpolation is periodic; boundary / Nyquist behavior matches the
% explicit real(ifft(...)) model and needs appropriate padding in applications.

    validateattributes(X, {'numeric'}, {'real', '2d', 'nonempty', 'finite'});
    validateattributes(theta, {'numeric'}, {'real', 'vector', 'numel', 3, 'finite'});
    validateattributes(gamma, {'numeric'}, {'real', 'scalar', 'finite', 'nonnegative'});
    validateattributes(grid_step, {'numeric'}, {'real', 'vector', 'numel', 2, 'positive', 'finite'});
    assert(isscalar(delay_sign) && any(delay_sign == [-1, 1]), 'delay_sign must be -1 or +1.');
    theta = theta(:);
    weights = cfg.symbol_power_weight(:);
    K = numel(weights);
    validateattributes(weights, {'numeric'}, {'real', 'nonempty', 'finite', 'nonnegative'});
    bins = cfg.migration_compensate_bin;
    assert(isequal(size(bins), [2, K]) && all(isfinite(bins(:))), 'Compensation must be finite 2-by-K.');
    assert(all(bins(:) == round(bins(:))), 'Fixed compensation must contain integer search-grid bins.');
    shifts = -bins.';
    [groups, ~, membership] = unique(shifts, 'rows', 'stable');

    [Mt, Mv] = size(X);
    rt = ifftshift(-floor(Mt/2):ceil(Mt/2)-1).';
    rv = ifftshift(-floor(Mv/2):ceil(Mv/2)-1);
    U = -1i*2*pi*rt/(Mt*grid_step(1));
    V = -1i*2*pi*rv/(Mv*grid_step(2));
    XF = fft(fft(X, [], 1), [], 2);
    rows = cfg.echo_rows;
    cols = cfg.echo_cols;
    response = zeros(numel(rows), numel(cols), 10);
    pairs = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];

    for g = 1:size(groups, 1)
        sums = complex(zeros(Mt, Mv, 10));
        members = find(membership == g);
        for kk = 1:numel(members)
            index = members(kk);
            k = index - 1;
            r = delay_sign*gamma*k;
            b = 0.5*delay_sign*gamma*k^2;
            tau_k = theta(1) + r*theta(2) + b*theta(3);
            nu_k = theta(2) + k*theta(3);
            phase = weights(index)*exp(U*tau_k + V*nu_k);
            B = cat(3, U + zeros(1, Mv), r*U + V, b*U + k*V);
            sums(:, :, 1) = sums(:, :, 1) + phase;
            for i = 1:3
                sums(:, :, 1+i) = sums(:, :, 1+i) + B(:, :, i).*phase;
            end
            for ij = 1:6
                i = pairs(ij, 1);
                j = pairs(ij, 2);
                sums(:, :, 4+ij) = sums(:, :, 4+ij) + B(:, :, i).*B(:, :, j).*phase;
            end
        end
        full_response = real(ifft(ifft(XF.*sums, [], 1), [], 2));
        search_response = full_response(rows, cols, :);
        response = response + circshift(search_response, [groups(g, :), 0]);
    end
    response = response/K;
    A = response(:, :, 1);
    dA = response(:, :, 2:4);
    ddA = zeros(numel(rows), numel(cols), 3, 3);
    for ij = 1:6
        i = pairs(ij, 1);
        j = pairs(ij, 2);
        ddA(:, :, i, j) = response(:, :, 4+ij);
        ddA(:, :, j, i) = response(:, :, 4+ij);
    end
end
