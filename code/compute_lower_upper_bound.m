function [qLower, qUpper] = compute_lower_upper_bound(lowerTol, upperTol, t, qDesired)

% Determine lower and upper curves
if abs(lowerTol) == abs(upperTol)
    % Compute parallel curves to desired trajectory
    % [tInner, qInner, tOuter, qOuter, ~, ~, ~, ~] = ...
    %     parallel_curve(t, qDesired, abs(lowerTol), 0, 0);
    tInner = t;
    tOuter = t;
    qInner = qDesired - abs(lowerTol);
    qOuter = qDesired + abs(lowerTol);
    % Determine lower and upper curves
    if qOuter(1) >= qDesired(1)
        qUpper = qOuter;
        tUpper = tOuter;
        qLower = qInner;
        tLower = tInner;
    else
        qUpper = qInner;
        tUpper = tInner;
        qLower = qOuter;
        tLower = tOuter;
    end
else
    % Compute parallel curves to desired trajectory
    % [tInner1, qInner1, tOuter1, qOuter1, ~, ~, ~, ~] = ...
    %     parallel_curve(t, qDesired, abs(lowerTol), 0, 0);
    % [tInner2, qInner2, tOuter2, qOuter2, ~, ~, ~, ~] = ...
    %    parallel_curve(t, qDesired, abs(upperTol), 0, 0);
    tInner1 = t;
    tOuter1 = t;
    qInner1 = qDesired - abs(lowerTol);
    qOuter1 = qDesired + abs(lowerTol);
    tInner2 = t;
    tOuter2 = t;
    qInner2 = qDesired - abs(upperTol);
    qOuter2 = qDesired + abs(upperTol);
    % Determine lower curve
    if qOuter1(1) >= qDesired(1)
        qLower = qInner1;
        tLower = tInner1;
    else
        qLower = qOuter1;
        tLower = tOuter1;
    end
    % Determine upper curve
    if qOuter2(1) >= qDesired(1)
        qUpper = qOuter2;
        tUpper = tOuter2;
    else
        qUpper = qInner2;
        tUpper = tInner2;
    end
end

% Resample curves to be consistent with simulation time vector
% tLower, tUpper must be a monotonically increasing column vector
% qLower, qUpper must be a column vector
% t must be a column vector
% qLower = interp1q(tLower,qLower,t);
% qUpper = interp1q(tUpper,qUpper,t);

end
