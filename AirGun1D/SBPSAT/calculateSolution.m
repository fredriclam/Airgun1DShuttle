% Calculates the solution of discretization for a given set of times.
%    discr      -- A discretisation object
%    timeUnit   -- A time interval used to specify which times are returned
%    multiples -- vector with integer multiples of timeUnit to return the solution for
function [q, bubble, T] = calculateSolution(discr, timeUnit, multiples)
    k_max = discr.getTimestep();
    [k,N] = alignedTimestep(k_max, timeUnit);

    opt.k = k;
    ts = discr.getTimestepper(opt);

    q = {};
    bubble = {};
    T = [];

    if multiples(1) == 0
        r = discr.getTimeSnapshot(0);
        T(1) = r.t;
        q{1} = r.q;
        bubble{1} = r.bubble;

        multiples(1) = [];
    end

    diffs = [multiples(1) diff(multiples)];
    for i = 1:length(diffs)
        ts.stepN(N*diffs(i),true);
        r = discr.getTimeSnapshot(ts);

        T(end + 1) = r.t;
        q{end + 1} = r.q;
        bubble{end + 1} = r.bubble;
    end
end