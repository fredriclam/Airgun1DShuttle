function sol = runEulerCode_initBubbleVol(nx,initBubbleVol,airgunPressure, airgunLength, airgunPortArea, airgunDepth)
%sol = runEulerCode(nx)
    d = DiscrAirgun_initBubbleVol(nx,3,initBubbleVol,airgunPressure, airgunLength, airgunPortArea, airgunDepth);

    q0 = d.q0;
    t0 = d.t0;
    bubble0 = d.bubble0;
    RHS = d.RHS;

    N = length(bubble0);
    function dy = odefun(t,y)
        bubble = y(1:N);
        q = y(N+1:end);

        [dq, dBubble] = RHS(q,t,bubble);
        dy = [dBubble; dq];
    end

    y0 = [bubble0; q0];
    tspan = [0; 1];
    options = odeset('RelTol',1e-6);
    sol = ode45(@odefun, tspan, y0,options);
end
