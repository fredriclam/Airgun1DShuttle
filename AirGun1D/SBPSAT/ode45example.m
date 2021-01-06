function sol = ode45example();
%[t, y] = ode45example()
    d = DiscrAirgun(200,3);

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

    sol = ode45(@odefun, tspan, y0);
    %[t, y] = ode45(@odefun, tspan, y0);

    %plot(t,y(:,1))
    %ylabel('R')
    %xlabel('t')

end
