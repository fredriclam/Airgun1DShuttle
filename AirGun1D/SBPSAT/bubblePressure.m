function p = bubblePressure(y, physConst)
    R = y(1);
    E = y(4);

    gamma = physConst.gamma;

    V = 4*pi*R^3/3;
    p = E*(gamma - 1)/V;
end