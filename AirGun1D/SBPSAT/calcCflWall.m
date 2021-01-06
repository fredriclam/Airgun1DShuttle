function calcCflWall()
    m = 200;
    order = [2,3,4,5,6,7,8,9];
    T = 1;
    threshold = 1e5;
    silent = true;

    for o = order
        noname.testCfl(DiscrWallTest(m,o,icDensityBell()),'rk4euler', T, [0 , 1],[],threshold,silent);
    end
end
