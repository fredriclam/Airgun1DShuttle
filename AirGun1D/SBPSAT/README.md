# README #

To make it run first make sure you have the sbplib toolbox installed.

By running
    `noname.animate(discr,speed)`
you'll see an animation of the problem described by discr at the given speed.

An example would be:
`noname.animate(DiscrAirgun(200,3),0.0005)`
The first parameter of DiscrAirgun decides the number of grid points per meter and the second is the order of accuracy.

To fiddle with the airgun+bubble model the main files of interest are(in no particular order):

* DiscrAirgun.m
* bubblePressure.m
* bubbleRHS.m
* bubbleTest.m
* configAirgun.m
* icAirgunShocktube.m
* setupEulerPlot.m

To change the initial conditions see DiscrAirgun.m.