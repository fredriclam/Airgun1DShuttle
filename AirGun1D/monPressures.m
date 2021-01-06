figure(77);
subplot(3,1,1)
plot(monitorStates(13,:), monitorStates(10,:))
hold on
plot(monitorStates(13,:), monitorStates(15,:))
plot(monitorStates(13,:), monitorStates(16,:))
hold off;
legend({'pFiring ->', 'pRear ->', 'pFront <-'})
xlabel 't [s]'
subplot(3,1,2)
plot(monitorStates(13,:), monitorStates(5,:))
xlabel 't [s]'
ylabel 'uPort [m/s]'
subplot(3,1,3)
plot(monitorStates(13,:), monitorStates(3,:), '.')
xlabel 't [s]'
ylabel 'State (#)'