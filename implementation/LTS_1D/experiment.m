ltslvl=[1 2 3 4 5 6];
caltime=[1.354 0.847 0.701 0.664 0.621 0.622];
improvement=[1 1.334/0.847 1.334/0.701 1.334/0.664 1.334/0.621 1.334/0.622];
plot(ltslvl, caltime,'-b*')
hold on
plot(ltslvl, improvement,'-rx')
legend('run time','speed up')
title('LTS vs GTS computational efficiency graph')
xlabel('LTS level: GTS-level = 1')
ylabel('blue line: computation time(sec), red line: speed up comparing to GTS')