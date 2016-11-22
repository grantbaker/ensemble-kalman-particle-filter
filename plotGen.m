x = EulerMaruyama(@GBWB, [0,1,0,1], 10, 10000);
figure
plot(x(:,1),x(:,2));
figure
plot(x(:,1),x(:,3));
figure
plot(x(:,1),x(:,4));