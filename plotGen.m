x = EulerMaruyama(@GBWB, [0,0.2,0,0], 10, 100000);
figure
subplot(2,2,1);
plot(x(:,1),x(:,2));
subplot(2,2,2);
plot(x(:,1),x(:,3));
subplot(2,2,3);
plot(x(:,1),x(:,4));