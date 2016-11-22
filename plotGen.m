x = EulerMaruyama(@GBWB, [0,1,0,2], 100, 100000);
figure
subplot(2,2,1);
plot(x(:,1),x(:,2));
subplot(2,2,2);
plot(x(:,1),x(:,3));
subplot(2,2,3);
plot(x(:,1),x(:,4));
