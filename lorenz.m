function dy = lorenz(t,y)
dy = zeros(3,1);
P = 10;
r = 28;
b = 8/3;
dy(1) = P*(y(2) - y(1));
dy(2) = -y(1)*y(3) + r*y(1) - y(2);
dy(3) = y(1)*y(2) - b*y(3);