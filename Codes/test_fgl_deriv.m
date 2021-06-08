
close all;

% Sample sine between 0 and 2pi with steps of .01 rad
h = 0.01; 
t = 0:h:2*pi;
y = sin(t);

% Compute fractional derivatives with different orders
order = 0:0.1:1;
yd    = zeros( length(order), length(t) );

for i=1:length(order)
    yd(i,:) = fgl_deriv( order(i), y, h );
end

% Plot the results
figure(1)

[X,Y] = meshgrid(t,order);
mesh(X,Y,yd);
axis equal;

title(sprintf('Fractional derivatives of sine with order between %.2f and %.2f', ...
    order(1), order(end)));
xlabel('t');ylabel('\alpha');zlabel('y');