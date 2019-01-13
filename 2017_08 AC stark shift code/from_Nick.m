close all
m1 = 10000;
m2 = 9000;
m3 = 8000;
k = 3e6;

odefun = @(t,x) [x(4) ;
                 x(5) ;
                 x(6) ;
                 k/m1*(x(2)-x(1)) - k/m1*x(1) ;
                 k/m2*(x(3)-x(2)) - k/m2*(x(2)-x(1)) ;
                 -k/m3*(x(3)-x(2)) ];
[t, x] = ode45(odefun, linspace(0,20,1000), [.05, .08, .1, 0, 0, 0]);

x1 = x(1,1);
x2 = x(1,2);
x3 = x(1,3);
bs = .5; % box size
v = bs*[-1 -1 0; 1 -1 0; 1 1 0; -1 1 0];
p0 = patch('vertices', v+[0  0 0], 'faces', 1:4, 'facecolor', [0 0 0]);
p1 = patch('vertices', v+[x1 1 0], 'faces', 1:4, 'facecolor', [1 0 0]);
p2 = patch('vertices', v+[x2 2 0], 'faces', 1:4, 'facecolor', [0 1 0]);
p3 = patch('vertices', v+[x3 3 0], 'faces', 1:4, 'facecolor', [0 0 1]);

xmax = max(x(:));
xlim([-xmax xmax])

for k = 2:length(x)
    % New positions:
    x1 = x(k,1);
    x2 = x(k,2);
    x3 = x(k,3);
    
    set(p1, 'vertices', v+[x1 1 0])
    set(p2, 'vertices', v+[x2 2 0])
    set(p3, 'vertices', v+[x3 3 0])
    xlim([-xmax xmax])
    
    % Draw and pause for a moment
    drawnow
    pause(eps)
end