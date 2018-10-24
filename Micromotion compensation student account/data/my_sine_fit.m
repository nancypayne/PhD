function [tt, my_fit, p] = my_sine_fit(time, data, a, f, p, b)

% check number of inputs
if nargin < 6
    error('Function requires 6 inputs: [time data amp freq phase bg]');
end

p  = [a, f, p, b];
% ub = [150, 2, 2*pi, 300];
% lb = [100, 0, 0, 100];

% define function you're fitting the data to
my_fit = @(p, t) p(1)*sin(p(2)*t + p(3)) + p(4);

p = lsqcurvefit(my_fit, p, time, data);
tt = linspace(time(1), time(end), 1001);

figure(1), clf, hold on
plot(time, data, '*')
plot(tt, my_fit(p, tt), 'Linewidth', 1)
hold off
ylim([0 1.1*max(data)])
shg

end 