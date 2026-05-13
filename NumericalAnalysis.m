clc
clear
close all

syms x
prompt = 'Enter your function: ';
str = input(prompt, "s");
f = str2func(str);

prompt = 'Enter your knot points: ';
knots = input(prompt);
n = length(knots) - 1;

y = f(knots);

%% Newton
disp('Newton Interpolation');

divided_diff = zeros(n + 1, n + 1); % Table for divided differences
divided_diff(:, 1) = y; % First column is y_values

% Compute divided differences table
for j = 2:n + 1
    for i = 1:(n-j+2)
        divided_diff(i, j) = (divided_diff(i+1, j-1) - divided_diff(i, j-1)) / ...
                             (knots(i+j-1) - knots(i));
    end
end

% Compute the polynomial P(x)
P = divided_diff(1, 1); % Start with y(1)
product_term = 1; % Initialize the product term
for i = 2:n + 1
    product_term = product_term * (x - knots(i-1));
    P = P + divided_diff(1, i) * product_term;
end
disp(P)

%% Natural Cubic Spline
disp('Natural Cubic Spline');

A = 4 * eye(n - 1) + diag(ones(1, n - 2), -1) + diag(ones(1, n - 2), 1);

B = [];
for i = 2: n
    B(i - 1) = 12 * DivDiff(knots(i - 1: i + 1), y(i - 1: i + 1));
end
B = B';

M = [A B];
R = rref(M);
m = [0; R(:, end); 0];
 
h = [];
for i = 1: n
    h(i) = knots(i + 1) - knots(i);
end
h = h';

p = [];
for i = 1: n
    p(i) = DivDiff(knots(i: i + 1), y(i: i + 1)) - (m(i + 1) - m(i)) * h(i) / 6;
end
p = p';

q = [];
for i = 1: n
    q(i) = y(i) - m(i) * h(i)^2 / 6;
end
q = q';

s1 = @(x, i) (((m(i + 1) * (x - knots(i))^3) / (6 * h(i))) + ((m(i) * (knots(i + 1) - x)^3) / (6 * h(i))) + (p(i) * (x - knots(i))) + q(i));

S1 = sym('s_%d',[1 n]);
for i = 1: n
    disp(S1(i))
    disp(s1(x, i))
end

%% Clamped Cubic Spline
disp('Clamped Cubic Spline');

a = y;
b = [];
c = [];
d = [];

df = diff(f, x);
O = subs(df, x, knots(1));
N = subs(df, x, knots(n+1));

h = [];
alpha = [];
l = [];
m = [];
z = [];

for i = 1: n
    h(i) = knots(i+1) - knots(i);
end

alpha(1) = 3*(a(2) - a(1))/h(1) - 3*O;
for i = 2: n
    alpha(i) = 3*(a(i+1) - a(i))/h(i) - 3*(a(i)-a(i-1))/h(i-1);
end
alpha(n+1) = 3*N - 3*(a(n+1) - a(n))/h(n);

l(1) = 2*h(1);
m(1) = 0.5;
z(1) = alpha(1)/l(1);
for i = 2: n
    l(i) = 2*(knots(i+1) - knots(i-1)) - h(i-1)*m(i-1);
    m(i) = h(i)/l(i);
    z(i) = (alpha(i) - h(i-1)*z(i-1))/l(i);
end
l(n+1) = h(n)*(2 - m(n));
z(n+1) = (alpha(n+1) - h(n)*z(n))/l(n+1);

c(n+1) = z(n+1);
for i = n: -1: 1
    c(i) = z(i) - m(i)*c(i+1);
    b(i) = (a(i+1) - a(i))/h(i) - h(i)*(c(i+1) + 2*c(i))/3;
    d(i) = (c(i+1) - c(i))/(3*h(i));
end

s2 = @(x, i) a(i) + b(i)*(x - knots(i)) + c(i)*(x - knots(i))^2 + d(i)*(x - knots(i))^3;

S2 = sym('s_%d',[1 n]);
for i = 1: n
    disp(S2(i))
    disp(s2(x, i))
end

%%
figure(1);

% Plot f(x) and store the handle
subplot(2, 1, 1);
h1 = fplot(f, [knots(1), knots(end)], 'k', 'Linewidth', 2);
hold on;

% Plot Newton's method and store the handle
h2 = fplot(P, [knots(1), knots(end)], '--', 'Color', 'r');
hold on;

% Plot the natural cubic spline and store the handles in an array
h3 = gobjects(1, n); % Preallocate array for handles
for i = 1:n
    subplot(2, 1, 1);
    h3(i) = fplot(s1(x, i), [knots(i), knots(i + 1)], '-.*', 'Color', 'm');
    hold on;
end

% Plot the clamped cubic spline and store the handles in an array
h4 = gobjects(1, n); % Preallocate array for handles
for i = 1:n
    subplot(2, 1, 1);
    h4(i) = fplot(s2(x, i), [knots(i), knots(i + 1)], '-o', 'Color', '[0, 0.5, 0]');
    hold on;
end

% Add legend with unique handles
legend([h1, h2, h3(1), h4(1)], 'f(x)', 'Newton', 'Natural Cubic Spline', 'Clamped Cubic Spline');

% Add title
title('Interpolations');


% Error
subplot(2, 1, 2);

% Plot the error for Newton's method and store the handle
h1 = fplot(abs(f - P), [knots(1), knots(end)], '--', 'Color', 'r');
hold on;

% Plot the error for the natural cubic spline and store the handles in an array
h2 = gobjects(1, n); % Preallocate array for handles
for i = 1:n
    h2(i) = fplot(abs(f - s1(x, i)), [knots(i), knots(i + 1)], '-.*', 'Color', 'm');
    hold on;
end

% Plot the error for the clamped cubic spline and store the handles in an array
h3 = gobjects(1, n); % Preallocate array for handles
for i = 1:n
    h3(i) = fplot(abs(f - s2(x, i)), [knots(i), knots(i + 1)], '-o', 'Color', '[0, 0.5, 0]');
    hold on;
end

% Add legend with unique handles
legend([h1, h2(1), h3(1)], 'Newton', 'Natural Cubic Spline', 'Clamped Cubic Spline');

% Add title
title('Errors');


%%
function div = DivDiff(x, y)
    difference = y;
    
    for i = 2: length(y)
        temp = difference;
        
        for j = i: length(y)
            difference(j) = (temp(j) - temp(j - 1)) / (x(j) - x(j - (i - 1)));
        end
    end
    
    div = difference(length(y));
end
