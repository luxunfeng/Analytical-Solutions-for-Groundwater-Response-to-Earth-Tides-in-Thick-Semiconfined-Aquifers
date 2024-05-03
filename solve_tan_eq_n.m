function x = solve_tan_eq_n(b, n, tol)
% Solve x*tan(x) = b for the first n values of x
% Inputs:
%   b - the value of b in the equation x*tan(x) = b
%   n - the number of values of x to find
%   tol - the tolerance of the solution
% Output:
%   x - an array containing the first n solutions of the equation x*tan(x) = b

if nargin < 3
    tol = 1e-10;
end

x = zeros(1, n);
i = 1;
f = @(x) x * tan(x) - b;
internal_min = 0;
internal_max = pi/2;
while i <= n
    % if b > 1e+9
    %     x = pi/2 + pi*[0:199];
    %     break;
    % else
        if i == 1
            % Use the bisection method with default interval
            x(i) = solve_tan_eq(b);
        else
            % Use the previous solution as the lower bound of the interval
            x(i) = solve_tan_eq(b, x(i-1)+internal_min, x(i-1)+internal_max, tol);
        end

        if b < 10000
            error = 0.1;
        elseif b < 500000
            error = 10;
        elseif b < 1e+6
            error = 30;
        elseif b < 2e+6
            error = 110;
        elseif b < 1e+8
            error = 1e+6;
        else
            error = 1e+12;
        end
        % f(x(i))
        if (abs(f(x(i))) < error & i == 1) | (abs(f(x(i))) < error & abs(x(i)-x(i-1))>pi/2)
            i = i + 1;
            internal_min = 0;
            internal_max = pi/2;
        else
            internal_min = internal_min + pi/2;
            internal_max = internal_max + pi/2;
        end
        % i
    % end
end
end

function x = solve_tan_eq(b, a, c, tol)
% Solve x*tan(x) = b using the bisection method
% Inputs:
%   b - the value of b in the equation x*tan(x) = b
%   a, c - the endpoints of the interval to search for the solution
%   tol - the tolerance of the solution
% Output:
%   x - the solution of the equation x*tan(x) = b

if nargin < 4
    tol = 1e-10;
end
if nargin < 3
    c = pi/2;
end
if nargin < 2
    a = 0;
end

f = @(x) x * tan(x) - b;
fa = f(a);
fc = f(c);
if abs(b) < tol
    % If b is very small, adjust the interval to focus on x=0
    x = b/(1/3 - b/45);
else
    while (c - a) > tol
        b = (a + c)/2;
        fb = f(b);
        if fa*fb < 0
            c = b;
            fc = fb;
        else
            a = b;
            fa = fb;
        end
    end
    x = (a + c)/2;
end
end

