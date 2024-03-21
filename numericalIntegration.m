function [trap, mid, simp] = numericalIntegration(f, a, b, r)
% The function numericalIntegration accepts an anonymous function f, left endpoint of an 
% interval, a, right endpoint b, and the number of sub-intervals r, then calculates the 
% approximate definite integral. The program sums the partitions of the interval using the 
% techniques of the composite trapezoidal rule, composite midpoint rule, and the composite 
% Simpson's rule. The function outputs 3 scalar values in the above order, i.e. the function 
% needs to be pass these numbers to 3 separate variables. A non-zero positive integer is
% expected for r, however, if a negative integer is sent to the function, it's sign will be
% changed. If r is zero or non-integer, then the program will quit and assign NaN to the
% output variables. Additionally, the user of this function should notice that Simpson's 
% rule is intended to be executed with an even number of sub-intervals. If an odd number is 
% used, an approximating is still calculated, however, it is less accurate in general. 
%
%          f = The anonymous function that the user wishes to integrate
%          a = The LEFT endpoint of the interval we want to integrate over (real number)
%          b = The RIGHT endpoint of the interval we want to integrate over (real number)
%          r = The number of equidistant points we want to partition the interval with 
%               (non-zero integer)

    % If r is 0 or a non-natural number, assign NaN to the ouput variables
    if (r == 0 || r - floor(r) ~= 0)
        disp('Error: r must be a natural number.');
        trap = NaN; 
        mid  = NaN; 
        simp = NaN;
        return
    end
    
    % If r is negative make it positive
    if (r < 0)
        disp('Sign changed for variable r.');
        r = (-1)*r;
    end

    % Set the length of the sub-interval
    h = (b - a)/r;

    % Create a vector of sub-interval points and sub-interval midpoints
    subPts = a:h:b;
    midPts = (subPts(:, 1:end - 1) + subPts(:, 2:end)) / 2;
    
    % Create a vector of quadrature abscissae for evaluating Simpson's method
    t2k   = subPts(3:2:end - 2);
    t2k_1 = subPts(2:2:end);

    % Trapezoidal, Midpoint, and Simpson's Rule
    trap = h/2*(f(a) + sum(2*f(subPts(2:end - 1))) + f(b));
    mid  = h  *(sum(f(midPts)));
    simp = h/3*(f(a) + 2*sum(f(t2k)) + 4*sum(f(t2k_1)) + f(b));
end