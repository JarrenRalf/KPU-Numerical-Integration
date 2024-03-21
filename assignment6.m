%% Question 3 (a) & (b)
clear; clc; close all;

syms x

% Set the integrands up as anonymous functions
f1 = @(x, y) 4./(1 + x.^2);
f2 = @(x) sqrt(x);

% Set the endpoints
a = 0;       b = 1;

% Choose the number of sub-intervals, r, and calculate sub-interval length, h
e = 1:5;     r = 2.^e;     h = (b - a)./r;

% Initialize the matrices of solutions
A1 = zeros(length(r), 3);
A2 = zeros(length(r), 3); 

% Execute the integral approximations and fill the matrices
for i = 1:length(r)
    [t1, m1, s1] = numericalIntegration(f1, a, b, r(i));
    [t2, m2, s2] = numericalIntegration(f2, a, b, r(i));
    A1(i, :) = [t1 m1 s1];
    A2(i, :) = [t2 m2 s2];
end

% Find the appropriate derivatives for calculating the error estimations
d2f1 = matlabFunction( (-1)*diff( diff( f1(x)))); % Second derivative for trap and midpoint
d2f2 = matlabFunction( (-1)*diff( diff( f2(x))));
d4f1 = matlabFunction( (-1)*diff( diff( diff( diff( f1(x)))))); % Fourth deriv for Simpson's
d4f2 = matlabFunction( (-1)*diff( diff( diff( diff( f2(x))))));

% Find the maximum of the derivatives in order to use the error estimation formulas
[x1, fval_d2_1] = fminbnd(d2f1, a, b); % Notice that this function finds the minimum value,
[x2, fval_d2_2] = fminbnd(d2f2, a, b); % so that is why the functions above are multiplied
[x3, fval_d4_1] = fminbnd(d4f1, a, b); % by -1. We minimize -f which is equivalent to
[x4, fval_d4_2] = fminbnd(d4f2, a, b); % finding the max.

% Use matlab's built in function to determine an accurate value of each integral
q1 = integral(f1, a, b);
q2 = integral(f2, a, b);

% Print the results of function 1
fprintf('For function 1:\nThe true value of the integral is %0.17f\n\n', q1);
fprintf('r \t trapezoidal \t\t midpoint \t\t\t Simpson''s\n');
for i = 1:length(r)
    fprintf('%g \t %0.12f \t %0.12f \t %0.12f \n', r(i), A1(i, 1), A1(i, 2), A1(i, 3));
end
fprintf('\nERROR ESTIMATES:\n');
fprintf('r \t trapezoidal \t\t midpoint \t\t\t Simpson''s\n');
for i = 1:length(r)
    fprintf('%g \t %0.12f \t %0.12f \t %0.12f \n', r(i), abs(fval_d2_1)/12*(b - a)*h(i)^2, ...
                       abs(fval_d2_1)/24*(b - a)*h(i)^2, abs(fval_d4_1)/180*(b - a)*h(i)^4);
end
fprintf('\nABSOLUTE ERROR:\n');
fprintf('r \t trapezoidal \t\t midpoint \t\t\t Simpson''s\n');
for i = 1:length(r)
    fprintf('%g \t %0.12f \t %0.12f \t %0.12f \n', r(i), ...
                                abs(q1 - A1(i, 1)), abs(q1 - A1(i, 2)), abs(q1 - A1(i, 3)));
end

% Print the results for function 2
fprintf('\n\nFor function 2:\nThe true value of the integral is  %0.16f\n\n', q2);
fprintf('r \t trapezoidal \t\t midpoint \t\t\t Simpson''s\n');
for i = 1:length(r)
    fprintf('%g \t %0.12f \t %0.12f \t %0.12f \n', r(i), A2(i, 1), A2(i, 2), A2(i, 3));
end
fprintf('\nERROR ESTIMATES:\n');
fprintf('r \t trapezoidal \t\t midpoint \t\t\t Simpson''s\n');
for i = 1:length(r)
    fprintf('%g \t %0.12f \t %0.12f \t %0.12f \n', r(i), abs(fval_d2_2)/12*(b - a)*h(i)^2, ...
                       abs(fval_d2_2)/24*(b - a)*h(i)^2, abs(fval_d4_2)/180*(b - a)*h(i)^4);
end
fprintf('\nABSOLUTE ERROR:\n');
fprintf('r \t trapezoidal \t\t midpoint \t\t\t Simpson''s\n');
for i = 1:length(r)
    fprintf('%g \t %0.12f \t %0.12f \t %0.12f \n', r(i), ...
                                abs(q2 - A2(i, 1)), abs(q2 - A2(i, 2)), abs(q2 - A2(i, 3)));
end

S = load('train');
sound(S.y, S.Fs)
%% Question 3
clear; clc; close all;

% Set the values of n
n = 0:9;

% Initialize a square matrix of NaNs
% The advantage of this is that many functions ignore NaNs and if zeros() was used
% for example, the zeros in the matrix may be mistaken for actual Gauss points
R = NaN(length(n));

for i = n
    
    % Find the coefficients of the legendre polynomial up to order i + 1
    p = legendrepol(i + 1);
    
    % Find the roots of the legendre polynomial using the coefficients
    R(i + 1, 1:i + 1) = roots(p(i + 2, :));

    % Graph the interpolation for each value of epsilon
    subplot(length(n) + 1, 1, i + 1);
    plot(R(i + 1, :), 0, 'ok', 'MarkerFaceColor', 'r'); % plot() ignores NaNs by default
    set(gca, 'YTick',      []); % Remove YTicks and YTickLabels in order
    set(gca, 'YTickLabel', []); % to clean up the display of the graph
    axis([-1 1 -0.00000000001 0.00000000001])
    title(['n = ', num2str(i)]);
end
%% Question 5
clear; clc; close all;

format long

% Set the integrands up as anonymous functions
f = @(x) 4./(1 + x.^2);

% Set the desired accuracy
s = 5;

% Set the endpoints
a = 0;
b = 1;

% Choose the number of sub-intervals
r = 1;

% Calculate the sub-interval length
h = (b - a)/r;

% Initialize the Lower Triangular Romberg Matrix (table)
R = NaN(s);

% Assigning the output of the numericalIntegration function to one variable
% will give us the trapezoidal rule
R(1, 1) = numericalIntegration(f, a, b, r);

for j = 1:s - 1

    % Set the element in the first column of the table
    % as you loop through starting with the second row
    h = h/2;
    R(j + 1, 1) = (1/2)*R(j, 1) + h*sum(f(a + h*(2*(1:r*2^(j - 1)) - 1)));

    % Calculate the rest of the elements in the row
    for k = 2:j + 1
        R(j + 1, k) = R(j + 1, k - 1)+(R(j + 1, k - 1) - R(j, k - 1))/(4^(k - 1) - 1);
    end
end

disp(R);
%% Tupper Formula (See second example below this one)
clear; clc; close all;

% use the symbolic toolbox to represent the big integer k
k =  sym(['960939379918958884971672962127852754715004339660129306651505'...
    '519271702802395266424689642842174350718121267153782770623355993237'...
    '280874144307891325963941337723487857735749823926629715517173716995'...
    '165232890538221612403238855866184013235585136048828693337902491454'...
    '229288667081096184496091705183454067827731551705405381627380967602'...
    '565625016981482083418783163849115590225610003652351370343874461848'...
    '378737238198224849863465033159410054974700593138339226497249461751'...
    '545728366702369745461014655997933798537483143786841806593422227898'...
    '388722980000748404719']);

[x, y] = meshgrid(0:1:106, 0:1:16);

% evaluate the tupper formula
tupper = rem(floor(floor((y + k)/17).*2.^(-17*x - rem((y + k), 17))), 2);

% convert from symbolic to Matlab's native double precision
tupper = double(tupper);

% display it!
image(fliplr((1 - tupper)*255));
colormap gray
axis equal

title('The plot of $\frac{1}{2} < \lfloor \mathrm{mod}\left(\lfloor\frac{y}{17}\rfloor \cdot 2^{-17 \lfloor x \rfloor - \mathrm{mod}(\lfloor y \rfloor, 17)} , 2\right) x \rfloor$ is itself!', 'interpreter', 'latex');
set(gca, 'XTick', [], 'YTick', []);
%%
clear; clc; close all;

% Now for a more personal example
k =  sym(['556014034746548166146205054173063297862594853701573742373872'...
    '175112122727492963917901056239208795869656877418598395459314703438'...
    '281795305077787595579262007704941397684296104625695528274221544262'...
    '718896215851682809114362015192708725022964672763023913982492824354'...
    '702841353267453853116957056003262024527982941886999364038567107321'...
    '469297746683565991506866930024289075081954135640598479778748355586'...
    '525111821247278019041923679168041069006326019221117402291878865862'...
    '250348784333799204064815265501581152867106173566302312167301258591'...
    '371883188021262574']);

[x, y] = meshgrid(0:1:106, 0:1:16);

% Evaluate the tupper formula
tupper = rem(floor(floor((y + k)/17).*2.^(-17*x - rem((y + k), 17))), 2);

tupper = double(tupper);

% Some background info
web('https://www.youtube.com/watch?v=_s5RFgd59ao');

% I made the message ussing the following website
%     http://keelyhill.github.io/tuppers-formula/

% Display the message
image(fliplr((1 - tupper)*255));
colormap gray
axis equal
set(gca, 'XTick', [], 'YTick', []);