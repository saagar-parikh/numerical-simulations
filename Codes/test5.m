%% Numerical simulations and stability analysis of a COVID-19 model using fractional derivatives
% 

tic                         % Start the clock to calculate the execution time

%% Initalize variables
clear;

size_k = 121;               % Size of the vectors

% A = Number of susceptible people
% B = Number of exposed people
% C = Number of infected people
% D = Number of asymptotically people (people showing no symptoms)
% E = Number of recovered people
% F = Number of people affected by outbreak at reservoir like market place 

% Initializing vectors for the subgroups used for the fractional order derivatives
A = zeros(size_k,1);
B = zeros(size_k,1);
C = zeros(size_k,1);
D = zeros(size_k,1);
E = zeros(size_k,1);
F = zeros(size_k,1);

% Initializing vectors for integer-order derivative
global Avec Bvec Cvec Dvec Evec Fvec

Avec = zeros(size_k,1);
Bvec = zeros(size_k,1);
Cvec = zeros(size_k,1);
Dvec = zeros(size_k,1);
Evec = zeros(size_k,1);
Fvec = zeros(size_k,1);

% Initializing the parameters
global h mu1 mu2 a1 b1 b2 sigma d f1 f2 e1 e2 g1 g2 N

% Initialize the values of these parameters
h = 0.01;                       % Step size
mu1 = 107644.22451;             % Natural birth rate
mu2 = 0.01;                     % Removing rate of the virus from the seafood market
a1 = 1/76.79;                   % Natural Death Rate
b1 = 0.05;                      % A and D are related by sigma b1 A D
b2 = 0.000001231;               % Disease spread coefficient
sigma = 0.02;                   % Transmissibility multiple of D to C
d = 0.1243;                     % Proportion of the asymptomatic infection
f1 = 0.00047876;                % Spread rate after completing the incubation period and becomes infected
f2 = 0.001;                     % Rate at which D contributes the virus to the seafood market
e1 = 0.005;                     % Spread rate joining the classes C and D
e2 = 0.000398;                  % Rate at which C contributes the virus to the seafood market
g1 = 0.09871;                   % Recovery rate
g2 = 0.854302;                  % Removal rate
L = 1e5;                        % Memory length

% Given Initial values
A(1) = 8065518;
B(1) = 200000;
C(1) = 282;
D(1) = 200;
E(1) = 0;
F(1) = 50000;
Avec(1) = 8065518;
Bvec(1) = 200000;
Cvec(1) = 282;
Dvec(1) = 200;
Evec(1) = 0;
Fvec(1) = 50000;

% Different fractional values of the time derivative
p1 = 0.9;
p2 = 0.9;
p3 = 0.9;
p4 = 0.9;
p5 = 0.9;
p6 = 0.9;

% p1 = [0.6 0.7 0.8 0.9 0.99];
% p2 = [0.6 0.7 0.8 0.9 0.99];
% p3 = [0.6 0.7 0.8 0.9 0.99];
% p4 = [0.6 0.7 0.8 0.9 0.99];
% p5 = [0.6 0.7 0.8 0.9 0.99];
% p6 = [0.6 0.7 0.8 0.9 0.99];


N = A(1) + B(1) + C(1) + D(1) + E(1);         % Total Number of people in the case study


% Storing the values in the defined functions
for k2 = 0:h:(size_k*h)-h
    k2;
    funcA(k2);
    funcB(k2);
    funcC(k2);
    funcD(k2);
    funcE(k2);
    funcF(k2);
end

%% Solving the fractional-order Grünwald–Letnikov derivative

for k = 2:(size_k)-1
    
    % GL derivative for A(k)
    sum1_A = 0;             % temp variable to store the summation
    j1 = 0;                 % summation variable
    while (k-j1)>0
        if j1==0
            frac = 1;       % cannot divide by zero, so set the value = 1
        else
            frac = p1/j1;   % fraction value used in the equation
        end
        sum1_A = sum1_A + (1/(h^p1))*((1)^j1)*frac*(Avec(k-j1));    % GL derivative eqn summation
        j1 = j1+1;          % increase counter
    end
    sum_A = 0;              % temp variable to store the summation
    for j = v(L,k,h):k
        sum_A = sum_A + c(p1,j)*A(k-j+1);        % Second summation part of the GL derivative solution
    end
    A(k) = sum1_A*(h^p1) - sum_A;                % Final solution of A(k)
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Repeat these steps for B(k),C(k),D(k),E(k),F(k)
    
    % GL derivative for B(k)
    sum_B = 0;
    for j = v(L,k,h):k
        sum_B = sum_B + c(p2,j)* B(k-j+1);        
    end
    sum1_B = 0;
    j1 = 0;
    beta = 2.0195e5 + B(1);          % correction constant
    while (k-j1)>0
        if j1==0
            frac = 1;
        else
            frac = p2/j1;
        end
        sum1_B = sum1_B + (1/(h^p2))*((-1))*(frac)*(Bvec(k-j1));
        j1 = j1+1;
    end
    B(k) = sum1_B*(h^p2) - sum_B + beta;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GL derivative for C(k)
    sum_C = 0;
    for j = v(L,k,h):k
        sum_C = sum_C + c(p3,j)*C(k-j+1);        
    end
    sum1_C = 0;
    j1 = 0;
    alpha = 100;
    beta = 2.7969e4 + C(1);          % correction constant
    while (k-j1)>0
        if j1==0
            frac = 1;
        else
            frac = p3/j1;
        end
        sum1_C = sum1_C + (1/(h^p3))*((-1))*(frac)*(Cvec(k-j1))*alpha;
        j1 = j1+1;
    end
    C(k) = sum1_C*(h^p3) - sum_C + beta;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GL derivative for D(k)
    sum_D = 0;
    for j = v(L,k,h):k
        sum_D = sum_D + c(p4,j)*D(k-j+1);        
    end
    sum1_D = 0;
    j1 = 0;
    alpha = 10;
    beta = 1.9741e3 + D(1);          % coreection constant
    while (k-j1)>0
        if j1==0
            frac = 1;
        else
            frac = p4/j1;
        end
        sum1_D = sum1_D + (1/(h^p4))*((-1))*(frac)*(Dvec(k-j1))*alpha;
        j1 = j1+1;
    end
    D(k) = sum1_D*(h^p4) - sum_D + beta;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GL derivative for E(k)
    sum_E = 0;
    for j = v(L,k,h):k
        sum_E = sum_E + c(p5,j)*E(k-j+1);        
    end
    sum1_E = 0;
    j1 = 0;
    alpha = 2e3;
    while (k-j1)>0
        if j1==0
            frac = 1;
        else
            frac = p5/j1;
        end
        sum1_E = sum1_E + (1/(h^p5))*((1)^j1)*(frac)*(Evec(k-j1))*alpha;
        j1 = j1+1;
    end
    E(k) = sum1_E*(h^p5) - sum_E;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GL derivative for F(k)
    sum_F = 0;
    for j = v(L,k,h):k
        sum_F = sum_F + c(p6,j)*F(k-j+1);        
    end
    sum1_F = 0;
    j1 = 0;
    while (k-j1)>0
        if j1==0
            frac = 1;
        else
            frac = p6/j1;
        end
        sum1_F = sum1_F + (1/(h^p6))*((1)^j1)*(frac)*(Fvec(k-j1));
        j1 = j1+1;
    end
    F(k) = sum1_F*(h^p6) - sum_F;


end


%% Plotting A, B, C, D, E, F with respect to time

t = 0:(size_k)-1;           % Time vector


figure;
plot(t, A(t+1), 'Color', 'r');
axis([0 120 0 9e6]);
xlabel('t');
ylabel("A(t)");

figure;
plot(t, B(t+1), 'Color', 'r');
axis([0 120 0 5e5]);
xlabel('t');
ylabel("B(t)")

figure;
plot(t, C(t+1), 'Color', 'r');
axis([0 120 0 4e4]);
xlabel('t');
ylabel("C(t)")

figure;
plot(t, D(t+1), 'Color', 'r');
axis([0 120 0 5500]);
xlabel('t');
ylabel("D(t)")

figure;
plot(t, E(t+1), 'Color', 'r');
axis([0 120 0 2e4]);
xlabel('t');
ylabel("E(t)")

figure;
plot(t, F(t+1), 'Color', 'r');
axis([0 120 0 5e4]);
xlabel('t');
ylabel("F(t)")


%% Stability analysis
% Define 2 equilibrium points: disease-free and endemic-equilibrium point
% Define the 2 Jacobian matrices at the two equilibrium points
% Obtain the 6 degree polynomial equations in eigenvalues(lambda)
% Find roots using the Newton-Raphson method 
% Check the stability of the model based on the eigenvalues obtained


syms J1 J2 lambda
syms func1(lambda) func2(lambda)

global func1 func2

% Identity matrix
I = [1, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1];

% Jacobian matrix at the disease-free equilibrium point
J1 = [ -0.0130 0 -0.0500 -0.0010 0 -10.1754 ;
       0 -0.0141 0.0500 0.0010 0 10.1754 ;
       0 0.0004 -0.1117 0 0 0  ; 
       0 0.0006 0 -0.8673 0 0 ; 
       0 0 0.0987 0.8543 -0.0130 0 ;
       0 0 0.0003 0.0010 0 -0.0100 ];
 
% Jacobian matrix at the endemic-equilibrium point
J2 = [ -0.0023 0 -0.2885 -0.0058 0 -58.7180 ;
       -0.0108 -0.0141 0.2885 0.0058 0 58.7180 ; 
       0 0.0004 -0.1117 0 0 0 ;
       0 0.0006 0 -0.8673 0 0 ;
       0 0 0.0987 0.8543 -0.0130 0 ;
       0 0 0.0003 0.0010 0 -0.0100 ];
   
% Functions for the characteistic equations to find the eigenvalues
func1 = @(lambda) mat_det(J1 - lambda*I);
func2 = @(lambda) mat_det(J2 - lambda*I);

% We will use the Newton Raphson method to find the roots of the above two
% 6 degree polynomial equations for lambda

% For func1, we use the following 6 initial estimates to find the roots
estimates1 = [-0.0140 -0.0120 -0.1110 -0.0060 -0.0170 -0.8670];

for i = 1:6
    % Eigenvalues corresponding to matrix J1
    root1(i) = newtonraphson1(estimates1(i)); 
end

% For func2, we use the following 6 initial estimates to find the roots
estimates2 = [-0.0140 -0.8670 -0.1110 -0.0200 -0.0125 -0.0061];

for i = 1:6
    % Eigenvalues corresponding to matrix J1
    root2(i) = newtonraphson2(estimates2(i)); 
end

root1
root2

% We see that the obtained eigenvalues corresponding to both the Jacobian
% matrices are all negative, which shows that the chosen equilibrium points
% are stable

% Plots of func1 and func2
% 
% fplot(lambda, func1)
% line(xlim, [0,0], 'Color', 'Black', 'LineWidth', 0.5);
% fplot(lambda, func2)
% line(xlim, [0,0], 'Color', 'Black', 'LineWidth', 0.5);


% Time taken for computation
timetaken = toc
%% Function definitions

% Solving integral function of the integer-order derivative using Trapezoidal method
% Difference between two data points (t2 - t1) = dt

function answer = funcA(t)
    global h mu1 a1 b1 sigma N
    global Avec
    
    if Avec(floor((t/h)+1)) ~= 0            % If value at particular t is already stored, no need to calculate again (Memoization)
        answer = Avec(floor((t/h)+1));
    elseif t == 0
        answer = 8065518;                   % Initial condition
        Avec(floor((t/h)+1)) = answer;
    else
        dt = h;
        answer = 0;
        it = 0;
        % Trapezoidal method for calculating the integral from 0 to t
        while (it+dt)<t
            answer = answer + (((mu1 - a1*funcA(it) - (b1*funcA(it)*(funcC(it) + (sigma*funcD(it))))/N) ...
                            + (mu1 - a1*funcA(it+dt) - (b1*funcA(it+dt)*(funcC(it+dt) + (sigma*funcD(it+dt))))/N ))/2)*dt;
            it = it + dt;
        end
        Avec(floor((t/h)+1)) = answer;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the same method as funcA to define funcB,funcC,funcD,funcE,funcF


function answer = funcB(t)
    global h a1 b1 b2 sigma d f1 e1 N
    global Bvec
    
    if Bvec(floor((t/h)+1)) ~= 0
        answer = Bvec(floor((t/h)+1));
    elseif t == 0
        answer = 200000;
        Bvec(floor((t/h)+1)) = answer;
    else
        dt = h;
        answer = 0;
        it = 0;
        while (it+dt)<t
            answer = answer + ((((b1*funcA(it)*(funcC(it) + (sigma*funcD(it))))/N + ...
                                b2 * funcA(it) * funcF(it) - ((1 - d)* f1 * funcB(it)) - ...
                                d * e1 * funcB(it) - a1 * funcB(it))...
                            + ((b1*funcA(it+dt)*(funcC(it+dt) + (sigma*funcD(it+dt))))/N + ...
                                b2 * funcA(it+dt) * funcF(it+dt) - ((1 - d)* f1 * funcB(it+dt)) - ...
                                d * e1 * funcB(it+dt) - a1 * funcB(it+dt)))/2)*dt;
            it = it + dt;
        end
        Bvec(floor((t/h)+1)) = answer;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function answer = funcC(t)
    global h a1 d f1 g1
    global Cvec
    
    if Cvec(floor(floor((t/h)+1))) ~= 0
        answer = Cvec(floor((t/h)+1));
    elseif t == 0
        answer = 282;
        Cvec(floor((t/h)+1)) = answer;
    else
        dt = h;
        answer = 0;
        it = 0;
        while (it+dt)<t
            answer = answer + ((((1-d)*f1*funcB(it) - (g1 + a1)*funcC(it))...
                            + ((1-d)*f1*funcB(it+dt) - (g1 + a1)*funcC(it+dt)))/2)*dt;
            it = it + dt;
        end
        Cvec(floor((t/h)+1)) = answer;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function answer = funcD(t)
    global h a1 d e1 g2
    global Dvec
    
    if Dvec(floor((t/h)+1)) ~= 0
        answer = Dvec(floor((t/h)+1));
    elseif t == 0
        answer = 200;
        Dvec(floor((t/h)+1)) = answer;
    else
        dt = h;
        answer = 0;
        it = 0;
        while (it+dt)<t
            answer = answer + (((d*e1*funcB(it) - (g2 + a1)*funcD(it))...
                            + (d*e1*funcB(it+dt) - (g2 + a1)*funcD(it+dt)))/2)*dt;
            it = it + dt;
        end
        Dvec(floor((t/h)+1)) = answer;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function answer = funcE(t)
    global h a1 g1 g2
    global Evec
    
    if Evec(floor((t/h)+1)) ~= 0
        answer = Evec(floor((t/h)+1));
    elseif t == 0
        answer = 0;
        Evec(floor((t/h)+1)) = answer;
    else
        dt = h;
        answer = 0;
        it = 0;
        while (it+dt)<t
            answer = answer + (((g1*funcC(it) + g2*funcD(it) - a1*funcE(it))...
                            + (g1*funcC(it+dt) + g2*funcD(it+dt) - a1*funcE(it+dt)))/2)*dt;
            it = it + dt;
        end
        Evec(floor((t/h)+1)) = answer;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function answer = funcF(t)
    global h mu2 f2 e2
    global Fvec
    
    if Fvec(floor((t/h)+1)) ~= 0
        answer = Fvec(floor((t/h)+1));
    elseif t == 0
        answer = 50000;
        Fvec(floor((t/h)+1)) = answer;
    else
        dt = h;
        answer = 0;
        it = 0;
        while (it+dt)<t
            answer = answer + (((e2*funcC(it) + f2*funcD(it) - mu2*funcF(it))...
                            + (e2*funcC(it+dt) + f2*funcD(it+dt) - mu2*funcF(it+dt)))/2)*dt;
            it = it + dt;
        end
        Fvec(floor((t/h)+1)) = answer;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculating coefficients of derivatives 
function output = c(p, j)
    output = 1;
    for i = 1:j
        output = (1 - ((1+p)/i))* output;
    end
end

% Calculating lower index by using short memory principle
function output = v(L,k,h)

    if (k <= (L/h))
        output = 1;
    else
        output = k - (L/h);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Newton Raphson method for finding the roots of func1 (6 degree polynomial eqn)
function root = newtonraphson1(x0)
    
    global func1

    h = 1e-10;                           % For differentiation
    x_old = x0;
    while 1
        
        g = (func1(x_old + h) - func1(x_old))/h;    % f'(x)
        x_new = x_old - func1(x_old)/g;             % Find the next best estimate

        % Condition for breaking the loop
        if (abs(func1(x_new))<1e-15)
            root = x_new;
            break 
        end
        
        % The new estimate becomes the old estimate for the next loop
        x_old = x_new;
    end
end

% Newton Raphson method for finding the roots of func1 (6 degree polynomial eqn)
function root = newtonraphson2(x0)
    
    global func2

    h = 1e-10;                           % For differentiation
    x_old = x0;
    while 1
        
        g = (func2(x_old + h) - func2(x_old))/h;    % f'(x)
        x_new = x_old - func2(x_old)/g;             % Find the next best estimate

        % Condition for breaking the loop
        if (abs(func2(x_new))<1e-12)
            root = x_new;
            break 
        end
        
        % The new estimate becomes the old estimate for the next loop
        x_old = x_new;
    end
end


% Function to find the determinant of a matrix A without using inbuilt function
function [sum] = mat_det(A)
    [m,n] = size(A);                          % size of matrix A 
    if(n == 1)                                % matrix with only one element
        sum = A(1,1);
    elseif(n == 2)                            % matrix with 4 elements
        sum = A(1,1)*A(2,2) - A(1,2)*A(2,1);
    else
        sum = 0;                              % temp variable to store the summation 
        for i = 1:m
           if(mod (i,2) == 0)                 % when i is even number
               p = 2*i - 1;
           else                               % when i is odd
               p = 2*i;
           end
            B = minor_gen(A,m,i);             % minor for coefficients  
            sum = sum+(((-1)^p)*A(1,i))*B;    % determinant summation
        end
    end
end

% Calculating the minor for coefficients
function [sum] = minor_gen(A,n,k)
    t1 = 0; B = [];                           % initializing B matrix
    for i = 2:n                               % element of ith row
        t1 = t1 + 1; t2 = 0;
        if(i == 1)
            i = i + 1;
        else
        for j = 1:n                           % element of jth column
            if(j == (k))
                j = j + 1;
            else
             t2 = t2 + 1;
            B(t1,t2) = A(i,j);
            end
        end
        end
    end
    f = size(B);                               % size of matrix B
     if(f(1) == 2)
         sum = B(1,1)*B(2,2) - B(1,2)*B(2,1);  % for matrix with 4 elements
     else
      sum = mat_det(B);
     end
end