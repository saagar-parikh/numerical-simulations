
syms J1 J2 lambda
syms func1(lambda) func2(lambda)

global func1 func2

I = [1, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1];

J1 = [ -0.0130 0 -0.0500 -0.0010 0 -10.1754 ;
       0 -0.0141 0.0500 0.0010 0 10.1754 ;
       0 0.0004 -0.1117 0 0 0  ; 
       0 0.0006 0 -0.8673 0 0 ; 
       0 0 0.0987 0.8543 -0.0130 0 ;
       0 0 0.0003 0.0010 0 -0.0100 ];
 
J2 = [ -0.0023 0 -0.2885 -0.0058 0 -58.7180 ;
       -0.0108 -0.0141 0.2885 0.0058 0 58.7180 ; 
       0 0.0004 -0.1117 0 0 0 ;
       0 0.0006 0 -0.8673 0 0 ;
       0 0 0.0987 0.8543 -0.0130 0 ;
       0 0 0.0003 0.0010 0 -0.0100 ];
   
   
func1 = @(lambda) mat_det(J1 - lambda*I);
func2 = @(lambda) mat_det(J2 - lambda*I);


estimates1 = [-0.0140 -0.0120 -0.1110 -0.0060 -0.0170 -0.8670];

for i = 1:6
   root1(i) = newtonraphson1(estimates1(i)); 
end

estimates2 = [-0.0140 -0.8670 -0.1110 -0.0200 -0.0125 -0.0061];

for i = 1:6
   root2(i) = newtonraphson2(estimates2(i)); 
end



fplot(lambda, func2)
line(xlim, [0,0], 'Color', 'Black', 'LineWidth', 0.5);


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
    if(m ~= n)                                
        disp('Determinant cannot be found');  % determinent is defined only for square mattrices 
    elseif(n == 1)                            % matrix with only one element
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