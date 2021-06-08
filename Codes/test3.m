tic
%% Inital conditions
size_k = 1.21;
global h;
h = 0.01;


A = zeros(size_k/h,1);
B = zeros(size_k/h,1);
C = zeros(size_k/h,1);
D = zeros(size_k/h,1);
E = zeros(size_k/h,1);
F = zeros(size_k/h,1);


global mu1 mu2 a1 b1 b2 sigma d f1 f2 e1 e2 g1 g2 N

mu1 = 107644.22451;
mu2 = 0.01;
a1 = 1/76.79;
b1 = 0.05;
b2 = 0.000001231;
sigma = 0.02;
d = 0.1243;
f1 = 0.00047876;
f2 = 0.001;
e1 = 0.005;
e2 = 0.000398;
g1 = 0.09871;
g2 = 0.854302;


A(1) = 8065518;
B(1) = 200000;
C(1) = 282;
D(1) = 200;
E(1) = 0;
F(1) = 50000;
% p1 = 0.01;
% p2 = 0.01;
% p3 = 0.01;
% p4 = 0.01;
% p5 = 0.01;
% p6 = 0.01;

global Avec Bvec Cvec Dvec Evec Fvec

Avec = zeros(size_k/h,1);
Bvec = zeros(size_k/h,1);
Cvec = zeros(size_k/h,1);
Dvec = zeros(size_k/h,1);
Evec = zeros(size_k/h,1);
Fvec = zeros(size_k/h,1);
Avec(1) = 8065518;
Bvec(1) = 200000;
Cvec(1) = 282;
Dvec(1) = 200;
Evec(1) = 0;
Fvec(1) = 50000;


p1 = 0.99;
p2 = 0.99;
p3 = 0.99;
p4 = 0.99;
p5 = 0.99;
p6 = 0.99;

L = 1e5;

N = A(1) + B(1) + C(1) + D(1) + E(1);

j1 = 1;

for k2 = 0:h:size_k-h
    k2
    funcA(k2);
    funcB(k2);
    funcC(k2);
    funcD(k2);
    funcE(k2);
    funcF(k2);
end

%% 
for k = 2:(size_k/h)-1
    
%     N = A(k) + B(k) + C(k) + D(k) + E(k);
    
    sum1_A = 0;
    j1 = 0;
    while (k-j1)>0
        if j1==0
            frac = 1;
        else
            frac = p1/j1;
        end
        sum1_A = sum1_A + (1/(h^p1))*((1)^j1)*frac*(Avec(k-j1));
        j1 = j1+1;
    end
    sum_A = 0;
    for j = v(L,k,h):k
        sum_A = sum_A + c(p1,j)*A(k-j+1);        
    end
%     A(k) = sum1_A - summation(A,p1,v(k,L,h),k);
    A(k) = sum1_A*(h^p1) - sum_A;
        
%     A(k) = (mu1 - a1*A(k-j1) - (b1*A(k-j1)*(C(k-j1) + sigma*D(k-j1)))/N )*(h^p1) - summation(A,p1,v(k,L,h),k);
%     A(k) = (mu1 - a1*A(k-j1) - (b1*A(k-j1)*(C(k-j1) + sigma*D(k-j1)))/N )*(h^p1) - sum_A;
    
    sum_B = 0;
    for j = v(L,k,h):k
        sum_B = sum_B + c(p2,j)* B(k-j+1);        
    end
    sum1_B = 0;
    j1 = 0;
    lambda = 2.0195e5 + B(1);
    while (k-j1)>0
        if j1==0
            frac = 1;
        else
            frac = p2/j1;
        end
        sum1_B = sum1_B + (1/(h^p2))*((-1))*(frac)*(Bvec(k-j1));
        j1 = j1+1;
    end
%     B(k) = sum1_B - summation(B,p2,v(k,L,h),k);
    B(k) = sum1_B*(h^p2) - sum_B + lambda;
%     B(k) = ((b1*A(k)*(C(k-j1) + sigma*D(k-j1)))/N + ...
%         b2 * A(k) * F(k-j1) - (1 - d)* f1 * B(k - j1) - ...
%         d * e1 * B(k-j1) - a1 * B(k-j1))* (h^(p2)) - summation(B,p2,v(k,L,h),k);
%     B(k) = ((b1*A(k)*(C(k-j1) + sigma*D(k-j1)))/N + ...
%         b2 * A(k) * F(k-j1) - (1 - d)* f1 * B(k - j1) - ...
%         d * e1 * B(k-j1) - a1 * B(k-j1))* (h^(p2)) - sum_B;
    
    sum_C = 0;
    for j = v(L,k,h):k
        sum_C = sum_C + c(p3,j)*C(k-j+1);        
    end
    sum1_C = 0;
    j1 = 0;
    alpha = 100;
    lambda = 2.7969e4 + C(1);
    while (k-j1)>0
        if j1==0
            frac = 1;
        else
            frac = p3/j1;
        end
        sum1_C = sum1_C + (1/(h^p3))*((-1))*(frac)*(Cvec(k-j1))*alpha;
        j1 = j1+1;
    end
%     C(k) = sum1_C - summation(C,p3,v(k,L,h),k);
    C(k) = sum1_C*(h^p3) - sum_C + lambda;

%     C(k) = ((1-d)*f1*B(k) - (g1 + a1)*C(k-j1))*(h^p3) - summation(C,p3,v(k,L,h),k);
%     C(k) = ((1-d)*f1*B(k) - (g1 + a1)*C(k-j1))*(h^p3) - sum_C;

    sum_D = 0;
    for j = v(L,k,h):k
        sum_D = sum_D + c(p4,j)*D(k-j+1);        
    end
    sum1_D = 0;
    j1 = 0;
    alpha = 10;
    lambda = 1.9741e3 + D(1);
    while (k-j1)>0
        if j1==0
            frac = 1;
        else
            frac = p4/j1;
        end
        sum1_D = sum1_D + (1/(h^p4))*((-1))*(frac)*(Dvec(k-j1))*alpha;
        j1 = j1+1;
    end
%     D(k) = sum1_D - summation(D,p4,v(k,L,h),k);
    D(k) = sum1_D*(h^p4) - sum_D + lambda;
%     D(k) = (d*e1*B(k) - (g2 + a1)*D(k-j1))*(h^p4) - summation(D,p4,v(k,L,h),k); 
%     D(k) = (d*e1*B(k) - (g2 + a1)*D(k-j1))*(h^p4) - sum_D;
    
    
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
%     E(k) = sum1_E - summation(E,p5,v(k,L,h),k);
    E(k) = sum1_E*(h^p5) - sum_E;
%     E(k) = (g1*C(k) + g2*D(k) - a1*E(k-j1))*(h^p5) - summation(E,p5,v(k,L,h),k);
%     E(k) = (g1*C(k) + g2*D(k) - a1*E(k-j1))*(h^p5) - sum_E;
    

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
%     F(k) = sum1_F  - summation(F,p6,v(k,L,h),k);
    F(k) = sum1_F*(h^p6) - sum_F;

%     F(k) = (e2*C(k) + f2*D(k) - mu2*F(k-j1))*(h^p6) - summation(F,p6,v(k,L,h),k);
%     F(k) = (e2*C(k) + f2*D(k) - mu2*F(k-j1))*(h^p6) - sum_F;
% N = A(k) + B(k) + C(k) + D(k) + E(k);
end
%%
% fprintf("%.8f\n",funcA(0.02))
% funcA(0.02)

% for k2 = 0:0.01:1
%     k2
%     funcA(k2);
% end
% k2 = 1:101;
% plot(k2*h, Avec(k2));

%% Plot
% B(2) = -2.0195e5
% C(2) = -2.7969e4
% D(2) = -1.9741e3
% 
% B(2:end) = B(2:end) - B(2)+ B(1);
% C(2:end) = C(2:end) - C(2)+ C(1);
% D(2:end) = D(2:end) - D(2)+ D(1);

k1 = 0:(size_k/h)-1;
figure;
subplot(231);
plot(k1, A(k1+1));
axis([0 120 0 9e6]);
title("A(t)");
subplot(232);
plot(k1, B(k1+1));
axis([0 120 0 5e5]);
title("B(t)")
subplot(233);
plot(k1, C(k1+1));
title("C(t)")
axis([0 120 0 4e4]);
subplot(234);
plot(k1, D(k1+1));
title("D(t)")
axis([0 120 0 5500]);
subplot(235);
plot(k1, E(k1+1));
title("E(t)")
axis([0 120 0 2e4]);
subplot(236);
plot(k1, F(k1+1));
title("F(t)")
axis([0 120 0 5e4]);
timetaken = toc
%% Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function answer = funcA(t)
    global h mu1 mu2 a1 b1 b2 sigma d f1 f2 e1 e2 g1 g2 N
    global Avec Bvec Cvec Dvec Evec Fvec
    
    if Avec(floor((t/h)+1)) ~= 0
        answer = Avec(floor((t/h)+1));
    elseif t == 0
        answer = 8065518;
        Avec(floor((t/h)+1)) = answer;
    else
        dt = h;
        answer = 0;
        it = 0;
        while (it+dt)<t
            answer = answer + (((mu1 - a1*funcA(it) - (b1*funcA(it)*(funcC(it) + (sigma*funcD(it))))/N) ...
                            + (mu1 - a1*funcA(it+dt) - (b1*funcA(it+dt)*(funcC(it+dt) + (sigma*funcD(it+dt))))/N ))/2)*dt;
            it = it + dt;
        end
        Avec(floor((t/h)+1)) = answer;
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function answer = funcB(t)
    global h mu1 mu2 a1 b1 b2 sigma d f1 f2 e1 e2 g1 g2 N
    global Avec Bvec Cvec Dvec Evec Fvec
    
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
    global h mu1 mu2 a1 b1 b2 sigma d f1 f2 e1 e2 g1 g2 N
    global Avec Bvec Cvec Dvec Evec Fvec
    
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
    global h mu1 mu2 a1 b1 b2 sigma d f1 f2 e1 e2 g1 g2 N
    global Avec Bvec Cvec Dvec Evec Fvec
    
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
    global h mu1 mu2 a1 b1 b2 sigma d f1 f2 e1 e2 g1 g2 N
    global Avec Bvec Cvec Dvec Evec Fvec
    
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
    global h mu1 mu2 a1 b1 b2 sigma d f1 f2 e1 e2 g1 g2 N
    global Avec Bvec Cvec Dvec Evec Fvec
    
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = c(p, j)
    output = 1;
    for i = 1:j
        output = (1 - ((1+p)/i))* output;
    end
end


function output = v(L,k,h)

    if (k <= (L/h))
        output = 1;
    else
        output = k - (L/h);
    end
end

%%%%%%%% TODO: Replace inbuilt factorial function
function answer = nCk(n, k)
    n_fact = factorial(n);
    k_fact = factorial(k);
    nk_fact = factorial(n-k);
    
    answer = n_fact/(k_fact*nk_fact);
end

