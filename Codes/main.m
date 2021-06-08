clear;

%% Inital conditions
size_k = 120;
h = 0.01;

A = zeros(size_k,1);
B = zeros(size_k,1);
C = zeros(size_k,1);
D = zeros(size_k,1);
E = zeros(size_k,1);
F = zeros(size_k,1);

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


p1 = 0.66;
p2 = 0.66;
p3 = 0.66;
p4 = 0.66;
p5 = 0.66;
p6 = 0.66;

L = 1e5;

N = A(1) + B(1) + C(1) + D(1) + E(1);

j1 = 1;


%% 
% for k = 2:size_k
%     
% %     N = A(k) + B(k) + C(k) + D(k) + E(k);
%     
%     sum_A = 0;
%     for j = v(L,k,h):k
%         sum_A = sum_A + c(p1,j)*A(k-j+1);        
%     end
%     sum1_A = 0;
%     for j1 = v(L,k,h)-1:k-1
%         sum1_A = sum1_A + (mu1 - a1*A(k-j1) - (b1*A(k-j1)*(C(k-j1) + sigma*D(k-j1)))/N )*(h^p1);
%     end
% %     A(k) = sum1_A - summation(A,p1,v(k,L,h),k);
%     A(k) = sum1_A - sum_A;
%         
% %     A(k) = (mu1 - a1*A(k-j1) - (b1*A(k-j1)*(C(k-j1) + sigma*D(k-j1)))/N )*(h^p1) - summation(A,p1,v(k,L,h),k);
% %     A(k) = (mu1 - a1*A(k-j1) - (b1*A(k-j1)*(C(k-j1) + sigma*D(k-j1)))/N )*(h^p1) - sum_A;
%     
%     sum_B = 0;
%     for j = v(L,k,h):k
%         sum_B = sum_B + c(p2,j)* B(k-j+1);        
%     end
%     sum1_B = 0;
%     for j1 = v(L,k,h)-1:k-1
%     sum1_B =  sum1_B + ((b1*A(k)*(C(k-j1) + sigma*D(k-j1)))/N + ...
%         b2 * A(k) * F(k-j1) - (1 - d)* f1 * B(k - j1) - ...
%         d * e1 * B(k-j1) - a1 * B(k-j1))* (h^(p2));
%     end
% %     B(k) = sum1_B - summation(B,p2,v(k,L,h),k);
%     B(k) = sum1_B - sum_B;
% %     B(k) = ((b1*A(k)*(C(k-j1) + sigma*D(k-j1)))/N + ...
% %         b2 * A(k) * F(k-j1) - (1 - d)* f1 * B(k - j1) - ...
% %         d * e1 * B(k-j1) - a1 * B(k-j1))* (h^(p2)) - summation(B,p2,v(k,L,h),k);
% %     B(k) = ((b1*A(k)*(C(k-j1) + sigma*D(k-j1)))/N + ...
% %         b2 * A(k) * F(k-j1) - (1 - d)* f1 * B(k - j1) - ...
% %         d * e1 * B(k-j1) - a1 * B(k-j1))* (h^(p2)) - sum_B;
%     
%     sum_C = 0;
%     for j = v(L,k,h):k
%         sum_C = sum_C + c(p3,j)*C(k-j+1);        
%     end
%     sum1_C = 0;
%     for j1 = v(L,k,h)-1:k-1
%     sum1_C =  sum1_C + ((1-d)*f1*B(k) - (g1 + a1)*C(k-j1))*(h^p3);
%     end
% %     C(k) = sum1_C - summation(C,p3,v(k,L,h),k);
%     C(k) = sum1_C - sum_C;
% 
% %     C(k) = ((1-d)*f1*B(k) - (g1 + a1)*C(k-j1))*(h^p3) - summation(C,p3,v(k,L,h),k);
% %     C(k) = ((1-d)*f1*B(k) - (g1 + a1)*C(k-j1))*(h^p3) - sum_C;
% 
%     sum_D = 0;
%     for j = v(L,k,h):k
%         sum_D = sum_D + c(p4,j)*D(k-j+1);        
%     end
%     sum1_D = 0;
%     for j1 = v(L,k,h)-1:k-1
%     sum1_D =  sum1_D + (d*e1*B(k) - (g2 + a1)*D(k-j1))*(h^p4);
%     end
% %     D(k) = sum1_D - summation(D,p4,v(k,L,h),k);
%     D(k) = sum1_D - sum_D;
% %     D(k) = (d*e1*B(k) - (g2 + a1)*D(k-j1))*(h^p4) - summation(D,p4,v(k,L,h),k); 
% %     D(k) = (d*e1*B(k) - (g2 + a1)*D(k-j1))*(h^p4) - sum_D;
%     
%     
%     sum_E = 0;
%     for j = v(L,k,h):k
%         sum_E = sum_E + c(p5,j)*E(k-j+1);        
%     end
%     sum1_E = 0;
%     for j1 = v(L,k,h)-1:k-1
%     sum1_E =  sum1_E + (g1*C(k) + g2*D(k) - a1*E(k-j1))*(h^p5);
%     end
% %     E(k) = sum1_E - summation(E,p5,v(k,L,h),k);
%     E(k) = sum1_E - sum_E;
% %     E(k) = (g1*C(k) + g2*D(k) - a1*E(k-j1))*(h^p5) - summation(E,p5,v(k,L,h),k);
% %     E(k) = (g1*C(k) + g2*D(k) - a1*E(k-j1))*(h^p5) - sum_E;
%     
% 
%     sum_F = 0;
%     for j = v(L,k,h):k
%         sum_F = sum_F + c(p6,j)*F(k-j+1);        
%     end
%     sum1_F = 0;
%     for j1 = v(L,k,h)-1:k-1
%     sum1_F =  sum1_F + (e2*C(k) + f2*D(k) - mu2*F(k-j1))*(h^p6);
%     end
% %     F(k) = sum1_F  - summation(F,p6,v(k,L,h),k);
%     F(k) = sum1_F - sum_F;
% 
% %     F(k) = (e2*C(k) + f2*D(k) - mu2*F(k-j1))*(h^p6) - summation(F,p6,v(k,L,h),k);
% %     F(k) = (e2*C(k) + f2*D(k) - mu2*F(k-j1))*(h^p6) - sum_F;
% % N = A(k) + B(k) + C(k) + D(k) + E(k);
% end

%%
% for tk = 2:(size_k/h)
%     k = tk*h
%     
%     
%     sum_A = 0;
%     for j = v(L,tk,h):tk
%         j;
%         sum_A = sum_A + c(p1,j)*A(tk-j+1);        
%     end
% %     sum1_A = 0;
% %     for j1 = v(L,k,h)-1:k-1
% %         sum1_A = sum1_A + (mu1 - a1*A(k-j1) - (b1*A(k-j1)*(C(k-j1) + sigma*D(k-j1)))/N )*(h^p1);
% %     end
% %     A(k) = sum1_A - summation(A,p1,v(k,L,h),k);
% %     A(k) = sum1_A - sum_A;
%         
% %     A(k) = (mu1 - a1*A(k-j1) - (b1*A(k-j1)*(C(k-j1) + sigma*D(k-j1)))/N )*(h^p1) - summation(A,p1,v(k,L,h),k);
%     A(tk) = (mu1 - a1*A(tk-j1) - (b1*A(tk-j1)*(C(tk-j1) + sigma*D(tk-j1)))/N )*(h^p1) - sum_A;
%     
%     sum_B = 0;
%     for j = v(L,tk,h):tk
%         sum_B = sum_B + c(p2,j)* B(tk-j+1);        
%     end
% %     sum1_B = 0;
% %     for j1 = v(L,k,h)-1:k-1
% %     sum1_B =  sum1_B + ((b1*A(k)*(C(k-j1) + sigma*D(k-j1)))/N + ...
% %         b2 * A(k) * F(k-j1) - (1 - d)* f1 * B(k - j1) - ...
% %         d * e1 * B(k-j1) - a1 * B(k-j1))* (h^(p2));
% %     end
% %     B(k) = sum1_B - summation(B,p2,v(k,L,h),k);
% %     B(k) = sum1_B - sum_B;
% %     B(k) = ((b1*A(k)*(C(k-j1) + sigma*D(k-j1)))/N + ...
% %         b2 * A(k) * F(k-j1) - (1 - d)* f1 * B(k - j1) - ...
% %         d * e1 * B(k-j1) - a1 * B(k-j1))* (h^(p2)) - summation(B,p2,v(k,L,h),k);
%     B(tk) = ((b1*A(tk)*(C(tk-j1) + sigma*D(tk-j1)))/N + ...
%         b2 * A(tk) * F(tk-j1) - (1 - d)* f1 * B(tk - j1) - ...
%         d * e1 * B(tk-j1) - a1 * B(tk-j1))* (h^(p2)) - sum_B;
%     
%     sum_C = 0;
%     for j = v(L,tk,h):tk
%         sum_C = sum_C + c(p3,j)*C(tk-j+1);        
%     end
% %     sum1_C = 0;
% %     for j1 = v(L,k,h)-1:k-1
% %     sum1_C =  sum1_C + ((1-d)*f1*B(k) - (g1 + a1)*C(k-j1))*(h^p3);
% %     end
% %     C(k) = sum1_C - summation(C,p3,v(k,L,h),k);
% %     C(k) = sum1_C - sum_C;
% 
% %     C(k) = ((1-d)*f1*B(k) - (g1 + a1)*C(k-j1))*(h^p3) - summation(C,p3,v(k,L,h),k);
%     C(tk) = ((1-d)*f1*B(tk) - (g1 + a1)*C(tk-j1))*(h^p3) - sum_C;
% 
%     sum_D = 0;
%     for j = v(L,tk,h):tk
%         sum_D = sum_D + c(p4,j)*D(tk-j+1);        
%     end
% %     sum1_D = 0;
% %     for j1 = v(L,k,h)-1:k-1
% %     sum1_D =  sum1_D + (d*e1*B(k) - (g2 + a1)*D(k-j1))*(h^p4);
% %     end
% %     D(k) = sum1_D - summation(D,p4,v(k,L,h),k);
% %     D(k) = sum1_D - sum_D;
% %     D(k) = (d*e1*B(k) - (g2 + a1)*D(k-j1))*(h^p4) - summation(D,p4,v(k,L,h),k); 
%     D(tk) = (d*e1*B(tk) - (g2 + a1)*D(tk-j1))*(h^p4) - sum_D;
%     
%     
%     sum_E = 0;
%     for j = v(L,tk,h):tk
%         sum_E = sum_E + c(p5,j)*E(tk-j+1);        
%     end
% %     sum1_E = 0;
% %     for j1 = v(L,k,h)-1:k-1
% %     sum1_E =  sum1_E + (g1*C(k) + g2*D(k) - a1*E(k-j1))*(h^p5);
% %     end
% %     E(k) = sum1_E - summation(E,p5,v(k,L,h),k);
% %     E(k) = sum1_E - sum_E;
% %     E(k) = (g1*C(k) + g2*D(k) - a1*E(k-j1))*(h^p5) - summation(E,p5,v(k,L,h),k);
%     E(tk) = (g1*C(tk) + g2*D(tk) - a1*E(tk-j1))*(h^p5) - sum_E;
%     
% 
%     sum_F = 0;
%     for j = v(L,tk,h):tk
%         sum_F = sum_F + c(p6,j)*F(tk-j+1);        
%     end
% %     sum1_F = 0;
% %     for j1 = v(L,k,h)-1:k-1
% %     sum1_F =  sum1_F + (e2*C(k) + f2*D(k) - mu2*F(k-j1))*(h^p6);
% %     end
% %     F(k) = sum1_F  - summation(F,p6,v(k,L,h),k);
% %     F(k) = sum1_F - sum_F;
% 
% %     F(k) = (e2*C(k) + f2*D(k) - mu2*F(k-j1))*(h^p6) - summation(F,p6,v(k,L,h),k);
%     F(tk) = (e2*C(tk) + f2*D(tk) - mu2*F(tk-j1))*(h^p6) - sum_F;
% % N = A(k) + B(k) + C(k) + D(k) + E(k);
% 
%     N = A(tk) + B(tk) + C(tk) + D(tk) + E(tk);
% end




%% 
k1 = 1:size_k;
subplot(231);
plot(k1*h, A(k1));
% axis([0 120 0 9e6]);
title("A(t)");
subplot(232);
plot(k1, B(k1));
% axis([0 120 0 6e6]);
title("B(t)")
subplot(233);
plot(k1, C(k1));
title("C(t)")
% axis([0 120 0 2.5e4]);
subplot(234);
plot(k1, D(k1));
title("D(t)")
% axis([0 120 0 4500]);
subplot(235);
plot(k1, E(k1));
title("E(t)")
% axis([0 120 0 3.5e5]);
subplot(236);
plot(k1, F(k1));
title("F(t)")
% axis([0 120 0 5e4]);


function sum = summation(func,p,v,k)

    sum = 0;
    for j = v:k
        sum = sum + c(p,j)*func(k-j+1);        
    end
    
    
end

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

