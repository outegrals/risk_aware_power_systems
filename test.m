clc; clear;

syms a2 real
syms beta1 beta2 beta3 real
syms x1o x2o x3o real  % forecasts x_i^o

% gammas
gamma12 = x2o/x1o; gamma13 = x3o/x1o;
gamma21 = x1o/x2o; gamma23 = x3o/x2o;
gamma31 = x1o/x3o; gamma32 = x2o/x3o;

% Build A and B per Eq. (11) specialized to n=3
A = sym(zeros(3,3));
B = sym(zeros(3,1));

betas  = [beta1; beta2; beta3];
xos    = [x1o; x2o; x3o];
G = [  0      gamma12 gamma13;
     gamma21    0     gamma23;
     gamma31  gamma32   0   ];

for i = 1:3
    A(i,i) = 2*(1 + 2*a2*betas(i));
    for j = 1:3
        if j ~= i
            A(i,j) = (1 + 4*a2*betas(i)) * G(i,j);
        end
    end
    B(i) = (1 + 4*a2*betas(i)) * sum([1, G(i,1:3)]);
end

% Solve for alpha*
alpha_star = simplify(A \ B);   % exact symbolic solution

% alpha_star = simplify(A\B);  % from your n=3 symbolic solve
alpha1_star = alpha_star(1);
alpha2_star = alpha_star(2);
alpha3_star = alpha_star(3);

syms gamma12 gamma13 gamma21 gamma23 gamma31 gamma32 real

% -------- Step 1: replace "sum ratios" (Type A patterns) --------
% These correspond to: (x1o+x2o+x3o)/xio = 1 + xjo/xio + xko/xio
S1 = (x1o + x2o + x3o)/x1o;
S2 = (x1o + x2o + x3o)/x2o;
S3 = (x1o + x2o + x3o)/x3o;

alpha1_g = subs(alpha1_star, S1, 1 + gamma12 + gamma13);
alpha2_g = subs(alpha2_star, S2, 1 + gamma21 + gamma23);
alpha3_g = subs(alpha3_star, S3, 1 + gamma31 + gamma32);

% -------- Step 2: replace pairwise ratios with gammas --------
% Define all ratio->gamma substitutions
ratio_list  = [ x2o/x1o, x3o/x1o, x1o/x2o, x3o/x2o, x1o/x3o, x2o/x3o ];
gamma_list  = [ gamma12, gamma13, gamma21, gamma23, gamma31, gamma32 ];

alpha1_g = subs(alpha1_g, ratio_list, gamma_list);
alpha2_g = subs(alpha2_g, ratio_list, gamma_list);
alpha3_g = subs(alpha3_g, ratio_list, gamma_list);

% -------- Step 3: simplify / factor / collect --------
alpha1_g = factor(simplify(alpha1_g));
alpha2_g = factor(simplify(alpha2_g));
alpha3_g = factor(simplify(alpha3_g));

alpha1_g
alpha2_g
alpha3_g
