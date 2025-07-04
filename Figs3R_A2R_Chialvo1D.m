%% Matlab Code for generating the right column of plots in  Figs 3 and A2 of the article
%% "Finite-time scaling on low-dimensional map bifurcations"
%% Works on MatlabR2018a
%% It follows most naming conventions and definitions from "Fig1_Logistic.m" and 
%% "Figs3L_A2L_Logistic.m"
close all
clear all
% Parameters
mu_start = 3-log(3)-0.002;
mu_end = 3-log(3)+0.002;

num_mu = 5000; % Number of mu values to test
muvec = linspace(mu_start, mu_end, num_mu);
mup=3-log(3);; %exact value of mu*
disp('Figure 101 shows results for the right columns in Fig 3 of the paper')
disp('Figure 201 shows results for the right columns in Fig A2 of the paper')
transient=1;

 step=0;
 tic
 pavec=10:1:18;
for pa=pavec
    step=step+1;
 iterations=2^pa;
 % Initialize result arrays
lyapunov_exponents = zeros(size(muvec));
susceptibilities = zeros(size(muvec));

for i = 1:length(muvec)
    mu = muvec(i);
    x = 0.8;  
     % Initial condition (in basin of attraction)
    dx_dmu = 0;     % Initial derivative of x with respect to mu
    
    % Discard transient iterations
    for t = 1:transient
       x = (x^2) * exp(mu - x);
       
    end
%     
    % Initialize accumulators
    lyapunov_sum = 0;
    chi_sum = 0;
    
    % Main calculation loop
    for t = 1:iterations
             % Current map derivative (for Lyapunov exponent)
        df_dx = (2*x - x^2) * exp(mu - x);
        lyapunov_sum = lyapunov_sum + log(abs(df_dx));
        
        % Accumulate susceptibility
        chi_sum = chi_sum + dx_dmu^2;
        
        % Update for next iteration
      x_new = (x^2) * exp(mu - x);
        dx_dmu_new = (2*x*exp(mu-x) - x^2*exp(mu-x)) * dx_dmu + x^2*exp(mu-x);
        
        % Update variables
        x = x_new;
        dx_dmu = dx_dmu_new;
    
    end
    
    % Store results
    lyapunov_exponents(i) = lyapunov_sum / iterations;
    susceptibilities(i) = chi_sum;
end  
 
 
  AllLya(step,:)=lyapunov_exponents; %Save all Lyapunov Exponents
  AllSus(step,:)=susceptibilities;  %Save all Susceptibility values
 
    

  L(step)=iterations; %Number of map iterations
 
[maxL,ind1]=max(lyapunov_exponents);
mupeak(step)=muvec(ind1);  %\mu at the Peak of lyapunov exponent for a given number of iterations
lyappeak(step)=abs(maxL);  % Peak of lyapunov exponent for a given number of iterations


[maxS,ind2]=max(susceptibilities);
mupeak2(step)=muvec(ind2); %\mu at the Peak of susceptibility for a given number of iterations
Susppeak(step)=(maxS); %Peak of susceptibility for a given number of iterations
 



% Plot the results
figure(100);

subplot(311)
semilogy(muvec, susceptibilities, 'g-', 'LineWidth', 1.5);
xlabel('\mu');
ylabel('\chi ');
 title('\chi for  Chialvo1D map');
grid on;
hold on;
% Highlight maximun Suscept
 plot(muvec(ind2), susceptibilities(ind2),'or')
 
 
 
subplot(312)
plot(muvec, lyapunov_exponents, 'b-', 'LineWidth', 1.5);
xlabel('\mu');
ylabel('\lambda Exponent');
title('\lambda for  Chialvo1D map');
grid on;
% Highlight maximun exponents  
hold on;
plot(muvec(ind1),lyapunov_exponents(ind1),'or')
 




 
end





subplot(313)
loglog(L,Susppeak)
hold on
plot([100:100:100000]*10,[100:100:100000].^(2),'--')
loglog(L,1./lyappeak)
plot([100:100:100000]*10,[100:100:100000]./log([100:100:100000]),'--')
legend('max \chi', '', 'max 1/|\lambda|')
xlabel('l: number of iterations')
ylabel('max \chi and 1/|\lambda| ')


%%%%%%%%%%%%%%%
figure(200)
subplot(211)
loglog(L,mupeak2-mup)
xlabel('l: number of iterations')
ylabel('\mu_l*-\mu*')
title('distance to mu*')
hold on
plot([100:100:100000],[100:100:100000].^(-1))

subplot(212)
for i=1:size(pavec,2)
    hold on
 plot((muvec-mup)*(2^pavec(i))^(1),AllSus(i,:)*(2^pavec(i))^(-2),'.')   
 
end
xlabel('(\mu-\mu*) l^{1/\nu}')
ylabel('\chi_l l^{\gamma/\nu}')
%gamma=2, nu=1
title('Susceptibility Collapse')
xlim([-5 10])
 