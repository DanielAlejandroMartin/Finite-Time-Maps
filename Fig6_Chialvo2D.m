%% Matlab Code for generating Fig 6 of the article
%% "Finite-time scaling on low-dimensional map bifurcations"
%% Works on MatlabR2018a
%% It follows most naming conventions and definitions from "Fig1_Logistic.m" and 
%% "Figs5_A3_Hopf.m"

clear all
close all

a=0.89; b=0.6; c=0.28 %Model Parameters
kvec = [0.0290: 0.000025: 0.030]  %values of parameter k.
icmax=30; %Number of initial conditions
%% First find the fixed points
kindex=0
for k=kvec
    kindex=kindex+1
    x(1)=0.05;  y(1)=0.05; %arbitrary initial condition
    if kindex>1; x(1)=xf(kindex-1); y(1)=yf(kindex-1); end
    for l=2:100000 %extremely long iteration
            y0=y(l-1);
            x0=x(l-1);
            x(l)=x0^2*exp(y0-x0)+k;
            y(l)=a*y0-b*x0+c;      
    end
        %save final position (our estimate of fixed point)
    xf(kindex)=x(end);
    yf(kindex)=y(end); 
    % Find Eigenvalues at the fixed point
    xx=xf(kindex)
    yy=yf(kindex)
            DF=[(2*xx-xx^2)*exp(yy-xx) xx^2*exp(yy-xx) ;-b a];
        [~,l]=eig(DF);
        lam1(kindex)=l(1,1);
        lam2(kindex)=l(2,2);
        allC(kindex)=sqrt(real(lam1(kindex)*lam2(kindex))); %square root of the product of eigenvalues
end
%plot(xf,yf,'o'); hold on; plot (xf(1:end-4),yf(1:end-4),'o') %Verify whether the last  points  converge to fixed point


%% RUN
kindex=0
for k=kvec
    kindex=kindex+1
    
    for ic=1:icmax %index of initial conditions
    %%Choose Initial Conditions on the basin of attraction and limited to a
    %%single 'curve'
        if ic==1
            x(1)=xf(kindex)+0.001;  y(1)=yf(kindex); %put a small perturbation on x so that it is still on the basin of atraction.
        else %Garantee that all initial conditions belong to the same curve. 
            x(1)=(allx(ic-1,kindex,5))/allC(kindex)^(4.)+xf(kindex);
            y(1)=(ally(ic-1,kindex,5))/allC(kindex)^(4.)+yf(kindex); 
            
        end
        

        
        %% Run
         for l=2:600 
            y0=y(l-1);
            x0=x(l-1);
            x(l)=x0^2*exp(y0-x0)+k;
            y(l)=a*y0-b*x0+c;      
    end%dynamics
        allx(ic,kindex,:)=x-xf(kindex); %save timeseries relative to the fixed point (Move to the center of mass)
        ally(ic,kindex,:)=y-yf(kindex);
        
    end %ic
    
    
    
    
end %k


%    Verify that initial conditions are comparable for all k values
   figure(2)
   hold on
for kindex=1:size(kvec,2)
   plot(allx(:,kindex,1),ally(:,kindex,1),'-o') 
end


%% Plot

figure(200)
colorVec = hsv(6)
subplot(121)
for tindex=1:6
time=fix(25*sqrt(2)^(tindex+1))
    for kindex=1:3:35
 
   plot3( ones(1,icmax)*kvec(kindex), squeeze(allx(:,kindex,time)),squeeze(ally(:,kindex,time)), '-o','linewidth',0.25,'Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 2 )
   hold on
end

end
xlabel({'$k $'},'Interpreter','latex')
ylabel({'$ x-p_x $'},'Interpreter','latex')
zlabel({'$ y-p_y $'},'Interpreter','latex')
grid on
view(31,10)



subplot(122)
for tindex=1:6
time=fix(25*sqrt(2)^(tindex+1))
    for kindex=1:3:35
 
 
   plot3( ones(1,icmax)*((squeeze(allC(kindex)))-1)*time, sqrt(time)*squeeze(allx(:,kindex,time)),sqrt(time)*squeeze(ally(:,kindex,time)), '-o','linewidth',0.25,'Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 2 )
   hold on
end

end
xlabel({'$z=l [\sqrt{\sigma_1 \sigma_2} -1 ] $'},'Interpreter','latex')
ylabel({'$\sqrt{l} (x-p_x) $'},'Interpreter','latex')
zlabel({'$\sqrt{l} (y-p_y) $'},'Interpreter','latex')
xlim([-5 -0])

grid on

view(31,10)