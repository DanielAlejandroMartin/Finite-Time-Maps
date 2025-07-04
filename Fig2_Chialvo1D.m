%% Matlab Code for generating Fig 2 of the article
%% "Finite-time scaling on low-dimensional map bifurcations"
%% Works on MatlabR2018a
%% It follows most naming conventions and definitions from "Fig1_Logistic.m"



clear all
close all
 mup=3-log(3) % \mu Parameter at the bifurcation
im=0  %Counter for \mu values
for mu=  cat(2,[1.88:0.002:1.9],[ 1.9:0.0005:1.92-0.0005], [1.92:0.002:1.94] ) 
    im=im+1
    muvec(im)=mu; %vector of used \mu values
    ic=0 %Counter for initial conditions
    for x0=  cat(2,[3.6:-0.1:3.4],[2.2:0.1:2.4]) % initial conditions
        ic=ic+1;
        x(1)=x0;
        for t=2:100000
            xx=x(t-1);
            x(t)=xx*xx*exp(mu-xx);
        end
        

        
        allx(im,ic,:)=x;
        s1=max(x(end-1:end));
       s2=min(x(end-1:end));
       pp=(s1+s2)/2;
       pm(im,ic)=pp; %Midpoint
       mvec(im)=(2*pp-pp*pp)*exp(mu-pp); %Derivative at the midpoint 
      
         distP(im,ic,:)=x-x(end); if mu>mup;  distP(im,ic,:)=x-pp;end
               %distance to the fixed point (for mu<mup=3-log(3)) or the midpoint midpoint
        %p_m
    end
end



colorVec = hsv(9)
figure(100)
subplot(121)
for tindex=1:6
    time=100*2^(tindex)
    for i=1:ic
       
    plot(muvec,squeeze(allx(:,i,time)),'-o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
hold on
    end     
end

xlabel('mu')
ylabel('x_l')
title('original data')



subplot(122)
for tindex=7:-1:2
     time=100*2^(tindex-1)
     for i=1:ic
   plot(-(mvec+1)*time,(time)^(1/2)*squeeze(distP(:,i,time)),'o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
    
hold on
     end
    
 
end
xlabel('z')
ylabel('sqrt(l)  (x_l-p_m)')
z=[-10:0.1:15]
 kz=z*2
 G=kz.*exp(kz)./(exp(kz)-1)
plot(z,sqrt(3*G), '-.','linewidth',2,'color','k')
plot(z,-sqrt(3*G), '-.','linewidth',2,'color','k')
xlim([-10 10])

