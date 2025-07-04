%% Matlab Code for generating  Figs 4 of the article
%% "Finite-time scaling on low-dimensional map bifurcations"
%% Works on MatlabR2018a
%% It follows most naming conventions and definitions from "Fig1_Logistic.m"


      
      
      
clear all
close all
im=0  %Counter for \mu values
for mu=-0.25:0.0000025:-0.24975  % \mu values
    im=im+1
    muvec(im)=mu;  % vector of \mu values

 %stable anu unstalbe solutions
     s1(im)= sqrt(1/2*(1+sqrt(1+4*mu)));
       s2(im)= sqrt(1/2*(1-sqrt(1+4*mu)));
       %derivative on stable solution
deriv(im)=1+mu+3* s1(im)^2-5* s1(im).^4;


    ic=0; %Counter for Initial conditions
   for x0=  cat(2, [1/sqrt(2)],[0.8 :0.05:0.95]) % Vector of initial conditions
        ic=ic+1;
        x(1)=x0;
        for t=2:7000 %Iterate the map
            xx=x(t-1);
            x(t)=(1+mu)*xx+xx^3-xx^5;
        end
             
        allx(im,ic,:)=x; %Save alll iterations
        distFPS(im,ic,:)=x- s1(im); %Save the distance to the atttractive solution
      
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colorVec = hsv(9)
figure(100)
subplot(221)
time=1
   i=1  % Show 1 initial condition only for clarity reasons
    plot(muvec,squeeze(allx(:,i,time)),'-o','Color', 'k', 'MarkerSize', 3 )
hold on

for tindex=1:7
    time=10*2^(tindex-1);
    for i=1:1 % Show 1 initical condition only
       
    plot(muvec,squeeze(allx(:,i,time)),'-o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
hold on
    end
       
end
plot(muvec,s1, '--','linewidth', 2, 'color','k')
plot(muvec,s2, '--','linewidth', 2, 'color','k')
xlabel('mu')
ylabel('x_l')
legend('l=0','l=10','l=20','l=40','l=80','l=160','l=320')
title('original data, x(0)=1/sqrt(2)')
%%%%%%%%%%%%%%%%%%%%%%%%
 
subplot(222)
time=1
   i=2  % Show 1 initial condition only for clarity reasons
    plot(muvec,squeeze(allx(:,i,time)),'-o','Color', 'k', 'MarkerSize', 3 )
hold on

for tindex=1:7
    time=10*2^(tindex-1);
    for i=2:2 % Show 1 initical condition only
       
    plot(muvec,squeeze(allx(:,i,time)),'-o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
hold on
    end
       
end
plot(muvec,s1, '--','linewidth', 2, 'color','k')
plot(muvec,s2, '--','linewidth', 2, 'color','k')
xlabel('mu')
ylabel('x_l')
legend('l=0','l=10','l=20','l=40','l=80','l=160','l=320')
title('original data, x(0)=0.8')
%%%%%%%%%%%%%%%%%%%%%%%%



subplot(223)
 for tindex=1:7 %-1:1
      time=100*2^(tindex-1)
      for i=1:1
     plot(-(deriv-1)*time,(time)^(1)*squeeze(distFPS(:,i,time)),'o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
   
 hold on
      end
     
 end
 
  xlabel('z')
 ylabel('|x_l-p| l')
 title('Collapse for x(0)=1/sqrt(2)')
 
z=[-10.05:0.001:0]
 kz=z*1
 G=kz.*exp(kz)./(exp(kz)-1)
plot(-z,-G*2/2.828, '-.','linewidth',1,'color','k')
plot(-z,G*2/2.828, '-.','linewidth',1,'color','k')
xlim([0 10])



 %%%%%%%%%%%%%%%%%
subplot(224)
 for tindex=1:7 %-1:1
      time=100*2^(tindex-1)
      for i=2:2
     plot(-(deriv-1)*time,(time)^(1)*squeeze(distFPS(:,i,time)),'o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
   
 hold on
      end
     
 end
 xlabel('z')
 ylabel('|x_l-p| l')
 title('Collapse for x(0)=0.8')
 
z=[-10.05:0.001:0]
 kz=z*1
 G=kz.*exp(kz)./(exp(kz)-1)
plot(-z,-G*2/2.828, '-.','linewidth',1,'color','k')
plot(-z,G*2/2.828, '-.','linewidth',1,'color','k')
xlim([0 10])




 
 
