%% Matlab Code for generating Fig A1 (Appendix) of the article
%% "Finite-time scaling on low-dimensional map bifurcations"
%% Works on MatlabR2018a
%% It follows most naming conventions and definitions from "Fig1_Logistic.m"

clear all
close all
im=0
%% 
%This code works for any value of Kappa, just changing the line below Take
%Kappa =3 and Kappa=8 to reproduce Figure A1
Kappa=2; 
%%
odd=Kappa-2*fix(Kappa/2); %If Kappa is odd, odd=1, else odd=0
k=Kappa
if odd==0
    k=2*Kappa-1
end

for mu=-0.03:0.0005:0.03 %mu values
    im=im+1;
    muvec(im)=mu; %vector of mu values
   
    ic=0
    for x0=cat(2,[0.7:-0.05:0.6],[-0.6:-0.05:-0.7]) %Initial conditions
        ic=ic+1;
        x(1)=x0;
        for l=2:7000 %dynamics
            xx=x(l-1);
            x(l)=(-1-mu)*xx+xx^Kappa;
        end
        
        
        allx(im,ic,:)=x;
    
        
    end
    p_m(im)=(x(end)+x(end-1))/2; %Middpoint
end
%% Plot Resutls

colorVec = hsv(9)
figure(100)
subplot(211)
for tindex=1:7
    time=100*2^(tindex-1)
    for i=1:ic
       
    plot(muvec,squeeze(allx(:,i,time)),'-o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
hold on
    end
    
    
end

xlabel('mu')
ylabel('x_l')
title('original data')
 
 
subplot(212)
for tindex=1:7
     time=100*2^(tindex-1)
     for i=1:ic
   plot((muvec)*time,(time)^(1/(k-1))*(squeeze(allx(:,i,time))-p_m'),'o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
    
  
hold on
     end
     
end
xlabel('mu l')
ylabel(' |x_l-p_m| l^{(1/k-1)}')


    z=[-20:0.05:20]
 kz=z*(k-1)
 G=kz.*exp(kz)./(exp(kz)-1)
if odd==1   
plot(z,(G/(k-1)).^(1/(k-1)), '-.','linewidth',2,'color','k')
plot(z,-(G/(k-1)).^(1/(k-1)), '-.','linewidth',2,'color','k')
end


  
