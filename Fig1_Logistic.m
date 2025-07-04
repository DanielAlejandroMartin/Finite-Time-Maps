%% Matlab Code for generating Fig 1 of the article
%% "Finite-time scaling on low-dimensional map bifurcations"
%% Works on MatlabR2018a


clear all
close all
im=0 %Counter for \mu values
for mu=3-0.03:0.001:3+0.03; % \mu values
    im=im+1
    muvec(im)=mu; %vector of used \mu values
    ic=0;%Counter for initial conditions
    for x0=   cat(2,[0.8:-0.02:0.7], [0.4:0.02:0.5]); %Initial conditions.
        ic=ic+1; 
        x(1)=x0;
        for t=2:7000 %Iterate the map
            xx=x(t-1);
            x(t)=mu*xx*(1-xx);
        end
        
    
        
        allx(im,ic,:)=x; % Matrix of all timeseries for all initial conditions an \mu values.
        s1=(mu+1+sqrt((mu-3)*(mu+1)))/2/mu; % solutions for \mu>3
        s2=(mu+1-sqrt((mu-3)*(mu+1)))/2/mu;
        distFPS(im,ic,:)=abs(x-(mu-1)/mu); if mu>3;  distFPS(im,ic,:)=min(abs(x-s1),abs(x-s2));end
        %distance to the closest solution (for mu<3 there is only one
        %solution)
        
        distP(im,ic,:)=x-(mu-1)/mu; if mu>3;  distP(im,ic,:)=x-(mu+1)/2/mu;end
        %distance to the fixed point (for mu<3) or the midpoint midpoint
        %p_m (for mu>3)
    end
end


figure(100)
subplot(221)
colorVec = hsv(9)
for tindex=1:6
    time=100*2^(tindex)
    for i=[1 ic]    
    plot(muvec,squeeze(allx(:,i,time)),'-o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
hold on
    end
     
end


%% Export variables as .txt
%time=100*2.^(1:4);
%TT=table([muvec' reshape(squeeze(allx(:,[1 ic] ,time)), [61 8])])
%writetable(TT,sprintf('Fig1A.txt'),'WriteVariableNames',false,'WriteRowNames',false,'Delimiter', ' ')

xlabel('mu')
ylabel('x_l')
title('original data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(222)
 for tindex=6:-1:1
      time=100*2^(tindex)
      for i=[1 ic]
     plot(muvec,squeeze(distFPS(:,i,time)),'o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )  
 hold on
      end
 end
 xlabel('mu')
ylabel('|x_l-p|')
title('Distance to the closest fixed point')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subplot(223)
for tindex=6:-1:1
     time=100*2^(tindex)
     for i=[1 ic]
   plot((muvec-3)*time,(time)^(1/2)*squeeze(distP(:,i,time)),'o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
    
hold on
     end
 

end

z=[-100:5:100]
 kz=z*2
 G=kz.*exp(kz)./(exp(kz)-1)
plot(z,sqrt(G/18), '-.','linewidth',1,'color','k')
plot(z,-sqrt(G/18), '-.','linewidth',1,'color','k')


xlabel('z')
ylabel('sqrt(l) (x_l-p_m)')

title({'Collapse to  fixed point (mu<3) ',' or  midpoint p_m (for mu>3)'})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subplot(224)
 for tindex=6:-1:1
      time=100*2^(tindex)
      for i =[1 ic]
     plot((muvec-3)*time,(time)^(1/2)*squeeze(distFPS(:,i,time)),'o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
   
 hold on
      end
      
 end
 xlabel('z')
 ylabel('|x-p_m|')
 %Now we add Function G
 z=[-10.05:0.001:0]
 kz=z*2
 G=kz.*exp(kz)./(exp(kz)-1)
plot(z,sqrt(G/18), '-.','linewidth',1,'color','k')
xlim([-5 5])


z=[0:0.001:10]
 kz=z*2
 G=kz.*exp(kz)./(exp(kz)-1)
plot(z,sqrt(G/18)-sqrt(z)/3, '-.','linewidth',1,'color','b')
xlim([-5 5])

 title('Collapse to closest solution')

 
 

