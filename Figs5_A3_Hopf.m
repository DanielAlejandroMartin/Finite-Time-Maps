%% Matlab Code for generating Fig 5 of the article
%% "Finite-time scaling on low-dimensional map bifurcations"
%% Works on MatlabR2018a
%% It follows most naming conventions and definitions from "Fig1_Logistic.m"

clear all
close all
%Equation parameters
dt=0.1
w=0.01
disp('Figure 100 shows results for Fig 4 of the paper')
disp('Figure 200 shows results  Fig A3 of the paper')

im=0 % counter for mu vaules
for mu=-0.03:0.0005:0.03
    im=im+1
    muvec(im)=mu; %vector of mu values
    %% Compute the product of the eigenvalues
    J=zeros(2); J(1,1)=1+dt*mu;  J(1,2)=-dt*w;  J(2,1)=dt*w;  J(2,2)=1+dt*mu;
  lam=eig(J); %Eigenvalues  of the Jacobian Matrix
  newm(im)=sqrt(lam(1)*lam(2)); 

  newm2(im)=newm(im);
 if mu>0  %Eigenvalues  of the Jacobian Matrix for r>0 (i.e., for mu>0)
        theta=0 %the actual value of t is irrelevant
    K=zeros(2); K(1,1)=1-2*dt*mu*cos(theta)^2;  K(1,2)=-dt*(w+2*mu*cos(theta)*sin(theta));  K(2,1)=-dt*(-w+2*mu*cos(theta)*sin(theta));  K(2,2)=1-2*dt*mu*sin(theta)^2;
  lam=eig(K)
  newm2(im)=sqrt(lam(1)*lam(2));
  end
  
    ic=0  % counter for initial conditions
    for theta=0:0.2:2*pi+0.05 %initial condition angles (initial radius is equal to 1)
        ic=ic+1;
        x(1)=cos(theta)      ;
        y(1)=sin(theta);
        for l=2:7000 %Iterate the map
            xx=x(l-1);
            yy=y(l-1);
            rr2=xx^2+yy^2;
            x(l)=xx+dt*((mu-rr2)*xx-w*yy); 
            y(l)=yy+dt*((mu-rr2)*yy+w*xx);
            r(l)=sqrt(x(l)*x(l)+y(l)*y(l));
            
        end
        

        allx(im,ic,:)=x; %save all values of x, im is the index for initial mu values, ic is the index for initial conditions 
        %and the third variable is l, the number of iterations
        ally(im,ic,:)=y; %save all values of y
        allr(im,ic,:)=r; %save all values of r=sqrt(x,y)
        distFPS(im,ic,:)=r; if mu>0;  distFPS(im,ic,:)=r-sqrt(mu);end %distance to the fixed point solution ( the fixed value r_p=\sqrt{\mu})
        
        
        distP(im,ic,:)=r; %distance to the origin 
    end
end

%% Plot Results

colorVec = hsv(9)
figure(100)
subplot(221)
for tindex=1:6
    time=100*2^(tindex-1)
    for i=1:ic
       
    plot(muvec,squeeze(allr(:,i,time)),'-o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
hold on
    end
         
end

xlabel('mu')
ylabel('r_l')
title('original data')


subplot(222)
 for tindex=6:-1:1
      time=100*2^(tindex-1)
      for i=1:ic
     plot(muvec,squeeze(distFPS(:,i,time)),'o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
    
    
 hold on
      end
   
 end
 xlabel('mu')
ylabel('|r_l-r_p|')
title('distance to the fixed value of r')




 
subplot(223) 
for tindex=6:-1:1
     time=100*2^(tindex-1)
     for i=1:ic
   plot((newm-1)*time,(time)^(1/2)*squeeze(distP(:,i,time)),'o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
hold on
     end
end
xlabel('z')
ylabel('sqrt(l) r_l')
title('Collapse of the distance to the origin')
z=[-20:0.2:50]


subplot(224)
 for tindex=7:-1:1
      time=100*2^(tindex-1)
      for i=1:ic
     plot((newm(1:60)-1)*time,(time)^(1/2)*squeeze(distFPS(1:60,i,time)),'o','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
 hold on 
 % for mu>0 when computing the collapse to the distance to fixed r, we evaluate the Jacobinan on any point on the fixed r
 plot(-(newm2(62:121)-1)*time,(time)^(1/2)*squeeze(distFPS(62:121,i,time)),'x','Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
      end
 end
 xlabel('z')
ylabel('sqrt(l) r_l')
title('Collapse of the distance the fixed r')
 





%%%%%%%%%%%%%%%%%
figure(200)
immin=1
immax=121
colorVec = hsv(6)
subplot(121)
for tindex=1:6
time=25*2^tindex
    for im=immin:8:immax
 
   plot3( ones(1,32)*muvec(im), squeeze(allx(im,:,time)),squeeze(ally(im,:,time)), '-o','linewidth',0.5,'Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
   hold on
end

end
xlabel({'$\mu $'},'Interpreter','latex')
ylabel({'$ x $'},'Interpreter','latex')
zlabel({'$ y $'},'Interpreter','latex')


xlim([-0.020 0.02])
grid on
view(31,10)

subplot(122)



for tindex=1:6
time=25*2^tindex
    for im=immin:8:immax
 
   plot3( ones(1,32)*muvec(im)*time, (time)^(1/2)*squeeze(allx(im,:,time)),(time)^(1/2)*squeeze(ally(im,:,time)), '-o','linewidth',0.5,'Color', colorVec(tindex,:),'MarkerfaceColor',colorVec(tindex,:), 'MarkerSize', 3 )
   hold on
    end
end
xlabel({'$l \mu $'},'Interpreter','latex')
ylabel({'$\sqrt{l} x $'},'Interpreter','latex')
zlabel({'$\sqrt{l} y $'},'Interpreter','latex')
xlim([-10 10])
grid on
view(31,10)

