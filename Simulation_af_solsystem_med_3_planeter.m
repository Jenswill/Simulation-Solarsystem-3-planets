clear all
close all
clc


G = -6.674*10^-11;

mearth = 5973.6 *10^24 ;

msol = 2e30;



%Jorden

positionJord = [149597871*1000 0]

AU = 149597871*1000;

n=900;

r = sqrt(positionJord(1)^2 + positionJord(2)^2);

Vjord = [0 sqrt(-G*msol/r)];
 

Fjord = [0 0];

A= [0 0];
T = 24*3600;

dt = T;

% Merkur

mmerkur = 3.285e23;

positionMerkur = [0.387*AU 0];

Rmerkur = sqrt(positionMerkur(1)^2 + positionMerkur(2)^2);
Fmerkur = [0 0];
Vmerkur = [0 sqrt(-G*msol/Rmerkur)];
Amerkur = [0 0];

%Mars

mMars = 6.39e23;
positionMars = [1.524*AU 0];
RMars = sqrt(positionMars(1)^2 + positionMars(2)^2);
FMars = [0 0];
VMars = [0 sqrt(-G*msol/RMars)];
AMars = [0 0];




for i = 1:n
   
    %Jorden
    
    Fjord(i+1,:) = (G*mearth*msol/r(i).^2) * positionJord(i,:)/r(i);
    
    A(i+1,:) = (Fjord(i+1,:)/mearth);
        
    Vjord(i+1,:) = Vjord(i,:) + A(i+1,:).*dt; 
    
    positionJord(i+1,:) = positionJord(i,:) + Vjord(i+1,:).*dt;
    
    r(i+1) = sqrt(positionJord(i,1).^2 + positionJord(i,2).^2);
    
    % Merkur
    
    Fmerkur(i+1,:) = (G*mmerkur*msol/Rmerkur(i).^2) * positionMerkur(i,:)/Rmerkur(i);
    
    Amerkur(i+1,:) = (Fmerkur(i+1,:)/mmerkur);
        
    Vmerkur(i+1,:) = Vmerkur(i,:) + Amerkur(i+1,:).*dt; 
    
    positionMerkur(i+1,:) = positionMerkur(i,:) + Vmerkur(i+1,:).*dt;
    
    Rmerkur(i+1) = sqrt(positionMerkur(i,1).^2 + positionMerkur(i,2).^2);
    
    
    %Mars
    
    FMars(i+1,:) = (G*mMars*msol/RMars(i).^2) * positionMars(i,:)/RMars(i);
    
    AMars(i+1,:) = (FMars(i+1,:)/mMars);
        
    VMars(i+1,:) = VMars(i,:) + AMars(i+1,:).*dt; 
    
    positionMars(i+1,:) = positionMars(i,:) + VMars(i+1,:).*dt;
    
    RMars(i+1) = sqrt(positionMars(i,1).^2 + positionMars(i,2).^2);
 
    
    pause(0.0001)
    
    plot(positionJord(i,1),positionJord(i,2),'O','color','blue')
    hold on
    plot(0.0,'O','color','yellow')
    axis([-0.5e12 0.5e12 -0.5e12 0.5e12]);
    hold on
    plot(positionMars(i,1),positionMars(i,2),'Og')
    hold on
    plot(positionMerkur(i,1),positionMerkur(i,2),'Or')
    hold on
    plot(positionJord(:,1),positionJord(:,2),'.','color','blue')
    hold on
    plot(positionMars(:,1),positionMars(:,2),'.g')
    hold on
    plot(positionMerkur(:,1),positionMerkur(:,2),'.r')
  set(legend(['Jorden'],['Solen'],['Mars'],['Merkur'],['Jordens bane'],['Mars bane'],['Merkurs bane']),'location','northoutside')

    hold off
    
    %axis([min(positionJord(i,2)) max(positionJord(i,2)) min(positionJord(i,1)) max(positionJord(i,1))])
    
end
title('Simulering af solsystemet')

xlabel('afstand [m]')
ylabel('afstand[m]')