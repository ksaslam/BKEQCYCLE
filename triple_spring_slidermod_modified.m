% Earthquake simulation 
% Triple spring-slider system with rate-state friction
% THREE spring-sliders in series all attached to the far field stress

% CONSTITUTIVE PARAMETERS
   d_c1 = 10^-3;                % critical displacement (m)
   d_c2 = 10^-3;                % critical displacement (m)
   d_c3 = 10^-3;                % critical displacement (m)
   A2 =   0.01;  %slider 2 is pulled at far field velocity and by slider 3
   B2 =   0.005;
   A1 =   0.01;  %slider 1 is pulled by slider 2 and far field velocity
   B1 =   0.015;
   A3 =   0.01;  %slider 3 is pulled at far vield velocity
   B3 =   0.015;
   mu_0 = 0.6;                % nominal friction coefficient

% BOUNDARY CONDITIONS
   sigma1 = 30;               % normal stress (MPa)
   sigma2 = 30;               % normal stress (MPa)
   sigma3 = 30;               % normal stress (MPa)
   
   %spring stiffness, k
   G=3*10^4;                  % Earth elastic modulus (MPa)
   L1=10^4;  k1=G*.9167*0.4/L1; ks1=G*.9167/L1;     % L1 is radius of circular slip patch in m
   L2=10^4;  k2=G*.9167*0.4/L2; ks2=G*.9167/L2;     % L2 is radius of circular slip patch in m
   L3=10^4;  ks3=G*.9167/L3;      % L2 is radius of circular slip patch in m
   
   v_inf  =   0.03/31536000.;           % load-point velocity in m/yr
   

   % SEISMIC RADIATION DAMPING TERM
   eta = 10^3*1.0e-7 * 31536000.0;       % units of MPa-yr/m

et = 30 * 31536000.0; %units of MPa-yr 
kw=3; %units of MPa
 
% INITIAL CONDITIONS

   v02 = v_inf;
   v01 = 10^-3/31536000.0;
   %v03 = 10^-3;
   %steady state
   tau02    =    sigma2*(mu_0 + (A2-B2)*log(v02));
   theta02 = d_c2/v02;
   tau01    =    sigma1*(mu_0 + (A1-B1)*log(v01));
   theta01 = d_c1/v01;
   %tau03    =    sigma3*(mu_0 + (A3-B3)*log(v03));
   %theta03 = d_c3/v03;
   x0      =   [v01; theta01; tau01; v02; theta02; tau02];

     
  
   
% START AND STOP TIME
   t0      =  0;
   tf= 977615001;
   %tf      =  977615000;% time in years
   %tspan= [t0 tf]
tspan= t0:100:tf; 
B= 977615001:.001: 977619000;
C= 977619001:100:80*31536000.0;
tspan= [tspan B C];
  % tspan= [0:100:300*31536000]
%  RUN ODE SOLVER
 
 const = [d_c1; d_c2; d_c3; A1; A2; A3; B1; B2; B3; sigma1; sigma2; sigma3; k1; k2; ks1; ks2; ks3; v_inf; eta; et; kw];
 opt = [];
  
  [t,x] = ode23s('triple_rate_statemod_modified', tspan,x0,opt,const);
  u1 = cumsum(x(:,1));  % displacement
  %u2 = cumsum(x(:,7));  % displacement
  delt=diff(t);
  delt(end+1)=delt(end);
  %
%PLOT OUTPUT
%
%plots for slider 2
figure
%phase plot (log velocity vs. stress)
    subplot(3,2,1); semilogx(x(:,4), x(:,6),'b'); title('phase plane plot - block2'); hold on
    %plot steady-state line
    vel=logspace(-5e-7,5e-5);
    tss=sigma2*(mu_0+(A2-B2)*log(vel));
    semilogx(vel, tss,'r:')
%state variable
    subplot(3,2,2); plot(t,x(:,5),'b'); title('theta2 vs time'); hold on
    %compare with steady-state value
    plot(t,d_c2./x(:,4),'r')
%velocity
    subplot(3,2,3); plot(t,x(:,4),'b'); title('v2 vs time'); hold on
%displacement
    subplot(3,2,4); plot(t,cumsum(delt.*x(:,4)),'b'); title('disp2 vs time'); hold on
%stress
    subplot(3,2,5:6); plot(t,x(:,6),'b'); title('stress2 vs time'); hold on

%plots for slider 1
figure
%phase plot (log velocity vs. stress)
    subplot(3,2,1); semilogx(x(:,1), x(:,3),'b'); title('phase plane plot - block1'); hold on
    %plot steady-state line
    vel=logspace(-5,5);
    tss=sigma1*(mu_0+(A1-B1)*log(vel));
    semilogx(vel, tss,'r:')
%state variable
    subplot(3,2,2); plot(t,x(:,2),'b'); title('theta1 vs time'); hold on
%velocity
    subplot(3,2,3); plot(t,x(:,1),'b'); title('vel vs time'); hold on
%displacement
    subplot(3,2,4); plot(t,cumsum(delt.*x(:,1)),'b'); title('disp1 vs time'); hold on
%stress
    subplot(3,2,5:6); plot(t,x(:,3),'b'); title('stress1 vs time'); hold on
    

% Stress Vs time for all  blocks
figure
   subplot(1,1,1); plot(t,x(:,6),'b'); title('stress vs time'); hold on
	          % plot(t,x(:,9),'r'); title('stress vs time'); hold on
                   plot(t,x(:,3),'g'); title('stress vs time'); hold on
                   
%displacement vs time for all blocks
figure
%displacement
    subplot(1,1,1); plot(t,cumsum(delt.*x(:,4)),'b'); title('displacement vs time'); hold on
                    plot(t,cumsum(delt.*x(:,1)),'g'); title('disp1 vs time'); hold on
                    %plot(t,cumsum(delt.*x(:,7)),'r'); title('disp3 vs time'); hold on
                    
%velocity v time for all blocks
figure
%velocity
    subplot(1,1,1); plot(t,x(:,4),'b'); title('velocity vs time'); hold on
                    plot(t,x(:,1),'g'); title('velocity vs time'); hold on
                    %plot(t,x(:,7),'r'); title('velocity vs time'); hold on
                    

figure
%displacements
	subplot(3,1,1); plot(t,cumsum(delt.*x(:,4)),'b'); title('displacement vs time'); hold on
        subplot(3,1,2); plot(t,cumsum(delt.*x(:,1)),'g'); title('disp1 vs time'); hold on
       % subplot(3,1,3); plot(t,cumsum(delt.*x(:,7)),'r'); title('disp3 vs time'); hold on

       % Stress zoomed
figure
   subplot(1,1,1); plot(t,x(:,6),'b'); title('stress vs time'); hold on
	         %  plot(t,x(:,9),'r'); title('stress vs time'); hold on
                   plot(t,x(:,3),'g'); title('stress vs time'); hold on
                 
% displacement zoomed
figure
%displacement
    subplot(1,1,1); plot(t/31536000.,cumsum(delt.*x(:,4)),'b'); title('displacement vs time'); hold on
                    plot(t/31536000.,cumsum(delt.*x(:,1)),'g'); title('disp1 vs time'); hold on
                   % plot(t,cumsum(delt.*x(:,7)),'r'); title('disp3 vs time'); hold on
                 
% velocity zoomed
figure
%velocity
    subplot(1,1,1); plot(t/31536000,x(:,4),'b'); title('velocity vs time'); hold on
                    plot(t,x(:,1),'g'); title('velocity vs time'); hold on
                   % plot(t,x(:,7),'r'); title('velocity vs time'); hold on
                   


