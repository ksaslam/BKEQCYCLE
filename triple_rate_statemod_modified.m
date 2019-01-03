function out = rate_state(t,x,flag,const)
%%       out = rate_state(t,x,flag,const)

%
% INPUTS:
%       t = times
%       x = model vector
%	x(1) is velocity for slider 1
%	x(2) is state  variable for slider 1
%	x(3) is shear stress for slider 1
%	x(4) is velocity for slider 2
%	x(5) is state  variable for slider 2
%	x(6) is shear stress for slider 2
%	x(7) is velocity for slider 3
%	x(8) is state  variable for slider 3
%	x(9) is shear stress for slider 3
%
% OUTPUTS: time derivatives of model vector
%   out(1) is dx(1)/dt (acceleration) for slider 1
%   out(2) is dx(2)/dt (d_theta/dt) for slider 1
%   out(3) is dx(3)/dt (stressing rate) for slider 1
%   out(4) is dx(1)/dt (acceleration) for slider 2
%   out(5) is dx(2)/dt (d_theta/dt) for slider 2
%   out(6) is dx(3)/dt (stressing rate) for slider 2
%   out(7) is dx(1)/dt (acceleration) for slider 3
%   out(8) is dx(2)/dt (d_theta/dt) for slider 3
%   out(9) is dx(3)/dt (stressing rate) for slider 3

%  unwrap constants
	d_c1 = const(1);
	d_c2 = const(2);
    d_c3 = const(3);
	A1 = const(4); 
    A2 = const(5);
    A3 = const(6);
	B1 = const(7); 
	B2 = const(8); 
	B3 = const(9); 
	sigma1 = const(10);
    sigma2 = const(11);
    sigma3 = const(12);
    
	k1 = const(13); 
	k2 = const(14);
    ks1 = const(15);
    ks2 = const(16);
    ks3 = const(17);
	v_inf = const(18);
	eta = const(19);
       et = const(20);
kw = const(21);

    v1     = x(1);
	theta1 = x(2);
	tau1   = x(3);
    v2     = x(4);
	theta2 = x(5);
	tau2   = x(6);
    %v3     = x(7);
	%theta3 = x(8);
	%tau3   = x(9);

     %return rates
     out(2) = 1-theta1*v1/d_c1; %theta1
     %out(3) = (ks1+(kw*et/(et+kw*t)))*(v_inf - v1)+k1*(v2 - v1);
     out(3) = (ks1+(kw*et/(et+kw*t)))*(v_inf - v1)+k1*(v2 - v1)+ 18e20* ( exp(- ( (t-31*31536000.).^2/1000000^2 ) )  ).*sin(20.*t);  %tau1
     out(1) =  (out(3)/sigma1 - B1*out(2)/theta1 ) /(eta/sigma1 + A1/v1); % v1
     
     out(5) = 1-theta2*v2/d_c2; %theta2
     out(6) = (ks2+(kw*et/(et+kw*t)))*(v_inf - v2)-k1*(v2-v1);  %tau2
     out(4) =  (out(6)/sigma2 - B2*out(5)/theta2 ) /(eta/sigma2 + A2/v2); % v2
     
     %out(8) = 1-theta3*v3/d_c3; %theta3
     %out(9) = (ks3+(kw*et/(et+kw*t)))*(v_inf - v3)-k2*(v3-v2);  %tau3
     %out(7) =  (out(9)/sigma3 - B3*out(8)/theta3 ) /(eta/sigma3 + A3/v3); % v3
     
     out = out';

