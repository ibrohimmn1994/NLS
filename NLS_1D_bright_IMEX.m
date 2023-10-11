clear all;
format short
video_on = 0;
%name = "dark_II";
%name = "dark_I";
name = "bright";

Control = "yes";

    % Number of gridpoints
m = 1601
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x_l=-4;x_r=12;  % dark_II
%x_l=-15;x_r=30;  % dark_I
%x_l=-15;x_r=20;  % bright origional
x_l=-3;x_r=3;
len = x_r - x_l;                   
h=(x_r-x_l)/(m-1);

    % End time origional bright
t_1 = 1.0;
    % maximum step size is given with 5*h^2
 
k=0.5*h^1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if video_on
    theAxes=[x_l x_r -2 2]; 
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SBP4;
L = [d_1; d_m];
I = speye(m);
P = sparse(I-HI*L'*((L*HI*L')\L) );
D2 = sparse(D2);
BC_ = sparse( HI*L'*(inv(L*HI*L') ) );


x = linspace(x_l,x_r,m)';	
V = sparse( u_exact(x,0,name) ) ;
t=0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = 0;
beta = -1; % bright
%beta = +1; % dark

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMEX Runge–Kutta schemes 
% (2) ESDIRK-ERK 4o6s
embedded_order = 3;
step_reduced = 0.1;          % failed step reduction factor
step_safety = 0.9;          % adaptivity safety factor
step_growth = 10;           % adaptivity growth bound
ONEMSM   = 1-sqrt(eps);     % coefficients to account for
ONEPSM   = 1+sqrt(eps);     %   floating-point roundoff
Tolerance   = 1.5;             % upper bound on allowed step error
rej = 0;
atol = 1e-8;
rtol = 1e-8;

ESDIRK_ERK;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iteration_Matrix = sparse(I-k*a_s(2,2)*0.5*1i*P*D2);
[l,u,p,q,d] = lu(iteration_Matrix);
%d = sparse((1./diag(d)).*I);
%pd = p*d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_1=0.1 ;error_2=0.1;
K = [0*V 0*V 0*V 0*V 0*V 0*V];
K_hat = [0*V 0*V 0*V 0*V 0*V 0*V];

t_start = cputime;

while t<t_1
    K(:,1) =  0.5*1i*D2*V;
    K_hat(:,1) = - beta* 1i*F(V);
    F_s_total = b_s(1)*K(:,1);  
    G_ns_total = b_ns(1)*K_hat(:,1);
    
    F_s_embedded = b_s_hat(1)*K(:,1);
    G_ns_embedded = b_ns_hat(1)*K_hat(:,1);

    for i=2:6
       

        BC = boundary(x,t+k*c_s(i),name, BC_);
       
        F_s = 0*V;
        G_ns = 0*V;

        for j=1:i-1

          
            F_s = F_s +  a_s(i,j)* K(:,j);
            G_ns = G_ns + a_ns(i,j)* K_hat(:,j);
        end        
           
     
        W_guess = q* (u\ (l\ (p*(d\ (P*(V  + k*F_s +k*G_ns)+BC)))) );
   
      
   
       
        W_guess = P*(W_guess)+BC;
        K(:,i) =  0.5*1i*D2*W_guess; 
        K_hat(:,i) = - beta* 1i*F(W_guess);

        F_s_total = F_s_total + b_s(i)*K(:,i) ;
        G_ns_total = G_ns_total +b_ns(i)*K_hat(:,i);

        F_s_embedded =  F_s_embedded +b_s_hat(i)*K(:,i);
        G_ns_embedded = G_ns_embedded + b_ns_hat(i)*K_hat(:,i);

       

    end
    %---------------------------------------------------
    BC = boundary(x,t+k,name, BC_);
  

 if Control == "Yes"
    
    V_new= P*(V + k*F_s_total + k*G_ns_total   )+BC  ;
    V_hat = P*(V+ k*F_s_embedded +k*G_ns_embedded) + BC;
    
    V_error = V_hat-V_new;
    step_error = max(norm(V_error./(rtol*V_new + atol),inf), eps);
    
    %step_error = norm(V_error);%./(V_new ),inf), eps);
    k_old = k;
    if (step_error > Tolerance*ONEPSM*0.1)
        rej = rej + 1;
        k = k*step_reduced;           
    else
        V = V_new; 
        t = t+k;
       %epss = 1e-1;
       epss = 0.1;
       %(PID-controller)
     
       k = step_safety*k*(epss/step_error)^(0.49/embedded_order)*(error_1/epss)^(0.34/embedded_order)*(epss/error_2)^(0.1/embedded_order);

       error_2 = error_1; 
       error_1 = step_error;
       
       %(I-controller)
       %k = step_safety * k * step_error^(-1.0/(embedded_order-1))
    end
   k = min(step_growth*k_old, k) ;
   iteration_Matrix = sparse(I-k*a_s(2,2)*0.5*1i*P*D2);
   [l,u,p,q,d] = lu(iteration_Matrix);

  
 else
    
     V_new = P*(V + k*F_s_total + k*G_ns_total)+BC ;
     V = V_new;
     t = t+k;

 end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_end =  cputime - t_start
error = sqrt( real(u_exact(x,t,name)) - real(V) )'*H*sqrt( real(u_exact(x,t,name))-real(V) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = u_exact(x,t,name)
switch name
    case "dark_II"
        a=3; v=2; k=2; x0=0; Q0=0; 
        u = ( a*tanh(a*(x-v*t-x0)) + 1i*(v-k)).*exp(1i*( k*x + Q0 - 1/2*(k^2+2*a^2+2*(v-k)^2)*t ));
    case "dark_I"
        a=2; v=1; x0=0; Q0=0;
        u = ( a*tanh(a*(x-v*t-x0)) + 1i*v )*exp( 1i*(-(a^2+v^2)*t+Q0) );
    case "bright"
        a=2; v=1; x0=0; Q0=0;
        u = a*sech( a*(x-v*t-x0) ).*exp( 1i*(v*x-1/2*(v^2-a^2)*t + Q0) );

    otherwise
        return;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BC = boundary(x,t,name, BC_)
switch name
    case "dark_II"
        a=3; v=2; k=2; x0=0; Q0=0;
        dudx_left = - k*exp(Q0*1i - t*(a^2 + k^2/2 + (k - v)^2)*1i + k*x(1)*1i)*(k*1i - v*1i + a*tanh(a*(x0 - x(1) + t*v)))*1i - a^2*exp(Q0*1i - t*(a^2 + k^2/2 + (k - v)^2)*1i + k*x(1)*1i)*(tanh(a*(x0 - x(1) + t*v))^2 - 1);
        dudx_right = - k*exp(Q0*1i - t*(a^2 + k^2/2 + (k - v)^2)*1i + k*x(end)*1i)*(k*1i - v*1i + a*tanh(a*(x0 - x(end) + t*v)))*1i - a^2*exp(Q0*1i - t*(a^2 + k^2/2 + (k - v)^2)*1i + k*x(end)*1i)*(tanh(a*(x0 - x(end) + t*v))^2 - 1);

    case "dark_I"
        a=2; v=1; x0=0; Q0=0;
        dudx_left = -a^2*exp(Q0*1i - t*(a^2 + v^2)*1i)*(tanh(a*(x0 - x(1) + t*v))^2 - 1);
        dudx_right = -a^2*exp(Q0*1i - t*(a^2 + v^2)*1i)*(tanh(a*(x0 - x(end) + t*v))^2 - 1);

    case "bright" % setting v=0 renders the soliton oscillatory with out moving 
        a=2; v=1; x0=0; Q0=0;
        dudx_left = (a^2*sinh(a*(x0 - x(1) + t*v))*exp(Q0*1i + v*x(1)*1i + t*(a^2/2 - v^2/2)*1i))/cosh(a*(x0 - x(1) + t*v))^2 + (a*v*exp(Q0*1i + v*x(1)*1i + t*(a^2/2 - v^2/2)*1i)*1i)/cosh(a*(x0 - x(1) + t*v));
        dudx_right = (a^2*sinh(a*(x0 - x(end) + t*v))*exp(Q0*1i + v*x(end)*1i + t*(a^2/2 - v^2/2)*1i))/cosh(a*(x0 - x(end) + t*v))^2 + (a*v*exp(Q0*1i + v*x(end)*1i + t*(a^2/2 - v^2/2)*1i)*1i)/cosh(a*(x0 - x(end) + t*v));
    otherwise
        return;
end
g = [dudx_left; dudx_right];
BC = BC_*g;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = F(u)
f =  (abs(u).^2).*u ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
