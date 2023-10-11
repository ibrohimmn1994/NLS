%{
NLS 2D with IMEX as an time integrator
%}
clear all;
video_on = 0;
Control = "No";
contol_tol = 1e-7;

t_1= 1;
m = 101%
x_l=-2.0;x_r=2.0;
len = x_r - x_l;                  
n=m*m;
h=(x_r-x_l)/(m-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if video_on
    theAxes=[x_l x_r x_l x_r -3 3]; 
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
   % vidObj = VideoWriter('System');
    %open(vidObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0.5*h^1;

Im=speye(m);
In=speye(n);

e1=[1 0];
e2=[0 1];

SBP4;

Dxx = kron(sparse(D2), Im); Dyy = kron(Im, sparse(D2));
HIx = kron(sparse(HI), Im); HIy = kron(Im, sparse(HI));
HI = HIx*HIy;

Hx = kron(sparse(H), Im); Hy = kron(Im, sparse(H));
H = Hx*Hy;

clear HIx HIy D1 M Q M_U Q_U S_U ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_E = kron(sparse(e_m), Im); e_W = kron(sparse(e_1), Im);
e_N = kron(Im, sparse(e_m)); e_S = kron(Im, sparse(e_1));

Lx = [ e_W'; e_E' ];
Ly = [ e_S'; e_N' ];

Px = In - HI*Lx'*((Lx*HI*Lx')\Lx);
Py = In - HI*Ly'*((Ly*HI*Ly')\Ly);
P = sparse(Px*Py);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCx = sparse( HI*Lx'*((Lx*HI*Lx')\eye(2*m)) ) ;
BCy = sparse( HI*Ly'*((Ly*HI*Ly')\eye(2*m)) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC = ones(n,1);
CC(1)=0.5; CC(m)=0.5; CC(end)=0.5; CC(end-m+1)=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=linspace(x_l,x_r,m);	% discrete x values
[X,Y] = meshgrid(x,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
V = sparse(reshape(u_exact(X,Y,0),n,1));
V_d = sparse(reshape(Vd(X,Y),n,1).*In);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMEX Runge–Kutta schemes for reaction–diffusion equations
% (2) ESDIRK-ERK 4o6s
embedded_order = 3;
step_reduced = 0.1;          % failed step reduction factor
step_safety = 0.9;          % adaptivity safety factor
step_growth = 10;           % adaptivity growth bound
ONEMSM   = 1-sqrt(eps);     % coefficients to account for
ONEPSM   = 1+sqrt(eps);     %   floating-point roundoff
Tolerance   = 1.5;             % upper bound on allowed step error
rej = 0;



ESDIRK_ERK;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = 0;
t=0.0;

iteration_Matrix = sparse(In-k*a_s(2,2) * (1i/2*P*(Dxx + Dyy)  +1i*P*V_d ) );
[l,u,p,q,d] = lu(iteration_Matrix);




error_1=0.1 ;error_2=0.1;
K = [0*V 0*V 0*V 0*V 0*V 0*V];
K_hat = [0*V 0*V 0*V 0*V 0*V 0*V];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_start = cputime;
while t< t_1
    k;
    K(:,1) =  1i/2*(Dxx + Dyy)*V + 1i*V_d*V;
    K_hat(:,1) =  1i*F(V) ;
    F_s_total = b_s(1)*K(:,1);  
    G_ns_total = b_ns(1)*K_hat(:,1);
    
    F_s_embedded = b_s_hat(1)*K(:,1);
    G_ns_embedded = b_ns_hat(1)*K_hat(:,1);

    for i=2:6
       
        BC =  boundary(x',t+k*c_s(i),BCx,BCy,CC);

       
        F_s = 0*V;
        G_ns = 0*V;

        for j=1:i-1

          
            F_s = F_s +  a_s(i,j)* K(:,j);
            G_ns = G_ns + a_ns(i,j)* K_hat(:,j);
        end        
  
     
        W_guess = q* (u\ (l\ (p*(d\ (P*(V  + k*F_s +k*G_ns)+BC)))) );
   
       
        W_guess = P*(W_guess)+BC;
        K(:,i) =  1i/2*(Dxx + Dyy)*W_guess + 1i*V_d*W_guess;  
        K_hat(:,i) = 1i*F(W_guess);% + 1i*V_d*W_guess;

        F_s_total = F_s_total + b_s(i)*K(:,i) ;
        G_ns_total = G_ns_total +b_ns(i)*K_hat(:,i);

        F_s_embedded =  F_s_embedded +b_s_hat(i)*K(:,i);
        G_ns_embedded = G_ns_embedded + b_ns_hat(i)*K_hat(:,i);

    end
    %---------------------------------------------------
    BC =  boundary(x',t+k,BCx,BCy,CC);
   
 
 if Control == "Yes"
    
    V_new= P*(V + k*F_s_total + k*G_ns_total   )+BC  ;
    V_hat = P*(V+ k*F_s_embedded +k*G_ns_embedded) +BC;
    
    V_error = (V_hat-V_new);
    step_error = max( (norm(V_error./(contol_tol*(V_new) + contol_tol),inf)), eps);
    %step_error = (max( (sum(real(V_error./(rtol*V_new + atol))))/length(V_error), eps));
   % step_error = (max(norm(V_error./(V_new),inf), eps));
    
    %step_error = norm(V_error);%./(V_new ),inf), eps);
    k_old = k;
    if (step_error > Tolerance)
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
   iteration_Matrix = sparse(In-k*a_s(2,2) * (1i/2*P*(Dxx + Dyy)  + 1i*P*V_d  ) );%+ diag(1i*V_d)
   [l,u,p,q,d] = lu(iteration_Matrix);
  % d = sparse((1./diag(d)).*In);
  % pd = p*d;
  
    freq = freq+1;


 else
   
     V_new = P*(V + k*F_s_total + k*G_ns_total)+BC ;
     V = V_new;
     
     t = t+k;

 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_end = cputime - t_start

U_approx = reshape(V,n,1);
U_exact = reshape(u_exact(X,Y,t),n,1);

err = sqrt(real(U_exact)-real(U_approx))'*H*sqrt(real(U_exact)-real(U_approx))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ui = u_exact(x,y,t)
xei = 1; k = 1;

ui = sqrt(xei) * exp(-k/2* (x.^2 + y.^2))*exp(-1i*k*t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = boundary(x,t,BCx,BCy,C)
xei = 1; k = 1;

West = sqrt(xei) * exp(-k/2* (x(1).^2 + x.^2))*exp(-1i*k*t); 
East = sqrt(xei) * exp(-k/2* (x(end).^2 + x.^2))*exp(-1i*k*t);


South = sqrt(xei) * exp(-k/2* (x.^2 + x(1).^2))*exp(-1i*k*t);
North = sqrt(xei) * exp(-k/2* (x.^2 + x(end).^2))*exp(-1i*k*t); 

gx =  [West ; East];
gy =  [South; North];

b = (BCx*gx + BCy*gy).*C;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = F(u)
f =  (abs(u).^2).*u ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = Vd(x,y)
xei = 1; k = 1;
v = -k/2 * (x.^2 + y.^2) - xei*exp(-k*(x.^2 + y.^2));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
