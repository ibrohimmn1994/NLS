%{
NLS 1D with SBP-SAT as an time integrator
%}


clear all;
format short
video_on = 0;


%name = "dark_II";
%name = "dark_I";
name = "bright";


quadrature = "SBP4_in_time";
%quadrature = "Gauss_Lobatto_4Nodes";
%quadrature = "Gauss_4Nodes";
%quadrature = "DI_Gauss_4Nodes"; % need larger step size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of gridpoints
m = 51
nr_blocks = 4; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x_l=-4;x_r=12;  % dark_II
%x_l=-15;x_r=30;  % dark_I
%x_l=-15;x_r=20;  % bright origional
x_l=-3;x_r=3;
len = x_r - x_l;                    
n=2*m;
h=(x_r-x_l)/(m-1);

                    % End time origional bright
t_1 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if video_on
    %theAxes=[x_l x_r -2 2]; % unifrom
    theAxes=[x_l x_r -2 2]; % Ma
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SBP4;
L = [d_1; d_m];
I = eye(m);
P=I-HI*L'*((L*HI*L')\L);
g = zeros(2,1);

x=linspace(x_l,x_r,m)';	
 V = sparse( u_exact(x,0,name) ) ;
t=0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = 0;
beta = -1; % bright
%beta = +1; % dark

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
NOTE:
SET mt = 9 for SBP4 in time and mt = 4 for the rest
%}

if quadrature == "SBP4_in_time"
    mt = 8;
else
    mt = 4;
end
block = t_1/nr_blocks;
dt = block/(mt-1);
seg = mt*nr_blocks - nr_blocks;
mx = m;
tol = 1e-12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Ix = speye(mx);
It = speye(mt);
I  = speye(mt*mx);

Quadratures_SAT

e_end = kron(e_mt', Ix);

D2_bar = sparse( kron(It, D2*P) );
D2 = sparse( kron(It, D2) );
P_bar= sparse( kron(It, P) );
BC_ = sparse( HI*L'*(inv(L*HI*L')) );
ONE = ones(mt,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch quadrature
    case "SBP4_in_time"
        T = 0:dt:block;
    case "Gauss_Lobatto_4Nodes"
        mm = (block)/2; c = (block)/2;
        T = mm*T_GL+c;
    case "Gauss_4Nodes"
        mm = (block)/2; c = (block)/2;
        T = mm*T_GL+c;
    case "DI_Gauss_4Nodes"
        T = (block)*T_GL ;
    otherwise
        disp("wrong name")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A =  sparse( Dt_bar - 1i/2*D2_bar  )  ;
[l,u,p,q,d] = lu(A);
d = sparse((1./diag(d)).*I);
pd = p*d;


j_tot = 0;
t_start = cputime;
while t< t_1

    norm_du = inf;
    u0 = kron(ONE,V);
    U0 = u0;
 
    
    BC = reshape(Boundary(x,BC_,T),mt*mx,1 );
    
    %LHS =  ( Dt_bar - 1i/2*D2_bar + 1i*beta*diag(J( u0)) ) \ eye(mt*mx) ;%
    j = 0;
    while norm_du > tol
       j = j+1;
       % u_ = P_bar*u0 + BC;
        
        b = -R*U0  + 1i/2*D2*BC - 1i*beta*F(u0) + (-Dt_bar + 1i/2*D2_bar)*(u0); % 

    
        du = q* (u\ (l\ (pd* (b)))) ;
       
      
        u0 =  (u0 + du);
        
        norm_du = norm(abs(du));
    end
    j_tot = j_tot +j;
    u0 = P_bar*u0 + BC;
    V = e_end*u0;

    t = t+block;
    T = T+block;
    if video_on && (mod(freq,1)==0)
        t ;    
        plot(x,real(V),'b','LineWidth',1);
        hold on
        plot(x,real(u_exact(x,t,name)),'r','LineWidth',1);

        title(['Numerical solution at t = ',num2str(t)]);
        axis(theAxes);
        grid;xlabel('x');
        legend('V')
        ax = gca;         
        ax.FontSize = 16;
        currFrame = getframe;
    end
    freq = freq + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j_tot
t_end = cputime - t_start

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
        disp("wrong name")
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BC = Boundary(x,BC_,t)


a=2; v=1; x0=0; Q0=0;
        
dudx_left = (a^2*sinh(a*(x0 - x(1) + t*v)).*exp(Q0*1i + v*x(1)*1i + t*(a^2/2 - v^2/2)*1i))./cosh(a*(x0 - x(1) + t*v)).^2 + (a*v*exp(Q0*1i + v*x(1)*1i + t*(a^2/2 - v^2/2)*1i)*1i)./cosh(a*(x0 - x(1) + t*v));
dudx_right = (a^2*sinh(a*(x0 - x(end) + t*v)).*exp(Q0*1i + v*x(end)*1i + t.*(a^2/2 - v^2/2)*1i))./cosh(a*(x0 - x(end) + t*v)).^2 + (a*v*exp(Q0*1i + v*x(end)*1i + t.*(a^2/2 - v^2/2)*1i)*1i)./cosh(a*(x0 - x(end) + t*v));
g =  [dudx_left; dudx_right];

BC = BC_*g;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = F(u)
f =  (abs(u).^2).*u ;
end

%{
function j = J( u)
j = abs(u).^2 + (u.*abs(u).*(u + conj(u)))./(u.*conj(u)).^(1/2);
end
%}
