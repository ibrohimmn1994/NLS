%{
NLS 2D with SBP-SAT as an time integrator
from:
Exact Solutions of two-Dimensional Nonlinear Schrödinger Equations with
External Potentials
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

video_on = 0;
%quadrature = "SBP4_in_time";
quadrature = "Gauss_Lobatto_4Nodes";
%quadrature = "Gauss_4Nodes";
%quadrature = "DI_Gauss_4Nodes"; % need larger step size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 51 % number of grid points
nr_blocks = 3;
tol = 1e-3;

t_1 = 1.0;  % end time
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
P = Px*Py;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCx = sparse( HI*Lx'*((Lx*HI*Lx')\eye(2*m)) );
BCy = sparse( HI*Ly'*((Ly*HI*Ly')\eye(2*m)) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Px Py Lx Ly HI In e_E e_N e_W e_S Hx Hy
C = ones(n,1);
C(1)=0.5; C(m)=0.5; C(end)=0.5; C(end-m+1)=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBP integrator part
mx = n;
if quadrature == "SBP4_in_time"
    mt = 9;
else
    mt = 4;
end
block = t_1/nr_blocks;
dt = block/(mt-1);
seg = mt*nr_blocks - nr_blocks;
mt_total = seg + 1;
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = kron( ones(mt,1),C);
D2 = sparse(Dxx + Dyy);

clear Dxx Dyy 
Ix = speye(mx);
It = speye(mt);
I = speye(mt*mx);
CC = CC.*I;
%BCx = kron(It,BCx);
%BCy = kron(It,BCy );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D2_bar = kron(It,sparse(D2*P));%<------
P_bar =  kron(It,sparse(P));%
D2 = kron(It,sparse(D2));%

x=linspace(x_l,x_r,m);	% discrete x values
[X,Y] = meshgrid(x,x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = sparse(reshape(u_exact(X,Y,0),n,1));

V_d = sparse(kron(ones(mt,1),reshape(Vd(X,Y),n,1)));
V_potential = sparse(V_d.*I);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Quadratures_SAT
e_end = kron(e_mt', Ix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A  = sparse(Dt_bar - 1i/2*D2_bar - 1i*V_potential );
[l,u,p,q,d] = lu(A);
d = sparse((1./diag(d)).*I);
pd = p*d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%n
%%%%%%%%%%%%%%%%

clear Ix It
switch quadrature
    case "SBP4_in_time"
        T = 0:dt:block;
    case "Gauss_Lobatto_4Nodes"
        m = (block)/2; c = (block)/2;
        T = m*T_GL+c;
    case "Gauss_4Nodes"
        m = (block)/2; c = (block)/2;
        T = m*T_GL+c;
    case "DI_Gauss_4Nodes"
        T = (block)*T_GL ;
    otherwise
        disp("wrong name")
end
[TT,XX] = meshgrid(T,x);
%TT= sparse(TT); XX = sparse(XX);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = 0;
t=0.0;
ONE = ones(mt,1);
t_start = cputime;

while t <  t_1-block/2
   t;
    BC =  Boundary(TT,x,XX,BCx,BCy,CC,mx*mt);

    u0 = kron(ONE,V);  
    U0 = u0; 
    norm_du = inf;
    %u0 = zeros(size(u0));

   j=0;
   
    while norm_du > tol
        j = j+1;

        
        u_ = P_bar*u0 +BC;
               
        b =  -R*U0 + 1i*F(u_)  -Dt_bar*u0 + (1i/2*D2 + 1i*V_potential)*(u_) ;

        %b = ( -R*U0 + 1i/2*D2*BC + 1i*F(u0) + (-Dt_bar + 1i/2*D2_bar + 1i*V_potential)*(u0) );
              
        du = q* (u\ (l\ (pd* (b)))) ;
        %du = q* (u\ (l\ (p*(d\ (b))))) ;
         
        u0 = u0 + du;
        norm_du = normest(abs(du));
    end
    j;
    u0 = P_bar*u0 +BC;
    V = e_end*u0;    

    t = t+block;
    TT = TT+block;
 
   
    if video_on &&(mod(freq,1)==0)
         
        surf(X,Y,reshape(real(V),m,m));  
        shading interp
       % view(90,0)
        colorbar
        title(['Numerical solution at t = ',num2str(t)]);
        %theAxes=[x_l x_r y_d y_u];
        axis(theAxes);
        grid;xlabel('x');
        legend('u')
        ax = gca;          % current axes
        ax.FontSize = 16;
        currFrame = getframe;
    end
    freq = freq + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_end = cputime - t_start

U_approx = reshape(V(1:n),n,1);
U_exact = reshape(u_exact(X,Y,t),n,1);
err = sqrt(real(U_exact)-real(U_approx))'*H*sqrt(real(U_exact)-real(U_approx))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BC = Boundary(t,x,X,BCx,BCy,C,len)
%

xei = 1; k = 1;      
East = sqrt(xei) * exp(-k/2* (x(1).^2 + X.^2)).*exp(-1i*k*t);
West = sqrt(xei) * exp(-k/2* (x(end).^2 + X.^2)).*exp(-1i*k*t);

South = sqrt(xei) * exp(-k/2* (X.^2 + x(1).^2)).*exp(-1i*k*t); 
North = sqrt(xei) * exp(-k/2* (X.^2 + x(end).^2)).*exp(-1i*k*t);

gx =  ([East ; West ]);
gy =  ([South ; North]);

BC = C*(reshape(BCx*gx + BCy*gy,len,1)) ;%.*C;
%BC = BCx*gx + BCy*gy;%.*C;
%size(BC)
%BC = [BC(:,1) ; BC(:,2); BC(:,3); BC(:,4)].*C;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ui = u_exact(x,y,t)
xei = 1; k = 1;

ui = sqrt(xei) * exp(-k/2* (x.^2 + y.^2))*exp(-1i*k.*t);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = F(u)
f =  ( (abs(u).^2).*u ) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = Vd(x,y)
xei = 1; k = 1;
v = -k/2 * (x.^2 + y.^2) - xei*exp(-k*(x.^2 + y.^2));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function j = J( u)
j = abs(u).^2 + (u.*abs(u).*(u + conj(u)))./(u.*conj(u)).^(1/2);
end