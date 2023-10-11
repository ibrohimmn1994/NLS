%{
NLS 2D saturated with SBP-projection as an time integrator and GMRE/FGMRES
iterative solver
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

m = 201
nr_blocks = 14 ;
tol = 1e-8;

restrt=2;
max_it=4;
tol_gmres = 1e-1; % GMRES tolerance

disp("-----")
t_1= 1;
x_l=-2.0;x_r=2.0;
len = x_r - x_l;                  
n=m*m;

h=(x_r-x_l)/(m-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if video_on
    theAxes=[x_l x_r x_l x_r -3 10]; 
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
   % vidObj = VideoWriter('System');
    %open(vidObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

clear HIx HIy M Q M_U Q_U S_U ;
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
    mt = 8;
else
    mt = 4;
end
block = t_1/nr_blocks;
dt = block/(mt-1);
seg = mt*nr_blocks - nr_blocks;
mt_total = seg + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = kron( ones(mt,1),C);
D2 = sparse(Dxx + Dyy);

clear Dxx Dyy C
Ix = speye(mx);
It = speye(mt);
I = speye(mt*mx);
CC = CC.*I;

D2_bar = kron(It,sparse(D2*P));%<------
P_bar =  kron(It,sparse(P));%
D2 = kron(It,sparse(D2));%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=linspace(x_l,x_r,m);	% discrete x values
[X,Y] = meshgrid(x,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
V = sparse(reshape(u_exact(X,Y,0),n,1));

V_d = sparse(kron(ones(mt,1),reshape(Vd(X,Y),n,1)));
V_potential = V_d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Quadratures_projection;
e_end = kron(e_mt', Ix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = sparse(I - JF*( 1i*D2_bar ));

tc = cputime;
[l,u,p,q,d] = lu(A);

%[l,u,p] = ilu(A);
%q=1;d=1;
%options.droptol = 1e-2;
%options.type = "crout";
%[l,u,p] = ilu(A,options);
%q=1;d=1;

tc = cputime - tc

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
[TT,XX] = meshgrid(T,x);
ONE = ones(mt,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0.0;
freq = 0;
du = sparse(zeros(mx*mt,1));

t_start = cputime;
while t< t_1 -block/2

  
    BC =  Boundary(TT,x,XX,BCx,BCy,CC,mt*mx);
     
    u0 = kron(ONE,V);
    U0 = u0; 
    norm_du = inf;
    j=0;
    

    while norm_du > tol
    j = j+1;
  
      
        u_ = P_bar*u0 +BC;

        A = (I - JF*( 1i*D2_bar - 1i*jac(u_,V_potential).*I ));
      
       

        b =  U0 - (u0) + JF*(  1i*D2*(u_)  - 1i*F(u_,V_potential)) ;


      %  [du,count] = GMRES_main( A, du, b, restrt, max_it, tol_gmres,l,u,p,q,d);%pqd

        [du,count] = F_GMRES( A, du, b, restrt, max_it, tol_gmres,l,u,p,q,d);
        %du = q* (u\ (l\ (p*(d\ (b))))) ;
        
        u0 = u0 + du;    
        
        norm_du = norm(abs(du));
    end
    j;
    u0 = P_bar*u0+BC;
    
    V = e_end*u0;          
    t = t + block;
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
       % writeVideo(vidObj,currFrame);
    end
    freq = freq + 1;

end

t_end = cputime - t_start

U_approx = V;
U_exact = reshape(u_exact(X,Y,t),n,1);
err = sqrt(real(U_exact)-real(U_approx))'*H*sqrt(real(U_exact)-real(U_approx))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BC = Boundary(t,x,X,BCx,BCy,C,len)

k = 1;xi=1;

West = sqrt(xi)*exp(-k/2* (x(1).^2 + X.^2 + 4*1i*t) ) ;
East = sqrt(xi)*exp(-k/2* (x(end).^2 + X.^2 + 4*1i*t) ) ;


South =sqrt(xi)*exp(-k/2* (X.^2 + x(1).^2 + 4*1i*t) ) ;
North =sqrt(xi)*exp(-k/2* (X.^2 + x(end).^2 + 4*1i*t) ) ;


gx =  ([East ; West ]);
gy =  ([South ; North]);

BC = C*(reshape(BCx*gx + BCy*gy,len,1)) ;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ui = u_exact(x,y,t)

k = 1; xi = 1;

ui = sqrt(xi)*exp(-k/2* (x.^2 + y.^2 + 4*1i*t) ) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = F(u,Vd)
E0=1;
f =  (E0*u)./(1+ Vd + abs(u).^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = Vd(x,y)
k = 1; xi=1; E0=1;
v = (  ( E0 - (k^2)*(x.^2 + y.^2) )./((k^2)* (x.^2 + y.^2))  ) ...
    - xi*exp(-k*(x.^2 + y.^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function j = jac( u,Vd)
E0=1;
j = ( E0./(abs(u).^2 + Vd + 1) - (E0*u.*abs(u).*(u + conj(u)))./((u.*conj(u)).^(1/2).*(abs(u).^2 + Vd + 1).^2) );
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GMRE
function [x,count] = GMRES_main( A, x, b, restrt, max_it, tol,L,U,P,Q,D)
iter = 0;                                         
count = 0;
bnrm2 = normest( b );

if  ( bnrm2 == 0.0 )
    disp("Change the GMRES parameters")
    bnrm2 = 1.0; 
end

n = length(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r =  Q* (U\ (L\ (P* (D\ (b-A*x)))) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_norm = normest(r);
error = r_norm / bnrm2;

if ( error < tol ) 
  
    disp("Not convergent")
    return 
end
                            
m = restrt;
V = zeros(n,m+1);
H = zeros(m+1,m);
cs = zeros(m,1);
sn = zeros(m,1);
e1  = zeros(n,1);
e1(1) = 1.0;

while iter <= max_it
  
     V(:,1) = r / r_norm;
     s = r_norm*e1;
     %%%%%%%%%%%%%%%
     for i = 1:m  
         count = count + 1;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         w = Q* (U\ (L\ (P* (D\ (A*V(:,i))))) );
 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
         for k = 1:i
             H(k,i)= w'*V(:,k);
             w = w - H(k,i)*V(:,k);
         end
         H(i+1,i) = normest( w );
         V(:,i+1) = w / H(i+1,i); %<------------------ % just a number
         for k = 1:i-1                              
             temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
             H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
             H(k,i)   = temp;
         end
         [cs(i),sn(i)] = givens_rotation( H(i,i), H(i+1,i) );
         temp   = cs(i)*s(i);                       
         s(i+1) = -sn(i)*s(i);
         s(i)   = temp;
         H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
         H(i+1,i) = 0.0;
         error  = abs(s(i+1)) / bnrm2;
         if ( error <= tol )                       
            y = H(1:i,1:i) \ s(1:i);     
            x = x + V(:,1:i)*y;
            break;
         end
     end
     %%%%%%%%%%%%%%%%%%
     iter = iter + m;

     if ( error <= tol )
       
         break
     end

     y = H(1:m,1:m) \ s(1:m); %<---------- Here, done every m times and it is of size m
     x = x + V(:,1:m)*y; 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     r =  Q* (U\ (L\ (P* (D\ (b-A*x)))) );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      
     r_norm = normest(r);
     s(i+1) = r_norm;
     error = s(i+1) / bnrm2; 
    
     if ( error <= tol )
        
         break
     end
end

if ( error > tol )
    disp("Not convergent")
end                
function [cs, sn] = givens_rotation(v1, v2)
  if (v1 == 0)
    cs = 0;
    sn = 1;
  else
    t = sqrt(v1^2 + v2^2);
    cs = v1 / t;  % see http://www.netlib.org/eispack/comqr.f
    sn = v2 / t;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%count  % total number of iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% FGMRES




function [x,count,count_t] = F_GMRES( A, x, b, restrt, max_it, tol,L,U,P,Q,D)
iter = 0;                                   
count = 0;
count_t = 0;
b_norm = normest( b );

if  ( b_norm == 0.0 )
    %disp("This might not be used")
    b_norm = 1.0; 
end

n = length(A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r =  ( b-A*x );  %\M   % not recorded and MVP is not counted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_norm = normest(r);
error = r_norm / b_norm;
                                % initialize workspace
m = restrt;
V = zeros(n,m+1);
H = zeros(m+1,m);
cs = zeros(m,1);
sn = zeros(m,1);
e1    = zeros(n,1);
e1(1) = 1.0;

Z = zeros(n,m+1);
z = zeros(size(Z(:,1)));
while iter <= max_it

     V(:,1) = r / r_norm;
     s = r_norm*e1;
     %%%%%%%%%%%%%%%
     for i = 1:m  % Arnoldi
         count = count + 1;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
         [z, count_] = S_GMRES( A, z, V(:,i), 4,4, 1e-2,L,U,P,Q,D);
         count_t = count_t + count_;
      
         Z(:,i) = z;
         w =  A*Z(:,i); 
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         for k = 1:i
             H(k,i)= w'*V(:,k);
             w = w - H(k,i)*V(:,k);
         end
         H(i+1,i) = normest( w );
         V(:,i+1) = w / H(i+1,i); %<------------------ % just a number
         for k = 1:i-1                              % apply Givens rotation
             temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
             H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
             H(k,i)   = temp;
         end
         [cs(i),sn(i)] = givens_rotation( H(i,i), H(i+1,i) ); % form i-th rotation matrix
         temp   = cs(i)*s(i);                        % approximate residual norm
         s(i+1) = -sn(i)*s(i);
         s(i)   = temp;
         H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
         H(i+1,i) = 0.0;
         error  = abs(s(i+1)) /b_norm;
         if ( error <= tol )                      
            y = H(1:i,1:i) \ s(1:i);      
    
            %x = x + V(:,1:i)*y;
            x = x + Z(:,1:i)*y;
           
            break;
         end
     end
     %%%%%%%%%%%%%%%%%%
     iter = iter + m;

     if ( error <= tol )
       
         break
     end

     y = H(1:m,1:m) \ s(1:m); %<---------- Here, done every m times and it is of size m
     %x = x + V(:,1:m)*y;
     x = x + Z(:,1:m)*y; 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
     r =  ( b-A*x );
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
     r_norm = normest(r);
     s(i+1) = r_norm;
     error = s(i+1) / b_norm;  
    
     if ( error <= tol )
        
         break
     end
end

function [cs, sn] = givens_rotation(v1, v2)
    t = sqrt(v1^2 + v2^2);
    cs = v1 / t;  
    sn = v2 / t;
end

end

function [x,count] = S_GMRES( A, x, b, restrt, max_it, tol,L,U,P,Q,D)
iter = 0;                                      
count = 0;
b_norm = normest( b );


if  ( b_norm == 0.0 )
    b_norm = 1.0; 
end
n = length(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r =  Q* (U\ (L\ (P*(D\ (b-A*x))) ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_norm = normest(r);
error = r_norm / b_norm;


                                  % initialize workspace
m = restrt;

V = zeros(n,m+1);
H = zeros(m+1,m);
cs = zeros(m,1);
sn = zeros(m,1);
e1  = zeros(n,1);
e1(1) = 1.0;

while iter <= max_it
     V(:,1) = r / r_norm;
     s = r_norm*e1;
     %%%%%%%%%%%%%%%
     for i = 1:m  
         count = count + 1;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         w = Q* (U\ (L\ (P*(D\ (A*V(:,i))))) ); 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         for k = 1:i
             H(k,i)= w'*V(:,k);
             w = w - H(k,i)*V(:,k);
         end
         H(i+1,i) = normest( w );
         V(:,i+1) = w / H(i+1,i); %<------------------ % just a number
         for k = 1:i-1                             
             temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
             H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
             H(k,i)   = temp;
         end
         [cs(i),sn(i)] = givens_rotation( H(i,i), H(i+1,i) ); 
         temp   = cs(i)*s(i);                        
         s(i+1) = -sn(i)*s(i);
         s(i)   = temp;
         H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
         H(i+1,i) = 0.0;
         error  = abs(s(i+1)) / b_norm;
         if ( error <= tol )                       
            y = H(1:i,1:i) \ s(1:i);      %<------   Here, done once and it is smaller than m
            x = x + V(:,1:i)*y;
            break;
         end
     end
     %%%%%%%%%%%%%%%%%%
     iter = iter + m;

     if ( error <= tol )
       
         break
     end
    
     y = H(1:m,1:m) \ s(1:m); %<---------- Here, done every m times and it is of size m
     x = x + V(:,1:m)*y; % update approximation
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
     r =  Q* (U\ (L\ (P*(D\ (b-A*x))) ));   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      
     r_norm = normest(r);
     s(i+1) = r_norm;
     error = s(i+1) / b_norm;  
    
     if ( error <= tol )
        
         break
     end
end

function [cs, sn] = givens_rotation(v1, v2)
  if (v1 == 0)
    cs = 0;
    sn = 1;
  else
    t = sqrt(v1^2 + v2^2);
    cs = v1 / t;  % see http://www.netlib.org/eispack/comqr.f
    sn = v2 / t;
  end
end
end

