
%{
two stationary solitons 1D  with RK4
from
Numerical comparison of mass-conservative schemes for the
Gross-Pitaevskii equation ∗
and

Exact Solutions to the Nonlinear Schrodinger Equation

This is strongly nonlinear equation
%}



clear all;
format short
video_on = 0;


%quadrature = "SBP4_in_time";
quadrature = "Gauss_Lobatto_4Nodes";
%quadrature = "Gauss_4Nodes";
%quadrature = "DI_Gauss_4Nodes"; % need larger step size

% Parameters
m = 401
nr_blocks = 200;
restrt=10;
max_it=20;
tol_gmres = 1e-2; % GMRES tolerance
tol = 1e-12;  % upper iteration tolarence         
%7.8598e-04


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_l=-2;x_r=2;
len = x_r - x_l;                   
n=2*m;
h=(x_r-x_l)/(m-1);

t_1 = 1;
CFL = 0.1;  
k=CFL*h^2; % this is to check how it will be with RK4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if video_on
    theAxes=[x_l x_r -4 5]; 
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SBP4;
L = [e_1'; e_m'];
I = eye(m);
P= I-HI*L'*((L*HI*L')\L) ;
g = zeros(2,1);

x=linspace(x_l,x_r,m)';	
V = sparse( u_exact(x,0) ) ;  
t=0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         %  Number of blocks
if quadrature == "SBP4_in_time"
    mt = 8;
else
    mt = 4;
end
block = t_1/nr_blocks;  % block length
block_ = round(block,10); % this to avoind irrational number
t_1 = block_*nr_blocks;
block = block_;

dt = block/(mt-1);      % time step in case of using SBP4 in time as an integrator
seg = mt*nr_blocks - nr_blocks; % number of time intervalls
mt_total = seg + 1;     % total number of nodes
mx = m;



%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ix = speye(mx);
It = speye(mt);
I = speye(mt*mx);


Quadratures_projection;
e_end = kron(e_mt', Ix);

D2_bar = kron(It,sparse(D2*P));
P_bar =  kron(It,sparse(P));
D2 = kron(It,sparse(D2));
BC_ = sparse(HI*L'*(inv(L*HI*L')) );
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
A =  sparse(I - JF*(1i*D2_bar  ));

%GMRES
[l,u,p,q,d] = lu(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time = zeros(nr_blocks+1,1);
%Space = x;
%Solution = zeros(nr_blocks+1,mx);
%Solution(1,:) = V';
%num = 1;
du = sparse(zeros(mx*mt,1));
t_start = cputime;
while t<  t_1-block/2
 
    norm_du = inf;
    u0 = kron(ONE,V);
    U0 = u0;

    BC = reshape(Boundary(x,BC_,T),mt*mx,1 );

   
   while norm_du > tol

      
        u_= P_bar*u0 + BC;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % A =  (I - JF*(1i*D2_bar +2i*diag(jac(u_))  ));
        %[l,u] = ilu(A);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        %
                
        b =  U0 - (u0) + JF*( 1i*D2*u_ + 1i*2*F(u_)  )  ; 
        %%%%%%%%% GMRES
       % [du,count] = GMRES_main( A, du, b, restrt, max_it, tol_gmres,l,u,p,q,d);%pqd

        %%%%%%%%% F-GMRES
       % [du,count] = F_GMRES( A, du, b, restrt, max_it, tol_gmres,l,u,p,q,d);
     

        %%%%%%%%% without iterative solver
        du = q* (u\ (l\ (p*(d\(b)))));

        u0 = ( u0 + du );
        norm_du = normest(du);

        
    end
  
    u0 = P_bar*u0 + BC;
    U0 = u0;
    V = e_end*u0;
    t = t+block;
    T = T+block;

    %{
    num = num+1
    Time(num,1) = t;
    Solution(num,:)=V';
    %}
    
    if video_on 
          
        plot(x,real(V),'b','LineWidth',1);
       % hold on
        %plot(x,real(u_exact(x,t)),'r','LineWidth',1);

        title(['Numerical solution at t = ',num2str(t)]);
        axis(theAxes);
        grid;xlabel('x');
        legend('V')
        ax = gca;          
        ax.FontSize = 16;
        currFrame = getframe;
    end

  


end

%count_tot
t_end = cputime - t_start

error = sqrt( real(u_exact(x,t)) - real(V) )' * H * sqrt( real(u_exact(x,t))-real(V) )






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = u_exact(x,t)
    
        u = ( 8*exp(4i*t) * (9*exp(-4*x) + 16*exp(4*x)) - 32*exp(16i*t) * (4*exp(-2*x) + 9*exp(2*x)) ) ...
            ./ ( -128*cos(12*t) + 4*exp(-6*x) + 16*exp(6*x) + 81*exp(-2*x) + 64*exp(2*x));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BC = Boundary(x,BC_,t)



    u_left = ( 8*exp(4i*t) * (9*exp(-4*x(1)) + 16*exp(4*x(1))) - 32*exp(16i*t) * (4*exp(-2*x(1)) + 9*exp(2*x(1))) ) ...
            ./ ( -128*cos(12*t) + 4*exp(-6*x(1)) + 16*exp(6*x(1)) + 81*exp(-2*x(1)) + 64*exp(2*x(1)));
    u_right = ( 8*exp(4i*t) * (9*exp(-4*x(end)) + 16*exp(4*x(end))) - 32*exp(16i*t) * (4*exp(-2*x(end)) + 9*exp(2*x(end))) ) ...
            ./ ( -128*cos(12*t) + 4*exp(-6*x(end)) + 16*exp(6*x(end)) + 81*exp(-2*x(end)) + 64*exp(2*x(end)));
 
g = ([u_left; u_right]);
BC = BC_*g;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = F(u)
f =  ((abs(u).^2).*u) ;
end

function j = jac(u)
j = abs(u).^2 + (u.*abs(u).*(u + conj(u)))./(u.*conj(u)).^(1/2);
j = sparse(j);
end
function j = jac2(u)
j = (u + conj(u)).^2./(2*conj(u)) +...
    (2*u.*abs(u))./(u.*conj(u)).^(1/2) + ...
    (2*abs(u).*(u + conj(u)))./(u.*conj(u)).^(1/2) ...
    - (u.*abs(u).*(u + conj(u)).^2)./(2*(u.*conj(u)).^(3/2));
 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GMRE
function [x,count] = GMRES_main( A, x, b, restrt, max_it, tol,L,U,P,Q,D)
iter = 0;                                         
count = 0;
bnrm2 = normest( b );

if  ( bnrm2 == 0.0 )
    disp("This might not be used")
    bnrm2 = 1.0; 
end

n = length(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r =  Q* (U\ (L\ (P* (D\ (b-A*x)))) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_norm = normest(r);
error = r_norm / bnrm2;

if ( error < tol ) 
  
    disp("This is statement might not be used")
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
   
         [z, count_] = S_GMRES( A, z, V(:,i),  4,4, 1e-3,L,U,P,Q,D);
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






