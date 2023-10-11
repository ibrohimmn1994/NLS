clear all;
format short
video_on = 0;

%quadrature = "SBP4_in_time";
%quadrature = "Gauss_Lobatto_4Nodes";
quadrature = "Gauss_4Nodes";
%quadrature = "DI_Gauss_4Nodes"; % need larger step size



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m =401
nr_blocks =20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_l=-10;x_r=10;
len = x_r - x_l;                    % Number of gridpoints
n=2*m;
h=(x_r-x_l)/(m-1);

                   % End time 
t_1 = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if video_on
 
    theAxes=[x_l x_r 0 10]; 
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SBP6_Higher_2018

L = [e_l';e_r';d1_r];

I = eye(m);
P = I-HI*L'*((L*HI*L')\L) ;



x=linspace(x_l,x_r,m)';	

mx = m;

if quadrature == "SBP4_in_time"
    mt = 8;
else
    mt = 4;
end
block = t_1/nr_blocks;
dt = block/(mt-1);
seq = mt*nr_blocks - nr_blocks;

tol = 1e-12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%      SBP-SAT     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Ix = speye(mx);
It = speye(mt);

Quadratures_projection;

e_end = kron(e_mt', Ix);


V = sparse( u_exact(x,0) ) ;
ONE = ones(mt,1);

D3_bar = sparse( kron(It, (D3)*P) );
D3 = sparse( kron(It, (D3)) );
D1 = sparse( kron(It, D1) );
P_bar= sparse( kron(It, P) );
BC = sparse( HI*L'*(inv(L*HI*L')) );


I = speye(size(D3));
A =  sparse(I - JF*(-D3_bar ));

[l,u,p,q,d] = lu(A);
pd = sparse(  p*( (1./diag(d)).*I )      );



t =0;
t_start = cputime;
while t< t_1 -block/2
    
    norm_du = inf;
   
    % BCT is the boundary data
    BCT = reshape(Boundary_GL(x,t,t+block,BC,T_GL,quadrature,dt),mt*mx,1 );
 


    u0 = kron(ONE,V); 
    du = ones(size(u0));
    U0 = u0; 

 
   
    while norm_du > tol
    
        u_ = P_bar*u0+BCT;
        
       % b = -R*U0  - Dt_bar*u0 - ( D3)*(u_) - 6*diag(u_)*(D1*u_);

        b =  U0 - (u0) + JF*( -D3*u_ - 6*diag(u_)*(D1*u_))  ;

        du = q* (u\ (l\ (pd* (b))) );
        u0 =  (u0 + du);
 
       norm_du = norm(abs(du));
      
    end

    u0 = P_bar*u0 +BCT;
    V = e_end*u0;

    t = t+block;


    if video_on 
   
        figure(1)
        plot(x,real(V),'b','LineWidth',1);
   
        title(['Numerical solution at t = ',num2str(t)]);
        axis(theAxes);
        grid;xlabel('x');
        legend('V')
        ax = gca;          % current axes
        ax.FontSize = 16;
        currFrame = getframe;
    end
   

end

t_end = cputime - t_start

U_exact = (u_exact(x,t));
error = sqrt(U_exact-V)' * H * sqrt(U_exact-V)









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = u_exact(x,t)
u = 1/2*( sech(1/2*(x-t)) ).^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = Boundary_GL(x,t_start, t_end,BC,T_GL,name,dt)


switch name
    case "SBP4_in_time"
        t = t_start:dt:t_end;
    case "Gauss_Lobatto_4Nodes"
        m = (t_end - t_start)/2; c = (t_end + t_start)/2;
        t = m*T_GL+c;
    case "Gauss_4Nodes"
         m = (t_end - t_start)/2; c = (t_end + t_start)/2;
         t = m*T_GL+c;
    case "DI_Gauss_4Nodes"
        t = (t_end-t_start)*T_GL + t_start;
    otherwise
        disp("wrong name")
end

Vl = 1/2*( sech(1/2*(x(1)-t)) ).^2;
Vr = 1/2*( sech(1/2*(x(end)-t)) ).^2;

Vr_x = sinh(t./2 - x(end)/2)./(2*cosh(t./2 - x(end)/2).^3);

g =  sparse([Vl; Vr ; Vr_x]);
b = BC*g;

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

