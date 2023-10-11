
clear all;
format short
video_on = 0;

%quadrature = "SBP4_in_time";
%quadrature = "Gauss_Lobatto_4Nodes";
quadrature = "Gauss_4Nodes";
%quadrature = "DI_Gauss_4Nodes"; % need larger step size



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 401 %2.2897e-04
nr_blocks =20 ;
tol = 1e-12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_l=-10;x_r=10;
len = x_r - x_l;                    % Number of gridpoints
n=2*m;
h=(x_r-x_l)/(m-1);

                   % End time 
t_1 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if video_on
    %theAxes=[x_l x_r -2 2]; % unifrom
    theAxes=[x_l x_r 0 10]; % Ma
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SBP6_Higher_2018



x=linspace(x_l,x_r,m)';	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mx = m;
if quadrature == "SBP4_in_time"
    mt = 8;
else
    mt = 4;
end
block = t_1/nr_blocks;
dt = block/(mt-1);
seq = mt*nr_blocks - nr_blocks;

        %%      SBP-SAT     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Ix = speye(mx);
It = speye(mt);

Quadratures_SAT;

e_end = kron(e_mt', Ix);
V = sparse( u_exact(x,0) ) ;
ONE_t = ones(mt,1);
ONE_x = ones(mx,1);


D3 = sparse( kron(It, (D3)) );
D1 = sparse( kron(It, (D1)) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HI_bar = kron(It, HI);

alpha1 = -1; alpha2 = 1; alpha3 = -0.5;





E_l =  sparse(repmat(e_l, mt, 1));
E_r = sparse(repmat(e_r, mt, 1));

D1_r = sparse(repmat(d1_r, 1, mt));




sat1 = alpha1 * HI*e_l   * (d2_l );
sat2 = alpha2 * HI*e_r   * (d2_r );

sat3 = alpha3 * HI*d1_r' * (d1_r );

satt = sat1 + sat2 + sat3;

SAT_L = kron(It,satt);
%SAT_L = sparse(sat1 + sat2 + sat3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aHe_l = sparse(diag(alpha1 * HI_bar*E_l))   ;
aHe_r = sparse(diag(alpha2 * HI_bar*E_r))   ;
aHd_r = sparse(diag(alpha3 * HI_bar*D1_r')) ;


E_lt = sparse(kron(It, e_l'));
E_rt = sparse(kron(It, e_r'));

D1_r = sparse(kron(It, d1_r));
D1_l = sparse(kron(It, d1_l));

D2_r = sparse(kron(It, d2_r));
D2_l = sparse(kron(It, d2_l));






%
A =  sparse( Dt_bar + D3 - SAT_L  ) ;
I = speye(size(A));
[l,u,p,q,d] = lu(A);
pd = sparse(  p*( (1./diag(d)).*I )      );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch quadrature
    case "SBP4_in_time"
        TT = 0:dt:block;
    case "Gauss_Lobatto_4Nodes"
        m = (block)/2; c = (block)/2;
        TT = m*T_GL+c;
    case "Gauss_4Nodes"
         m = (block)/2; c = (block)/2;
         TT = m*T_GL+c;
    case "DI_Gauss_4Nodes"
        TT = (block)*T_GL ;
    otherwise
        disp("wrong name")
end
TT = TT';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=0.0;
t_start = cputime;
while t<  t_1 -block/2
    
    norm_du = inf;

    u0 = kron(ONE_t,V); 
    du = ones(size(u0));
    U0 = u0; 
 
    [g1, g2, g3] = Boundary_data(x,TT);
   
       
     
 
    while norm_du > tol
    
        sat1 = D2_l*u0  + 2*((E_lt*u0)+abs(E_lt*u0)).*(E_lt*u0) - g1;
        sat2 = D2_r*u0  + 2*((E_rt*u0)-abs(E_rt*u0)).*(E_rt*u0) - g2;
        sat3 = D1_r*u0 - g3;
       % Ur = kron(E_rt*u0,ONE_x);

        SAT1 = aHe_l*  kron( sat1   ,ONE_x  );
        SAT2 = aHe_r*  kron( sat2   ,ONE_x  );
        SAT3 = aHd_r*  kron( sat3   ,ONE_x  );

       
        b = -R*U0  - (Dt_bar + D3)*(u0 ) - 6*diag(u0)*(D1*(u0)) + SAT1 + SAT2 + SAT3 ;
  
        du = q* (u\ (l\ (pd* (b))) );
        u0 =  (u0 + du) ;
       
        norm_du = abs(normest(du));
    
      
    end

   
    V = e_end*u0;
    t = t+block;
    TT = TT + block;

    if video_on 

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
%
function  [g1, g2, g3] = Boundary_data(x,t)

Vl = 1/2*( sech(1/2*(x(1)-t)) ).^2;
Vr = 1/2*( sech(1/2*(x(end)-t)) ).^2;
Vr_x = sinh(t./2 - x(end)/2)./(2*cosh(t./2 - x(end)/2).^3);
Vl_xx = (3*sinh(t./2 - x(1)/2).^2)./(4*cosh(t./2 - x(1)/2).^4) - 1./(4*cosh(t./2 - x(1)/2).^2);
Vr_xx = (3*sinh(t./2 - x(end)/2).^2)./(4*cosh(t./2 - x(end)/2).^4) - 1./(4*cosh(t./2 - x(end)/2).^2);

g1 = sparse(Vl_xx +2*((Vl)+abs(Vl)).*(Vl) );
g2 = sparse(Vr_xx +2*((Vr)-abs(Vr)).*(Vr) );
g3 = sparse(Vr_x);



end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
