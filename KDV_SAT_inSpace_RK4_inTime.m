%{
On the numerical solution of Korteweg–de Vries equation
by the iterative splitting method
%}


clear all;
format short
video_on = 0;

m = 401;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_l=-10;x_r=10;
len = x_r - x_l;                    % Number of gridpoints
n=2*m;
h=(x_r-x_l)/(m-1);

                   % End time 
t_1 = 1;
CFL=0.1;  
k=CFL*h^3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if video_on
    %theAxes=[x_l x_r -2 2]; % unifrom
    theAxes=[x_l x_r 0 10]; % Ma
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBP4;
SBP6_Higher_2018



x=linspace(x_l,x_r,m)';	
V = sparse( u_exact(x,0) ) ;  
t=0.0;

        %%      RK4     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D3 = sparse(D3);
D1 = sparse(D1);
HI = sparse(HI);
e_l = sparse(e_l); e_r = sparse(e_r); d2_l=sparse(d2_l);d2_r = sparse(d2_r);
d1_r = sparse(d1_r);

alpha1 = -1; alpha2 = 1; alpha3 = -1/2;

freq = 0;

t_Start = cputime;
while t < t_1
    %_____________ Stage one   
    [Vl,Vr, Vr_x,Vl_xx,Vr_xx] = Boundary_data(x,t);
    W = V ;

    SAT1 = alpha1 * HI*e_l   * (d2_l*W  + 2*((e_l'*W)+abs(e_l'*W))*(e_l'*W)...
                                          - (Vl_xx +2*((Vl)+abs(Vl))*(Vl) ));
    SAT2 = alpha2 * HI*e_r   * (d2_r*W  + 2*((e_r'*W)-abs(e_r'*W))*(e_r'*W)...
                                          - (Vr_xx +2*((Vr)-abs(Vr))*(Vr) ));
    SAT3 = alpha3 * HI*d1_r' * (d1_r*W  - Vr_x);
    SAT = SAT1 + SAT2 + SAT3;

    w1 = -D3*W - 6*diag(W)*(D1*W) + SAT ;
   
    %_____________ Stage two
    [Vl,Vr, Vr_x,Vl_xx,Vr_xx] = Boundary_data(x,t+k/2);
    W = V + k/2*w1 ;

    SAT1 = alpha1 * HI*e_l   * (d2_l*W  + 2*((e_l'*W)+abs(e_l'*W))*(e_l'*W) ...
                                        - (Vl_xx +2*((Vl)+abs(Vl))*(Vl) ) );
    SAT2 = alpha2 * HI*e_r   * (d2_r*W  + 2*((e_r'*W)-abs(e_r'*W))*(e_r'*W) ...
                                        - (Vr_xx +2*((Vr)-abs(Vr))*(Vr) ));
    SAT3 = alpha3 * HI*d1_r' * (d1_r*W  - Vr_x);
    SAT = SAT1 + SAT2 + SAT3;

    w2 = -D3*W - 6*diag(W)*(D1*W) + SAT ;
   %_____________Stage three
    [Vl,Vr, Vr_x,Vl_xx,Vr_xx] = Boundary_data(x,t+k/2);
    W = V + k/2*w2 ;

    SAT1 = alpha1 * HI*e_l   * (d2_l*W  + 2*((e_l'*W)+abs(e_l'*W))*(e_l'*W) ...
                                        - (Vl_xx +2*((Vl)+abs(Vl))*(Vl) ) );
    SAT2 = alpha2 * HI*e_r   * (d2_r*W  + 2*((e_r'*W)-abs(e_r'*W))*(e_r'*W) ...
                                        - (Vr_xx +2*((Vr)-abs(Vr))*(Vr) ) );
    SAT3 = alpha3 * HI*d1_r' * (d1_r*W  - Vr_x);
    SAT = SAT1 + SAT2 + SAT3;

    w3 = -D3*W - 6*diag(W)*(D1*W) + SAT ;
    %_____________Stage four
    [Vl,Vr, Vr_x,Vl_xx,Vr_xx] = Boundary_data(x,t+k);
    W = V +k*w3 ;

    SAT1 = alpha1 * HI*e_l   * (d2_l*W  + 2*((e_l'*W)+abs(e_l'*W))*(e_l'*W) ...
                                        - (Vl_xx +2*((Vl)+abs(Vl))*(Vl) ) );
    SAT2 = alpha2 * HI*e_r   * (d2_r*W  + 2*((e_r'*W)-abs(e_r'*W))*(e_r'*W) ...
                                        - (Vr_xx +2*((Vr)-abs(Vr))*(Vr) ) );
    SAT3 = alpha3 * HI*d1_r' * (d1_r*W  - Vr_x);
    SAT = SAT1 + SAT2 + SAT3;

    w4 = -D3*W - 6*diag(W)*(D1*W) + SAT ;
    %_______________________________________
  
 
    V =  V + k/6*(w1 + 2*w2 + 2*w3 + w4)  ; 
    t = t+k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if video_on && (t>=0.0)&&(mod(freq,100)==0)
 
        plot(x,real(V),'b',x,real(u_exact(x,t)),'r','LineWidth',1);
        title(['Numerical solution at t = ',num2str(t)]);
        axis(theAxes);
        grid;xlabel('x');
        legend('V')
        ax = gca;          % current axes
        ax.FontSize = 16;
        currFrame = getframe;
    end
    freq = freq + 1;
end

t_End = cputime - t_Start
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U_exact = (u_exact(x,t));
error = sqrt(U_exact-real(V))' *H* sqrt(U_exact-real(V))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = u_exact(x,t)
u = 1/2*( sech(1/2*(x-t)) ).^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [Vl, Vr, Vr_x, Vl_xx, Vr_xx] = Boundary_data(x,t)

Vl = 1/2*( sech(1/2*(x(1)-t)) )^2;
Vr = 1/2*( sech(1/2*(x(end)-t)) )^2;
Vr_x = sinh(t/2 - x(end)/2)/(2*cosh(t/2 - x(end)/2)^3);
Vl_xx = (3*sinh(t/2 - x(1)/2)^2)/(4*cosh(t/2 - x(1)/2)^4) - 1/(4*cosh(t/2 - x(1)/2)^2);
Vr_xx = (3*sinh(t/2 - x(end)/2)^2)/(4*cosh(t/2 - x(end)/2)^4) - 1/(4*cosh(t/2 - x(end)/2)^2);
end
