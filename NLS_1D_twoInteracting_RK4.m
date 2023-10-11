%{
two stationary solitons 1D  with RK4
from
Numerical comparison of mass-conservative schemes for the
Gross-Pitaevskii equation ∗
and

Exact Solutions to the Nonlinear Schrodinger Equation
%}



clear all;
format short
video_on = 0;


 % Number of gridpoints
m = 801

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_l=-2;x_r=2;
len = x_r - x_l;                   
n=2*m;
h=(x_r-x_l)/(m-1);

                   % End time origional bright
t_1 = 1;
CFL = 0.1;  
k=CFL*h^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if video_on
    %theAxes=[x_l x_r -2 2]; % unifrom
    theAxes=[x_l x_r -2 2]; % Ma
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SBP4;
L = [e_1'; e_m'];
I = speye(m);
P = sparse( I-HI*L'*((L*HI*L')\L) );
BC_ = sparse( HI*L'*(inv(L*HI*L')) );
D2 = sparse(D2);

x = linspace(x_l,x_r,m)';	
V = sparse( u_exact(x,0) ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0.0;
freq = 0;


        %%      RK4     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BC = 0;
t_start = cputime;
while t < t_1
    %_____________stage one   
    W = V ;
    w1 =  1i*D2*W + 1i*2*F(W)    ;
    %_____________stage two
    BC = boundary(x,t+k/2, BC_);
    W = P*(V+k/2*w1)+BC;
    w2 =  1i*D2*W + 1i*2*F(W)    ;  
    %____________stage three
    BC = boundary(x,t+k/2, BC_);
    W = P*(V+k/2*w2)+BC;
    w3 =  1i*D2*W + 1i*2*F(W)    ;
    %_____________stage four
    BC = boundary(x,t+k, BC_);
    
    W = P*(V+k*w3)+BC;
    w4 =    1i*D2*W + 1i*2*F(W)    ; 

    V=P*(V + k/6*(w1+2*w2+2*w3+w4)   ) +BC ;  %+ k*source(x,t+k,eps) 
    t=t+k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if video_on && (mod(freq,100)==0)
        t ;    
       % plot(x,real(V),'b','LineWidth',1);
        plot(x,real(u_exact(x,t)),'r','LineWidth',1);
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_end = cputime - t_start

error = sqrt( real(u_exact(x,t)) - real(V) )'*H*sqrt( real(u_exact(x,t))-real(V) )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = u_exact(x,t)
    
        u = ( 8*exp(4i*t) * (9*exp(-4*x) + 16*exp(4*x)) - 32*exp(16i*t) * (4*exp(-2*x) + 9*exp(2*x)) ) ...
            ./ ( -128*cos(12*t) + 4*exp(-6*x) + 16*exp(6*x) + 81*exp(-2*x) + 64*exp(2*x));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BC = boundary(x,t, BC_)

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













