%{
Bright soliton 1D SBP4 with RK4 time integrato
%}

clear all;
format short
video_on = 0;
%name = "dark_II";
%name = "dark_I";
name = "bright";


 % Number of gridpoints
m = 101
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x_l=-4;x_r=12;  % dark_II
%x_l=-15;x_r=30;  % dark_I
%x_l=-15;x_r=20;  % bright origional
x_l=-3;x_r=3;
len = x_r - x_l;                   
n=2*m;
h=(x_r-x_l)/(m-1);

%t_1= 5.0;                    % End time origional bright
t_1 = 1;
CFL=0.5;  
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
L = [d_1; d_m];
I = speye(m);
P = sparse( I-HI*L'*((L*HI*L')\L) );
BC_ = sparse( HI*L'*(inv(L*HI*L')) );
D2 = sparse(D2);

x = linspace(x_l,x_r,m)';	
V = sparse( u_exact(x,0,name) ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0.0;
freq = 0;
beta = -1; % bright
%beta = +1; % dark
        %%      RK4     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BC = 0;
t_start = cputime;
while t < t_1
    %_____________stage one   
    W = V ;
    w1 =  0.5*1i*D2*W - beta* 1i*F(W)  ;
    %_____________stage two
    BC = boundary(x,t+k/2,name, BC_);
    W = P*(V+k/2*w1)+BC;
    w2 =   0.5*1i*D2*W - beta* 1i*F(W) ;
    %____________stage three
    BC = boundary(x,t+k/2,name, BC_);
    W = P*(V+k/2*w2)+BC;
    w3 =  0.5*1i*D2*W - beta* 1i*F(W)  ;
    %_____________stage four
    BC = boundary(x,t+k,name, BC_);
    
    W = P*(V+k*w3)+BC;
    w4 =    0.5*1i*D2*W - beta* 1i*F(W)  ;
 
    V=P*(V+k/6*(w1+2*w2+2*w3+w4) ) +BC ; 
    t=t+k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if video_on && (mod(freq,10)==0)
        t ;    
        plot(x,real(V),'b','LineWidth',1);
       % plot(x,real(u_exact(x,t,name)),'r','LineWidth',1);
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
        u = a*sec\h( a*(x-v*t-x0) ).*exp( 1i*(v*x-1/2*(v^2-a^2)*t + Q0) );

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
