clear all;
format short
video_on = 0;

m = 51;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_l=-10;x_r=10;
len = x_r - x_l;                    % Number of gridpoints
n=2*m;
h=(x_r-x_l)/(m-1);

                   % End time 
t_1 = 1.0;
CFL=0.1;  
k=CFL*h^3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if video_on   
    theAxes=[x_l x_r 0 10];
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SBP6_Higher_2018


L = [e_l';e_r'; d1_r];

I = eye(m);
P= I-HI*L'*((L*HI*L')\L) ;

x=linspace(x_l,x_r,m)';	
V = sparse( u_exact(x,0) ) ;  
t=0.0;

        %%      RK4     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = sparse(P);
D3 = sparse(D3);
D1 = sparse(D1);
BC_ = sparse(HI*L'*inv(L*HI*L'));

freq = 0;

t_Start = cputime;
while t < t_1
    %_____________Stage one   
  
   
  
    W = V ;
    w1 =   -D3*W - 6*diag(W)*(D1*W) ;  
    %_____________Stage two

    BC = Boundary(BC_, x, t+k/2);
  
    W = P*(V +k/2*w1)+BC;
    w2 =    -D3*W - 6*diag(W)*(D1*W) ;  
    %_____________Stage three
    BC = Boundary(BC_, x, t+k/2);
  
    W = P*( V +k/2*w2 )+BC;
    w3 =    -D3*W - 6*diag(W)*(D1*W)  ;  
     %_____________Stage four
    BC = Boundary(BC_, x, t+k);
  
    W = P*( V + k*w3 )+BC ;
    w4 =   -D3*W - 6*diag(W)*(D1*W) ;  
 
 
    V = P*(V + k/6*(w1 + 2*w2 + 2*w3 + w4) ) +BC ; 
    t = t+k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if video_on &&(mod(freq,100)==0)
 
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

function BC = Boundary(BC_,x,t)
Vl = 1/2*( sech(1/2*(x(1)-t)) ).^2;
Vr = 1/2*( sech(1/2*(x(end)-t)) ).^2;
Vr_x = sinh(t/2 - x(end)/2)/(2*cosh(t/2 - x(end)/2)^3);


g = sparse([ Vl; Vr ; Vr_x]);
BC = BC_*g;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
