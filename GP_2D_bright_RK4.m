%{
NLS 2D with RK4 as an time integrator
%}
clear all;
video_on = 0;


t_1= 1;
m = 25%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CFL =0.1;
k=CFL*h^2;

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
P = sparse(Px*Py);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BCx = sparse( HI*Lx'*(inv(Lx*HI*Lx')) );
BCy = sparse( HI*Ly'*(inv(Ly*HI*Ly')) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC = ones(n,1);
CC(1)=0.5; CC(m)=0.5; CC(end)=0.5; CC(end-m+1)=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=linspace(x_l,x_r,m);	% discrete x values
[X,Y] = meshgrid(x,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = zeros(n,1);        
V(1:n) = reshape(u_exact(X,Y,0),n,1);

V = sparse(V);
V_d = sparse(reshape(Vd(X,Y),n,1).*In);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = 0;
t=0.0;
t_start = cputime;
BC = 0;
while t< t_1 
    % stage1
    W = V ;
  
    w1 =  1i/2*(Dxx + Dyy)*W + 1i*F(W) + 1i*V_d*W;   
    %______________________________________________________________________
    %stage2
    BC =  boundary(x',t+k/2,BCx,BCy,CC);

    W = P*(V+k/2*w1) +BC;
    w2 =  1i/2*(Dxx + Dyy)*W + 1i*F(W) + 1i*V_d*W;     
    %______________________________________________________________________
    %stage3
    BC =  boundary(x',t+k/2,BCx,BCy,CC);
 
    W = P*(V+k/2*w2) +BC;
    w3 =  1i/2*(Dxx + Dyy)*W + 1i*F(W) + 1i*V_d*W;   
   
    %______________________________________________________________________
    %stage 4
    BC =  boundary(x',t+k,BCx,BCy,CC);
    
    W = P*(V+k*w3) +BC;
    w4 =   1i/2*(Dxx + Dyy)*W + 1i*F(W) + 1i*V_d*W;      
    %_____________________________________________________________________
    V= P*(V+k/6*(w1+2*w2+2*w3+w4)) +BC;
        
    t=t+k;
    if video_on && (t>=0)&&(mod(freq,20)==0)
    
        surf(X,Y,reshape(real(V),m,m)); 
      %  surf(X,Y,reshape(u_exact(X,Y,t),m,m) )
        shading interp
       % view(90,0)
        colorbar
        title(['Numerical solution at t = ',num2str(t)]);
        axis(theAxes);
        grid;xlabel('x');
        legend('u')
        ax = gca;         
        ax.FontSize = 16;
        currFrame = getframe;
       % writeVideo(vidObj,currFrame);
    end
    freq = freq + 1;

end
t_end = cputime - t_start

U_approx = reshape(V,n,1);
U_exact = reshape(u_exact(X,Y,t),n,1);

err = sqrt(real(U_exact)-real(U_approx))'*H*sqrt(real(U_exact)-real(U_approx))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ui = u_exact(x,y,t)
xei = 1; k = 1;

ui = sqrt(xei) * exp(-k/2* (x.^2 + y.^2))*exp(-1i*k*t);
end

function b = boundary(x,t,BCx,BCy,C)
xei = 1; k = 1;

West = sqrt(xei) * exp(-k/2* (x(1).^2 + x.^2))*exp(-1i*k*t); 
East = sqrt(xei) * exp(-k/2* (x(end).^2 + x.^2))*exp(-1i*k*t);


South = sqrt(xei) * exp(-k/2* (x.^2 + x(1).^2))*exp(-1i*k*t);
North = sqrt(xei) * exp(-k/2* (x.^2 + x(end).^2))*exp(-1i*k*t); 

gx =  [West ; East];
gy =  [South; North];
size(BCx)
size(gx)
size(C)

b = (BCx*gx + BCy*gy).*C;
end

function f = F(u)
f =  (abs(u).^2).*u ;
end

function v = Vd(x,y)
xei = 1; k = 1;
v = -k/2 * (x.^2 + y.^2) - xei*exp(-k*(x.^2 + y.^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
