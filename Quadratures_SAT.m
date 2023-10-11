%{

Quadrature for SBPSAT in time
This is the main version to be used
Does not include the non quadrature version
From the paper 
    <<  HIGH-ORDER IMPLICIT TIME-MARCHING METHODS BASED ON
GENERALIZED SUMMATION-BY-PARTS OPERATORS∗>>

Date: 3/03/2023

mt: number of temporal nodes, here it hsould be 4 for all of them
Ix: should be prefefiend which is an identity matrix of size mx*mx
and mx is the number of spacial nodes
block: should be predefined and it is the length of the time domain to
be integrated

%}


switch quadrature
    case "SBP4_in_time"
        T_GL = 0;

        Ht=diag(ones(mt,1),0);
        Ht(1:4,1:4)=diag([17/48 59/48 43/48 49/48]);
        Ht(mt-3:mt,mt-3:mt)=fliplr(flipud(diag([17/48 59/48 43/48 49/48])));
        Ht=Ht*dt;
        HIt=inv(Ht);

   
        Qt=(-1/12*diag(ones(mt-2,1),2)+8/12*diag(ones(mt-1,1),1)-8/12*diag(ones(mt-1,1),-1)+1/12*diag(ones(mt-2,1),-2));
        Q_Ut = [0 0.59e2 / 0.96e2 -0.1e1 / 0.12e2 -0.1e1 / 0.32e2; -0.59e2 / 0.96e2 0 0.59e2 / 0.96e2 0; 0.1e1 / 0.12e2 -0.59e2 / 0.96e2 0 0.59e2 / 0.96e2; 0.1e1 / 0.32e2 0 -0.59e2 / 0.96e2 0;];
        Qt(1:4,1:4)=Q_Ut;
        Qt(mt-3:mt,mt-3:mt)=flipud( fliplr(-Q_Ut(1:4,1:4) ) );

        e_1t=zeros(mt,1);e_1t(1)=1;
        e_mt=zeros(mt,1);e_mt(mt)=1;

        D1t=HIt*(Qt -1/2*e_1t*e_1t'+1/2*e_mt*e_mt') ;

        alpha  = -1;

        Dt_bar = D1t - alpha*HIt*(e_1t*e_1t');
        Dt_bar = sparse( kron(Dt_bar, Ix) );
        R = sparse( kron(diag(alpha*HIt*e_1t), Ix ));
   
    case "Gauss_Lobatto_4Nodes" % [-1, +1]
        e_1t=zeros(mt,1);e_1t(1) = 1;
        e_mt=zeros(mt,1);e_mt(mt) = 1;

        T_GL = [-1 -1/5*sqrt(5) 1/5*sqrt(5) 1];

        Ht = eye(mt)*5; 
        Ht(1,1) = 1; Ht(end,end) = 1; 
        Ht = (0.5*block)/6 * Ht;
        HIt = inv(Ht);

        D1t  =1/(0.5*block)*...
                 [-3                   -(5*sqrt(5))/(sqrt(5)-5) -(5*sqrt(5))/(sqrt(5)+5)   1/2;...
                (sqrt(5))/(sqrt(5)-5)   0                      sqrt(5)/2               -(sqrt(5))/(sqrt(5)+5);...
                (sqrt(5))/(sqrt(5)+5)  -sqrt(5)/2               0                   -(sqrt(5))/(sqrt(5)-5);...
                -1/2                 (5*sqrt(5))/(sqrt(5)+5) (5*sqrt(5))/(sqrt(5)-5)            3];


        Qt=Ht*D1t ;

        alpha  = -1;

        E0 = diag(e_1t);
        Dt_bar = D1t - alpha*HIt*E0; 
        
      
        Dt_bar = sparse( kron(Dt_bar, Ix) );

        R = sparse( kron(alpha*HIt*E0, Ix) );

    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case "DI_Gauss_4Nodes" % [0, +1]
       
        e_1t = [0.8808689243587871 0.9884420520048577 -0.6011474168414327 -0.2681635595222120]';
        e_mt = [0.9928785357819795 -0.4986129934126102 0.4691078563418350 0.03662660128879568]';

        T_GL = [0.5975501145870646 0.1236947892666459 0.9813648784844768 0.2188347157850838];
        
        
        Ht = (block)*diag([0.5263633266867775 0.3002573924935185 0.1447678514141155 0.0286114294055885]);
        HIt = inv(Ht);

        D1t  =1/(block)*...
                [0.1993658318073258 -1.654157580888287 1.006020084619771 0.4487716644611903;...
                -1.648792506689303 -1.212963928918776 1.978966716941006 0.8827897186670728;...
                3.217338082860363 -1.615712813301921 -0.4880781006041668 -1.113547168954275;...
                1.271022350640990 -0.6382938457303877 0.6005231745715582 -1.233251679482160];


        Qt=Ht*D1t ; 
        %D1t=HIt*(Qt -1/2*e_1t*e_1t'+1/2*e_mt*e_mt') ;

        alpha  = -1;

        Dt_bar = D1t - alpha*HIt*(e_1t*e_1t');
        Dt_bar = sparse( kron(Dt_bar,Ix) );
        R = sparse( kron(diag(alpha*HIt*e_1t), Ix ));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case "Gauss_4Nodes" % [-1, +1]
      
        e_1t = [1.5267881254572668 -0.8136324494869273 0.4007615203116504 -0.1139171962819899]';
       
        e_mt = [-0.1139171962819899 0.4007615203116504 -0.8136324494869273 1.5267881254572668]';
      

        T_GL = [-0.8611363115940526 -0.3399810435848563 0.3399810435848563 0.8611363115940526];

        Ht = (0.5*block)*diag([0.3478548451374539 0.6521451548625461 0.6521451548625461 0.3478548451374539]);
        HIt = inv(Ht);

        D1t  =1/(0.5*block)*...
                [-3.3320002363522817  4.8601544156851962   -2.1087823484951789   0.5806281691622644;...
                -0.7575576147992339   -0.3844143922232086   1.4706702312807167   -0.3286982242582743;...
                0.3286982242582743   -1.4706702312807167    0.3844143922232086   0.7575576147992339;...
                -0.5806281691622644   2.1087823484951789    -4.8601544156851962   3.3320002363522817];

                 

        Qt=Ht*D1t ;
    
        alpha  = -1;
     
      
        Dt_bar = D1t - alpha*HIt*(e_1t*e_1t');
        Dt_bar = sparse( kron(Dt_bar, Ix) );
        R = sparse( kron(diag(alpha*HIt*e_1t), Ix ));
 
    otherwise
        disp("Wrong name")
end