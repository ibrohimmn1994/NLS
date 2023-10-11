switch quadrature


    case "SBP4_in_time"
        T_GL = 0;

        Ht=diag(ones(mt,1),0);
        Ht(1:4,1:4)=diag([17/48 59/48 43/48 49/48]);
        Ht(mt-3:mt,mt-3:mt)=fliplr(flipud(diag([17/48 59/48 43/48 49/48])));
        Ht=Ht*dt;
        HIt=inv(Ht);

  
        Qt=(-1/12*diag(ones(mt-2,1),2)+8/12*diag(ones(mt-1,1),1)-8/12*diag(ones(mt-1,1),-1)+1/12*diag(ones(mt-2,1),-2));
        Q_U = [0 0.59e2 / 0.96e2 -0.1e1 / 0.12e2 -0.1e1 / 0.32e2; -0.59e2 / 0.96e2 0 0.59e2 / 0.96e2 0; 0.1e1 / 0.12e2 -0.59e2 / 0.96e2 0 0.59e2 / 0.96e2; 0.1e1 / 0.32e2 0 -0.59e2 / 0.96e2 0;];
        Qt(1:4,1:4)=Q_U;
        Qt(mt-3:mt,mt-3:mt)=flipud( fliplr(-Q_U(1:4,1:4) ) );

        e_1t = zeros(mt,1); e_1t(1) = 1;
        e_mt = zeros(mt,1); e_mt(mt) = 1;
        r = dt;
        dt = 1;
      
        QQt = Qt - 1/2*e_1t*e_1t' + 1/2*e_mt*e_mt';
        D1t = HIt*(QQt) ;

        D_adj = 1/dt*HIt*(dt*D1t')*Ht*1/dt;

        o =  null((D_adj)); % no sign
      %  check0 = o' *Ht*D1t  ;     % <<--- verify that o is correct, should be all almost zero vector 
        o_adj = (o'*Ht*1/dt);
        o2 = o'*Ht*o*1/dt ; % should equal to dt
        F = (It- (o*o_adj)/(o2));
        %{
        alpha=-1;
        %D1t = D1t - alpha*HIt*(e_1t*e_1t');
        J = inv(D1t);
        %}
        J = pinv(D1t);
        A = J*F;
        
        for j=1:mt
           A(:,j) = A(:,j)-e_1t'*A(:,j);
        end
        JF = kron(sparse(A),Ix);
        dt = r;
        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        dt = 0.5*block;
        D_adj = 1/dt*HIt*(dt*D1t')*Ht*1/dt;

        o =  null((D_adj)); % no sign
        check0 = o' *Ht*D1t  ;     % <<--- verify that o is correct, should beall zero vector 
        o_adj = (o'*Ht*1/dt);
        o2 = o'*Ht*o*1/dt ; % should equal to dt
        F = (It- (o*o_adj)/(o2));

        %alpha=-1;
        %D1t = D1t - alpha*HIt*(e_1t*e_1t');
        %J = inv(D1t);

       
        J = pinv(D1t);

       A = J*F;
    
       for j=1:mt
           A(:,j) = A(:,j)-e_1t'*A(:,j);
       end
       JF = kron(sparse(A),Ix);
    
    
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
        dt = block; %<<------------------------------------
        D_adj = 1/dt*HIt*(dt*D1t')*Ht*1/dt;

        o =  null((D_adj));
        check0 = o' *Ht*D1t  ;    
        o_adj = (o'*Ht*1/dt);
        o2 = o'*Ht*o*1/dt ; 
        F = (It- (o*o_adj)/(o2));

       J = pinv(D1t);

       A = J*F;
       for j=1:mt
           A(:,j) = A(:,j)-e_1t'*A(:,j);
       end
       JF = kron(sparse(A),Ix);
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

        dt = 0.5*block;
        D_adj = 1/dt*HIt*(dt*D1t')*Ht*1/dt;

        o =  null((D_adj)); 
        check0 = o' *Ht*D1t  ;    
        o_adj = (o'*Ht*1/dt);
        o2 = o'*Ht*o*1/dt ;
        F = (It- (o*o_adj)/(o2));

        J =  pinv(D1t);


        A = J*F;
        for j=1:mt
            A(:,j) = A(:,j)-e_1t'*A(:,j);
        end
       
        JF = kron(sparse(A),Ix);

       

    otherwise
        disp("Wrong name")
end


