A = [0 1 0 0; -2.1 0 -0.69 0;0 0 0 1;-0.92 0 -4.308 0];
B = [0 0; -1 1;0 0;1.467 2.533];
f = [0; 1;0;12/9];
x0 = [0.1; -0.1; 0; 0.33];
    t0 = 0;
    tf = 15;
    g_d = [-1;-1;-1;-1];
    g_u = [1;1;1;1];
    H = eye(4);
    N=750;
    h = (tf-t0)/N;
    Th = 0:h:tf-h;
    w =@(t) 0.01*sin(2*t);
    F = @(t) expm(A*t); %F(tf,t)
  
    l = @(t) H*F(tf-t)*f;
    
    dn = zeros(N,4);
    dv = ones(N,4);
     d = @(s) H*integral(@(t)F(tf-t)*B, s,s+h,'ArrayValued', true);
          Clp=h*ones(4,N);
     Alp = [];
     for s = Th
         Alp = [Alp, d(s), -d(s)];
     end
      
    gama = @(s) w(s)*integral(@(t)abs(l(t)),t0,tf,'ArrayValued', true);
    blp_n = g_d + gama(t0) -H*F(tf-t0)*x0;
    blp_v = g_u - gama(t0) -H*F(tf-t0)*x0;
    u0= linprog(Clp,[Alp;-Alp],[blp_v;-blp_n],[],[],dn,dv);
    u=reshape(u0, 4, N);
    u0=u(1,:)+u(2,:)+u(3,:)+u(4,:);
    u0=[u0'; u0(end)];
% ---
   display(u0) 
    
  figure(1);
  subplot(3,1,1);  grid on;hold on;  
  stairs(0:h:tf, u0);
  ylim([-1.1, 1.1])

