H = [1, 0, 1, 0; 0, 1, 0, 1; -1, 0, -1, 0; 0, -1, 0, -1];
G = [9; 9; -4; -4];
G0 = [3; 3; 1; 1];
A = [0 1 0 0; -2.1 0 -0.69 0;0 0 0 1;-0.92 0 -4.308 0];
B = [0 0; -1 1;0 0;1.467 2.533];
f = [0; 1;0;12/9];
H0 = [1 1 1 1];
t0 = 0;
tf = 15;
N = 750;
h = (tf - t0)/N;
g_d = [-1;-1;-1;-1];
g_u = [1;1;1;1];
lb = -1*ones(N,N);

l = @(t) H*F(tf-t)*f;
x_h_values = zeros(4, N+1);
c_h_values = zeros(2, N);
d_h_values = zeros(4, N);
F=@(t)expm(A*t);
x0 = [0.1; -0.1; 0; 0.33];
w =@(t)0.01*sin(2*t);
options = optimoptions('fmincon','Algorithm','active-set');

funForC = @(t)H0*F(tf-t)*B;
funForD = @(t)H*F(tf-t)*B;

c_h = @(s)integral(funForC,s,s+h,'ArrayValued', true);
d_h = @(s)integral(funForD,s,s+h,'ArrayValued', true);
t=d_h(t0);
d_h_values(:,1) =t(:,1);
for i = 1:N
    c_h_values(:,i) = c_h(t0+h*i-h);
    d_h_values(:,2*i) =t(:,2);
    t=d_h(t0+h*i);
    if(i~=N)
    d_h_values(:,2*i+1) =t(:,1);
    end
end
 gama = @(s) w(s)*integral(@(t)abs(l(t)),t0,tf,'ArrayValued', true);
 sz=size(gama(t0));
  blp_n = g_d + gama(t0) - H*F(tf-t0)*x0;
  blp_v = g_u - gama(t0) -H*F(tf-t0)*x0;
 %lb',-1*lb'
u0 = linprog(-c_h_values,[d_h_values;-d_h_values],[blp_v;-blp_n],zeros(2*N,2*N),zeros(N,2),lb',-1*lb');
u=reshape(u0, 2, N);
    u0=u(1,:)+u(2,:);
    u0=[u0'; u0(end)];
   figure(1);
  subplot(3,1,1);  grid on;hold on;  
  stairs(0:h:tf, u0);
  ylim([-1.1, 1.1])
  x_h_values(:,1)=x0;  
    for i = 1:N
        t1 = i*h;
        interForKoshi = @(t)(F(t1-t)*B);
        v=u0(i)+w(t1);
        buff=v*integral(interForKoshi,t1-h,t1,'ArrayValued', true);
        x_h_values(:,i+1) = F(h)*x_h_values(:,i)+buff(:,1)+buff(:,2);
    end
%display(x_h_values);
subplot(3,1,2); 
hold on;
plot(x_h_values(1,:),x_h_values(2,:));
subplot(3,1,3); 
hold on;
plot(x_h_values(3,:),x_h_values(4,:));

