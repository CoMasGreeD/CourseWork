H = [1, 0; 0, 1; -1, 0; 0, -1];
G = [9; 9; -4; -4];
G0 = [3; 3; 1; 1];
A = [0 1; -1 0];
B = [0; 1];
H0 = [1 1];
t0 = 0;
tf = 20;
N = 200;
h = (tf - t0)/N;

lb = -1*ones(1,N);
x_h_values = zeros(2, N+1);
c_h_values = zeros(1, N);
d_h_values = zeros(4, N);
Gi0 = zeros(4, 1);
F=@(t)expm(A*t);

options = optimoptions('fmincon','Algorithm','active-set');

for i=1:4
    FforMax = @(x)(-H(i,:)*F(tf-t0)*x);
    [x1,fval1] = fmincon(FforMax,[-G0(3),-G0(4)]',[],[],[],[],[-G0(3),-G0(4)]',[G0(1),G0(2)]',[],options);
    fval1 = -fval1;
    Gi0(i) = G(i) - fval1;
end

display(Gi0);

funForC = @(t)H0*F(tf-t)*B;
funForD = @(t)H*F(tf-t)*B;

c_h = @(s)integral(funForC,s,s+h,'ArrayValued', true);
d_h = @(s)integral(funForD,s,s+h,'ArrayValued', true);

for i = 1:N
    c_h_values(i) = c_h(t0+h*i-h);
    d_h_values(:,i) = d_h(t0+h*i-h);
end

u = linprog(-c_h_values,d_h_values,Gi0,zeros(N,N),zeros(N,1),lb',-1*lb');

for j=1:6
    x_h_values(:,1)=[0, (j-1)*0.4]';
    for i = 1:N
        t1 = i*h;
        interForKoshi = @(t)(F(t1-t)*B);
        x_h_values(:,i+1) = F(h)*x_h_values(:,i)+u(i)*integral(interForKoshi,t1-h,t1,'ArrayValued', true);
    end
    
    display(x_h_values);
    hold on;
    plot(x_h_values(1,:),x_h_values(2,:));
end

rectangle('Position',[-G0(3),-G0(4),G0(1)+G0(3),G0(2)+G0(4)]);
rectangle('Position',[-G(3),-G(4),G(1)+G(3),G(2)+G(4)]);
