%% PentaNet: a graph on five vertices conjured up by me to test the 
%  directional residual generator of Speyer and White.
%  Sam Nazari
%  03-Sept-2015
clear,
clc

%% Graph Structure
Ag = [0 1 1 1 0 0 0;
     1 0 1 1 1 0 0;
     1 1 0 0 0 1 0;
     1 1 0 0 0 0 1;
     0 1 0 0 0 0 1;
     0 0 1 0 0 0 1;
     0 0 0 1 1 1 0];
 
Deg = eye(7);
for i = 1:7
    Deg(i,i) = sum(Ag(i,:));
end

C = eye(7)

L = Deg-Ag

x0 = [0 1 1 0 0 0 0]

A = -L

% Condition: Distinct Eigenvalues
eig(A)

%% Construct Fault Vectors
f1 = [1 0 0 0 0 0 0]'  % Vertex one is the intruder
f2 = [0 1 0 0 0 0 0]'  % Vertex two is the intruder
f3 = [0 0 1 0 0 0 0]'  % Vertex three is the intruder
f4 = [0 0 0 1 0 0 0]'  % Vertex four is the intruder
f5 = [0 0 0 0 1 0 0]'  % Vertex five is the intruder
f6 = [0 0 0 0 0 1 0]'  % Vertex six is the intruder
f7 = [0 0 0 0 0 0 1]'  % Vertex seven is the intruder

%E = [f1 f2 f3 f4 f5 f6 f7]
%E = [f2 f3 f4 f5 f6 f7]
E = [f2 f3 f4]
% Choose the agent to be attacked
flt1  = 0
flt2  = 1
flt3  = 1
flt4  = 1
flt5  = 0
flt6  = 0
flt7  = 0

% Choose a magnitude for the attack
f1Val = 10
f2Val = 10
f3Val = 10
f4Val = 10
f5Val = 10
f6Val = 10
f7Val = 10

% Chose the attack time
tf1   = 12
tf2   = 2
tf3   = 2
tf4   = 2
tf5   = 4
tf6   = 5
tf7   = 7

%% Step 1: Obtain M1..M5

% Check to see if the faults are output seperable rank(C*E)=dim(A)
rank(C*E)

c1 = (eye(7)-C*f1*pinv(C*f1)*C)
c2 = (eye(7)-C*f2*pinv(C*f2)*C)
c3 = (eye(7)-C*f3*pinv(C*f3)*C)
c4 = (eye(7)-C*f4*pinv(C*f4)*C)
c5 = (eye(7)-C*f5*pinv(C*f5)*C)
c6 = (eye(7)-C*f6*pinv(C*f6)*C)
c7 = (eye(7)-C*f7*pinv(C*f7)*C)

k1 = A*(eye(7)-f1*pinv(C*f1)*C)
k2 = A*(eye(7)-f2*pinv(C*f2)*C)
k3 = A*(eye(7)-f3*pinv(C*f3)*C)
k4 = A*(eye(7)-f4*pinv(C*f4)*C)
k5 = A*(eye(7)-f5*pinv(C*f5)*C)
k6 = A*(eye(7)-f6*pinv(C*f6)*C)
k7 = A*(eye(7)-f7*pinv(C*f7)*C)

m1 = [c1;c1*k1;c1*k1^2;c1*k1^3;c1*k1^4;c1*k1^5;c1*k1^6]
m2 = [c2;c2*k2;c2*k2^2;c2*k2^3;c2*k2^4;c2*k2^5;c2*k2^6]
m3 = [c3;c3*k3;c3*k3^2;c3*k3^3;c3*k3^4;c3*k3^5;c3*k3^6]
m4 = [c4;c4*k4;c4*k4^2;c4*k4^3;c4*k4^4;c4*k4^5;c4*k4^6]
m5 = [c5;c5*k5;c5*k5^2;c5*k5^3;c5*k5^4;c5*k5^5;c5*k5^6]
m6 = [c6;c6*k6;c6*k6^2;c6*k6^3;c6*k6^4;c6*k6^5;c6*k6^6]
m7 = [c7;c7*k7;c7*k7^2;c7*k7^3;c7*k7^4;c7*k7^5;c7*k7^6]

v1 = 7-rank(m1)
v2 = 7-rank(m2)
v3 = 7-rank(m3)
v4 = 7-rank(m4)
v5 = 7-rank(m5)
v6 = 7-rank(m6)
v7 = 7-rank(m7)

v  = v1+v2+v3+v4+v5+v6+v7

% D will be the same as K1.  D will be computed by hand..
% note that the subpace for each fault actually corresponds to the fault
% vector in this example, ie: v1 = f1 = [1 0 0 0 0]' etc..
% D = 2*eye(5)
D = [0 1 1 1 0 0 0;
    1 -1 1 1 1 0 0;
    1  1 0 0 0 1 0;
    1  1 0 0 0 0 1;
    0  1 0 0 1 0 1;
    0  0 1 0 0 1 1;
    0  0 0 1 1 1 0]


% K1 = D
% F  = A1-K1*C
% K2 = F*H
% 
% K = K1+K2

%% Compute projections

CV = C*E
pCV = CV*inv(CV'*CV)*CV'

%% Run model
sim('intDetEx3p4p1BJDFMDL')

%% plot
%  Note that only faults in vertices 2,3,4 (the vertices in the
%  neighborhood of vertex 1) will be sensitive to failures/intruders. 
%  Intruders at other vertices in the graph will not trigger an alarm.
figure,

subplot(711)
plot(tout,resf1),grid on,ylim([0 30]),title('Residual Signals')
ylabel('f_1')

subplot(712)
plot(tout,resf2),grid on,ylim([0 30])
ylabel('f_2')

subplot(713)
plot(tout,resf3),grid on,ylim([0 30])
ylabel('f_3')

subplot(714)
plot(tout,resf4),grid on,ylim([0 30])
ylabel('f_4')

subplot(715)
plot(tout,resf5),grid on,ylim([0 30])
ylabel('f_5')

subplot(716)
plot(tout,resf6),grid on,ylim([0 30])
ylabel('f_6')

subplot(717)
plot(tout,resf7),grid on,ylim([0 30])
ylabel('f_7')
xlabel('Time (sec)')
