clc
clear all

% Labels of Parameters
syms l1 l2 l3 real;
% l1=sym('l1','real');
% l2=sym('l2','real');
% l3=sym('l3','real');

syms m1 m2 m3 real;
% m1=sym('m1','real');         %Masse  of link1
% m2=sym('m2','real');         %Masse  of link2
% m3=sym('m3','real');         %Masse  of link3
m_s = [m1;m2;m3];

syms q1 q2 q3 qp1 qp2 qp3 real;
% q1=sym('q1','real');         %Coordinate angle 1
% q2=sym('q2','real');         %Coordinate angle 2
% q3=sym('q3','real');         %Coordinate angle 3
% qp1=sym('qp1','real');       %Coordinate angle 1, velocity
% qp2=sym('qp2','real');       %Coordinate angle 2, velocity
% qp3=sym('qp3','real');       %Coordinate angle 3, velocity

q=[q1;q2;q3];
qp=[qp1;qp2;qp3];     %Vector of velocity of angle

syms I111 I112 I113 I122 I123 I133 I211 I212 I213 I222 I223 I233 I311 I312 I313 I322 I323 I333 real;
% I111=sym('I111','real');      %Componentes of the inertia tensor
% I112=sym('I112','real');
% I113=sym('I113','real');
% I122=sym('I122','real');
% I123=sym('I123','real');
% I133=sym('I133','real');
% 
% I211=sym('I211','real');
% I212=sym('I212','real');
% I213=sym('I213','real');
% I222=sym('I222','real');
% I223=sym('I223','real');
% I233=sym('I233','real');
% 
% I311=sym('I311','real');
% I312=sym('I312','real');
% I313=sym('I313','real');
% I322=sym('I322','real');
% I323=sym('I323','real');
% I333=sym('I333','real');

I1=[I111 I112 I113;
    I112 I122 I123;
    I113 I123 I133];

I2=[I211 I212 I213;
    I211 I222 I223;
    I213 I223 I233];

I3=[I311 I312 I313;
    I312 I322 I323;
    I313 I323 I333];
I_s = [{I1} {I2} {I3}];
g=sym('g','real');           %gravity

syms pi,'real';     %% this is important, if it not do like this, cos(pi/2) will not be =0;


%compute Transformation Matrix
HT1_0 = dh2tran(q1, l1, 0, pi/2);

HT2_1 = dh2tran(q2, 0, l2, 0);
HT2_0 = HT1_0*HT2_1;

HT3_2 = dh2tran(q3, 0 , l3, 0);
HT3_0 = HT2_0*HT3_2;

HTcm1_0 = dh2tran(q1, l1, 0, pi/2);

HTcm2_1 = dh2tran(q2, 0, l2, 0);
HTcm2_0 = HT1_0*HTcm2_1;

HTcm3_1 = dh2tran(q3, 0 , l3, 0);
HTcm3_0 = HT2_0*HTcm3_1;


HTi0_s(1) = {HT1_0};
HTi0_s(2) = {HT2_0};
HTi0_s(3) = {HT3_0};
HTcmi0_s(1) = {HTcm1_0};
HTcmi0_s(2) = {HTcm2_0};
HTcmi0_s(3) = {HTcm3_0};

numOfJoint = 3;
isRevolute = [1 1 1]; 

syms g_vec  P  'real';
g_vec = [0;0;-g];

%Compute Jacobians (Jcm_s)
for idxOfCM = 1:numOfJoint
    curJcm = sym(zeros(6,numOfJoint));
    curHTcm = HTcmi0_s{idxOfCM};
    
    z = [0;0;1];
    tcm  = curHTcm(1:3,4);
    
    if(isRevolute(1))
        curJcm(1:3,1) = cross(z, tcm);
        curJcm(4:6,1) = z;
    else
        curJcm(1:3,1) = z;
        curJcm(4:6,1) = [0;0;0];
    end
    
    for idxOfJoint = 2:idxOfCM
        curHT = HTi0_s{idxOfJoint - 1};
        
        z = curHT(1:3,3);
        tlink = curHT(1:3,4);
        
        if(isRevolute(idxOfJoint))
            t = cross(z, (tcm - tlink));
            curJcm(1:3,idxOfJoint) = cross(z, (tcm - tlink));
            curJcm(4:6,idxOfJoint) = z;
        else
            curJcm(1:3,idxOfJoint) = z;
            curJcm(4:6,idxOfJoint) = [0;0;0];
        end
    end
    Jcm_s(idxOfCM) = {simple(curJcm)};
end

%Compute M
M = sym(zeros(3,3));
for idxOfJoint = 1:numOfJoint
    curJcm = Jcm_s{idxOfJoint};
    Jv = curJcm(1:3,1:3);
    Jw = curJcm(4:6,1:3);
    
    curHT = HTi0_s{idxOfJoint};
    R = curHT(1:3,1:3);
    
    M = simple(M + m_s(idxOfJoint)*Jv'*Jv + Jw'*R*I_s{idxOfJoint}*R'*Jw);
end
M = simple(M);

%Compute C
C = sym(zeros(3,3));
for k = 1:3
    for j = 1:3
        for i = 1:3
            C(k,j) = simple(C(k,j) + (1/2)*(diff(M(k,j), q(i)) + diff(M(k,i), q(j)) - diff(M(i,j), q(k)))*qp(i));
        end
    end
end
C = simple(C);

%compute G
%add up potential energy
for idxOfCM = 1:numOfJoint
    P = simple(P + m_s(idxOfCM)*g_vec'*(HTcmi0_s{idxOfCM}(1:3,4)));
end
for idxOfCM = 1:numOfJoint
    G(idxOfCM,1) = simple(diff(P, q(idxOfCM)));
end
G = simple(G);

M
C
G
syms qp1r qp2r qp3r real
syms qpp1r qpp2r qpp3r real
qpr = [qp1r;qp2r;qp3r];
qppr = [qpp1r;qpp2r;qpp3r];
tau = M*qppr + C*qpr + G

Theta(1,1)=I122;
Theta(2,1)=I322;
Theta(3,1)=I222;
Theta(4,1)=I211;
Theta(5,1)=I311;
Theta(6,1)=I212;
Theta(7,1)=l2^2*m2;
Theta(8,1)=l2^2*m3;
Theta(9,1)=l3^2*m3;
Theta(10,1)=I312;
Theta(11,1)=l2*l3*m3;
Theta(12,1)=I323;
Theta(13,1)=I213;
Theta(14,1)=I313;
Theta(15,1)=I223;
Theta(16,1)=I233;
Theta(17,1)=I333;
Theta(18,1)=l3*m3;
Theta(19,1)=l2*m3;
Theta(20,1)=l2*m2;

numOfPara = 0;
expandedTau = expand(tau);
for jointIdx = 1:length(expandedTau)
    asStr = char(expandedTau(jointIdx));
    termsAsStr = strsplit(str,delimiter);
end

Yr = sym(zeros(3,20));
Yr(1,1)= mycoeff(tau(1),I122);
Yr(1,2)= mycoeff(tau(1),I322);
Yr(1,3)= mycoeff(tau(1),I222);
Yr(1,4)= mycoeff(tau(1),I211);
Yr(1,5)= mycoeff(tau(1),I311);
Yr(1,6)= mycoeff(tau(1),I212);
Yr(1,7)= mycoeff2(mycoeff(tau(1),m2),l2);
Yr(1,8)= mycoeff2(mycoeff(tau(1),m3),l2);
Yr(1,9)= mycoeff2(mycoeff(tau(1),m3),l3);
Yr(1,10)=mycoeff(tau(1),I312);
Yr(1,11)=mycoeff(mycoeff(mycoeff(tau(1),m3),l2),l3);
Yr(1,12)=mycoeff(tau(1),I323);
Yr(1,13)=mycoeff(tau(1),I213);
Yr(1,14)=mycoeff(tau(1),I313);
Yr(1,15)=mycoeff(tau(1),I223);
Yr(1,16)=mycoeff(tau(1),I233);
Yr(1,17)=mycoeff(tau(1),I333);
Yr(1,18)=mycoeff0(mycoeff(mycoeff(tau(1),m3),l3),l2);
Yr(1,19)=mycoeff0(mycoeff(mycoeff(tau(1),m3),l2),l3);
Yr(1,20)=mycoeff0(mycoeff(mycoeff(tau(1),m2),l2),l3);
Yr(2,1)= mycoeff(tau(2),I122);
Yr(2,2)= mycoeff(tau(2),I322);
Yr(2,3)= mycoeff(tau(2),I222);
Yr(2,4)= mycoeff(tau(2),I211);
Yr(2,5)= mycoeff(tau(2),I311);
Yr(2,6)= mycoeff(tau(2),I212);
Yr(2,7)= mycoeff2(mycoeff(tau(2),m2),l2);
Yr(2,8)= mycoeff2(mycoeff(tau(2),m3),l2);
Yr(2,9)= mycoeff2(mycoeff(tau(2),m3),l3);
Yr(2,10)=mycoeff(tau(2),I312);
Yr(2,11)=mycoeff(mycoeff(mycoeff(tau(2),m3),l2),l3);
Yr(2,12)=mycoeff(tau(2),I323);
Yr(2,13)=mycoeff(tau(2),I213);
Yr(2,14)=mycoeff(tau(2),I313);
Yr(2,15)=mycoeff(tau(2),I223);
Yr(2,16)=mycoeff(tau(2),I233);
Yr(2,17)=mycoeff(tau(2),I333);
Yr(2,18)=mycoeff0(mycoeff(mycoeff(tau(2),m3),l3),l2);
Yr(2,19)=mycoeff0(mycoeff(mycoeff(tau(2),m3),l2),l3);
Yr(2,20)=mycoeff0(mycoeff(mycoeff(tau(2),m2),l2),l3);
Yr(3,1)= mycoeff(tau(3),I122);
Yr(3,2)= mycoeff(tau(3),I322);
Yr(3,3)= mycoeff(tau(3),I222);
Yr(3,4)= mycoeff(tau(3),I211);
Yr(3,5)= mycoeff(tau(3),I311);
Yr(3,6)= mycoeff(tau(3),I212);
Yr(3,7)= mycoeff2(mycoeff(tau(3),m2),l2);
Yr(3,8)= mycoeff2(mycoeff(tau(3),m3),l2);
Yr(3,9)= mycoeff2(mycoeff(tau(3),m3),l3);
Yr(3,10)=mycoeff(tau(3),I312);
Yr(3,11)=mycoeff(mycoeff(mycoeff(tau(3),m3),l2),l3);
Yr(3,12)=mycoeff(tau(3),I323);
Yr(3,13)=mycoeff(tau(3),I213);
Yr(3,14)=mycoeff(tau(3),I313);
Yr(3,15)=mycoeff(tau(3),I223);
Yr(3,16)=mycoeff(tau(3),I233);
Yr(3,17)=mycoeff(tau(3),I333);
Yr(3,18)=mycoeff0(mycoeff(mycoeff(tau(3),m3),l3),l2);
Yr(3,19)=mycoeff0(mycoeff(mycoeff(tau(3),m3),l2),l3);
Yr(3,20)=mycoeff0(mycoeff(mycoeff(tau(3),m2),l2),l3);

for i = 1:3
    for j = 1:20
        s(i,j) = {sprintf('Yr(%d,%d) = %s;\n', i ,j , char(Yr(i,j)))};
    end
end
ts = '';
for i = 1:3
    for j = 1:20
        ts = sprintf('%s%s',ts,s{i,j});
    end
end
    