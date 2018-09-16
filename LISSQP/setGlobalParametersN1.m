function setGlobalParameters(sigSq_t, DMrat, Vt, altIn, altErr)
%Function to set global parameters for the control problem, not the model

%global P nt Int K sig0sq Kact M T ;
%global altitude V B dist2 Xact D alpha k dt Vact;


%%%%%%%%%%%%%%   PHYSICAL GLOBAL SCALARS   %%%%%%%%%%%%%%
global P M D T B nt dt V alt sig0sq;

%Max instantaneous transmit power
P=100;

%Bandwidth in Hz
B=1e5;

%Max memory of node (bits)
M=1e9;

%initial data on node
D=DMrat*M;

%Simulation time, minutes
T=40;

%time length of each iteration (seconds)
dt = (T*60)/nt;

%Expected Node velocity m/s
V=Vt;

%Recieve Noise
sig0sq = sigSq_t;

%minimum distance from base station
alt=1000;
if(isempty(altIn)==0)
    alt = altIn*alt;
end
if(isempty(altErr))
    altAct=alt;
else
    altAct=altErr*alt;
end
d=0;%700;

%%%%%%%%%%%%%%   PHYSICAL GLOBAL SCALARS   %%%%%%%%%%%%%%
global Vact Xact Xest Int norm_e norm_a k K Ka alpha;

%path loss exponent > 1 corresponds to slightly obscured environment
alpha=1.5*ones(nt+1,1);

%Speed Uncertainty (%) - zero mean
su = 0;
SU = (su/2) - rand(nt+1,1)*su;

%Actual Node Velocity (time varying)
Vact = V*(1+SU);

%antenna gain - constant isotropic
C = ones(nt+1,1);

%Channel Coefficient
k=C./sig0sq;

%Discrete integration operator
Int=ones(nt+1, 1);
Int(1)=0.5;
Int(nt+1)=0.5;
Int=Int*dt; %asociate with time

%AP is at coordinates (0,0,0) -> define node trajectory w.r.t AP
x_init = -(V*T*60)/2; %Fly past AP at time t=T/2
Xact = ones(nt+1,1).*x_init + tril(ones(nt+1), -1)*(dt*Vact);
Xest = ones(nt+1,1).*x_init + tril(ones(nt+1), -1)*ones(nt+1,1)*(dt*V);

%Distances to AP
norm_e = ones(nt+1,1)*(alt^2+d^2);
norm_a = ones(nt+1,1)*(altAct^2+d^2);
for i=1:nt+1
    t=dt*(i-1);
    %ideal
    norm_e(i) = norm_e(i)+(V*t-0.5*V*T*60)^2; %symetry of trajectory
    %actual (used for fading horizon)
    norm_a(i) = (Xact(i))^2 + norm_a(i);
end

%Term to be used in calculating the SNR
% ricMag = 1 + sqrt(randn(nt+1,1).^2 + randn(nt+1,1).^2);
% K = ricMag./K;
K = k./(norm_e.^alpha);
Ka = k./(norm_a.^alpha);

end




%Noise power on reciever(W)
%sig0sq = 0.5e-10*rand(1,1)*ones(nt+1,1);


