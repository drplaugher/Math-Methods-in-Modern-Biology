% Code for Book Chapter
% By Daniel Plaugher 05/25/2024

%% 1. Add required paths

clear
clc
close all

addpath('...Matlab_dependencies/');


%% 2. Load or create model 

% ********************************************************
% Load user's model here, or use a pre-loaded model below
% ********************************************************

n = XXXX; % number of nodes/genes
p = 2; % number of states

%%% for SDDS_Build
syms x1 x2      %encode variables for nodes

f= []; % define functions for BN as polynomials, one fnx per line

[varF,nv,F]=SDDS_Build(syms,f,p); % automatically builds vars. for SDDS


%% -- Simple 2-cycle
n = 2; 
p = 2;

%%% for SDDS_Build
syms x1 x2

f= [x2;
    x1]; 
[varF,nv,F]=SDDS_Build(syms,f,p);


%% -- Simple 3-cycle
n = 3; 
p = 2;

%%% for SDDS_Build
syms x1 x2 x3

f= [x3
    x1;
    x2]; 
[varF,nv,F]=SDDS_Build(syms,f,p);

%% -- Repressillator

n = 3; 
p = 2;

%%% for SDDS_Build
syms x1 x2 x3

f= [x3+1
    x1+1;
    x2+1]; 
[varF,nv,F]=SDDS_Build(syms,f,p);

%% -- Reduced PC model

n = 22; 
p = 2;

syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22;

f=[x1; %f1 
    x2 ; %f2
    x3 ; %f3
    x4 ; %f4
    x6*x9*x14*x15*x16+x6*x9*x14*x15+x6*x9*x14*x16+x6*x9*x15*x16+x6*x14*x15*x16+x9*x14*x15*x16+x6*x9*x14+x6*x9*x15+x6*x14*x15+x9*x14*x15+x6*x9*x16+x6*x14*x16+x9*x14*x16+x6*x15*x16+x9*x15*x16+x14*x15*x16+x6*x9+x6*x14+x9*x14+x6*x15+x9*x15+x14*x15+x6*x16+x9*x16+x14*x16+x15*x16+x6+x9+x14+x15+x16; %f5 
    x1*x7+x1+x7 ; %f6
    x3*x5*x8+x3*x5+x3*x8+x5*x8+x3+x5 ; %f7
    x5*x7*x8+x5*x7+x5*x8+x5 ; %f8
    x3*x5+x3+x5 ; %f9
    x8*x9+x9 ; %f10
    x7*x9*x12 ; %f11
    x2*x3*x4*x6*x9+x2*x3*x4*x6+x2*x3*x4*x9+x3*x6*x9+x3*x6+x3*x9 ; %f12
    x13; %f13
    x9*x13*x14+x9*x13+x9*x14+x13*x14+x9+x13+x14 ; %f14
    x13; %f15
    x13*x14*x19+x13*x14+x13*x19+x14*x19+x13+x14; %f16
    x3*x15*x19+x3*x15+x3*x19+x15*x19+x3+x19 ; %f17
    x14*x15*x16*x19+x14*x15*x16+x14*x15*x19+x14*x16*x19+x15*x16*x19+x14*x15+x14*x16+x15*x16+x14*x19+x15*x19+x16*x19+x14+x15+x16 ; %f18
    x16*x17+x16*x19+x19+1 ; %f19
    x16*x18*x19+x16*x18+x16*x19+x18*x19+x18+x19 ; %f20
    x16*x18*x20+x14*x16+x16*x18+x16*x20+x18*x20+x16+x18+x20+1 ; %f21
    x14*x16*x17+x14*x16 %f22
    ];

[varF,nv,F]=SDDS_Build(syms,f,p);


%% 3. Set up for Markov Chains - for smaller models

% -- Finding attractors and fixed points
[~,Avec] = bnAsparse(F,varF,nv); % shows the atractors and their respective basins
[ab,dd] = bnAttractor(Avec); % dd~ steps from attr.
attrs =  unique(ab(ab<0)); % finds attractors

disp(attrs)
a1= find(ab==-1); % finds attractor 
% aXXXX= find(ab==-XXXX); % finds attractor 

% Finding fixed points only
for i= 1:length(attrs)
    attr_size = NaN;
    attr_size = sum(ab(:)==attrs(i));
    if attr_size == 1
        disp(attrs(i))
    end
end

ab_size = length(find(ab==-1)); % gives size of basin for attractor

% find binary representation of attractors
att=a1;
for i=1:length(att)
    x = dec2multistate(att(i)-1,p,n); % binary represention of attractor
end


% -- Create Markov Chain
c = 0.1*ones(2,n); % normal propensities for SDDS
TM=multistateA(F,varF,nv,c,p); % transistion matrix--probability of moving from one node to another 

%plotting 
TM_mc=dtmc(TM);
figure;
graphplot(TM_mc,'ColorEdges',true);

%% -- Approximate a stationary distribution

% --- Google Matrix 
K=(1/2^n)*ones(2^n,2^n); % K matrix for noise
g=0.9;
G=g*TM+(1-g)*K; % google matrix, adds noise to make regular for Perron-Frobenius Thm
G_mc=dtmc(G); % markov chain for google matrix
G_pwr=G^10000; % high powers make columns converge to stationary distribution
G_dist=asymptotics(G_mc); % stationary distribution of google matrix--time spent at node
[B, I]=sort(G_dist,'descend'); % sorts in descending order and oututs associated index
Ranks=[I' B']; % table with index and dist.

%% 4. Inductions and Control

% --- Option 1: Induce mutations and controls directly in F
f(xi)=1; % GOF mutation // Target agonist (knock-in)
f(xj)=0; % LOF mutation // Target inhibition (knock-out)

% f(9) = 0; % ERKs
% f(16) = 0; % PIP3c


% --- Option 2: Alter truth table for nodes 
F1 = TruthTable_del_n_temp(F,nv,varF,p, NODE, ACTION); % knock-in (1) or knock-out (0)

% --- Option 2: Alter truth table edges
F1 = TruthTable_del_a_temp(F,nv,varF,p,TAIL,HEAD, ACTION); 


% rebuild SDDS with the new conditions
[varF,nv,F]=SDDS_Build(syms,f,p);

%% 5. Simulation
nins = 1000; % number of initializations
nsteps=100; % number of steps for SDDS
g=0.01; % noise (optional - needs SDDS_simNoise)
c = 0.9*ones(2,n); % normal propensities for SDDS
 

%[Y,My]=SDDS_simNoise(g,F,varF,nv,p,c,n, nsteps,nins); 
[Y,My]=SDDS_sim(F,varF,nv,p,c,n, nsteps,nins); % simulation w/o noise
Ylast=Y(:,end); % long-term trajectories


%% 6. Graphing (customize)
X = 0:1:nsteps; % time steps

% adjust plots as needed
figure('Name', 'Simulation')
plot(X,Y)
legend()%
xlabel('Time Steps')
ylabel('Average Frequencies')
title('Descriptive title HERE')




% ----- For Reduced PCC model
% figure('Name', '')
% plot(X,Y(10,:),X,Y(11,:),X,Y(12,:),X,Y(8,:), X,Y(21,:),'k',X,Y(20,:),'b', X,Y(22,:),'r','LineWidth',1.5,'MarkerSize',10)
% legend('Prols','Migs', 'Acts','Apops','Autc','Apoc','Proc')
% xlabel('Time Steps')
% ylabel('Average Frequencies')
% title('Pancreatic Cancer Model')


%% -- used for zoom plots
zp = BaseZoom();
zp.plot;

