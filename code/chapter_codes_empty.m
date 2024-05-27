% Code for Book Chapter
% By Daniel Plaugher 05/25/2024

%% 1. Add required paths

% home 
clear
clc
close all

addpath('...Matlab_dependencies/');

%% 2. Load or create model 


n = XXXX; % number of nodes/genes
p = 2; % number of states

%%% for SDDS_Build
syms x1 x2   ...   %encode variables for nodes

f= []; % define functions for BN as polynomials, one fnx per line

[varF,nv,F]=SDDS_Build(syms,f,p); % automatically builds vars. for SDDS


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


% --- Option 2a: Alter truth table for nodes 
F1 = TruthTable_del_n_temp(F,nv,varF,p, NODE, ACTION); % knock-in (1) or knock-out (0)

% --- Option 2b: Alter truth table edges
F1 = TruthTable_del_a_temp(F,nv,varF,p,TAIL,HEAD, ACTION); % knock-in (1) or knock-out (0)


% rebuild SDDS with the new conditions
[varF,nv,F]=SDDS_Build(syms,f,p);

%% 5. Simulation
nins = 1000; % number of initializations
nsteps=100; % number of steps for SDDS
g=0.01; % noise (optional - needs SDDS_simNoise)
c = 0.9*ones(2,n); % normal propensities for SDDS
 

%[Y,My]=SDDS_simNoise(g,F,varF,nv,p,c,n, nsteps,nins); % simulation w noise
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



%% -- used for zoom plots
zp = BaseZoom();
zp.plot;

