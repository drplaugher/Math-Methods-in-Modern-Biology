function [traj, y] = SDDSRun(x,F,varF,nv,p,c,nsteps)

% y = delayRun(x,F,varF,nv,nsteps) - run a Boolean Network
% This function runs a BN starting from the initial starting state x (a
% binary row vector) for nsteps steps. It essentially uses delayNextState in
% an iterative fashion. The output is a matrix Y (nsteps+1-by-length(x))
% containing the history of the process. If x is not a vector, but the word
% 'rand', then a random starting state is used
% M : Maximum expressions

% Functions used: bnNextState.m

% Ilya Shmulevich; Aug. 14, 2001
% Modified May 22, 2003 by HL.
% Modified Feb 24, 2011 by Boris Aguilar

n = size(varF,2); % number of genes
if ischar(x) % if we want a random start
    x = rand(1,n)>0.5;
end

traj = zeros(1,nsteps+1);
pow = p.^(n-1:-1:0)';

decx = x*pow + 1;                           % convert state to decimal
traj(1) = decx; 

y = zeros(nsteps+1,n); % initialize Y
y(1,:) = x; % first one is initial state

for step = 2 : nsteps+1
   x = SDDSNextState(x,F,varF,nv,p,c);
   y(step,:) = x; % update history
   decx = x*pow + 1;                           % convert state to decimal
   traj(step) = decx;
end
