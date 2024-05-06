function [Y,My]=SDDS_sim(F,varF,nv,p,c,n, nsteps,nins)
% Simulation of an SDDS system, output allows you to graph the frequency of
% node expressions given a number of desired initializations
% Written by Daniel Plaugher 2019


% F : transition table 
% varF : variables for each node
% nv : number of variables per node
% p: number of possible states i.e. Boolean=2
% c : probabilites of activation (first row) and degradation (second row)
% nins: number of initializations
% nsteps: number of steps for SDDSRun


% Uses: SDDSRun 

My = cell(nins,1);


for i=1:nins
    x=randi([0 1], 1,n); % random binary vector
    [~, y] = SDDSRun(x,F,varF,nv,p,c,nsteps); % Runs SDDS for desired # of steps
    My{i}=y; % stores each output of y in a separate "cell"
end


freq =0;
for m=1:nins
    freq = freq + My{m,1}; % gives matrix for frequencies of node expressions
end 
Y=(1/nins).*freq'; % scales and transposes frequency matrix for plotting



end