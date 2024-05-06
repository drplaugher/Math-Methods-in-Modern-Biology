function [varF,nv,F] = SDDS_Build(syms,f,p)
% This code will build the variables needed for SDDS: varF, nv, F

% Inputs required: function variables as symbols (syms), 
%                   functions in polynomial form (f)
%                   number of possible states (p)
% Outside code required: dec2multistate(y,p,n)

% outputs: varF, nv, F

%written by Daniel Plaugher 3/25/22
p;
syms;

nf=size(f,1); %number of functions

nv = NaN(1,nf);
% stores nv 
for k=1:nf
    vars = sort(symvar(f(k))); % lists sorted variables of function
        %varF(:,k) = vars; % stores variables in varF
    nv(k) = size(vars,2); % number of variables in equation
end
MAXnv=max(nv); % maximum number of variables in one function

% stores F for SDDS
F = -1*ones(p^MAXnv,nf); 
for k=1:nf
    vars = sort(symvar(f(k))); % lists sorted variables of function
    numvar = size(vars,2); % number of variables in equation
    sizeSS = p^numvar; % size of state space for numvars
    decstates = 1:sizeSS; % lists decimal states

    B = NaN(sizeSS,numvar);
    %builds lexicographic truth table
    for i = 1:sizeSS
        B(i,:) = dec2multistate(i-1,p,numvar); 
    end
    
    % solves functions for tt values, stores as cols in F
    for j = 1:sizeSS
        F_eval = subs(f(k), vars, B(j,:)); 
        F(j,k) = mod(F_eval,p);
    end
end


% stores varF for SDDS
varF = -1*ones(MAXnv,nf);
for k=1:nf
    vars = sort(symvar(f(k))); % lists sorted variables of function
    numvar = size(vars,2); % number of variables in equation
    str = string(vars);% converts variables to string 
    newstr = erase(str,"x"); % makes new string without x's
    varname = str2double(newstr);
    for j=1:numvar
        varF(j,k) = varname(j); % stores variables in varF
    end
end





end