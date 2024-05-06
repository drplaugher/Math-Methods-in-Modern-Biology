function decx = multistate2dec(x,p,n)

% x : the state in multistate representation
% p : number of states
% n : number of nodes

% for conversion to decimal
pow = p.^(n-1:-1:0)';  
% convert state to decimal
decx = x*pow + 1;                         