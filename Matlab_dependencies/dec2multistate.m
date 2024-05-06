function B = dec2multistate(y,p,n)

% y : the state in decimal representation
% b : number of states
% n : number of nodes

c = y;
B = zeros(1,n);
P = [];
while c>0

    z=round(p*(c/p-fix(c/p)));
    c=fix(c/p);
    
    P = [z,P];

end

B(end - length(P) +1 :end) =  P;

