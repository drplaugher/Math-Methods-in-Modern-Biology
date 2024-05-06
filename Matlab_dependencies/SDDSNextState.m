function y = SDDSNextState(x,F,varF,nv,p,c)


% This function performs one step of SDDS. The current state is x
% (binary vector of 1 x n, where n is # of genes) For other parameters, see
% e.g.

% x : current state
% F : transition table 
% varF : variables for each node
% nv : number of variables per node
% c : probabilites of activation (first row) and degradation (second row)

% Boris Aguilar ; Feb 24 , 2011
% Modified by David Murrugarra; Mar 24, 2012

n = length(x);          % number of genes
y = zeros(size(x));     % initialize next state

b = p.^(size(varF,1)-1:-1:0)';

z = zeros(1,n);

for i = 1 : n 
    z(i) = F(x(varF(1:nv(i),i))*b(end-nv(i)+1:end)+1,i);     
end

% Update each node separately.
for i = 1 : n
    if ( x(i) > z(i) ) %degradation
        r = rand;
        if r < c(2,i)
            y(i) = z(i);
        else
            y(i) = x(i);
        end
        
    elseif ( x(i) < z(i) ) % activation
     r = rand;
       if r < c(1,i)
            y(i) = z(i);
        else
            y(i) = x(i);
        end
    else
        y(i) = x(i);
    end
        
end 