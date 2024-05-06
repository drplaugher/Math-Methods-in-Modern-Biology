function A = multistateA(F,varF,nv,c,p)

% [A,v] = pbnA(F,varF,nf,nv,cij,p) - state transition matrix of a PBN
% This function creates the 2^n x 2^n transition matrix A corresponding to
% a PBN. Another optional output argument of the function is the stationary
% distribution v.
% INPUT:
% p     - number of states
% Other inputs are defined e.g. in pbnRnd.m

% Functions used: pbnAij.m

% Ilya Shmulevich; Aug. 14, 2001
% Modified May 13, 2003 by HL.
% Modified Nov 9, 2010 by Boris Aguilar
% modified Sept 11, 2013 David Murrugarra

n = length(nv); % number of genes
A = zeros(p^n,p^n);

% Set each element of A.
for i = 1:p^n
    A(i,:) = multistateAi(F,varF,nv,c,i,p);
%     for j = 1:p^n,
%         A(i,j) = multistateAij(F,varF,nf,nv,c,i,j,p);
%     end
end
% if nargout == 2
%     [v,d] = eig(A');
%     [~,idxeig] = max(diag(d));             % largest eigenvalue
%     v = v(:,idxeig);
%     v = v/sum(v);
% end
