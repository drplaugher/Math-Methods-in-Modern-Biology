 function FF = TruthTable_del_n_temp(F,nv,varF,p, node,v)
 % node to delete
 
 % David Murrugarra; May 15, 2014

n = length(nv); 
FF = F;

k = 0;
for i = 1 : n
    for j = 1 : nv(i) 
        if varF(j,i) == node
            k = k+1;
            TT = TruthTable_del_a_temp(F,nv,varF,p,node,i,v);
            FF(:,i) = TT(:,i);
            break
        end
    end
end

FF(1:p^(nv(node)),node) = v*ones(p^(nv(node)),1);

 end 
 
 