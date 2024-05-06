function FF = TruthTable_del_a_temp(F,nv,varF,p,tail,head,v)
% arrow to delete (tail,head)
 
% David Murrugarra; Jul 23, 2018


if ismember(tail, varF(:,head))
   [~,loc] = ismember(tail, varF(:,head));
end
FF = F;

for i = 1 : p^(nv(head))
    x = dec2multistate(i-1,p,nv(head));
    if x(loc) ~= v
        x(loc) = v;
        decx = multistate2dec(x,p,nv(head));
        FF(i,head) = F(decx,head);
    end
end

