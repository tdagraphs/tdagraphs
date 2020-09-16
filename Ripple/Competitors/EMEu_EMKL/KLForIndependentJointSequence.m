function kl = KLForIndependentJointSequence(p,q,NGroup)
% p,q are row vectors
kl = zeros(NGroup,1);;
if (~isequal(size(p),size(q)))
  return;
end;
NVariable = floor(length(p)/NGroup);
if (NVariable>25)
  return;
end;  
ZERO_TH = 1E-7;
p(find(p>1-ZERO_TH)) = (1-ZERO_TH);
p(find(p<ZERO_TH)) = (ZERO_TH);
q(find(q>1-ZERO_TH)) = (1-ZERO_TH);
q(find(q<ZERO_TH)) = (ZERO_TH);
NConfiguration = 2^(NVariable);
PP = [p;1-p];
QQ = [q;1-q];
%{
for i_config = 0:(NConfiguration-1)
  selector = de2bi(i_config,NVariable)+1;
  inx = sub2ind(size(P),selector,1:NVariable);
  kl = kl+prod(P(inx))*sum(log(P(inx))-log(Q(inx)));
end;
%}

i_config = (0:(NConfiguration-1))';
selector = de2bi(i_config,NVariable)+1;
start_i = 0;
for group_i=1:NGroup
  P = PP(:,[group_i:NGroup:end]);
  %Q = QQ(:,[start_i+1:start_i+NVariable]);
  Q = QQ(:,[group_i:NGroup:end]);
  inx = sub2ind(size(P),selector,repmat(1:NVariable,NConfiguration,1));
  try
    kl(group_i) = sum(prod(P(inx),2).*sum(log(P(inx))-log(Q(inx)),2),1);
  catch
    kl(group_i) = 0;
  end;
end;
clear i_config;
end
