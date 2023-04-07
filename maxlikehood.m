function Opt = maxlikehood(RxSymbol,Hreal,cons)

val= repmat({cons'},size(Hreal,2),1);
[x{1 : numel(val)}] = ndgrid(val{:});
ret = reshape(cat(numel(val), x{:}), [], numel(val))';
[~,choice]=min(sum((RxSymbol-Hreal*ret).^2,1));
Opt=ret(:,choice);

end

