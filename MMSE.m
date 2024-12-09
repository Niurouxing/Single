function symest = MMSE(TxAntNum,RxSymbol,H,Nv,sym)
S_mmse = (H'*H+Nv*eye(2*TxAntNum)) \ H'*RxSymbol;
[~,Indice]=min(abs(S_mmse-sym),[],2);
symest = sym(Indice);


end

