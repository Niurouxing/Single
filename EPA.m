function  symest= EPA(TxAntNum,RxAntNum,slen,RxSymbol,H,Nv,sym,delta,iterNum)

delta=0.1;

Sigma=inv(H'*H/Nv+eye(2*TxAntNum));
Mu=Sigma*H'*RxSymbol/Nv;

diagSigma=diag(Sigma);

h2=diagSigma./(1-diagSigma);
t=Mu./(1-diagSigma);

V=1./h2;
rho=t./h2;

temp=V*Nv;

for i=1:iterNum
        
    disExp=exp(-0.5*V*sym'.^2+rho*sym');

    eta=sum(sym'.*disExp,2)./sum(disExp,2);

    m=RxSymbol-H*eta;

    rho_new=H'*m/Nv+V.*eta;

    rho=delta*rho_new+(1-delta)*rho;


end
[~,Indice]=min(abs(eta-sym'),[],2);
symest = sym(Indice);
 
end