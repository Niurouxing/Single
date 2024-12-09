function [Opt, CC_mnsd, average_K] = MNSD(TxAntNum, RxSymbol, Hreal, Nv, sym,TxSymbol_real)


Kmax=10;
num=10;

symest_mmse =(Hreal'*Hreal+Nv*eye(2*TxAntNum)) \ Hreal'*RxSymbol;

PED_thre = zeros(2*TxAntNum + 1,1);

[bQ, bR] = qr(Hreal);
R = bR(1:2*TxAntNum, :);
Q = bQ(:, 1:2*TxAntNum);
z = Q' * RxSymbol;

for i=2*TxAntNum:-1:1
    PED_thre(i) = PED_thre(i+1) + (z(i) - R(i,i:end)*symest_mmse(i:end))^2;
end
K_D = ones(1, 2*TxAntNum);
conSize = length(sym);
currentSurvive=[];
for step=1:length(z)
     expandedPath=[repmat(1:conSize,1,max(size(currentSurvive,2),1));kron(currentSurvive,ones(1,conSize))];
     PED=sum((z(end-step+1:end)- R(end-step+1:end,end-step+1:end)*sym(expandedPath)).^2,1);

     [sorted_PED,choice]=sort(PED, 'ascend'); 
     index = find(sorted_PED>PED_thre(2*TxAntNum+1-step), 1, 'first');
     if(isempty(index))
         index=1;
     end

     K_D(2*TxAntNum-step+1) = min(min(Kmax, index+num), length(PED));
     currentSurvive=expandedPath(:,choice(1:K_D(2*TxAntNum-step+1)));    
end
Opt= sym(currentSurvive(:,1))';

K_D(1) = 1;
average_K = mean(K_D);
K_D = [K_D 1];
RA_D = (  sum(  K_D((1:2*TxAntNum)+1).*( (2*TxAntNum-(1:2*TxAntNum)) +  3*min(K_D(1:2*TxAntNum), repmat(length(sym),[1,2*TxAntNum])) ) )  )   *32;
RM_D = (  sum(  K_D((1:2*TxAntNum)+1).*( (2*TxAntNum-(1:2*TxAntNum)) +  2*min(K_D(1:2*TxAntNum), repmat(length(sym),[1,2*TxAntNum]))  ) )   ) * 32^2;
CP_D =  sum( K_D((1:2*TxAntNum)+1).*  min(repmat(length(sym),[1,2*TxAntNum]), K_D(1:2*TxAntNum) ) .* log2(K_D(1:2*TxAntNum) +1)) * 32;
% CC_mnsd = RA_D + RM_D + CP_D;
CC_mnsd = RM_D;



end