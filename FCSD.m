function [symest, CC] = FCSD(TxAntNum, HMat, RxSymbol, symQAM, fulexplayers)
%Fixed-complexity sphere decoding
%Reference: 
% "Fixing the Complexity of the Sphere Decoder for MIMO Detection", TWC, 2008.

%Parameters:
% HMat: Complex channel matrix;
% RxSymbol: Complex received signal vector;
% symQAM: Complex modulation constellation;
% fulexplayers = Fully expansion layers

Ns = [zeros(1, TxAntNum- fulexplayers), ones(1, fulexplayers)];



[Q,R] = qr(HMat);
z = Q'*RxSymbol;
Mc = length(symQAM);
paths = 0;
PED = 0;
num_path = 1;
RVMs = 0;
RVAs = 0;

for layer = TxAntNum:-1:(TxAntNum-fulexplayers+1)
    paths_child_slides = cell(num_path,1);
    PED_slides = cell(num_path,1);
    fatherpath_slides = cell(num_path,1);

    for cnt = 1:num_path
        if layer<TxAntNum
            temp = R(layer,(layer+1):end)*paths(cnt, (layer+1):end).';
        else
            temp = 0;
        end
        PED_child_curpath = zeros(Mc,1);
        paths_child_curpath = zeros(Mc,1);
        for m = 1:Mc
            PED_child_curpath(m) = PED(cnt) + abs(z(layer)- R(layer, layer) * symQAM(m) - temp)^2;
            paths_child_curpath(m) = symQAM(m);
        end
        PED_slides{cnt} = PED_child_curpath;
        paths_child_slides{cnt} = paths_child_curpath;
        fatherpath_slides{cnt} = cnt*ones(Mc,1);

        RVMs = RVMs + num_path*(4*(TxAntNum-layer)+6*length(symQAM));
        RVAs = RVAs + num_path*(4*(TxAntNum-layer)+8*length(symQAM));
    end
    paths_child = cell2mat(paths_child_slides);
    PED_child = cell2mat(PED_slides);
    father = cell2mat(fatherpath_slides);

    num_path = num_path * Mc;

    temp_paths = zeros(num_path, TxAntNum);

    temp_paths(:,layer) = paths_child;
    temp_paths(:,layer+1:TxAntNum) = paths(father,layer+1:TxAntNum);

    paths = temp_paths;
    PED = PED_child;
end


parfor cnt = 1:num_path
    path = paths(cnt,:);
    for layer = (TxAntNum-fulexplayers): -1 :1
        [ED,i] = min(abs(z(layer)- R(layer, layer) * symQAM - R(layer,layer+1:end)*path(layer+1:end).').^2);
        path(layer) = symQAM(i);
        PED(cnt) = PED(cnt) + ED;
    end
    paths(cnt,:) = path;
end

[~,i] = min(PED);
symest = paths(i,:)';

CC = RVMs * 32^2 + RVAs * 32;

%norm(z-R*paths(1,:).',2)/sqrt(10)
ED = PED(1,1);
Tx_sig_esti = paths(1,:);
end