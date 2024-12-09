function [symest, visited_nodes, flops] = BFBB(TxAntNum, RxAntNum, Hreal, RxSymbol, sym, norm_f)
% Breadth-first branch and bound algorithm, cited from Low "Complexity Detection Algorithms
% in Large-Scale MIMO Systems", TWC, 2016
% Paramters:
% L is predefined the maximum depth;
% M is the predefined maximum width.

 L=16;
 M=2;
 alpha=0.0001;

visited_nodes = 0; %Number of visited nodes in the tree each detection process
flops = 0;   %Number of operations

%###################### Problem formulation ###############################
Q =  Hreal'*Hreal;                                                            %  flops = flops + 8*TxAntNum^2*RxAntNum;
b = -Hreal'*(RxSymbol/norm_f + (length(sym)-1)*Hreal*ones(2*TxAntNum,1))/2;   %  flops = flops + 8*TxAntNum*RxAntNum;

% Initial box-constraint condition
lb0 = zeros(size(b)); 
ub0 = (length(sym)-1) * ones(size(b));

% ##################### Find the root solution ############################
visited_nodes = visited_nodes + 1;
% Solve BQP problem through interior-point method.
options = optimset('Display', 'off',...                                         
                    'Algorithm', 'interior-point-convex',...
                    'MaxIter', 200);
z_star0     = quadprog(Q, b, [], [], [], [],lb0, ub0, [], options);             %flops未知
z_star_quan = max(min(round(z_star0),ub0),lb0);    %Quantilize the solution

i=find(abs(z_star_quan - z_star0)>0.01, 1,'first'); %The branching variable
if isempty(i) %If the first solution is best, detection over.
    symest = (2*z_star_quan - (length(sym)-1))*norm_f;                          
    return;
end

lb_temp             = lb0;
ub_temp             = ub0;
lb_temp(i)          = ceil(z_star0(i));
ub_temp(i)          = floor(z_star0(i));
Solutions_container = []; %each node cantains z, lb, ub, and its branching variable of the subproblem
Solutions_container = [Solutions_container; {{}, lb_temp, ub0, {}}];
Solutions_container = [Solutions_container; {{}, lb0, ub_temp, {};}];

fup = inf; %Initialize upper bound of the fvalue of the best solution
l   = 0; %BrF search start from the root (layer 0)
while(l<L)
    l = l + 1;
    new_solcontner = cell(0,4);
    ml = size(Solutions_container,1);
    visited_nodes = visited_nodes + ml;
    fval_set = [];
    remain_idx = [];
    for m = 1:ml
        lb = cell2mat(Solutions_container(m,2));
        ub = cell2mat(Solutions_container(m,3));
    
        [z_star, fval_star] = quadprog(Q, b, [], [], [], [],lb, ub, [], options);
        if fval_star > fup
            continue;
        end
        z_star_quan = max(min(round(z_star),ub),0);
        fval_star_quan = 0.5*z_star_quan' * Q * z_star_quan + b' * z_star_quan; flops = flops + 4*TxAntNum^2 + 4*TxAntNum;
        
    
        i=find(abs(z_star_quan - z_star)>0.01, 1,'first'); %The first solution is best
        if isempty(i) || abs(fval_star - fval_star_quan) <= alpha * abs(fval_star)
            symest_best = (2*z_star_quan - (length(sym)-1))*norm_f;
            fup = fval_star;
        else
            fval_set = [fval_set; fval_star];
            remain_idx = [remain_idx; m];
            Solutions_container(m,1) = {z_star};
            Solutions_container(m,4) = {[]};
        end
    end

    if isempty(fval_set)
        break;
    elseif l==L
        [~,idx] = min(fval_set);
        z_star_quan = max(min(round(cell2mat(Solutions_container(remain_idx(idx),1))),cell2mat(Solutions_container(remain_idx(idx),3))),cell2mat(Solutions_container(remain_idx(idx),2)));
        symest_best = (2*z_star_quan - (length(sym)-1))*norm_f;
    elseif(2*length(fval_set)<=M)
        for j=1:length(remain_idx)
            lb = cell2mat(Solutions_container(remain_idx(j),2));
            ub = cell2mat(Solutions_container(remain_idx(j),3));
            i = cell2mat(Solutions_container(remain_idx(j),4));
            lb_temp = lb;
            ub_temp = ub;
            lb_temp(i) = min(ceil(z_star(i)), ub(i));
            ub_temp(i) = max(floor(z_star(i)), lb(i));
            new_solcontner = [new_solcontner; {{}, lb_temp, ub0, {}}];
            new_solcontner = [new_solcontner; {{}, lb0, ub_temp, {}}];
        end
    else
        [~, idx] = sort(fval_set);
        while(size(new_solcontner, 1)<M)
        for j=1:length(remain_idx)
            lb = cell2mat(Solutions_container(remain_idx(j),2));
            ub = cell2mat(Solutions_container(remain_idx(j),3));
            i = cell2mat(Solutions_container(remain_idx(j),4));
            lb_temp = lb;
            ub_temp = ub;
            lb_temp(i) = min(ceil(z_star(i)), ub(i));
            ub_temp(i) = max(floor(z_star(i)), lb(i));
            new_solcontner = [new_solcontner; {{}, lb_temp, ub0, {}}];
            new_solcontner = [new_solcontner; {{}, lb0, ub_temp, {}}];
        end
        end
    end
    Solutions_container = new_solcontner;

end
flops = flops + 4*TxAntNum^3*visited_nodes;
symest = symest_best;

end

