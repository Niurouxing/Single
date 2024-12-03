function Opt = GML(R, sym, z) 

    % conSize = length(sym);
    % currentSurvive=[];
    % for step=1:length(z)
    %  expandedPath=[repmat(1:conSize,1,max(size(currentSurvive,2),1));kron(currentSurvive,ones(1,conSize))];
    %  PED=sum((z(end-step+1:end)- R(end-step+1:end,end-step+1:end)*sym(expandedPath)).^2,1);
    %  [~,choice]=mink(PED,k);
    %  currentSurvive=expandedPath(:,choice);
    % end
    % Opt= sym(currentSurvive(:,1))';

    
    L = 5;  
    ConSize = length(sym);

    Tx = size(z,1);
    stage = Tx / L;

    accumulated_path = [];

    % 生成所有可能组合
    combination = permute(sym(fullfact(repmat(ConSize, 1, L)))',[3,1,2]);

    for i = 1:stage

        % choose the target R area
        target_area = R(Tx - i * L + 1: Tx - (i - 1) * L, Tx - i * L + 1: Tx - (i - 1) * L);
        
        % reserve area for all accumulated path
        reserve_area = R(   Tx - i * L + 1: Tx - (i - 1) * L,Tx - (i - 1) * L + 1: Tx);

        % calculate the dot product of the target area and combination
        
        new_PED = squeeze(sum( target_area .* combination,2));

       

        % calculate the dot product of the reserve area and accumulated path
        accumulated_PED = sum(reserve_area * accumulated_path,2);

        % calculate the total PED
        total_PED =sum(( z(Tx - i * L + 1: Tx - (i - 1) * L) - new_PED - accumulated_PED).^2,1);

        % find the minimum PED
        [~,choice] = min(total_PED);

        % update the accumulated path
        accumulated_path = [ combination(1,:,choice)'; accumulated_path];

    end


    % 如果未能整除的情况
    remain = mod(Tx, L);
    if remain ~= 0
        % 生成新的combination
        combination = permute(sym(fullfact(repmat(ConSize, 1, remain)))',[3,1,2]);

        % choose the target R area
        target_area = R(1:remain, 1:remain);

        % reserve area for all accumulated path
        reserve_area = R(1:remain, remain + 1: Tx);

        % calculate the dot product of the target area and combination
        new_PED = squeeze(sum( target_area .* combination,2));

        % calculate the dot product of the reserve area and accumulated path
        accumulated_PED = sum(reserve_area * accumulated_path,2);

        % calculate the total PED
        total_PED =sum(( z(1:remain) - new_PED - accumulated_PED).^2,1);

        % find the minimum PED
        [~,choice] = min(total_PED);

        % update the accumulated path

        accumulated_path = [ combination(1,:,choice)'; accumulated_path];

    end

    
    Opt = accumulated_path;









end