function [tree,E_optimal] = pruneBPTminimizeEnergy(tree,param,opt)

for i=1:length(tree)
    tree(i).pruning = 0;
end

Lambda = [tree.energy];
Lambda_plus = [Lambda.lambdaplus];
Lambda_moins = [Lambda.lambdamoins];
RootPartialEnergy = cell2mat(tree(end).energy.partialenergy);
LambdaVal = RootPartialEnergy(1:3:end);
GofFVal = RootPartialEnergy(3:3:end);
ReguVal = RootPartialEnergy(2:3:end);


if strcmp(opt,'unconstrained')
    
    lambda = param;
    fprintf('Unconstrained energy minimization with lambda = %g \n',lambda);
    optimalcut = find(Lambda_plus<=lambda&Lambda_moins>lambda);
    regstring = 'regions';
    if length(optimalcut)==1
        regstring = 'region';
    end
    fprintf('%d %s in the optimal partition.\n',length(optimalcut),regstring)
    lambda_min = LambdaVal(find(LambdaVal<=lambda,1,'last'));
    lambda_max = LambdaVal(find(LambdaVal>lambda,1,'first'));
    if isempty(lambda_max)
        lambda_max = Inf;
    end
    fprintf('Partition is optimal for interval [%g:%g]\n',lambda_min,lambda_max);
    %     [optimalcut,tree] = findoptimalcut(tree,lambda);
    
elseif strcmp(opt,'constrained')
    
    CostOpti = param;
    fprintf('Cost constrained energy minimization with cost = %g\n',CostOpti);
    fprintf('Minimal achievable cost: %g\n',ReguVal(end))
    if CostOpti<ReguVal(end)
        fprintf('Required cost is below minimal achievable cost ')
        fprintf(' => required cost is replaced by minimal achievable cost.\n')
        CostOpti = ReguVal(end);
    end
    CostIdx = find(ReguVal<=CostOpti,1,'first');
    Cost_val = ReguVal(CostIdx);
    GofF_val = GofFVal(CostIdx);
    lambda_min = LambdaVal(CostIdx);
    if CostIdx == length(ReguVal)
        lambda_max = Inf;
    else
        lambda_max = LambdaVal(CostIdx+1);
    end
    optimalcut = find(Lambda_plus<=lambda_min&Lambda_moins>lambda_min);
    
    regstring = 'regions';
    if length(optimalcut)==1
        regstring = 'region';
    end
    fprintf('Optimal partition achieves cost = %g with minimal fitting term of %g\n',Cost_val,GofF_val);
    fprintf('Optimal partition has %d %s',length(optimalcut),regstring);
    fprintf(' and is optimal for lambda in [%g:%g]\n',lambda_min,lambda_max);
    
elseif strcmp(opt,'regions')
    
    numReg = param;
    regstring = 'regions';
    if numReg==1
        regstring = 'region';
    end
    fprintf('Region constrained energy minimization with %d %s expected.\n',numReg,regstring);
    
    lambda = LambdaVal(end);
    optimalcut = find(Lambda_plus<=lambda&Lambda_moins>lambda);
    cmpt = 0;
    while length(optimalcut)<numReg
        cmpt = cmpt+1;
        lambda = LambdaVal(end-cmpt);
        optimalcut = find(Lambda_plus<=lambda&Lambda_moins>lambda);
    end
    % when exit while loop, optimalcut has >= than numReg regions
    
    if length(optimalcut)==numReg % if == numReg
        lambda_min = LambdaVal(find(LambdaVal<=lambda,1,'last'));
        
        lambda_max = LambdaVal(find(LambdaVal>lambda,1,'first'));
        if isempty(lambda_max)
            lambda_max = Inf;
        end
        fprintf('An optimal partition was found with %d regions.\n',numReg);
        fprintf('Corresponding optimal lambda values are in [%g:%g]\n',lambda_min,lambda_max);
        
    else % > numReg regions
        
        lambda_less = LambdaVal(end-cmpt+1);
        optimalcut_less = find(Lambda_plus<=lambda_less&Lambda_moins>lambda_less);
        numReg_less = length(optimalcut_less);
        
        if numReg-numReg_less<length(optimalcut)-numReg
            % Partition with less regions is closer to optimal than
            % partition with more regions
            lambda = lambda_less;
            optimalcut = optimalcut_less;
        end
        
        lambda_min = LambdaVal(find(LambdaVal<=lambda,1,'last'));
        lambda_max = LambdaVal(find(LambdaVal>lambda,1,'first'));
        if isempty(lambda_max)
            lambda_max = Inf;
        end
        fprintf('No optimal partition was found with %d regions.\n',numReg);
        fprintf('Closest optimal partition found has %d regions.\n',length(optimalcut))
        fprintf('Corresponding optimal lambda values are in [%g:%g]\n',lambda_min,lambda_max);
        
    end
    
else
    error('Wrong input string argument.')
end

fprintf('\n')
E_optimal = zeros(length(optimalcut),3);
for i=1:length(optimalcut)
    tree(optimalcut(i)).pruning = 1;
    E_optimal(i,1) = optimalcut(i);
    E_optimal(i,2:3) = (tree(optimalcut(i)).energy.properenergy)';
end

% Guigues' algorithm for energy minimization
% ------------------------------------------
%
% function [optimalcut,tree] = findoptimalcut(tree,lambda)
%
% Nbnodes = length(tree);
% for ii=1:Nbnodes
%
%
%     % step 1 of Guigues' algorithm has to be done already (each node has
%     % its proper energy computed)
%     properEnergy = tree(ii).energy.properenergy;
%     GofF = properEnergy(1);
%     Regu = properEnergy(2);
%
%     % step 2 of Guigues' algorithm
%     tree(ii).energy.optimalenergy = GofF+lambda*Regu;
%     tree(ii).energy.optimalcut = ii;
%     tree(ii).energy.childrenenergy = +Inf;
%
%     % step 3 of Guigues' algorithm
%     if ~isempty(tree(ii).nodeinfo.children)
%         c = tree(ii).nodeinfo.children;
%         Ec1 = tree(c(1)).energy.optimalenergy;
%         Ec2 = tree(c(2)).energy.optimalenergy;
%         childrenEnergy = Ec1 + Ec2;
% %         childrenEnergy = max(Ec1,Ec2); % uncomment if energy is max-separable
%         tree(ii).energy.childrenenergy = childrenEnergy;
%
%         % step 4 of Guigues' algorithm
%         if childrenEnergy < tree(ii).energy.optimalenergy
%             tree(ii).energy.optimalenergy = childrenEnergy;
%             optimalcut1 = tree(c(1)).energy.optimalcut;
%             optimalcut2 = tree(c(2)).energy.optimalcut;
%             tree(ii).energy.optimalcut = cat(2,optimalcut1,optimalcut2);
%         end
%     end
%
% end
%
% optimalcut = tree(end).energy.optimalcut;
