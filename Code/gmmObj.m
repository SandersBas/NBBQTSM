function objValue = gmmObj(weight, psi, X, theta, W)

    result = arrayfun(@(i) weight(i).*psi(X(i,:)', theta), 1:size(X,1), 'UniformOutput', false);
    result = cell2mat(result)';
    
    if size(result,1)==1
        objValue = result*W*result';
    else
        objValue = sum(result)*W*sum(result)';
    end
end