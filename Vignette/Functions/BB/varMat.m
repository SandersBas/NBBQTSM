function varcovar = varMat(weight, psi, X, theta)

    result = arrayfun(@(i) weight(i).*psi(X(i,:)', theta), 1:size(X,1), 'UniformOutput', false);
    result = cell2mat(result)';
        
    dev = result - kron(ones(size(X,1),1),sum(result));
    
    varcovar = 0;
    for i = 1:size(X,1)
        varcovar = varcovar + weight(i).*dev(i,:)'*dev(i,:);
    end
end
