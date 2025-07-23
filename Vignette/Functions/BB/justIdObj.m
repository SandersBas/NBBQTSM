function objValue = justIdObj(weight,psi,X,theta)

    result = arrayfun(@(i) weight(i).*psi(X(i,:)', theta), 1:size(X,1), 'UniformOutput', false);
    result = cell2mat(result)';

    objValue = sum(result)';
end