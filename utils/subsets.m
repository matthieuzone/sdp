function powerSet = subsets(bucket)
%SUBSETS Generate the power set of a given bucket.

    % Convert the bucket to a row vector to ensure consistent subset orientation
    bucket = bucket(:).';
    n = numel(bucket);
    numSubsets = 2^n;
    powerSet = cell(numSubsets, 1);
    
    for i = 0:numSubsets - 1
        % Convert the current index to a binary string of length n
        binStr = dec2bin(i, n);
        % Use the binary string to select elements from the bucket
        subset = bucket(binStr == '1');
        powerSet{i + 1} = subset;
    end
end