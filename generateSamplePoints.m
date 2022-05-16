% sampling points specifically for squar shape
% n: points on each side (total 4 sides x n points/side  = 4*n points)
% s: length of each side

function [p] = generateSamplePoints(n,s)
    p = zeros(4*n,2);
    dr = s/n;
    r = s/2;
    count = 1;
    for j = 1:n
        p(count,1) = r;
        p(count,2) = -r+dr*(j-1);
        count = count+1;
    end
    
    for j = 1:n
        p(count,1) = r-dr*(j-1);
        p(count,2) = r;
        count = count+1;
    end
    
    for j = 1:n
        p(count,1) = -r;
        p(count,2) = r-dr*(j-1);
        count = count+1;
    end
    
    for j = 1:n
        p(count,1) = -r+dr*(j-1);
        p(count,2) = -r;
        count = count+1;
    end
end