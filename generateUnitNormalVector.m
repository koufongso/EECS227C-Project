function [u] = generateUnitNormalVector(p)
    x = p(1);
    y = p(2);
    if(abs(x)==abs(y))
        u = [sign(x)/sqrt(2),sign(y)/sqrt(2)];
    elseif(abs(x)>abs(y))
        u = [sign(x);0];
    else
        u = [0;sign(y)];
    end
end