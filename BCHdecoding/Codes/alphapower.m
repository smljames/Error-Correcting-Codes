function y = alphapower(x, gf)
    y = zeros(1,length(x));
    for i = 1:length(x)
        if x(i) == 0
            y(i) = -1;
        else
            y(i) = find(gf==x(i))-2;
        end
    end
end