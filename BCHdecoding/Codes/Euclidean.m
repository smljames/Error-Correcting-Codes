function [r, h] = Euclidean(a, b, f, g)
    while(b(1)==0)
        b = b(2:end);
    end
    [q, r] = deconv(a, b);
    tmp = conv(q, g);
    f = [zeros(1, length(tmp)-length(f)), f];
    h = f - conv(q, g);
    
    if degree(r) >= degree(h)
        [r, h] = Euclidean(b, r, g, h);
    end
    while(r(1)==0)
        r = r(2:end);
    end
    while(h(1)==0)
        h = h(2:end);
    end
end