function n = degree(a)
    n = length(a);
    i=1;
    while(a(i)==0)
        i = i+1;
    end
    n = n - i;
end