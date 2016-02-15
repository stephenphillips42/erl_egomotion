function value = definedOrDefault_long(name,default,args_in)
    ind = (find(strcmp(args_in,name)));
    if isempty(ind)
        value = default;
    else
        value = args_in{ind+1};
    end
end

