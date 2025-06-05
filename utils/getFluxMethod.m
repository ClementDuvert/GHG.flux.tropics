function out=getFluxMethod(x)
    x=x(~isundefined(x));
    if isempty(x)
        out=categorical("undefined",["measured","modelled","undefined"]);
    elseif any(contains(string(x),'measured'))
        out=categorical("measured",["measured","modelled","undefined"]);
    else
        out=categorical("modelled",["measured","modelled","undefined"]);
    end
end
