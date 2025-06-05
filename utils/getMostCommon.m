function mostCommon=getMostCommon(x)
    if iscell(x)
        u=unique(x);
        counts=zeros(size(u));
        for k=1:length(u)
            counts(k)=sum(strcmp(x,u{k}));
        end
    else
        u=unique(x);
        counts=zeros(size(u));
        for k=1:length(u)
            counts(k)=sum(x==u(k));
        end
    end
    [~,idx]=max(counts);
    mostCommon=u(idx);
end