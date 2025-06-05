function out=getHighestFreq(freqs,rank_map)
    freqs=unique(freqs); 
    ranks=arrayfun(@(f) rank_map(f),freqs); 
    [~,idx]=min(ranks);
    out=freqs(idx);
end
