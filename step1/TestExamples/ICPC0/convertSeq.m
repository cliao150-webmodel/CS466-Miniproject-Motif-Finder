function SeqStr = convertSeq(numArray)

SeqStr = cell(1,length(numArray));
for i=1:length(numArray)
    switch numArray(i)
        case 1
            SeqStr{i} = 'A';
        case 2
            SeqStr{i} = 'C';
        case 3
            SeqStr{i} = 'G';
        case 4
            SeqStr{i} = 'T';
    end
end

end

