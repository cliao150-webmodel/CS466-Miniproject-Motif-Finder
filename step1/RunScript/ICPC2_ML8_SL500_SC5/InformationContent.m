function dx = InformationContent(x,ICPC)

dx = sum(abs(x).*log2(abs(x)))+(1-sum(abs(x))).*log2(1-sum(abs(x)))+2-ICPC;
%dx = sum(x(1:end-1).*log2(x(1:end-1)))+(1-sum(x(1:end-1)))*log2(1-sum(x(1:end-1)))+2-ICPC;

end

