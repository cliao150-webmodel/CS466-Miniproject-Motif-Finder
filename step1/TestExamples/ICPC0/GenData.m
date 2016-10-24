function GenData(ICPC,ML,SL,SC,FileNo)

%   ICPC: information content per column
%   ML: motif length
%   SL: sequence length
%   SC: sequence count

assert(ICPC > -1e-10 && ICPC < 2.0 + 1e-10);
assert(ML > 0);
assert(SL >= ML);
assert(SC > 0);

%%  Generate SC random sequences (with uniform nucleotide frequencies). Each random sequence has length SL.
%   A(1), C(2), G(3), T(4)
Seq = zeros(SC,SL);
for i=1:SC
    Seq(i,:) = ceil(rand(1,SL)*4);
end

%%  Generate a random motif (position weight matrix) of length ML, with total information content being ICPC*ML
options = optimset('Display','off','Algorithm','Levenberg-Marquardt');
while 1
    x0 = [rand/4,rand/4,rand/4];
    [x,fval,exitflag] = fsolve(@InformationContent,x0,options,ICPC);
    
    if (exitflag > 0 && abs(fval) < 1e-6)
        break;
    end
    
end
profile = [abs(x),1-sum(abs(x))];  %   probability of base (A, C, G, T) with information content per column ICPC

%%  Generate SC binding sites by sampling from this random motif
bsites = zeros(SC,ML);
freq = zeros(ML,4);
for i=1:SC
    [bsites(i,:),count] = RandomNucleotide(rand(1,ML),profile);
    freq = freq + count;
end

%%  Plant one sampled site at a random location in each random sequence generated in step 2
startpos = zeros(1,SC);
for i=1:SC
    startpos(i) = ceil(rand*(SL-ML+1));
    Seq(i,startpos(i):startpos(i)+ML-1) = bsites(i,:);
end

%%  Write out the SC sequences into a FASTA format file called sequences.fa
for i=1:SC
    data(i).Sequence = char(convertSeq(Seq(i,:)))';
    data(i).Header   = strcat('seq',num2str(i));
end
fastawrite('sequences.fa', data);
%type('sequences.fa')

%%  In a separate text file (called "sites.txt") write down the location of the planted site in each sequence
fileID = fopen('sites.txt','w');
for i=1:SC
    fprintf(fileID,'seq%s\t\t%d\n',num2str(i),startpos(i));
end
fclose(fileID);

%%  In a separate text file (called "motif.txt") write down the motif that was generated in step 3.
fileID = fopen('motif.txt','w');
fprintf(fileID,'>MOTIF%d\t%d\n',FileNo,ML);
fclose(fileID);
% dlmwrite('motif.txt',freq,'delimiter','\t','-append');
dlmwrite('motif.txt',freq,'-append','delimiter','\t');
fileID = fopen('motif.txt','a');
fprintf(fileID,'<');
fclose(fileID);

%%  In a separate test file (called "motiflength.txt") write down the motif length
fileID = fopen('motiflength.txt','w');
fprintf(fileID,'%d',ML);
fclose(fileID);

end
