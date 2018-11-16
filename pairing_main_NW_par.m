function [dG,dH,dS,bp,loop,dangle,data] = pairing_main_NW_par(seq1,seq2,NN_dG,NN_dH,NN_dS,gap)
% pairing prediction between input sequences seq1 and seq2 based on NN
% parameters (basepairs and internal mismatches) and gap parameters

NN_dG = double(NN_dG)/100;
NN_dH = double(NN_dH)/10;
NN_dS = double(NN_dS)/10;

% transform input sequences to integers
intseq1 = nt2int(seq1);
intseq1(intseq1(:) > 4) = 5; 
intseq2 = nt2int(seq2);
intseq2(intseq2(:) > 4) = 5;
m = length(seq1) 
n = length(seq2) % sequences of equal length

% call function to obtain scoring matrix and optimal pairing path
[score,path,F,pointer] = pairing_scoring_traceback_NW(intseq1,m,intseq2,n,NN_dG,gap);
    
% pairing of reversed sequences
[score2,path2,F2,pointer2] = pairing_scoring_traceback_NW(intseq1(end:-1:1),m,intseq2,n,NN_dG,gap);

if score2 < score
    score = score2;
    path = path2;
    F = F2;
    pointer = pointer2;
end

% create & visualize dimer in the area of sequence overlap
sp = size(path,1);
dimer = repmat(('- -')',1,sp);
path = flipud(path);

dimer(1,path(:,1)>0) = seq1(path(path(:,1)>0,1)); % fill where not 0 = gap
dimer(3,path(:,2)>0) = seq2(path(path(:,2)>0,2));

for i = 1:sp 
    if path(i,1) && path(i,2) > 0
        if intseq1(path(i,1)) == -double(intseq2(path(i,2)))+5
           dimer(2,i) = '|';
        end   
    end
end

% offset on both sides of dimer
offset(1) = abs(path(find(path(:,1)>0,1,'first'),1)-path(find(path(:,2)>0,1,'first'),2));
offset(2) = abs((m-path(find(path(:,1)>0,1,'last'),1))-(n-path(find(path(:,2)>0,1,'last'),2)));

% calculate entalpy and entropy of dimer
% melting temperature Tm = 1000*dH/(dS-18.3028)-273.15
intdim(1,:) = nt2int(dimer(1,:)); 
intdim(1,intdim(1,:) > 4) = 5;    % 5 = gap
intdim(2,:) = nt2int(dimer(3,:));
intdim(2,intdim(2,:) > 4) = 5;
dH = NN_dH(1,1); % initialize
dS = NN_dS(1,1);
% Watson-Crick basepairs and internal mismatches
for i = 1:sp-1 
   dH = dH+NN_dH(10*intdim(1,i)+intdim(1,i+1),10*intdim(2,i+1)+intdim(2,i));
   dS = dS+NN_dS(10*intdim(1,i)+intdim(1,i+1),10*intdim(2,i+1)+intdim(2,i));
end
% terminal AT penalties
if offset(1) == 0 && ((intdim(1,1) == 1 && intdim(2,1) == 4) || (intdim(1,1) == 4 && intdim(2,1) == 1))
    dH = dH+NN_dH(1,2); % terminal AT penalty
    dS = dS+NN_dS(1,2);
end    
if offset(2) == 0 && ((intdim(1,end) == 1 && intdim(2,end) == 4) || (intdim(1,end) == 4 && intdim(2,end) == 1))
    dH = dH+NN_dH(1,2);
    dS = dS+NN_dS(1,2);
end
dG = -score;

% count loops
bindim = zeros(size(dimer));
bindim(1,:) = intdim(1,:) < 5;
bindim(3,:) = intdim(2,:) < 5;
bindim(2,:) = dimer(2,:) == '|';
dimloop = sum(bindim(:,find(bindim(2,:)==1,1,'first'):find(bindim(2,:)==1,1,'last')),1);
dimloop(dimloop==3) = 0;

% calculate ratio of bp in percent
tmp1 = 0.5*(m+n);
tmp2 = 0.5*(offset(1)+offset(2));
bp = sum(bindim(2,:))/tmp1;

% calculate ratio of loops in percent
if isempty(dimloop) == 0
    loop = 0.5*sum(dimloop)/(tmp1-tmp2);
else
    loop = 0;
end

% calculate ratio of dangling ends = 1-overlap
dangle = tmp2/tmp1;

% inner data
data.path = path;
data.F = F;
data.pointer = pointer;
data.dimer = dimer;
end