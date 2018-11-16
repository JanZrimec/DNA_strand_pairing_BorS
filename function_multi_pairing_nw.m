function [output] = function_multi_pairing_nw(seq,gap)
% parallel function to call pairing_main with multiple sequences

tic

lsiz = size(seq,1);  

output = zeros(lsiz,lsiz,6);
load('NN_dG.mat')
load('NN_dH.mat')
load('NN_dS.mat')
G = int16(NN_dG*100); % pass files to parrallel workers
H = int16(NN_dH*10);
S = int16(NN_dS*10);

% main loop
for i=1:lsiz
    i
    slice = seq(i,:);
    
    parfor j=i:lsiz % upper triangular
        [output(i,j,1),dH(j),dS(j),bp(j),loop(j),dangle(j)] = pairing_main_NW_par(slice,seq(j,:),G,H,S,gap);
    end
    
    for j = i:lsiz % edit sliced output outside of parfor
        output(i,j,2)=dH(j);
        output(i,j,3)=dS(j);
        output(i,j,4)=bp(j);
        output(i,j,5)=loop(j);
        output(i,j,6)=dangle(j);
    end
end

toc

% trasnform over diagonal
for i=1:lsiz
    for j=i:lsiz
        
        output(j,i,:)=output(i,j,:);
        
    end
end

end
