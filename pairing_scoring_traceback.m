function [score,path,F,pointer] = pairing_scoring_traceback(intseq1,m,intseq2,n,NN_dG,gap)
% dynamic algorithm for scoring and traceback of optimal pairing - local

%% initialize
% scoring matrix
init = -NN_dG(1,1);
F = zeros(n,m); 
F(:,1) = repmat(init,n,1);
F(1,:) = repmat(init,1,m);
currentFColumn = F(:,1);

% terminal AT penalty
if (intseq1(1) == 1 && intseq2(1) == 4) || (intseq1(1) == 4 && intseq2(1) == 1)
    F(1,1) = F(1,1)-NN_dG(1,2);
end

% back tracing pointer matrix
pointer(n,m) = uint8(0); 
ptr(n,1) = uint8(0); % column vector
path = [];
nnseq2 = [];

% column vector NN intseq2 
for i = 1:n-1 
    nnseq2 = [nnseq2; intseq2(i+1)*10+intseq2(i)];
end    

%% build scoring matrix F
% main loop runs through the matrix looking for maximal scores
for outer = 2:m
    
    % score current column
    scoredMatchColumn = -NN_dG(nnseq2,intseq1(outer-1)*10+intseq1(outer));
    
    % grab data from the matrices and initialize values
    lastFColumn = currentFColumn;
    currentFColumn = F(:,outer);
    best = F(1,outer);
           
    for inner = 2:n
        % scoring of three options
        up = best + gap; % the way it is used here, true value of gap is negative value of this
        left = lastFColumn(inner) + gap;   
        diagonal = lastFColumn(inner-1) + scoredMatchColumn(inner-1);

        % find highest non-negative score
        best = init;
        pos = 0;
         
        if up > init && up > left
            best = up ;
            pos = 2; % UP
        elseif left > init
            best = left;
            pos = 3; % LEFT
        end
        if diagonal >= best
            best = diagonal;
            ptr(inner) = 1; % DIAGONAL (is preferred)
        else
            ptr(inner) = pos;
        end

        currentFColumn(inner) = best;
        
    end % inner
    
    % save columns of data 
    F(:,outer) = currentFColumn;
    pointer(:,outer) = ptr;
end

% terminal AT penalty
if (intseq1(end) == 1 && intseq2(end) == 4) || (intseq1(end) == 4 && intseq2(end) == 1)
    F(end,end) = F(end,end)-NN_dG(1,2); 
end

%% traceback
% find highest scores in scoring matrix
score = max(max(F));

% return one pairing - best local alignment
[nn,mm] = find(F == score,1,'last');

% path from pointers
path = [mm, nn];
while 1 %F(nn,mm) > init
    switch pointer(nn,mm);
        case 0
            break
        case 1 % DIAG
            mm = mm-1;
            nn = nn-1;
            path = [path; mm,nn];
        case 2 % UP  
            nn = nn-1;
            path = [path; 0,nn]; 
        case 3 % LEFT
            mm = mm-1;
            path = [path; mm,0];
    end
end
path = [path; path(end,1)-1,path(end,2)-1]; % account for N in first NN

if isempty(path) == 1
   path = [0 0];
end

end