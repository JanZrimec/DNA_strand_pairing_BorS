function [score,path,F,pointer] = pairing_scoring_traceback_NW(intseq1,m,intseq2,n,NN_dG,gap)
% dynamic algorithm for scoring and traceback of optimal pairing - global

%% initialize
% scoring matrix
init = -NN_dG(1,1);
F = zeros(n,m); 
F(1,1) = init;
F(2:end,1) = gap*(1:n-1)'+init;
F(1,2:end) = gap*(1:m-1)+init;
currentFColumn = F(:,1);

% terminal AT penalty
if (intseq1(1) == 1 && intseq2(1) == 4) || (intseq1(1) == 4 && intseq2(1) == 1)
    F(1,1) = F(1,1)-NN_dG(1,2);
end

% back tracing pointer matrix
pointer= repmat(uint8(4),n,m);
pointer(:,1) = 2;  % up
pointer(1,1) = 1;
ptr = pointer(:,2); % ptr(1) is always 4
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
    best = currentFColumn(1);
           
    for inner = 2:n
        % scoring of three options
        up = best + gap; % the way it is used here, true value of gap is negative value of this
        left = lastFColumn(inner) + gap;   
        diagonal = lastFColumn(inner-1) + scoredMatchColumn(inner-1);

        % find highest non-negative score      
        if up > left
            best = up;
            pos = 2;
        else
            best = left;
            pos = 4;
        end
        if diagonal >= best % diagonal path preffered in pointer
            best = diagonal;
            ptr(inner) = 1;
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
% score at end of scoring matrix
score=F(end,end);
mm=m;
nn=n;

% path from pointers
path = [mm, nn];
while (mm-1 & nn-1) > 0
    switch pointer(nn,mm); % diagonal path is preferred
        case 1 % DIAG
            mm = mm-1;
            nn = nn-1;
            path = [path; mm,nn];
        case 2 % UP  
            nn = nn-1;
            path = [path; 0,nn]; 
        case 4 % LEFT
            mm = mm-1;
            path = [path; mm,0];
    end
end
path = [path; path(end,1)-1,path(end,2)-1]; % account for N in first NN

if isempty(path) == 1
   path = [0 0];
end

end