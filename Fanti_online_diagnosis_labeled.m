function [times, counter] = Fanti_online_diagnosis_labeled(Pre, Post, M0, Tu, Tf, labelfun, w)
% SETTING: The GUROBI solver has to be installed correctly on your computer
%          and the connection of the GUROBI solver with the MATLAB should
%          be corretly configured.

%OUTPUT:
%   times --- a row vector in which each element denotes the running time
%             of making a diagnsis for each corresponding event.
%   counter --- a row vector in which each element denotes the number of
%             times to solve the ILP problem for the diagnsos of the corresponding event.
%   Note that the diagnosis result will be printed on the screen.

%INPUT:
%   Pre, Post and M0 --- the pre-incidence, post-incidence and initial
%                        marking of a net, respectively.
%   Tu --- the set of unobservable transitions (including faulty transitions)
%          e.g., Tu=[2, 3, 5];
%   Tf --- the set of faulty transitions, e.g., Tf=[3,5];
%   labelfun --- the label function--A mapping (containers.Map). If we have 
%                                          a->{'t1','t3'},b->{'t2','t4'}, then 
%         labelfun = contains.Map; labelfun('a') = {'t1','t3'}; labelfun('b') = {'t2','t4'};
%   w --- the observed word, a cell array, e.g., w = {'a', 'b', 'c'}

nf = length(Tf);
seqLen = length(w);
A = {{}};
for i = 1:seqLen
    tic;
    my = 0;
    lastDiagResult = zeros(1, nf+1);
    lastDiagUseFirst = 1;
    ANew = {};
    label = w{i};
    trans = labelfun(label);
    transLen = length(trans);
    ALen = length(A);
    for j = 1:transLen
        tran = trans{j};
        for k=1:ALen
            tranSeq = A{k};
            tranSeq = [tranSeq, tran];
            [ diagResult, stop, c] = faultdiagnosis_ILP_forFanti(Pre, Post, M0, Tu, Tf, tranSeq);
            my = my + c;
            if stop == 0
                temp = {tranSeq};
                ANew = [ANew, temp];
                if lastDiagUseFirst
                    lastDiagResult = diagResult;
                    lastDiagUseFirst = 0;
                else
                    lastDiagResult = andOperation(lastDiagResult, diagResult);
                end
            end
        end
    end
    
    A = ANew;
    times(i) = toc;
    counter(i) = my;
    %disp('ResultResult:');
    disp(lastDiagResult);
end

end

%***************************************************************************************************
function [diagResult, stop, ILPCounter] = faultdiagnosis_ILP_forFanti(Pre, Post, M0, Tu, Tf, w)
%OUTPUT:
% diagResult --- a row vector [df1 df2 df3 dF], where dfi denotes the diagnosis of fi and dF the diagnosis of overall system status
% stop --- if the ILP problem has no solution, then set stop = 1; otherwise stop = 0. See "Fault detection by labeled Petri nets in centralized and distributed approaches" for the details.
% ILPCounter --- the number of times to solve the ILP problem

% INPUT:
%   Pre, Post and M0 --- the pre-incidence, post-incidence and initial
%                        marking of a net, respectively.
%   Tu --- the set of unobservable transitions (including faulty transitions)
%          e.g., Tu=[2, 3, 5];
%   Tf --- the set of faulty transitions, e.g., Tf=[3,5];
%   w --- the observed transition sequence, a cell array, e.g., w = {'t1', 't2', 't5'}

[m,n] = size(Pre);
nu = length(Tu);
nf = length(Tf);
wlen = length(w);

C = Post - Pre;
Cu = C(:, Tu);

stop = 0;
ILPCounter = 0;

b = zeros(m, 1); 
rh = []; 
lastTNum = 0;
zeroMat = zeros(m,nu);
%%%%%%%%%%%%% the first observable transition %%%%%%%%%%%%%%
t = w{1};
tnum = extracttnum(t);
lastTNum = tnum;

b = b - M0;
rh = b + Pre(:, tnum);

nVars = nu * 1;

model.A = sparse(Cu); % see the MATLAB documentation of GURIBO solver for the details.
model.vtype = repmat('I', nVars, 1);
model.sense = '>';
model.lb = zeros(nVars, 1);
model.ub = Inf(nVars, 1);
model.rhs = rh; 

for i = 2:wlen
    b = b - C(:, lastTNum);
    t = w{i};
    tnum = extracttnum(t);
    lastTNum = tnum;
    
    rh = [rh;zeros(m,1)];
    rh(m*(i-1)+1:end) = b + Pre(:, tnum); 
    
    nVars = nu * i;
    model.A = [model.A, sparse(repmat(zeroMat, i-1, 1))];
    model.A = [model.A;sparse(repmat(Cu, 1, i))];
    model.vtype = repmat('I', nVars, 1);
    model.sense = '>';
    model.lb = zeros(nVars, 1);
    model.ub = Inf(nVars, 1);
    model.rhs = rh; 
end
    
params.outputflag = 0;
diagResult = zeros(1,nf+1);
TfAt = zeros(1,nf);
for j = 1:nf
    model.modelsense = 'max';
    objective = zeros(1, nVars);

    ind = 1:nu:i*nu;
    at = find(Tu==Tf(j));
    TfAt(j) = at;
    ind = ind + (at-1);

    objective(ind) = 1;
    model.obj = objective;
    ILPCounter = ILPCounter + 1;
    result = gurobi(model, params);
    if ~strcmp(result.status, 'OPTIMAL')
        %fprintf('Error: fixed model is not optimal\n');
        stop = 1;
        return;
    end
    val = result.objval;
    if val == 0 %f must not occur
        diagResult(j) = 0;
    else % f may occur and further test is needed.
        model.modelsense = 'min';
        ILPCounter = ILPCounter + 1;
        result = gurobi(model, params);
        if ~strcmp(result.status, 'OPTIMAL')
            %fprintf('Error: fixed model is not optimal\n');
            stop = 1;
            return;
        end
        val = result.objval;
        if val == 0 % f may occur.
            diagResult(j) = 1;
        else % f must occur
            diagResult(j) = 2; 
        end
    end
end
%===============the diagnosis of overall system status =============
diagF = diagResult(1:1:nf);
if all(diagF == 0)
    diagResult(nf+1) = 0;
elseif any(diagF == 2)
    diagResult(nf+1) = 2;
elseif (sum(diagF == 1)) == 1 % there is only one fault transition whose diagnosis is 1
    diagResult(nf+1) = 1;
else
    model.modelsense = 'min';
    objective = zeros(1, nVars);
    objective(TfAt) = 1;
    model.obj = objective;
    ILPCounter = ILPCounter + 1;
    result = gurobi(model, params);
    if ~strcmp(result.status, 'OPTIMAL')
        %fprintf('Error: fixed model is not optimal\n');
        stop = 1;
        return;
    end
    val = result.objval;
    if val == 0
        diagResult(nf+1) = 1;
    else
        diagResult(nf+1) = 2;
    end
end

%============================
%disp(diagResult);
end






%***************************************************************
function [andResult] = andOperation(diagRet1, diagRet2)
    andTab = [0,1,1;1,1,1;1,1,2];
    len = length(diagRet1);
    andResult = zeros(1, len);
    for i = 1:len
        andResult(i) = andTab(diagRet1(i)+1, diagRet2(i) + 1);
    end
end

