function [] = faultdiagnosis_BRM_labeled(Pre, Post, M0, Tu, Tf, labelfun, w)
%   Pre, Post and M0 --- the pre-incidence, post-incidence and initial
%                        marking of a net, respectively.
%   Tu --- the set of unobservable transitions (including faulty transitions)
%          e.g., Tu=[2, 3, 5];
%   Tf --- the set of faulty transitions, e.g., Tf=[3,5];
%   labelfun --- the label function-A mapping (containers.Map) such as
%                                          a->{'t1','t3'},b->{'t2','t4'} 
%       for example, labelfun = containers.Map; labelfun('a') = {'t1','t3'};
%   w --- the observed sequence, a cell array, e.g., w = {'a', 'b', 'c'}

[m] = size(Pre, 1);
nu = length(Tu);
nf = length(Tf);
wlen = length(w);

C = Post - Pre;
Cu = C(:, Tu);

lastTotalBRMs = {struct('M', M0, 'y', zeros(1, nu))};

for i = 1:wlen
    lab = w{i};
    trans = labelfun(lab);
    transLen = length(trans);
    structArrLen = length(lastTotalBRMs);
    newTotalBRMs = {};
    lastDiagResult = zeros(1, nf+1);
    lastDiagUseFirst = 1;
    for p = 1:structArrLen
        diagResult = zeros(1, nf+1);
        lastBRMs = lastTotalBRMs{p};
        for q = 1:transLen
            t = trans{q};
            tnum = extracttnum(t);
            newBRMs = struct('M', {}, 'y', {});
            brmNum = 0;
            oldlen = length(lastBRMs);
            quitCurr = 0;
            for j = 1:oldlen
                M = lastBRMs(j).M;
                y = lastBRMs(j).y;
                ymin = getYminMt(Pre, Post, Tu, M, t);
                if(isempty(ymin))
                    quitCurr = 1;
                end
                yrow = size(ymin, 1);
                for k = 1:yrow
                    ele = ymin(k,:);
                    brmNum = brmNum + 1;
                    newBRMs(brmNum).M = M + Cu * ele';
                    newBRMs(brmNum).y = y + ele;
                end
            end
            if quitCurr
                continue;
            end
            % 去除 newBRMs 中大的元素
            newlen = length(newBRMs);
            cflag = zeros(1, newlen);
            zeroEle = find(cflag == 0, 1);
            while ~isempty(zeroEle)
                cflag(zeroEle) = 1;
                delInd = [];
                yitem = newBRMs(zeroEle).y;
                for j = 1:newlen
                    if j == zeroEle
                        continue;
                    end
                    if all(newBRMs(j).y >= yitem)
                        delInd = [delInd, j];
                    end
                end
                %删除元素
                for j = delInd
                    newBRMs(j) = [];
                    cflag(j) = [];
                end
                %重置变量
                newlen = length(newBRMs);
                zeroEle = find(cflag == 0, 1);
            end
            %开始做诊断
            model.A = sparse(Cu);
            model.rhs = Pre(:,tnum) - M;
            model.sense = '>';
            model.vtype = 'I';

            params.outputflag = 0;
            TfAt = zeros(1,nf);
            for j = 1:nf
                at = find(Tu==Tf(j));
                TfAt(j) = at;
                fvalue = zeros(1, newlen);
                for k = 1:newlen
                    tmpY = newBRMs(k).y;
                    fvalue(k) = tmpY(at);
                end
                if all(fvalue > 0)
                    diagResult(j) = 2; % 2表示 F
                elseif any(fvalue > 0)
                    diagResult(j) = 1; % 1表示 U
                else % fvalue的每一项都是0
                    mod = model;
                    addrow = zeros(1, nu);
                    addrow(at) = 1;
                    mod.A = [mod.A; sparse(addrow)];
                    mod.rhs = [mod.rhs; 1];
                    result = gurobi(mod, params);
                    if ~strcmp(result.status, 'OPTIMAL') %无解
                       diagResult(j) = 0; % 0 表示 N
                    else
                        diagResult(j) = 1; % 1表示 U
                    end
                end
            end
            %做全局诊断
            diagF = diagResult(1:1:nf);
            if all(diagF == 0)
                diagResult(nf+1) = 0;
            elseif any(diagF == 2)
                diagResult(nf+1) = 2;
            elseif (sum(diagF == 1)) == 1 % 有且只有一个 fU
                diagResult(nf+1) = 1;
            else
                allZero = [];
                for j = 1:newlen
                    tmpY = newBRMS(j).y;
                    if all(tmpY == 0)
                        allZero = [allZero, j];
                    end
                end
                if isempty(allZero)
                    diagResult(nf+1) = 2;
                else
                    rowf = zeros(1, nu);
                    rowf(TfAt) = 1;
                    model2.A = sparse([Cu;rowf]);
                    model2.vtype = 'I';
                    for j = allZero
                        tmpM = newBRMS(j).M;
                        model2.rhs = [Pre(:,tnum) - tmpM; 0];
                        model2.sense = [repmat('>', m, 1); '='];
                        result = gurobi(model2, params);
                        if ~strcmp(result.status, 'OPTIMAL') %无解
                           diagResult(nf+1) = 2; % 0 表示 N
                        else
                            diagResult(nf+1) = 1; % 1表示 U
                            break; % 只要找到一个解就表示总结果有可能正常
                        end
                    end
                end
            end
            if lastDiagUseFirst
                lastDiagResult = diagResult;
                lastDiagUseFirst = 0;
            else
                lastDiagResult = andOperation(lastDiagResult, diagResult);
            end
            for j = 1:newlen
                newBRMs(j).M = newBRMs(j).M + C(:,tnum);
            end
            newTotalBRMs = [newTotalBRMs, newBRMs];
        end
    end
   
    lastTotalBRMs = newTotalBRMs;
    disp(lastDiagResult);
end

end

function [andResult] = andOperation(diagRet1, diagRet2)
    andTab = [0,1,1;1,1,1;1,1,2];
    len = length(diagRet1);
    andResult = zeros(1, len);
    for i = 1:len
        andResult(i) = andTab(diagRet1(i)+1, diagRet2(i) + 1);
    end
end

