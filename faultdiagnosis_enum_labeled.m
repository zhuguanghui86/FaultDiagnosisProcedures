function [] = faultdiagnosis_enum_labeled(Pre, Post, M0, Tu, Tf, labelfun, w)
%   Pre, Post and M0 --- the pre-incidence, post-incidence and initial
%                        marking of a net, respectively.
%   Tu --- the set of unobservable transitions (including faulty transitions)
%          e.g., Tu=[2, 3, 5];
%   Tf --- the set of faulty transitions, e.g., Tf=[3,5];
%   labelfun --- the label function-A mapping (containers.Map) such as
%                                          a->{'t1','t3'},b->{'t2','t4'} 
%       for example, labelfun = contains.Map; labelfun('a') = {'t1','t3'};
%   w --- the observed sequence, a cell array, e.g., w = {'a', 'b', 'c'}

m = size(Pre, 1);
nu = length(Tu);
nf = length(Tf);
wlen = length(w);

C = Post - Pre;
Cu = C(:, Tu);

KSize = nf + 2;
KZero = zeros(1, KSize);

KM0 = KZero;
KM0(nf+1) = 1;
KM0(nf+2) = 1;
arrNodes = struct('M', M0, 'k', KM0);
%arrNodes = struct('M', {[0 1 0 1 0 0 0]',[0 0 0 1 1 0 0]',[0 0 0 1 0 1 0]',[0 0 0 1 0 0 1]',[0 0 1 1 0 0 0]'}, 'k', {[1 1 0],[0 3 3],[3 3 0],[0 2 2],[0 1 1]});
diagResult = zeros(1, nf+1);
for i = 1:wlen
    lab = w{i}; % the observed label
    
    G = digraph;
    nodenum = 0;
    
    arrLen = length(arrNodes);
    for j = 1:arrLen
        nodenum = nodenum + 1;
        G = addnode(G, num2str(nodenum));
    end
    
    j = 1;
    while j <= arrLen
        M = arrNodes(j).M;
        for tu = Tu
            if all(M >= Pre(:,tu)) % Ttu在M下可触发
                M1 = M + C(:,tu);
                ind = 0;
                len_1 = length(arrNodes);
                for k = 1:len_1 % 最耗时循环
                    if all(M1 == arrNodes(k).M)
                        ind = k;
                        break;
                    end
                end
                if ind % 已经存在结点
                    G = addedge(G, j, ind, tu);
                else
                    nodenum = nodenum + 1;
                    arrNodes(nodenum) = struct('M', M1, 'k', KZero);
                    G = addedge(G, j, nodenum, tu);
                end
            end
        end
        
        j = j + 1;
        arrLen = length(arrNodes);
    end
    
%    plot(G, 'EdgeLabel',G.Edges.Weight); %画图进行调试
    order = toposort(G);
    for j = order
        k1 = arrNodes(j).k;
        sucs = successors(G,j);
        sucs = sucs'; % 因为 for suc = sucs 这种循环方式只针对行向量，因此进行转置 
        for suc = sucs
            k = k1;
            tu = G.Edges.Weight(findedge(G,j,suc));
            at = find(Tf == tu, 1);
            mem = ismember(tu, Tf);
            
            k2 = arrNodes(suc).k;
            if all(k2 == 0) %还没有被计算
                if mem
                    k(at) = k(nf+1);
                    k(nf+2) = 0;
                    arrNodes(suc).k = k;
                else
                    arrNodes(suc).k = k;
                end
            else % 已经被计算过
                if mem
                    inds = 1:1:nf+1;
                    inds = setdiff(inds, at);
                    k2(inds) = k2(inds) + k(inds);
                    k2(at) = k2(at) + k(nf+1);
                    arrNodes(suc).k = k2;
                else
                    k2 = k2 + k;
                    arrNodes(suc).k = k2;
                end
            end
        end
    end
    
    newNodes = struct('M',{},'k',{});
    newNodesNum = 0;
    for j = 1:nodenum
        M = arrNodes(j).M;
        k = arrNodes(j).k;
        trans = labelfun(lab);
        transLen = length(trans);
        for q = 1:transLen
            t = trans{q};
            tnum = extracttnum(t);
            if all(M >= Pre(:,tnum))
                M1 = M + C(:,tnum);

                ind = 0;
                for s = 1:newNodesNum
                    if all(M1 == newNodes(s).M)
                        ind = s;
                        break;
                    end
                end
                if ind % 已经存在结点
                    newNodes(ind).k = newNodes(ind).k + k;
                else
                    newNodesNum = newNodesNum + 1;
                    newNodes(newNodesNum) = struct('M', M1, 'k', k);
                end
            end
        end
    end
    % 做诊断
    diagK = zeros(1,nf+2);
    for j = 1:newNodesNum
        diagK = diagK + newNodes(j).k;
    end
    total = diagK(nf+1);
    for j = 1:nf
        fv = diagK(j);
        if fv == 0
            diagResult(j) = 0;
        elseif fv == total
            diagResult(j) = 2;
        else
            diagResult(j) = 1;
        end
    end
    sv = diagK(nf+2);
    if sv == 0
        diagResult(nf+1) = 2;
    elseif sv == total
        diagResult(nf+1) = 0;
    else
        diagResult(nf+1) = 1;
    end
    
    disp(diagResult);
    arrNodes = newNodes;
end
end

