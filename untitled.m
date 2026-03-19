%% EEL 6683 - Homework #5 - Problem 4
% Distributed Graph Analysis: Neighbor Structure + SCCs
% Digraph edges from the figure:
%   1->2, 1->3
%   2->3, 2->4
%   3->7
%   4->8
%   5->7
%   6->1, 6->5
%   7->6, 7->8
%   8->5

clear; clc;

%% ============================================================
%  PART 0: Build the adjacency matrix
% =============================================================
n = 8;  % number of nodes

% Adjacency matrix A(i,j) = 1 means edge i->j
A = zeros(n, n);
edges = [1 2; 1 3; 2 3; 2 4; 3 7; 4 8; 5 7; 6 1; 6 5; 7 6; 7 8; 8 5];
for e = 1:size(edges, 1)
    A(edges(e,1), edges(e,2)) = 1;
end

fprintf('Adjacency Matrix A:\n');
disp(A);

%% ============================================================
%  PART 1: Neighboring Structure
%  Out-neighbors, In-neighbors, and Shortest Distances
% =============================================================
fprintf('\n========== PART 1: Neighboring Structure ==========\n');

for i = 1:n
    out_nb = find(A(i, :));        % i -> j  (row i)
    in_nb  = find(A(:, i))';       % j -> i  (col i)
    fprintf('Node %d | Out-neighbors: [%s] | In-neighbors: [%s]\n', ...
        i, num2str(out_nb), num2str(in_nb));
end

% ----- Distributed BFS / Shortest Distances via repeated matrix multiply -----
% D(i,j) = shortest directed distance from i to j
% Method: power the reachability matrix step by step
% R_k(i,j) = 1 if there is a path of length <= k from i to j

fprintf('\nShortest path distances (Inf = unreachable):\n');
Dist = inf(n, n);
for i = 1:n
    Dist(i, i) = 0;
end
% Initialize with direct edges
for i = 1:n
    for j = 1:n
        if A(i,j) == 1
            Dist(i,j) = 1;
        end
    end
end

% Relax up to n-1 times (Bellman-Ford / Floyd-Warshall style)
for k = 1:n-1
    for i = 1:n
        for j = 1:n
            for m = 1:n
                if Dist(i,m) + Dist(m,j) < Dist(i,j)
                    Dist(i,j) = Dist(i,m) + Dist(m,j);
                end
            end
        end
    end
end

fprintf('Distance matrix (rows=source, cols=dest, Inf=unreachable):\n');
% Print nicely
header = sprintf('%6s', '');
for j = 1:n; header = [header sprintf('%6d', j)]; end
fprintf('%s\n', header);
for i = 1:n
    row_str = sprintf('Node%2d', i);
    for j = 1:n
        if isinf(Dist(i,j))
            row_str = [row_str sprintf('%6s', 'Inf')];
        else
            row_str = [row_str sprintf('%6.0f', Dist(i,j))];
        end
    end
    fprintf('%s\n', row_str);
end

%% ============================================================
%  PART 2: Strongly Connected Components (Kosaraju's Algorithm)
%  Implemented in a distributed / explicit-step manner
% =============================================================
fprintf('\n========== PART 2: Strongly Connected Components ==========\n');

% --- Step 1: DFS on original graph to get finish order ---
visited  = false(1, n);
finish_order = [];  % nodes ordered by finish time

for start = 1:n
    if ~visited(start)
        [visited, finish_order] = dfs_iterative(A, start, visited, finish_order);
    end
end

% --- Step 2: Transpose the graph ---
AT = A';

% --- Step 3: DFS on transposed graph in reverse finish order ---
visited2 = false(1, n);
scc_label = zeros(1, n);
scc_id = 0;

for k = length(finish_order):-1:1
    v = finish_order(k);
    if ~visited2(v)
        scc_id = scc_id + 1;
        [visited2, scc_label] = dfs_label(AT, v, visited2, scc_label, scc_id);
    end
end

% --- Report SCCs ---
fprintf('SCC assignments per node:\n');
for i = 1:n
    fprintf('  Node %d -> SCC %d\n', i, scc_label(i));
end

fprintf('\nSCC groups:\n');
for s = 1:scc_id
    members = find(scc_label == s);
    fprintf('  SCC %d: { %s }\n', s, num2str(members));
end

%% ============================================================
%  Visualization
% =============================================================
figure('Name', 'Digraph - Problem 4');
G = digraph(A);
node_colors = scc_label';  % use SCC id as color index

% Map SCC labels to colors
cmap = lines(scc_id);
node_rgb = cmap(node_colors, :);

h = plot(G, 'Layout', 'force', 'ArrowSize', 12, 'LineWidth', 1.5);
for i = 1:n
    highlight(h, i, 'NodeColor', node_rgb(i,:), 'MarkerSize', 10);
end
title('Digraph with SCC coloring (same color = same SCC)');
labelnode(h, 1:n, arrayfun(@(x) sprintf('%d (SCC%d)', x, scc_label(x)), 1:n, 'UniformOutput', false));

%% ============================================================
%  Helper: Iterative DFS (finish-order)
% =============================================================
function [visited, finish_order] = dfs_iterative(A, start, visited, finish_order)
    n = size(A, 1);
    stack = {start};        % cell stack of [node, iterator_index]
    in_stack = false(1,n);
    iter_ptr = ones(1, n);  % next neighbor to visit
    visited(start) = true;
    dfs_stack = [start];    % explicit stack

    while ~isempty(dfs_stack)
        v = dfs_stack(end);
        neighbors = find(A(v, :));
        advanced = false;
        while iter_ptr(v) <= length(neighbors)
            nb = neighbors(iter_ptr(v));
            iter_ptr(v) = iter_ptr(v) + 1;
            if ~visited(nb)
                visited(nb) = true;
                dfs_stack(end+1) = nb;
                advanced = true;
                break;
            end
        end
        if ~advanced
            finish_order(end+1) = v;
            dfs_stack(end) = [];
        end
    end
end

%% ============================================================
%  Helper: DFS for labeling SCC
% =============================================================
function [visited, scc_label] = dfs_label(A, start, visited, scc_label, id)
    stack = [start];
    visited(start) = true;
    scc_label(start) = id;
    while ~isempty(stack)
        v = stack(end);
        stack(end) = [];
        neighbors = find(A(v, :));
        for nb = neighbors
            if ~visited(nb)
                visited(nb) = true;
                scc_label(nb) = id;
                stack(end+1) = nb;
            end
        end
    end
end