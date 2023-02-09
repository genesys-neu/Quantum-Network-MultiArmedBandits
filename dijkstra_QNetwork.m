function [shortestPath, totalCost] = dijkstra_QNetwork(netCostMatrix, s, d)
%==============================================================
% shortestPath: the list of nodes in the shortestPath from source to destination;
% totalCost: the total cost of the  shortestPath;
% farthestNode: the farthest node to reach for each node after performing the routing;
% n: the number of nodes in the network;
% s: source node index;
% d: destination node index;
%==============================================================
%  Code by:
% ++by Xiaodong Wang
% ++23 Jul 2004 (Updated 29 Jul 2004)
% ++http://www.mathworks.com/matlabcentral/fileexchange/5550-dijkstra-shortest-path-routing
% Modifications (simplifications) by Meral Shirazipour 9 Dec 2009
%==============================================================
% Metric modified by Vini Chaudhary (GENESYS, Northeastern University, Boston) to account for quantum physics-informed metric
%==============================================================
n = size(netCostMatrix,1);
for i = 1:n
    % initialize the farthest node to be itself;
    farthestPrevHop(i) = i; % used to compute the RTS/CTS range;
    farthestNextHop(i) = i;
end

% all the nodes are un-visited;
visited(1:n) = false;

distance(1:n) = inf;    % it stores the shortest distance between each node and the source node;    
                        %Initially, set 1/fidelity = 1/0 such that entangled pairs are non reachable to that node.
parent(1:n) = 0;

distance(s) = 1/1;      %fidelity_init = 1 (max possible), distance = 1/fidelity_init = 1/1 (min possible)

for i = 1:(n-1)
    temp = [];
    for h = 1:n
         if ~visited(h)  % in the tree;
             temp=[temp distance(h)];
         else
             temp=[temp inf];
         end
    end
     [t, u] = min(temp);            % it starts from node with the shortest distance to the source;
     visited(u) = true;             % mark it as visited;
     for v = 1:n                    % for each neighbors of node u;
         Metric = 1/CascadedFidelity(1/netCostMatrix(u,v),1/distance(u));
         if ( Metric < distance(v) )
             distance(v) = Metric;  % update the shortest distance when a shorter shortestPath is found;
             parent(v) = u;         % update its parent;
         end;             
     end;
end;

shortestPath = [];
if parent(d) ~= 0   % if there is a shortestPath!
    t = d;
    shortestPath = [d];
    while t ~= s
        p = parent(t);
        shortestPath = [p shortestPath];
        
        if netCostMatrix(t, farthestPrevHop(t)) < netCostMatrix(t, p)
            farthestPrevHop(t) = p;
        end;
        if netCostMatrix(p, farthestNextHop(p)) < netCostMatrix(p, t)
            farthestNextHop(p) = t;
        end;

        t = p;      
    end;
end;

totalCost = distance(d);

%return;
end
