function [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb)
[edges, element2edges, Db2edges, Nb2edges] = edge_index(n4e, Db, Nb);
%create new Elements via red refine
newNodes = size(c4n,1) + (1:size(edges,1))';
newIndices = newNodes(element2edges);
n4e = [n4e(:,1), newIndices(:,[3 2]), n4e(:,2), newIndices(:,[1 3]), ...
       n4e(:,3), newIndices(:,[2 1]), newIndices];
n4e = reshape(n4e', 3, [])';
%get coordinates for new Nodes
newCoord = .5 * (c4n( edges(:,1), : ) + c4n(edges(:,2), :));
c4n = [c4n; newCoord];
%refine boundary elements

if(size(Db,1) > 0 )
newNb = newNodes(Nb2edges);
Db = [Db(:,1), newDb, newDb, Db(:,2)];
Db = reshape(Db', 2, [])';
end

if( size(Nb,1) > 0 )
   newDb = newNodes(Db2edges);
   Nb = [Nb(:,1), newNb, newNb, Nb(:,2)];
   Nb = reshape(Nb', 2, [])';
end
end

function [edges,element2edges,Db2edges,Nb2edges] = edge_index(n4e,Db,Nb)
%find edge numbering
edges = n4e(:, [2 3, 1 3, 1 2] );% ecken zu kanten in jedem dreieck
edges = reshape(edges', 2, []); % umordnen der informationen, 2 Anzahl der Zeilen, [] so dass es passt
edges = sort(edges, 1);     % Kante soll von kleinerem zu groesserem Knoten gehen
[edges, ~, edgeNumbering] = unique(edges', 'rows', 'first'); % Eliminieren der doppelten Eintraege und Indizierung
nEdges = size(edges, 1); % Anzahl der Zeilen bzw. Kanten
element2edges = reshape(edgeNumbering', 3, [])'; %gibt dann Kanten der Elemente aus (Kanten pro Dreieck...)
%find boundary edge indices
[~, ~, bdy2edges] = unique([edges;sort(Db,2);sort(Nb,2)], 'rows', 'first'); % vergleicht die Kanten mit den Randkanten und speichert Indices
Db2edges = bdy2edges(nEdges + (1:size(Db,1))); %speichert Kantennummer des Dirichletrandes
Nb2edges = bdy2edges(nEdges + size(Db,1) + (1:size(Nb,1))); %speichert Kantennummer des Neumannrandes
end