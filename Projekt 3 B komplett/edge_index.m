function [element2edgeSign,edgeNumbering,edges,element2edges,Db2edges,Nb2edges] = edge_index(n4e,Db,Nb)
%find edge numbering
edges = n4e(:, [2 3, 1 3, 1 2] );% ecken zu kanten in jedem dreieck
edges = reshape(edges', 2, []); % umordnen der informationen, 2 Anzahl der Zeilen, [] so dass es passt
edges = sort(edges, 1);     % Kante soll von kleinerem zu groesserem Knoten gehen
[edges, idx_first, edgeNumbering] = unique(edges', 'rows', 'first'); % Eliminieren der doppelten Eintraege und Indizierung
nEdges = size(edges, 1); % Anzahl der Zeilen bzw. Kanten
element2edges = reshape(edgeNumbering', 3, [])'; %gibt dann Kanten der Elemente aus (Kanten pro Dreieck...)
%find boundary edge indices
[~, ~, bdy2edges] = unique([edges;sort(Db,2);sort(Nb,2)], 'rows', 'first'); % vergleicht die Kanten mit den Randkanten und speichert Indices
Db2edges = bdy2edges(nEdges + (1:size(Db,1))); %speichert Kantennummern(Zeilennummer in edges) des Dirichletrandes
Nb2edges = bdy2edges(nEdges + size(Db,1) + (1:size(Nb,1))); %speichert Kantennummern des Neumannrandes
%
element2edgeSign = -1*ones(3*size(n4e,1),1);
element2edgeSign(idx_first) = 1;
element2edgeSign = reshape(element2edgeSign', 3, [])';


% 
% element2edgeSign = zeros( size(element2edges,1),3);
% for e = 1:size(edges,1)
%     l = find(element2edges == e);
%     c = size(element2edges,1);
%     j= rem(l,c);
%     if j(1) < 1
%         j(1)= j(1)+c;
%     end
%     if size(l,1) < 2
%         k= find(element2edges(j(1),:) == e);
%         element2edgeSign(j(1),k) = 1;
%     else
%         k= find(element2edges(j(1),:) == e);
%         if j(2) < 1
%             j(2)= j(2)+c;
%         end
%         r= find(element2edges(j(2),:) == e);
%         element2edgeSign(j(1),k) = 1;
%         element2edgeSign(j(2),r) = -1;
%     end
%     
% end

end
