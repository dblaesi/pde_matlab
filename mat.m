function [A,B,C,x] = mat(n4e,Db,Nb,c4n)
for red = 1 : 1
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end
[element2edgeSign,edgeNumbering,edges,element2edges,Db2edges,Nb2edges] = edge_index(n4e,Db,Nb);
edg = edges(edgeNumbering,:);
nEdges = size(edges,1);
A = sparse(nEdges,nEdges);
B = sparse(nEdges,size(n4e,1));
C = sparse(nEdges,size(Nb,1));
G = zeros(size(Nb,1));
n_E = zeros(size(Nb,1),2);


% Berechnung von G
for k = 1 : size(Nb,1)
     mid = .5 * (c4n(Nb(k,1),:)+c4n(Nb(k,2),:));
     G(k) = g(mid);
     C(Nb2edges(k),k) = 1;
end



% Berechnung von A
for k = 1 : size(n4e,1)
    area_T = det([1,1,1;c4n(n4e(k,:),:)'])/2;    
    for i = 1 : 3
        for j = 1 : 3
            s=0;
            for l = 1:3
                for r = 1:3
                    if r == l
                        delta = 1;
                    else
                        delta = 0;
                    end
                    s = s + (c4n(n4e(k,l),:) - c4n(n4e(k,i),:))* ...
                        (c4n(n4e(k,r),:) - c4n(n4e(k,j),:))'* ...
                        (1 + delta)/12 ;
                    
                end
            end
            A(element2edges(k,i),element2edges(k,j)) = A(element2edges(k,i),element2edges(k,j)) + ...
                element2edgeSign(k,i) * element2edgeSign(k,j) * 1/(4 * area_T) * ...
                norm(c4n(edg(element2edges(k,i),1),:) - c4n(edg(element2edges(k,i),2),:)) * ...
                norm(c4n(edg(element2edges(k,j),1),:) - c4n(edg(element2edges(k,j),2),:)) * s;
        end
    end
end
% Berechnung von B
for k = 1 : size(n4e,1)
    for i = 1 : 3
        B(element2edges(k,i),k) = 0.5 * element2edgeSign(k,i) * ...
            norm(c4n(edg(element2edges(k,i),1),:) - c4n(edg(element2edges(k,i),2),:));
    end
end

S = [A,B,C];
%geaendert
S = [A,B,C;  B', zeros(size(n4e,1),size(S,2)-nEdges);C',zeros(size(Nb,1),size(S,2)-nEdges)];

rs = zeros(size(S,1),1);
for k = 1 :size(n4e,1)  %L^0(omega) hat soviele freiheitsgrade wie dreiecke
    rs(nEdges + k) = area_T;
end
for l = 1 : size(G)
    rs(l + size(S,1)-size(G))= G(l);
end

% Q = koeff(n4e, c4n, edges, element2edges, element2edgeSign);
% Q' * A * Q
x = S \ rs ;
u = x(1:size(A,1))
p = x(size(A,1) + (1:size(B',1)) )

show(c4n,n4e,u,p,element2edges,element2edgeSign);

end


% function Q = koeff(n4e, c4n, edges, element2edges, element2edgeSign)
% idx1 = [2 3 1]; idx2 = [3 1 2];
% Q = zeros(size(edges,1),1);
% for el = 1:size(n4e,1)
%     for k = 1:3
%         if( element2edgeSign(el,k) > 0 )
%             e = element2edges(el, k);
%             z_1 = n4e(el, idx1(k));
%             z_2 = n4e(el, idx2(k));
%             n = (c4n(z_1,:) - c4n(z_2,:)) * [0,1 ;-1,0 ]';
%             n = n / norm(n);
%             Q(e) = ([1,0] * n');
%         end
%     end
% end
% end

function val = g(x)
    val = -1;
end

