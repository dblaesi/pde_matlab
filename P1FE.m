function P1FE
load c4n.dat; load n4e.dat; load Db.dat; load Nb.dat;
for red = 1 : 5
    [c4n,n4e,Db,Nb] = red_refine(c4n,n4e,Db,Nb);
end

[element2edgeSign,edgeNumbering,edges,element2edges,Db2edges,Nb2edges] = edge_index(n4e,Db,Nb);
nEdges = size(edges,1);


[S,rs,A] = sysmat(n4e,Db,Nb,c4n);
Q = koeff(n4e, c4n, edges, element2edges, element2edgeSign);
Q' * A * Q;
x = S \ rs ;
u = x(1 : nEdges);
p = x(nEdges + (1:size(n4e,1)));

show(c4n,n4e,u,p,element2edges,element2edgeSign);

end


function Q = koeff(n4e, c4n, edges, element2edges, element2edgeSign)
idx1 = [2 3 1]; idx2 = [3 1 2];
Q = zeros(size(edges,1),1);
for el = 1:size(n4e,1)
    for k = 1:3
        if( element2edgeSign(el,k) > 0 )
            e = element2edges(el, k);
            z_1 = n4e(el, idx1(k));
            z_2 = n4e(el, idx2(k));
            n = (c4n(z_1,:) - c4n(z_2,:)) * [0,1 ;-1,0 ]';
            n = n / norm(n);
            Q(e) = ([1,0] * n');
        end
    end
end
end



