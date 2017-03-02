function show(c4n,n4e,p,u,element2edges,element2edgeSign)
nE = size(n4e,1);
p_T = zeros(nE,2);
mp_T = zeros(nE,2);
shift1 = [2 3 1]; shift2 = [3 1 2];
for j = 1 : nE
    area_T = det([1 1 1;c4n(n4e(j,:),:)'])/2;
    mp_T(j,:) = sum(c4n(n4e(j,:),:))/3;
    for m = 1 : 3
        side_m = element2edges(j,m);
        length_m = norm(c4n(n4e(j,shift1(m)),:)-c4n(n4e(j,shift2(m)),:));
        signum_m = element2edgeSign(j,m);
        p_T(j,:) = p_T(j,:)+p(side_m)*signum_m*length_m*(1/(2*area_T))...
            *(mp_T(j,:)-c4n(n4e(j,m),:));
    end
end
E = reshape(1:3*nE,3,nE)'; n4e_t = n4e';
X = c4n(n4e_t(:),1); Y = c4n(n4e_t(:),2); U = repmat(u,1,3)';
figure(1); clf; subplot(1,2,1); trisurf(E,X,Y,U(:));
subplot(1,2,2); triplot(n4e,c4n(:,1),c4n(:,2),':k'); hold on;
quiver(mp_T(:,1),mp_T(:,2),p_T(:,1),p_T(:,2));


