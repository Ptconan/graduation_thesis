function amm = mk_matrix(dm,ir,nm,nd)

%%%% make full-dimensional matrices

ir = ir';
amm = zeros(nd^2,nm);
for i=1:nm
    am(:,:) = zeros(nd+1,nd+1);
    % 
    for j=1:4
        for k=1:4
            am(ir(k,i),ir(j,i)) = dm(k,j,i);
            amm(:,i) = reshape(am(1:nd,1:nd), nd^2, 1);
        end
    end
end
amm = sparse(amm);

