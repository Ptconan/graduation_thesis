function dm = elem_matrix(xx,irr,dll,rho,nm)

%%%% make coordinate-transformation matrix : mat_T 
matT = zeros(4,4,nm);
for i=1:nm
    dx = xx(irr(i,2),:) - xx(irr(i,1),:);
    %
    T0 = [dx(1)/dll(i), dx(2)/dll(i);...
        -dx(2)/dll(i), dx(1)/dll(i)];
    % 
    matT(:,:,i) = [T0, zeros(2,2);...
            zeros(2,2), T0];
end

%%% make member mass matrix : dm(:,:,i)
dm = zeros(4,4,nm);
for i=1:nm
    mat_elm = zeros(4,4);
    for j=1:4
        mat_elm(j,j) = 1/3;
    end
    for j=1:2
        mat_elm(j,j+2) = 1/6;
        mat_elm(j+2,j) = 1/6;
    end
    mat_elm = rho * dll(i) * mat_elm;
    % 
    dm(:,:,i) = matT(:,:,i) * mat_elm * (matT(:,:,i)');
end

