function AspectRatios = CalculateMajorAndMinorAxesOfLattice(lattice)

% Aspect ratio is the ratio of the two Eigenvectors for the covariance matrix of
% the lattice points
    N = size(lattice.MainHex.Points,1);
    AspectRatios = zeros(N,1);
    
    for n = 1:N
        x = lattice.Points( lattice.MainHex.Points(n,:), 1);
        y = lattice.Points( lattice.MainHex.Points(n,:), 2);
        x = x-mean(x);
        y = y-mean(y);
        
        covx = cov([x,y]);
        [~,latent,~] = pcacov(covx);
        AspectRatios(n,1) = latent(1)/latent(2);
        
%         figure(1); clf
%         plot([0;latent(1)*coeff(1,1)],[0;latent(1)*coeff(2,1)],'r.-' ); 
%         hold on
%         plot([0;latent(2)*coeff(1,2)],[0;latent(2)*coeff(2,2)],'r.-' );
%         axis equal
%         DT = delaunayTriangulation(x,y);
%         C = convexHull(DT);
%         plot(DT.Points(C,1),DT.Points(C,2),'b') 
%         hold off
%         pause
    end
    