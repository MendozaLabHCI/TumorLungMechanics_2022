% This code creates a 100 x 100 lattice made of hexagons, and selects tumor hexagons (only for non-randomized lattice)
% 
% Create non-randomized matrix ---------------------------------------------------------------
% Lattice = Create_Hexagonal_Lattice_2D_Stretch_Variation3( 54, 190, 1, 0, 0, false, false); 
%
% Create randomized matrix -------------------------------------------------------------------
alpha = 0.1; % scaling constant for add
beta  = 4.24;
Lattice = Create_Hexagonal_Lattice_2D_Stretch_Variation3( 54, 190, 1, 0.1, 4.24, false, false); % (nHx, nHy, R, alpha, beta, PlotLattice, PlotCenters)
std( CalculateMajorAndMinorAxesOfLattice(Lattice))
MaxVal = max([max(Lattice.Points(:,1)); max(Lattice.Points(:,2))]);


Lattice.Points = Lattice.Points/MaxVal; % Normalize node coordinates
Lattice.MainHex.Centers = Lattice.MainHex.Centers/MaxVal;
Lattice.Points = (10^4)*Lattice.Points;

figure(1); clf
    %plot(Lattice.Points(:,1),Lattice.Points(:,2),'.b');  axis equal tight; 
    
   

     for m = 1:size(Lattice.AllSegments,1)
       hold on
       plot( [Lattice.Points(Lattice.AllSegments(m,1),1); Lattice.Points(Lattice.AllSegments(m,2),1)],...
             [Lattice.Points(Lattice.AllSegments(m,1),2); Lattice.Points(Lattice.AllSegments(m,2),2)],'-b','LineWidth',0.5)
     end
     hold on
     %plot(Lattice.MainHex.Centers(:,1), Lattice.MainHex.Centers(:,2), '.r')
     hold off
     
     MaxVal = max([max(Lattice.Points(:,1)); max(Lattice.Points(:,2))]);
     
     M = ceil(MaxVal);
     
     axis equal
     axis([0,M,0,M])
     nHex = length(Lattice.MainHex.Centers(:,1));
     display(nHex/(100^2))
     
     
     % THIS PART ONLY WORKS FOR NON-RANDOMIZED LATTICE. SO CALCULATE SEGMENTS FIRST THEN GO BACK AND RANDOMIZE
     % Check if width at y = 0.5 is even or odd
        idx = find(Lattice.MainHex.Centers(:,2) > 0.4999 & Lattice.MainHex.Centers(:,2) < 0.5001);
     % Center of tumor 
        n = idx(27);
        HC = Lattice.MainHex.Centers; % Hexagon centers
        TC = Lattice.MainHex.Centers(n,:); % Center of tumor coordinates
        D = sqrt( (TC(1,1)-HC(:,1)).^2 + (TC(1,2)-HC(:,2)).^2 );
        DS = [D,(1:length(D))'];
        DS = sortrows(DS,1);
        Tidx = DS(1:19,2); % Tumor indices are the first closest 19 hexagons to the center
        hold on
%         for k = 1:length(Tidx)
%             hold on
%             plot(HC(Tidx(k),1),HC(Tidx(k),2),'.r','MarkerSize',20)
%         end
        hold off
        
        THP = Lattice.MainHex.Points(Tidx,:);
        % Find which segments are in tumor
        Tsegs = [];
        for m = 1:size(THP,1)
            for p = 1:6
                for q = 1:6
                    
                        seg = [THP(m,p),THP(m,q)];
                        idx = find(  ismember(Lattice.AllSegments,seg,'rows')  );
                        if ~isempty(idx)
                            Tsegs = [Tsegs;idx];
                        end
                    
                end
            end
        end
        Tsegs = unique(Tsegs);
        
        TumorSegments = Lattice.AllSegments(Tsegs,:);
        
            X1 = Lattice.Points(TumorSegments(:,1),1);
            Y1 = Lattice.Points(TumorSegments(:,1),2);
            X2 = Lattice.Points(TumorSegments(:,2),1);
            Y2 = Lattice.Points(TumorSegments(:,2),2);
       hold on
       plot([X1';X2'] , [Y1';Y2'] , 'r-','LineWidth',2)
       hold off
        
        