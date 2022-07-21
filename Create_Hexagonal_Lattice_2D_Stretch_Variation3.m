function Lattice = Create_Hexagonal_Lattice_2D_Stretch_Variation3(nHx,nHy,R,alpha,beta,PlotLattice,PlotCenters)
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Hexogonal Lattice that is stretched in the X and Y direction only
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Lattice = Create_Hexagonal_Lattice_2D_Stretch(7,26,1,0.1,10,true,false);
   
    % Test values
%         nHz = 7;
%         nHy = 26;
%         R = 1;
%         epsilon = 0.05;
%         beta = 0.1;

    %Build Inital repetative 1-D block------------------------------------------------------
    hexrow = [ [0;2*R],[0;0] ]; % [[xpts],[ypts]]
    center = [R, 0];
   
    % Create hexagon section of lattice --------------------------------------------------------   
    centers1 = [];
    hexCent = [];
    pts1 = [];
    Points = [];
    
    % Build first Row ---------------------------------------------
    for k = 0:nHx-1
        pts1 = [pts1; [k*3*R+hexrow(:,1), hexrow(:,2) ] ];
        centers1 = [centers1; [k*3*R+center(1,1),center(1,2)] ];
    end
    
    % Now add rows, but shifting every other one  ----------------------------------------------
    for q = 1:(nHy-1)
        yshift = q*sqrt(3)/2*R;
        if rem(q,2) == 0
            xshift = 0;
            Points = [Points; [ pts1(:,1), yshift+pts1(:,2)]   ];
            hexCent = [hexCent; [xshift+centers1(:,1), yshift+centers1(:,2)]  ];
        else
            xshift = 3*R/2;
            Points = [Points;  [  [R/2; 3*R/2 + pts1(1:end-1,1)], [yshift+pts1(1,2); yshift+pts1(1:end-1,2)]  ] ];
            hexCent = [hexCent; [xshift+centers1(1:end-1,1),  yshift+centers1(1:end-1,2)]  ];
        end
    end
    
    %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    % Since I'm edited my code for a slightly different layout (nHy-1 at 26 and nHy-1 at 32) I have to shift everything down R/2 in the y direction  
        Points(:,2) = Points(:,2) - sqrt(3)/2*R;
        hexCent(:,2) = hexCent(:,2) - sqrt(3)/2*R;
     
    % Now find association of points related to hexagon centers (used in random rotation) -----------------------------
        BotHex.Centers  = [];
        BotHex.Points   = [];
        TopHex.Centers  = [];
        TopHex.Points   = [];    
        MainHex.Centers = [];
        MainHex.Points  = [];

        for n = 1:size(hexCent,1)
            d = sqrt( ( hexCent(n,1) - Points(:,1) ).^2 + ( hexCent(n,2) - Points(:,2) ).^2 );
            idx = find(d <= R+0.00001)';

            if length(idx) < 6 && isequal(hexCent(n,2),0)
                BotHex.Centers = [BotHex.Centers; hexCent(n,:)];
                BotHex.Points  = [BotHex.Points;  idx];
            elseif length(idx) < 6 && ~isequal(hexCent(n,2),0)
                TopHex.Centers = [TopHex.Centers; hexCent(n,:)];
                TopHex.Points  = [TopHex.Points;  idx];
            else% (hexCent(n,2) < sqrt(3)*(nHy-1)*R/2) && (hexCent(n,2) > sqrt(3)*R/2)
                MainHex.Centers = [MainHex.Centers; hexCent(n,:)];
                MainHex.Points  = [MainHex.Points;  idx];
            end
        end
        
%         MainHex.Centers = [MainHex.Centers; BotHex.Centers];
%         MainHex.Points = [MainHex.Points; BotHex.Points];
        
    % Now find all unique line segments ------------------------------------------------------
        com = nchoosek(1:size(Points,1),2); % Create all possible line connections without repeats then elimate based on distance
        D = sqrt((Points(com(:,1),1) - Points(com(:,2),1)).^2 + (Points(com(:,1),2) - Points(com(:,2),2)).^2);
        AllSegments = com(D <= R+0.00001,:);
        % Remove flat segments in first and top row  
%             idx1 = find(  (Points(AllSegments(:,1),2) == 0) & (Points(AllSegments(:,2),2) == 0) );
%             idx2 = find(  (Points(AllSegments(:,1),2) == max(Points(:,2)) ) & (Points(AllSegments(:,2),2) == max(Points(:,2)) ) );
%        AllSegments([idx1;idx2],:) = [];
        % Sort out points by locations: bottom row, top row, and main points
            idx3y  = find( Points(:,2) == 0 );
            idx4y  = find( Points(:,2) == max(Points(:,2)) );
            idx3x  = find( Points(:,1) == 0 ); 
                idx3x  = setdiff( setdiff(idx3x,idx3y), idx4y); % remove points it has in common with idx3y and idx4y
            idx4x  = find( Points(:,1) == max(Points(:,1)) );
                idx4x  = setdiff( setdiff(idx4x,idx3y), idx4y); % remove points it has in common with idx3y and idx4y
        BotPoints = idx3y;
        TopPoints = idx4y;
        RightPoints = idx3x;
        LeftPoints  = idx4x;
        MainPoints = (1:size(Points,1))';
        MainPoints( [BotPoints;TopPoints;RightPoints;LeftPoints] ) = []; % Remove points connected to the edge
        % Sort Segments by location: bottom connectors, top connectors, and main segments
            idx5 = find(  (Points(AllSegments(:,1),2) == 0) | (Points(AllSegments(:,2),2) == 0)  );
        BotSegments = AllSegments(idx5,:);
            idx6 = find(  (Points(AllSegments(:,1),2) == max(Points(:,2)) ) | (Points(AllSegments(:,2),2) == max(Points(:,2)) ) );
        TopSegments = AllSegments(idx6,:);
        MainSegments = AllSegments;
        MainSegments([idx5;idx6],:) = [];
        
    % Now find spokes for each MainPoint ----------------------------------------------------------    
        nMP = size(MainPoints,1);
        MainSpokes = cell(nMP,1);
        MainSpokesAllSegmentsReference = cell(nMP,1);
        for m = 1:nMP
            d = sqrt( (   Points(MainPoints(m,1),1) - Points(:,1)    ).^2 + ( Points(MainPoints(m,1),2) - Points(:,2) ).^2 );
            idx = find( (d <= R+0.001) & (d~=0) )';    
            MainSpokes{m,1} = idx;
            
            % Now find the Segments that correspond to the spokes for this
            % node in the AllSegments list
            idxn = [];
            for mm = 1:length(idx)
                [idxn1,~] = find( idx(mm) == AllSegments);
                [idxn2,~] = find( MainPoints(m,1) == AllSegments);
                idxn = [idxn; intersect(idxn1,idxn2)];
            end
            MainSpokesAllSegmentsReference{m,1} = idxn';
        end
        
    % Add Noise to points (excluding top and bottom rows and rotate hexagons slightly
        Points(MainPoints,:) = Points(MainPoints,:) + alpha*randn(size(Points(MainPoints,:)));
        
    % Add random rotation to main hexagons     
        for n = 1:size(MainHex.Centers,1)
            C = MainHex.Centers(n,:); % Hexagon center point
            V = Points(MainHex.Points(n,:),:); % Hexagon points
            thetad = beta*randn(1); % Random degree rotation
            R  = [cosd(thetad) -sind(thetad); sind(thetad) cosd(thetad)];
            V = (R*(V - C)')' + C;
            Points(MainHex.Points(n,:),:) = V;
        end
        
   % Assign to stuctured output variable ---------------------------------------- 
       Lattice.Points = Points;
       Lattice.TopPoints = TopPoints;
       Lattice.BottomPoints = BotPoints;
       Lattice.RightPoints  = RightPoints;
       Lattice.LeftPoints   = LeftPoints;
       Lattice.AllBorderPoints = [TopPoints; BotPoints; RightPoints; LeftPoints];
       Lattice.MainPoints = MainPoints;
       Lattice.MainSpokes = MainSpokes;
       Lattice.MainSpokesAllSegmentsReference = MainSpokesAllSegmentsReference;
       Lattice.TopHex = TopHex;
       Lattice.MainHex = MainHex;
       Lattice.BottomHex = BotHex;
       Lattice.AllSegments = AllSegments;
       Lattice.BottomSegments = BotSegments;
       Lattice.MainSegments = MainSegments;
       Lattice.TopSegments = TopSegments;    
        
        
    % Test plotting ---------------------------------------------------------------  
 
    
    if PlotLattice
        figure(1); clf
        plot(Points(:,1),Points(:,2),'.b');  axis equal tight; 
        
        if PlotCenters
            hold on
            plot(hexCent(:,1),hexCent(:,2),'or')
            hold off
        end
        
         for m = 1:size(AllSegments,1)
           hold on
           plot( [Points(AllSegments(m,1),1); Points(AllSegments(m,2),1)],...
                 [Points(AllSegments(m,1),2); Points(AllSegments(m,2),2)],'-b')
           hold off
        end
        
        for m = 1:size(MainSegments,1)
           hold on
           plot( [Points(MainSegments(m,1),1); Points(MainSegments(m,2),1)],...
                 [Points(MainSegments(m,1),2); Points(MainSegments(m,2),2)],'-b')
           hold off
        end
        for m = 1:size(BotSegments,1)
           hold on
           plot( [Points(BotSegments(m,1),1);Points(BotSegments(m,2),1)], [Points(BotSegments(m,1),2);Points(BotSegments(m,2),2)],'-g')
           hold off
        end
        for m = 1:size(TopSegments,1)
           hold on
           plot( [Points(TopSegments(m,1),1);Points(TopSegments(m,2),1)], [Points(TopSegments(m,1),2);Points(TopSegments(m,2),2)],'-m')
           hold off
        end
    end
