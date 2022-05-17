clear, clc

%% Extract map from .csv

[file_name, pathname] = uigetfile({'*.csv';'*.*'},...
   'Get Map from CSV','MultiSelect','off'); 

points = csvread(file_name);

Latitude = points(:,1)';
Longitude = points(:,2)';


%% Generate the map structure
Map.Longitude    = Longitude;
Map.Latitude     = Latitude;
Map.RefLongitude = Longitude(1);
Map.RefLatitude  = Latitude(1);

    % Calculate the local coordinates                       
    nPtos  = length(Map.Longitude);
    xEast  = zeros(1,nPtos);
    yNorth = zeros(1,nPtos);
    zUp    = zeros(1,nPtos);                             
    earth = referenceSphere('Earth'); %Spheroid object
    
    [xEast,yNorth,zUp] = geodetic2enu(Map.Latitude,Map.Longitude,0,Map.RefLatitude,Map.RefLongitude,0,earth);
                                        
Map.x            = xEast;
Map.y            = yNorth;
Map.distWP       = 1;
Map.nEH          = 100;
Map.nPtos        = nPtos;
Map.SpeedLimit   = uint8(40 * ones(1, Map.nPtos));
Map.id           = linspace(1, Map.nPtos, Map.nPtos);
Map.nLanes       = 2 *ones(1, Map.nPtos);
Map.WidthLanes   = 3.5 *ones(1, Map.nPtos);
Map.idprev       = zeros(1, Map.nPtos); 
Map.idpost       = zeros(1, Map.nPtos);

    % Generate points index
    for j = 1:Map.nPtos
      
        if (j == 1)
            Map.idprev(j) = Map.nPtos;
        else
            Map.idprev(j) = j-1;
        end
        
        if (j == Map.nPtos)
            Map.idpost (j) = 1;
        else
            Map.idpost(j) = j+1;
        end
    
    end

%% Save Map and clear variables
varlist = {'earth','j','Latitude','Longitude','nPtos','points','xEast','yNorth','zUp','file_name','pathname'};

clear(varlist{:}); clear varlist;
save Map.mat

clear, clc
disp('Map generation finished');
