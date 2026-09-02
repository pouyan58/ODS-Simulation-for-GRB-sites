function sites = listGrbOdsSpotBeamCases(catalogFile)
%listGrbOdsSpotBeamCases Display the selectable public cases.

arguments
    catalogFile (1,1) string {mustBeFile} = fullfile( ...
        fileparts(mfilename("fullpath")),"data","GRB_Receiver_Sites.csv")
end

sites = loadGrbReceiverSites(catalogFile);
disp(sites(:,["Index","Name","LatitudeDeg","LongitudeDeg", ...
    "Satellite","AntennaDiameterM","FeedHeightM","DefaultZaRadiusKm"]));
end
