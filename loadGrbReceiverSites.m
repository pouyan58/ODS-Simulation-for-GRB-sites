function sites = loadGrbReceiverSites(catalogFile)
%loadGrbReceiverSites Import and validate the public GRB site catalog.

arguments
    catalogFile (1,1) string {mustBeFile} = fullfile( ...
        fileparts(mfilename("fullpath")),"data","GRB_Receiver_Sites.csv")
end

requiredVariables = ["Site","LatitudeDeg","LongitudeDeg", ...
    "Satellite","FeedHeightM"];
raw = readtable(catalogFile,TextType="string", ...
    VariableNamingRule="preserve");
assert(all(ismember(requiredVariables,string(raw.Properties.VariableNames))), ...
    "ODS:SiteCatalogColumns", ...
    "The site catalog does not contain all required columns.");

name = cleanText(raw.Site);
latitudeDeg = parseNumber(raw.LatitudeDeg,"LatitudeDeg");
longitudeDeg = parseNumber(raw.LongitudeDeg,"LongitudeDeg");
satellite = cleanText(raw.Satellite);
feedHeightM = parseNumber(raw.FeedHeightM,"FeedHeightM");

isEast = strcmpi(satellite,"GOES East");
isWest = strcmpi(satellite,"GOES West");
isBackup = strcmpi(satellite,"GOES Backup");
assert(all(isEast | isWest | isBackup), ...
    "ODS:UnsupportedGoesAssignment", ...
    "Each row must use GOES East, GOES West, or GOES Backup.");
boresightLongitudeDeg = nan(height(raw),1);
boresightLongitudeDeg(isEast) = -75.2;
boresightLongitudeDeg(isWest) = -137.0;
boresightLongitudeDeg(isBackup) = -104.7;
boresightLongitudeSource = strings(height(raw),1);
boresightLongitudeSource(isEast | isWest) = ...
    "NOAA GOES-R fleet: GOES-East 75.2 W; GOES-West 137.0 W";
boresightLongitudeSource(isBackup) = ...
    "NOAA GOES-R fleet: backup orbital slot 104.7 W";

receiverId = strings(height(raw),1);
siteId = strings(height(raw),1);
for rowIndex = 1:height(raw)
    receiverId(rowIndex) = makeIdentifier(name(rowIndex));
    siteId(rowIndex) = receiverId(rowIndex) + "_" + ...
        makeIdentifier(satellite(rowIndex));
end

assert(all(latitudeDeg >= -90 & latitudeDeg <= 90), ...
    "ODS:SiteLatitude","Catalog latitudes must lie within [-90, 90].");
assert(all(longitudeDeg >= -180 & longitudeDeg <= 180), ...
    "ODS:SiteLongitude","Catalog longitudes must lie within [-180, 180].");
assert(all(feedHeightM > 0),"ODS:SiteFeedHeight", ...
    "Catalog feed heights must be positive.");
assert(numel(unique(siteId)) == height(raw),"ODS:DuplicateSiteId", ...
    "Each site and GOES assignment combination must be unique.");

defaultZaRadiusKm = 50*ones(height(raw),1);
defaultZaRadiusKm(ismember(receiverId,["fairmont_wv","wallops_va"])) = 30;
index = (1:height(raw)).';
antennaDiameterM = 4.5*ones(height(raw),1);
sourceCatalog = repmat(string(catalogFile),height(raw),1);
sourceRow = index + 1;
sites = table(index,siteId,receiverId,name,latitudeDeg,longitudeDeg, ...
    satellite,antennaDiameterM,feedHeightM,defaultZaRadiusKm, ...
    boresightLongitudeDeg,boresightLongitudeSource,sourceCatalog,sourceRow, ...
    VariableNames=["Index","SiteId","ReceiverId","Name", ...
    "LatitudeDeg","LongitudeDeg","Satellite","AntennaDiameterM", ...
    "FeedHeightM","DefaultZaRadiusKm","BoresightLongitudeDeg", ...
    "BoresightLongitudeSource","SourceCatalog","SourceRow"]);
end

function values = parseNumber(rawValues,variableName)
values = str2double(cleanText(rawValues));
assert(all(isfinite(values)),"ODS:SiteCatalogNumber", ...
    "Column '%s' contains a missing or nonnumeric value.",variableName);
end

function values = cleanText(rawValues)
values = strip(regexprep(string(rawValues),"\s+"," "));
assert(all(strlength(values) > 0),"ODS:SiteCatalogText", ...
    "The receiver catalog contains a blank required text value.");
end

function identifier = makeIdentifier(label)
identifier = lower(regexprep(label,"[^A-Za-z0-9]+","_"));
identifier = regexprep(identifier,"^_+|_+$","");
end
