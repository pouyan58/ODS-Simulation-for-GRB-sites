function results = runGrbOdsSpotBeamSimulator(siteSelection,options)
%runGrbOdsSpotBeamSimulator Select a GRB case and calculate ODS angles.
%   RESULTS = runGrbOdsSpotBeamSimulator() displays the public site and
%   GOES-boresight choices. A menu number, SiteId, or unique combined label
%   can be supplied for noninteractive operation.

arguments
    siteSelection = []
    options.CatalogFile (1,1) string {mustBeFile} = fullfile( ...
        fileparts(mfilename("fullpath")),"data","GRB_Receiver_Sites.csv")
    options.OutputRoot (1,1) string = string(fileparts( ...
        mfilename("fullpath")))
    options.GeodeticHeightM (1,1) double = NaN
    options.ZaRadiusKm (1,1) double = NaN
    options.DurationSeconds (1,1) double {mustBePositive} = 7200
    options.TimeStepSec (1,1) double {mustBePositive} = 2
    options.MonteCarloRuns (1,1) double {mustBeInteger,mustBePositive} = 4
    options.WriteDetailedTimeSeries (1,1) logical = true
    options.RegenerateAtmosphericLookup (1,1) logical = false
    options.ItuDigitalMapsFolder (1,1) string = ""
    options.DryRun (1,1) logical = false
end

projectFolder = string(fileparts(mfilename("fullpath")));
if ~isfolder(options.OutputRoot)
    mkdir(options.OutputRoot);
end
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(projectFolder,fullfile(projectFolder,"config"), ...
    fullfile(projectFolder,"src"));

sites = loadGrbReceiverSites(options.CatalogFile);
selectedSite = selectGrbReceiverSite(sites,siteSelection);
cfg = defaultGrbOdsSpotBeamConfig(selectedSite,projectFolder, ...
    options.OutputRoot,options.ZaRadiusKm);
cfg.simulation.duration = seconds(options.DurationSeconds);
cfg.simulation.timeStepSec = options.TimeStepSec;
cfg.simulation.monteCarloRuns = options.MonteCarloRuns;
cfg.output.writeDetailedTimeSeries = options.WriteDetailedTimeSeries;
if isfinite(options.GeodeticHeightM)
    cfg.site.heightM = options.GeodeticHeightM;
    cfg.site.heightSourceStatus = ...
        "User-supplied geodetic antenna-reference height.";
end

if options.DryRun
    validateOdsConfig(cfg);
    results.config = cfg;
    results.siteCatalog = sites;
    results.selectedSite = selectedSite;
    displayConfiguration(cfg);
    clear pathCleanup
    return
end

if options.RegenerateAtmosphericLookup || ...
        ~isfile(cfg.propagation.atmosphericLookupFile)
    assert(strlength(options.ItuDigitalMapsFolder) > 0 && ...
        isfolder(options.ItuDigitalMapsFolder),"ODS:MissingItuMaps", ...
        "Set ItuDigitalMapsFolder to an installed MathWorks ITU map folder.");
    fprintf("Generating site-specific GRB P.618 lookup...\n");
    generateP618AtmosphericLookup(cfg,options.ItuDigitalMapsFolder, ...
        string(cfg.propagation.atmosphericLookupFile));
end

validateOdsConfig(cfg);
results = runOdsAngleStudy(cfg);
writeGrbOdsSummary(results,selectedSite);

fprintf("\nGRB ODS Spot Beam Simulator result for %s | %s\n", ...
    cfg.site.name,cfg.site.boresight.satelliteAssignment);
fprintf("  Primary outer angle from boresight: %.2f deg\n", ...
    results.practical.outerAngleDeg);
fprintf("  Primary inner angle from boresight: %.2f deg\n", ...
    results.practical.innerAngleDeg);
fprintf("  Robust modeled availability: %.6f %%\n", ...
    results.practical.robustAvailabilityPercent);
fprintf("  Geographic ZA radius: %.2f km\n", ...
    results.config.za.protectedCoreRadiusKm);
fprintf("  Results folder: %s\n",cfg.output.folder);
clear pathCleanup
end

function displayConfiguration(cfg)
configuration = table(string(cfg.site.name),cfg.site.latitudeDeg, ...
    cfg.site.longitudeDeg,string(cfg.site.boresight.satelliteAssignment), ...
    cfg.site.boresight.geoLongitudeDeg,cfg.antenna.diameterM, ...
    cfg.receiver.centerFrequencyHz/1e6,cfg.receiver.bandwidthHz/1e6, ...
    cfg.propagation.receiverAntennaHeightAglM,cfg.site.heightM, ...
    cfg.za.protectedCoreRadiusKm,cfg.beams.numberPerSatellite, ...
    cfg.rf.maximumBeamEirpDensityDbwPerMHz, ...
    VariableNames=["Site","LatitudeDeg","LongitudeDeg", ...
    "GoesAssignment","BoresightLongitudeDeg","AntennaDiameterM", ...
    "CenterFrequencyMHz","ReceiverBandwidthMHz","FeedHeightAglM", ...
    "GeodeticHeightM","GeographicZaRadiusKm","BeamsPerSatellite", ...
    "MaximumBeamEirpDensityDbwPerMHz"]);
disp(configuration);
end
