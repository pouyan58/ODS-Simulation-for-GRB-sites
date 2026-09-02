function cfg = configureGrbOdsZaRadius(cfg, zaRadiusKm)
%configureGrbOdsZaRadius Apply one consistent geographic ZA radius.
%   CFG = configureGrbOdsZaRadius(CFG, ZARADIUSKM) updates every model
%   field that depends on the no-user/active-beam-center radius. Use this
%   function instead of changing individual cfg.za fields separately.

arguments
    cfg (1,1) struct
    zaRadiusKm (1,1) double {mustBeFinite,mustBeNonnegative}
end

requiredGroups = ["meta","site","za","ods","calibration","output"];
for groupName = requiredGroups
    assert(isfield(cfg,groupName),"ODS:ZaMissingConfigGroup", ...
        "Cannot configure ZA because cfg.%s is missing.",groupName);
end
assert(zaRadiusKm < cfg.calibration.fullCorrectionDistanceKm, ...
    "ODS:ZaCalibrationRange", ...
    "ZA radius must be less than the %.6g km full-calibration distance.", ...
    cfg.calibration.fullCorrectionDistanceKm);
assert(zaRadiusKm < cfg.ods.outerMinimumBeamCenterDistanceKm, ...
    "ODS:ZaRetaskRange", ...
    "ZA radius must be less than the %.6g km outer-retask distance.", ...
    cfg.ods.outerMinimumBeamCenterDistanceKm);

zaLabel = string(sprintf("%.6g",zaRadiusKm));
zaToken = replace(zaLabel,".","p");

cfg.za.enabled = true;
cfg.za.protectedCoreRadiusKm = zaRadiusKm;
cfg.za.minimumActiveBeamCenterDistanceKm = zaRadiusKm;
cfg.za.definition = compose( ...
    "No D2D served user/active beam center inside %s km of the protected receiver", ...
    zaLabel);
cfg.za.sourceStatus = compose( ...
    "User-configured %s-km geographic avoidance radius; not a coordination result.", ...
    zaLabel);

cfg.calibration.rampStartDistanceKm = zaRadiusKm;
cfg.calibration.sourceStatus = compose( ...
    "Provisional cross-service engineering calibration; receiver- and " + ...
    "band-specific measurements supersede it. The distance ramp begins " + ...
    "at the configured %s-km ZA boundary.",zaLabel);

if any(string(cfg.site.receiverId) == ["fairmont_wv","wallops_va"])
    defaultZaRadiusKm = 30;
else
    defaultZaRadiusKm = 50;
end
cfg.za.publicDefaultRadiusKm = defaultZaRadiusKm;
cfg.za.publicDefaultMapping = ...
    "30 km at Fairmont and Wallops; 50 km at Sioux Falls and Suitland";
if zaRadiusKm == defaultZaRadiusKm
    cfg.meta.baselineId = "GRB-10-SPOT-60-FOR-30DBW-MHZ-SITE-ZA";
    cfg.output.folder = fullfile(fileparts(cfg.output.folder), ...
        "results_" + string(cfg.site.siteId));
else
    cfg.meta.baselineId = "GRB-10-SPOT-60-FOR-30DBW-MHZ-ZA" + ...
        upper(zaToken);
    cfg.output.folder = fullfile(fileparts(cfg.output.folder), ...
        "results_" + string(cfg.site.siteId) + "_za" + zaToken + "km");
end
end
