function results = runOdsAngleStudy(cfg, geometry)
%runOdsAngleStudy Calculate minimum robust ODS outer and inner angles.

arguments
    cfg (1,1) struct
    geometry (1,1) struct = struct()
end

validateOdsConfig(cfg);
outputFolder = string(cfg.output.folder);
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

rngState = rng;
cleanup = onCleanup(@() rng(rngState));
rng(cfg.simulation.randomSeed, "twister");

if isempty(fieldnames(geometry))
    fprintf("Building geometry for %s...\n", cfg.meta.caseName);
    geometry = buildConstellationGeometry(cfg);
else
    validateSharedGeometry(cfg, geometry);
    fprintf("Reusing paired protected-site geometry for %s...\n", ...
        cfg.meta.caseName);
end
fprintf("Building per-beam link and traffic realizations...\n");
rng(cfg.simulation.randomSeed, "twister");
realizations = buildOdsRealizations(cfg, geometry);
if isfield(cfg.search, "adaptiveRefinementEnabled") && ...
        cfg.search.adaptiveRefinementEnabled
    fprintf("Running coarse ODS search over %d Monte Carlo runs...\n", ...
        cfg.simulation.monteCarloRuns);
    [coarseTable, coarseBest, coarsePractical, baseline] = searchAvoidanceAngles( ...
        cfg, geometry, realizations);
    coarseTable.SearchStage = repmat("coarse", height(coarseTable), 1);

    refinementCfg = centeredSearchConfig(cfg, ...
        [coarseBest; coarsePractical], ...
        cfg.search.refinementHalfWidthDeg, cfg.search.refinementStepDeg);
    fprintf("Refining strict %.1f/%.1f and practical %.1f/%.1f degrees...\n", ...
        coarseBest.OuterAngleDeg, coarseBest.InnerAngleDeg, ...
        coarsePractical.OuterAngleDeg, coarsePractical.InnerAngleDeg);
    [refinementTable, refinementBest, refinementPractical] = searchAvoidanceAngles( ...
        refinementCfg, geometry, realizations);
    refinementTable.SearchStage = repmat("refinement", ...
        height(refinementTable), 1);

    finalCfg = centeredSearchConfig(cfg, ...
        [refinementBest; refinementPractical], ...
        cfg.search.finalHalfWidthDeg, cfg.search.finalStepDeg);
    fprintf("Final half-degree search around strict and practical candidates...\n");
    [finalTable, bestRow, practicalBestRow] = searchAvoidanceAngles( ...
        finalCfg, geometry, realizations);
    finalTable.SearchStage = repmat("final", height(finalTable), 1);
    searchTable = [coarseTable; refinementTable; finalTable];
else
    fprintf("Searching %d angle pairs over %d Monte Carlo runs...\n", ...
        countAnglePairs(cfg, geometry.summary.visibleSkyMaximumSeparationDeg), ...
        cfg.simulation.monteCarloRuns);
    [searchTable, bestRow, practicalBestRow, baseline] = searchAvoidanceAngles( ...
        cfg, geometry, realizations);
end
results = assembleResults(cfg, geometry, realizations, searchTable, ...
    bestRow, practicalBestRow, baseline);
writeOdsResults(results, outputFolder);
plotOdsResults(results, outputFolder);

fprintf("Strict diagnostic outer/inner: %.2f/%.2f deg\n", ...
    results.selected.outerAngleDeg,results.selected.innerAngleDeg);
fprintf("Strict diagnostic worst-case exceedance: %.3f %% (limit %.3f %%)\n", ...
    results.selected.worstExceedancePercent, ...
    cfg.reporting.strictAllowedExceedancePercent);
fprintf("Primary 99.9%% outer/inner: %.2f/%.2f deg; robust availability %.3f %%\n", ...
    results.practical.outerAngleDeg, results.practical.innerAngleDeg, ...
    results.practical.robustAvailabilityPercent);
fprintf("Protected-carrier shutdown duty: %.3f %%; physical retask duty: %.3f %%\n", ...
    results.selected.shutdownDutyPercent, results.selected.retaskDutyPercent);
fprintf("%s snapshot steering-lattice hit diagnostic: %.3f %%; mean effective RFI contributors: %.3f\n", ...
    cfg.site.name, results.geometry.siteCellCoveragePercent, ...
    mean(results.selected.effectiveContributorCount));
clear cleanup
end

function validateSharedGeometry(cfg, geometry)
expectedTimeCount = round(seconds(cfg.simulation.duration) / ...
    cfg.simulation.timeStepSec) + 1;
expectedSatelliteCount = sum([cfg.constellation.shells.satellites]);
tolerance = 1e-10;
assert(numel(geometry.time) == expectedTimeCount, ...
    "ODS:IncompatibleSharedGeometry", ...
    "Shared geometry time grid is incompatible with this case.");
assert(size(geometry.angleDeg, 2) == expectedSatelliteCount, ...
    "ODS:IncompatibleSharedGeometry", ...
    "Shared geometry satellite count is incompatible with this case.");
assert(abs(geometry.summary.beamConeHalfAngleDeg - ...
    cfg.beams.coneHalfAngleDeg) <= tolerance, ...
    "ODS:IncompatibleSharedGeometry", ...
    "Shared geometry beam boundary is incompatible with this case.");
assert(geometry.summary.beamsPerSatellite == cfg.beams.numberPerSatellite, ...
    "ODS:IncompatibleSharedGeometry", ...
    "Shared geometry beam count is incompatible with this case.");
assert(abs(geometry.summary.siteLatitudeDeg - cfg.site.latitudeDeg) <= ...
    tolerance && abs(geometry.summary.siteLongitudeDeg - ...
    cfg.site.longitudeDeg) <= tolerance, ...
    "ODS:IncompatibleSharedGeometry", ...
    "Shared geometry protected-site coordinates are incompatible.");
end

function stageCfg = centeredSearchConfig(cfg, bestRows, halfWidthDeg, stepDeg)
stageCfg = cfg;
outerGrid = [];
innerGrid = [];
for rowIdx = 1:height(bestRows)
    outerGrid = [outerGrid centeredGrid(bestRows.OuterAngleDeg(rowIdx), ...
        halfWidthDeg, stepDeg)]; %#ok<AGROW>
    innerGrid = [innerGrid centeredGrid(bestRows.InnerAngleDeg(rowIdx), ...
        halfWidthDeg, stepDeg)]; %#ok<AGROW>
end
stageCfg.search.outerAngleRangeDeg = unique(outerGrid);
stageCfg.search.innerAngleRangeDeg = unique(innerGrid);
stageCfg.search.adaptiveRefinementEnabled = false;
end

function valuesDeg = centeredGrid(centerDeg, halfWidthDeg, stepDeg)
lowerDeg = max(0, centerDeg - halfWidthDeg);
upperDeg = centerDeg + halfWidthDeg;
valuesDeg = unique([lowerDeg:stepDeg:upperDeg centerDeg]);
end

function count = countAnglePairs(cfg, visibleSkyLimitDeg)
count = 0;
outerCandidates = cfg.search.outerAngleRangeDeg( ...
    cfg.search.outerAngleRangeDeg <= visibleSkyLimitDeg);
innerCandidates = cfg.search.innerAngleRangeDeg( ...
    cfg.search.innerAngleRangeDeg <= visibleSkyLimitDeg);
outerCandidates = unique([outerCandidates visibleSkyLimitDeg]);
innerCandidates = unique([innerCandidates visibleSkyLimitDeg]);
for outerAngle = outerCandidates
    count = count + nnz(innerCandidates <= outerAngle);
end
end

function results = assembleResults(cfg, geometry, realizations, searchTable, ...
    bestRow, practicalBestRow, baseline)
results.config = cfg;
results.geometry = geometry.summary;
results.traffic = realizations.summary;
results.searchTable = searchTable;
results.baseline = baseline;
results.selected.outerAngleDeg = bestRow.OuterAngleDeg;
results.selected.innerAngleDeg = bestRow.InnerAngleDeg;
results.selected.worstExceedancePercent = bestRow.WorstExceedancePercent;
results.selected.robustExceedancePercent = bestRow.RobustExceedancePercent;
results.selected.medianExceedancePercent = bestRow.MedianExceedancePercent;
results.selected.maximumIOverNDb = bestRow.MaximumIOverNDb;
results.selected.p95IOverNDb = bestRow.P95IOverNDb;
results.selected.practicalPercentileIOverNDb = ...
    bestRow.PracticalPercentileIOverNDb;
results.selected.robustAvailabilityPercent = ...
    bestRow.RobustAvailabilityPercent;
results.selected.retaskDutyPercent = bestRow.RetaskDutyPercent;
results.selected.shutdownDutyPercent = bestRow.ShutdownDutyPercent;
results.selected.taperDutyPercent = bestRow.RetaskDutyPercent;
results.selected.deepNullDutyPercent = bestRow.ShutdownDutyPercent;
results.selected.feasible = bestRow.Feasible;
results.selected.score = bestRow.Score;
results.selected.time = geometry.time;
results.selected.aggregateInterferenceDbm = bestRow.AggregateInterferenceDbm{1};
results.selected.iOverNDb = bestRow.IOverNDb{1};
results.selected.dominantSatelliteIndex = bestRow.DominantSatelliteIndex{1};
results.selected.visibleSatelliteCount = geometry.visibleCount;
results.selected.candidateSatelliteCount = geometry.contributingCount;
if string(cfg.aggregation.powerCombinationMode) == "dominantOnly"
    results.selected.effectiveContributorCount = double( ...
        results.selected.dominantSatelliteIndex > 0);
else
    results.selected.effectiveContributorCount = sum(double( ...
        realizations.activeScheduledLink(:, :, 1)), 2);
end
results.selected.eligibleContributorCount = geometry.eligibleContributorCount;

results.practical.outerAngleDeg = practicalBestRow.OuterAngleDeg;
results.practical.innerAngleDeg = practicalBestRow.InnerAngleDeg;
results.practical.worstExceedancePercent = ...
    practicalBestRow.WorstExceedancePercent;
results.practical.robustExceedancePercent = ...
    practicalBestRow.RobustExceedancePercent;
results.practical.maximumIOverNDb = practicalBestRow.MaximumIOverNDb;
results.practical.practicalPercentileIOverNDb = ...
    practicalBestRow.PracticalPercentileIOverNDb;
results.practical.robustAvailabilityPercent = ...
    practicalBestRow.RobustAvailabilityPercent;
results.practical.retaskDutyPercent = practicalBestRow.RetaskDutyPercent;
results.practical.shutdownDutyPercent = practicalBestRow.ShutdownDutyPercent;
results.practical.feasible = practicalBestRow.PracticalFeasible;
results.practical.score = practicalBestRow.Score;
results.practical.aggregateInterferenceDbm = ...
    practicalBestRow.AggregateInterferenceDbm{1};
results.practical.iOverNDb = practicalBestRow.IOverNDb{1};
end
