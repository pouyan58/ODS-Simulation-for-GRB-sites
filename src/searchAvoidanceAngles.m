function [searchTable, bestRow, practicalBestRow, baseline] = ...
    searchAvoidanceAngles(cfg, geometry, realizations)
%searchAvoidanceAngles Evaluate physical ODS candidate angle pairs.

noiseDbm = receiverNoiseDbm(cfg);
criterionDbm = noiseDbm + cfg.receiver.protectionCriterionDb;
visibleSkyLimitDeg = geometry.summary.visibleSkyMaximumSeparationDeg;
outerCandidatesDeg = cfg.search.outerAngleRangeDeg( ...
    cfg.search.outerAngleRangeDeg <= visibleSkyLimitDeg);
if isfield(cfg.ods, "minimumOperationalOuterAngleDeg")
    outerCandidatesDeg = outerCandidatesDeg(outerCandidatesDeg >= ...
        cfg.ods.minimumOperationalOuterAngleDeg - 1e-10);
end
innerCandidatesDeg = cfg.search.innerAngleRangeDeg( ...
    cfg.search.innerAngleRangeDeg <= visibleSkyLimitDeg);
outerCandidatesDeg = unique([outerCandidatesDeg visibleSkyLimitDeg]);
innerCandidatesDeg = unique([innerCandidatesDeg visibleSkyLimitDeg]);

candidateContext = buildCandidateContext(cfg, geometry, realizations, ...
    outerCandidatesDeg, innerCandidatesDeg, noiseDbm);
baselineMetrics = evaluatePairCached(cfg, candidateContext, 0, 0);
baseline.noiseDbm = noiseDbm;
baseline.criterionDbm = criterionDbm;
baseline.worstExceedancePercent = baselineMetrics.worstExceedancePercent;
baseline.robustExceedancePercent = baselineMetrics.robustExceedancePercent;
baseline.maximumIOverNDb = baselineMetrics.maximumIOverNDb;
baseline.practicalPercentileIOverNDb = ...
    baselineMetrics.practicalPercentileIOverNDb;

rows = cell(countPairs(outerCandidatesDeg, innerCandidatesDeg), 14);
rowIdx = 0;
for outerAngleDeg = outerCandidatesDeg
    for innerAngleDeg = innerCandidatesDeg
        if innerAngleDeg > outerAngleDeg
            continue
        end
        rowIdx = rowIdx + 1;
        metrics = evaluatePairCached(cfg, candidateContext, outerAngleDeg, ...
            innerAngleDeg);
        strictFeasible = metrics.worstExceedancePercent <= ...
            cfg.reporting.strictAllowedExceedancePercent;
        practicalFeasible = metrics.robustExceedancePercent <= ...
            cfg.reporting.practicalAllowedExceedancePercent && ...
            metrics.practicalPercentileIOverNDb <= ...
            cfg.receiver.protectionCriterionDb;
        score = outerAngleDeg + innerAngleDeg + ...
            cfg.search.serviceImpactWeight * ( ...
            metrics.retaskDutyPercent + 5 * metrics.shutdownDutyPercent);
        rows(rowIdx, :) = {outerAngleDeg, innerAngleDeg, ...
            metrics.worstExceedancePercent, metrics.robustExceedancePercent, ...
            metrics.medianExceedancePercent, metrics.maximumIOverNDb, ...
            metrics.p95IOverNDb, metrics.practicalPercentileIOverNDb, ...
            metrics.robustAvailabilityPercent, metrics.retaskDutyPercent, ...
            metrics.shutdownDutyPercent, strictFeasible, practicalFeasible, score};
    end
end

searchTable = cell2table(rows, VariableNames=["OuterAngleDeg", "InnerAngleDeg", ...
    "WorstExceedancePercent", "RobustExceedancePercent", ...
    "MedianExceedancePercent", "MaximumIOverNDb", "P95IOverNDb", ...
    "PracticalPercentileIOverNDb", "RobustAvailabilityPercent", ...
    "RetaskDutyPercent", "ShutdownDutyPercent", "Feasible", ...
    "PracticalFeasible", "Score"]);
bestSummary = selectCandidate(searchTable, "Feasible");
practicalSummary = selectCandidate(searchTable, "PracticalFeasible");
bestMetrics = evaluatePairCached(cfg, candidateContext, ...
    bestSummary.OuterAngleDeg, bestSummary.InnerAngleDeg);
practicalMetrics = evaluatePairCached(cfg, candidateContext, ...
    practicalSummary.OuterAngleDeg, practicalSummary.InnerAngleDeg);
bestRow = attachTimeSeries(bestSummary, bestMetrics);
practicalBestRow = attachTimeSeries(practicalSummary, practicalMetrics);
end

function selected = selectCandidate(searchTable, feasibilityVariable)
feasible = searchTable.(feasibilityVariable);
if any(feasible)
    ranked = sortrows(searchTable(feasible, :), ...
        ["Score", "OuterAngleDeg", "InnerAngleDeg"]);
else
    if feasibilityVariable == "Feasible"
        warning("ODS:NoStrictFeasibleAngles", ...
            "No candidate met strict zero exceedance; returning the least-exceeding pair.");
    else
        warning("ODS:NoPracticalFeasibleAngles", ...
            "No candidate met the practical availability criterion; returning the least-exceeding pair.");
    end
    ranked = sortrows(searchTable, ["WorstExceedancePercent", ...
        "RobustExceedancePercent", "MaximumIOverNDb", "Score", ...
        "OuterAngleDeg", "InnerAngleDeg"]);
end
selected = ranked(1, :);
end

function row = attachTimeSeries(row, metrics)
row.AggregateInterferenceDbm = {metrics.aggregateInterferenceDbm};
row.IOverNDb = {metrics.iOverNDb};
row.DominantSatelliteIndex = {metrics.dominantSatelliteIndex};
end

function count = countPairs(outerCandidatesDeg, innerCandidatesDeg)
count = 0;
for outerAngleDeg = outerCandidatesDeg
    count = count + nnz(innerCandidatesDeg <= outerAngleDeg);
end
end

function noiseDbm = receiverNoiseDbm(cfg)
boltzmann = physconst("Boltzmann");
noiseW = boltzmann * cfg.receiver.noiseTemperatureK * cfg.receiver.bandwidthHz;
noiseDbm = 10 * log10(noiseW) + 30;
end

function context = buildCandidateContext(cfg, geometry, realizations, ...
    outerCandidatesDeg, innerCandidatesDeg, noiseDbm)
% Cache normal power removed by outer action and physically retasked power
% inserted for the same links. Inner shutdown removes the retasked subset.

assert(string(cfg.aggregation.powerCombinationMode) == "linearSum", ...
    "ODS:CachedSearchLinearSum", ...
    "The operational cached search requires linear-sum aggregation.");
timeCount = numel(geometry.time);
runCount = cfg.simulation.monteCarloRuns;
outerCandidatesDeg = unique([0 outerCandidatesDeg]);
innerCandidatesDeg = unique([0 innerCandidatesDeg]);
outerCount = numel(outerCandidatesDeg);
innerCount = numel(innerCandidatesDeg);
context.outerCandidatesDeg = outerCandidatesDeg;
context.innerCandidatesDeg = innerCandidatesDeg;
context.noiseDbm = noiseDbm;
context.totalNormalPowerMw = zeros(timeCount, runCount);
context.outerNormalPowerMw = zeros(timeCount, outerCount, runCount);
context.outerRetaskedPowerMw = zeros(timeCount, outerCount, runCount);
context.innerRetaskedPowerMw = zeros(timeCount, innerCount, runCount);
context.outerActiveLinkCount = zeros(outerCount, runCount);
context.innerActiveLinkCount = zeros(innerCount, runCount);
context.activeLinkDenominator = zeros(runCount, 1);

for runIdx = 1:runCount
    normalPowerMw = 10.^(realizations.preMitigationPowerDbm(:, :, runIdx) / 10);
    retaskedPowerMw = 10.^(realizations.outerRetaskedPowerDbm(:, :, runIdx) / 10);
    activeScheduled = realizations.activeScheduledLink(:, :, runIdx);
    context.totalNormalPowerMw(:, runIdx) = sum(normalPowerMw, 2);
    context.activeLinkDenominator(runIdx) = max(1, nnz(activeScheduled));
    for candidateIdx = 1:outerCount
        protected = hystereticConeMask( ...
            realizations.effectiveAngleDeg(:, :, runIdx), ...
            realizations.predictedEntryAngleDeg(:, :, runIdx), ...
            realizations.protectionRelevant, ...
            outerCandidatesDeg(candidateIdx), cfg.ods.hysteresisDeg);
        context.outerNormalPowerMw(:, candidateIdx, runIdx) = ...
            sum(normalPowerMw .* protected, 2);
        context.outerRetaskedPowerMw(:, candidateIdx, runIdx) = ...
            sum(retaskedPowerMw .* protected, 2);
        context.outerActiveLinkCount(candidateIdx, runIdx) = ...
            nnz(protected & activeScheduled);
    end
    for candidateIdx = 1:innerCount
        protected = hystereticConeMask( ...
            realizations.effectiveAngleDeg(:, :, runIdx), ...
            realizations.predictedEntryAngleDeg(:, :, runIdx), ...
            realizations.protectionRelevant, ...
            innerCandidatesDeg(candidateIdx), cfg.ods.hysteresisDeg);
        context.innerRetaskedPowerMw(:, candidateIdx, runIdx) = ...
            sum(retaskedPowerMw .* protected, 2);
        context.innerActiveLinkCount(candidateIdx, runIdx) = ...
            nnz(protected & activeScheduled);
    end
end
end

function metrics = evaluatePairCached(cfg, context, outerAngleDeg, innerAngleDeg)
runCount = cfg.simulation.monteCarloRuns;
timeCount = size(context.totalNormalPowerMw, 1);
exceedancePercent = zeros(runCount, 1);
maximumIOverNDb = zeros(runCount, 1);
p95IOverNDb = zeros(runCount, 1);
practicalPercentileIOverNDb = zeros(runCount, 1);
retaskDutyPercent = zeros(runCount, 1);
shutdownDutyPercent = zeros(runCount, 1);
representativeAggregateDbm = [];
representativeIOverNDb = [];
representativeDominantSatelliteIndex = [];
outerIdx = find(abs(context.outerCandidatesDeg - outerAngleDeg) < 1e-10, 1);
innerIdx = find(abs(context.innerCandidatesDeg - innerAngleDeg) < 1e-10, 1);
assert(~isempty(outerIdx) && ~isempty(innerIdx), ...
    "ODS:CachedCandidateMissing", ...
    "Requested angle pair is absent from the cached candidate grids.");

for runIdx = 1:runCount
    outerNormalMw = context.outerNormalPowerMw(:, outerIdx, runIdx);
    outerRetaskedMw = context.outerRetaskedPowerMw(:, outerIdx, runIdx);
    innerRetaskedMw = context.innerRetaskedPowerMw(:, innerIdx, runIdx);
    unprotectedPowerMw = max(0, ...
        context.totalNormalPowerMw(:, runIdx) - outerNormalMw);
    outerOnlyPowerMw = max(0, outerRetaskedMw - innerRetaskedMw);
    % The inner protected-carrier chain is shut down exactly; it contributes
    % zero rather than a finite assumed null floor.
    aggregateMw = unprotectedPowerMw + outerOnlyPowerMw;
    aggregateDbm = 10 * log10(max(aggregateMw, realmin("double")));
    dominantSatelliteIndex = zeros(timeCount, 1);
    iOverNDb = aggregateDbm - context.noiseDbm;
    exceedancePercent(runIdx) = 100 * mean(iOverNDb(:) > ...
        cfg.receiver.protectionCriterionDb);
    maximumIOverNDb(runIdx) = max(iOverNDb(:));
    p95IOverNDb(runIdx) = prctile(iOverNDb(:), 95);
    practicalPercentileIOverNDb(runIdx) = prctile(iOverNDb(:), ...
        cfg.reporting.practicalIOverNPercentile);
    innerActiveCount = context.innerActiveLinkCount(innerIdx, runIdx);
    outerActiveCount = context.outerActiveLinkCount(outerIdx, runIdx);
    retaskDutyPercent(runIdx) = 100 * max(0, ...
        outerActiveCount - innerActiveCount) / ...
        context.activeLinkDenominator(runIdx);
    shutdownDutyPercent(runIdx) = 100 * innerActiveCount / ...
        context.activeLinkDenominator(runIdx);
    if runIdx == 1
        representativeAggregateDbm = aggregateDbm;
        representativeIOverNDb = iOverNDb;
        representativeDominantSatelliteIndex = dominantSatelliteIndex;
    end
end

metrics.worstExceedancePercent = max(exceedancePercent);
metrics.robustExceedancePercent = prctile(exceedancePercent, ...
    cfg.uncertainty.robustPercentile);
metrics.medianExceedancePercent = median(exceedancePercent);
metrics.maximumIOverNDb = max(maximumIOverNDb);
metrics.p95IOverNDb = prctile(p95IOverNDb, ...
    cfg.uncertainty.robustPercentile);
metrics.practicalPercentileIOverNDb = prctile( ...
    practicalPercentileIOverNDb, cfg.uncertainty.robustPercentile);
metrics.robustAvailabilityPercent = 100 - ...
    metrics.robustExceedancePercent;
metrics.retaskDutyPercent = prctile(retaskDutyPercent, ...
    cfg.uncertainty.robustPercentile);
metrics.shutdownDutyPercent = prctile(shutdownDutyPercent, ...
    cfg.uncertainty.robustPercentile);
metrics.aggregateInterferenceDbm = representativeAggregateDbm;
metrics.iOverNDb = representativeIOverNDb;
metrics.dominantSatelliteIndex = representativeDominantSatelliteIndex;
assert(numel(representativeIOverNDb) == timeCount, "ODS:TimeSeriesSize", ...
    "Representative time-series size is inconsistent.");
end

function protected = hystereticConeMask(angleDeg, entryAngleDeg, ...
    protectionRelevant, thresholdDeg, hysteresisDeg)
% Explicit pre-entry prediction controls entry; measured angle controls exit.

protected = false(size(angleDeg));
priorProtected = false(1, size(angleDeg, 2));
for timeIdx = 1:size(angleDeg, 1)
    currentProtected = protectionRelevant(timeIdx, :) & ( ...
        entryAngleDeg(timeIdx, :) <= thresholdDeg | ...
        (priorProtected & angleDeg(timeIdx, :) <= ...
        thresholdDeg + hysteresisDeg));
    protected(timeIdx, :) = currentProtected;
    priorProtected = currentProtected;
end
end
