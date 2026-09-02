function realizations = buildOdsRealizations(cfg, geometry)
%buildOdsRealizations Build operational traffic and link-power realizations.
%   The protected 0.4-MHz receiver channel lies in one 5-MHz carrier. A
%   burst scheduler activates at most one ZA-permitted cochannel beam on
%   each satellite. Normal and physically retasked receive powers are built
%   from the same schedule, load, and uncertainty realization.

timeCount = numel(geometry.time);
satelliteCount = size(geometry.angleDeg, 2);
beamCount = cfg.beams.numberPerSatellite;
runCount = cfg.simulation.monteCarloRuns;

realizations.angleErrorDeg = cfg.uncertainty.ephemerisAngleSigmaDeg * ...
    randn(1, satelliteCount, runCount, "single");
realizations.pointingErrorDeg = cfg.uncertainty.pointingSigmaDeg * ...
    randn(timeCount, 1, runCount, "single");
realizations.effectiveAngleDeg = zeros(timeCount, satelliteCount, runCount);
realizations.predictedEntryAngleDeg = zeros(timeCount, satelliteCount, runCount);
realizations.preMitigationPowerDbm = -Inf(timeCount, satelliteCount, runCount);
realizations.outerRetaskedPowerDbm = -Inf(timeCount, satelliteCount, runCount);
realizations.activeScheduledLink = false(timeCount, satelliteCount, runCount);
realizations.dominantPreMitigationPowerDbm = -Inf(timeCount, runCount);
realizations.dominantSatelliteIndex = zeros(timeCount, runCount);
realizations.scheduleDutyPercent = zeros(runCount, 1);
realizations.meanScheduledLoad = zeros(runCount, 1);

[predictionTimeSec, predictedVisible, protectionRelevant] = ...
    preparePredictionGeometry(cfg, geometry);
realizations.protectionRelevant = protectionRelevant;
for runIdx = 1:runCount
    effectiveAngleDeg = min( ...
        geometry.summary.visibleSkyMaximumSeparationDeg, max(0, ...
        double(geometry.angleDeg) + ...
        double(realizations.angleErrorDeg(:, :, runIdx)) + ...
        double(realizations.pointingErrorDeg(:, :, runIdx)) + ...
        cfg.antenna.pointingBiasDeg));
    predictedAngleDeg = interp1((0:(timeCount - 1)).' * ...
        cfg.simulation.timeStepSec, effectiveAngleDeg, predictionTimeSec, ...
        "linear", "extrap");
    entryAngleDeg = effectiveAngleDeg;
    entryAngleDeg(predictedVisible) = min( ...
        effectiveAngleDeg(predictedVisible), ...
        predictedAngleDeg(predictedVisible));

    schedule = generateOperationalSchedule(cfg, geometry);
    perBeamPowerOffsetDb = cfg.uncertainty.powerSigmaDb * ...
        randn(1, satelliteCount, beamCount, "single");
    [preMitigationPowerDbm, outerRetaskedPowerDbm, activeScheduledLink] = ...
        calculateReceivePower(cfg, geometry, effectiveAngleDeg, ...
        double(perBeamPowerOffsetDb), schedule);
    [dominantPowerDbm, dominantSatelliteIndex] = max( ...
        preMitigationPowerDbm, [], 2);
    noEligibleInterferer = ~isfinite(dominantPowerDbm);
    dominantSatelliteIndex(noEligibleInterferer) = 0;

    realizations.effectiveAngleDeg(:, :, runIdx) = effectiveAngleDeg;
    realizations.predictedEntryAngleDeg(:, :, runIdx) = entryAngleDeg;
    realizations.preMitigationPowerDbm(:, :, runIdx) = preMitigationPowerDbm;
    realizations.outerRetaskedPowerDbm(:, :, runIdx) = outerRetaskedPowerDbm;
    realizations.activeScheduledLink(:, :, runIdx) = activeScheduledLink;
    realizations.dominantPreMitigationPowerDbm(:, runIdx) = dominantPowerDbm;
    realizations.dominantSatelliteIndex(:, runIdx) = dominantSatelliteIndex;
    visibleSchedule = schedule.activeBeamMask & geometry.contributing;
    realizations.scheduleDutyPercent(runIdx) = 100 * ...
        nnz(visibleSchedule) / max(1, nnz(geometry.contributing) * ...
        cfg.resources.maximumSimultaneousCochannelBeamsPerSatellite);
    if any(visibleSchedule, "all")
        realizations.meanScheduledLoad(runIdx) = ...
            mean(double(schedule.loadByBeam(visibleSchedule)));
    end
end
realizations.summary.robustScheduleDutyPercent = prctile( ...
    realizations.scheduleDutyPercent, cfg.uncertainty.robustPercentile);
realizations.summary.meanScheduledLoad = mean(realizations.meanScheduledLoad);
realizations.summary.maximumCochannelBeamsPerSatellite = ...
    cfg.resources.maximumSimultaneousCochannelBeamsPerSatellite;
realizations.summary.calibrationCorrectionDb = ...
    cfg.calibration.sidelobeCorrectionDb;
end

function [receivePowerDbm, retaskedPowerDbm, activeScheduledLink] = ...
    calculateReceivePower(cfg, geometry, angleDeg, perBeamPowerOffsetDb, schedule)
% Calculate scheduled normal and outer-retasked power for every satellite.

timeCount = size(angleDeg, 1);
satelliteCount = size(angleDeg, 2);
gainDbi = receiverAntennaGain(angleDeg, cfg);
wavelengthM = physconst("LightSpeed") / cfg.receiver.centerFrequencyHz;
effectiveAreaDbM2 = gainDbi + 20 * log10(wavelengthM) - ...
    10 * log10(4 * pi);
slantRangeM = double(geometry.slantRangeM);
spreadingLossDbM2 = 10 * log10(4 * pi * slantRangeM.^2);
commonReceiveDb = effectiveAreaDbM2 + 30 - ...
    cfg.rf.polarizationLossDb - cfg.rf.atmosphericLossDb - ...
    cfg.rf.implementationLossDb - ...
    earthSpaceExcessLoss(cfg, double(geometry.elevationDeg));

aggregateBeamPowerMw = zeros(timeCount, satelliteCount);
aggregateRetaskedPowerMw = zeros(timeCount, satelliteCount);
for beamIdx = cfg.resources.cochannelBeamIndices
    activeBeam = schedule.activeBeamMask(:, :, beamIdx) & ...
        geometry.contributing & geometry.zaAllowedBeam(:, :, beamIdx);
    interferenceBandwidthHz = ...
        cfg.resources.beamInterferenceBandwidthHz(beamIdx);
    bandwidthFactorDb = 10 * log10(interferenceBandwidthHz / 1e6);
    offAxisAngleDeg = double(geometry.beamOffAxisAngleDeg(:, :, beamIdx));
    halfPowerBeamwidthDeg = double( ...
        geometry.beamHalfPowerWidthDeg(:, beamIdx)).';
    scanAngleDeg = double(geometry.beamScanAngleDeg(:, beamIdx)).';
    [~, relativeGainDb] = satelliteTransmitGain(offAxisAngleDeg, ...
        halfPowerBeamwidthDeg, scanAngleDeg, cfg);
    beamCenterDistanceKm = double( ...
        geometry.beamCenterGroundSeparationKm(:, :, beamIdx));
    relativeGainDb = relativeGainDb + empiricalSidelobeCorrection( ...
        beamCenterDistanceKm, cfg);
    linkEirpDensityDbwPerMHz = calculateBeamEirpDensity( ...
        relativeGainDb, perBeamPowerOffsetDb(:, :, beamIdx), cfg);

    needsRetask = beamCenterDistanceKm < ...
        cfg.ods.outerMinimumBeamCenterDistanceKm & ...
        geometry.outerRetaskFeasible;
    retaskedOffAxisDeg = offAxisAngleDeg;
    retaskedOffAxisDeg(needsRetask) = double( ...
        geometry.outerRetaskedBeamOffAxisAngleDeg(needsRetask));
    retaskedScanDeg = scanAngleDeg + zeros(timeCount, satelliteCount);
    retaskedScanDeg(needsRetask) = double( ...
        geometry.outerRetaskedBeamScanAngleDeg(needsRetask));
    [~, retaskedRelativeGainDb] = satelliteTransmitGain( ...
        retaskedOffAxisDeg, halfPowerBeamwidthDeg, retaskedScanDeg, cfg);
    retaskedDistanceKm = beamCenterDistanceKm;
    retaskedDistanceKm(needsRetask) = ...
        cfg.ods.outerMinimumBeamCenterDistanceKm;
    retaskedRelativeGainDb = retaskedRelativeGainDb + ...
        empiricalSidelobeCorrection(retaskedDistanceKm, cfg);
    retaskedEirpDensityDbwPerMHz = calculateBeamEirpDensity( ...
        retaskedRelativeGainDb, perBeamPowerOffsetDb(:, :, beamIdx), cfg);

    pfdDensityDbwM2MHz = linkEirpDensityDbwPerMHz - spreadingLossDbM2;
    retaskedPfdDensityDbwM2MHz = retaskedEirpDensityDbwPerMHz - ...
        spreadingLossDbM2;
    if cfg.rf.maximumPfdCapEnabled
        pfdDensityDbwM2MHz = min(pfdDensityDbwM2MHz, ...
            cfg.rf.maximumPfdCapDbwM2MHz);
        retaskedPfdDensityDbwM2MHz = min(retaskedPfdDensityDbwM2MHz, ...
            cfg.rf.maximumPfdCapDbwM2MHz);
    end

    loadingDb = 10 * log10(max(double( ...
        schedule.loadByBeam(:, :, beamIdx)), ...
        realmin("double")));
    beamPowerDbm = pfdDensityDbwM2MHz + commonReceiveDb + ...
        bandwidthFactorDb + loadingDb;
    beamRetaskedPowerDbm = retaskedPfdDensityDbwM2MHz + ...
        commonReceiveDb + bandwidthFactorDb + loadingDb;
    beamPowerDbm(~activeBeam) = -Inf;
    beamRetaskedPowerDbm(~activeBeam) = -Inf;
    aggregateBeamPowerMw = aggregateBeamPowerMw + 10.^(beamPowerDbm / 10);
    aggregateRetaskedPowerMw = aggregateRetaskedPowerMw + ...
        10.^(beamRetaskedPowerDbm / 10);
end

activeScheduledLink = any(schedule.activeBeamMask,3) & ...
    geometry.contributing;
receivePowerDbm = 10 * log10(max(aggregateBeamPowerMw, realmin("double")));
retaskedPowerDbm = 10 * log10(max(aggregateRetaskedPowerMw, ...
    realmin("double")));
receivePowerDbm(~activeScheduledLink) = -Inf;
retaskedPowerDbm(~activeScheduledLink) = -Inf;
end

function correctionDb = empiricalSidelobeCorrection(distanceKm, cfg)
% Ramp the measured-envelope correction outside the geographic ZA.

if ~cfg.calibration.enabled
    correctionDb = zeros(size(distanceKm));
    return
end
fraction = (distanceKm - cfg.calibration.rampStartDistanceKm) ./ ...
    (cfg.calibration.fullCorrectionDistanceKm - ...
    cfg.calibration.rampStartDistanceKm);
fraction = min(1, max(0, fraction));
correctionDb = cfg.calibration.sidelobeCorrectionDb .* fraction;
end

function schedule = generateOperationalSchedule(cfg, geometry)
% Independent two-state schedulers, one per overlapping protected carrier.

timeCount = numel(geometry.time);
satelliteCount = size(geometry.angleDeg, 2);
beamCount = cfg.beams.numberPerSatellite;
protectedCarriers = cfg.resources.protectedCarrierIndices(:).';
activeBeamMask = false(timeCount,satelliteCount,beamCount);
loadByBeam = zeros(timeCount,satelliteCount,beamCount,"single");
selectedBeamByCarrier = zeros(timeCount,satelliteCount, ...
    numel(protectedCarriers),"uint8");

for carrierOffset = 1:numel(protectedCarriers)
    carrierIndex = protectedCarriers(carrierOffset);
    candidateBeams = find(cfg.resources.beamCarrierIndex == carrierIndex & ...
        ismember(1:beamCount,cfg.resources.cochannelBeamIndices));
    assert(~isempty(candidateBeams),"ODS:EmptyProtectedCarrier", ...
        "Protected carrier %d has no candidate beams.",carrierIndex);
    [selectedBeam,scheduledLoad] = generateCarrierSchedule( ...
        cfg,geometry,candidateBeams);
    selectedBeamByCarrier(:,:,carrierOffset) = selectedBeam;
    for beamIdx = candidateBeams
        selected = selectedBeam == beamIdx;
        activeBeamMask(:,:,beamIdx) = selected;
        beamLoad = zeros(timeCount,satelliteCount,"single");
        beamLoad(selected) = scheduledLoad(selected);
        loadByBeam(:,:,beamIdx) = beamLoad;
    end
end

schedule.activeBeamMask = activeBeamMask;
schedule.loadByBeam = loadByBeam;
schedule.selectedBeamByCarrier = selectedBeamByCarrier;
end

function [selectedBeam,scheduledLoad] = generateCarrierSchedule( ...
    cfg,geometry,candidateBeams)
% One ZA-aware scheduled beam on a single overlapping D2D carrier.

timeCount = numel(geometry.time);
satelliteCount = size(geometry.angleDeg, 2);
candidateBeams = candidateBeams(:).';
candidateCount = numel(candidateBeams);
activeState = false(timeCount, satelliteCount);
selectedBeam = zeros(timeCount, satelliteCount, "uint8");
preferredRank = uint8(randi(candidateCount, 1, satelliteCount));
activeState(1, :) = rand(1, satelliteCount) < ...
    cfg.resources.meanSatelliteCarrierDutyCycle;
pTurnOff = min(1, cfg.simulation.timeStepSec / ...
    cfg.resources.meanOnDurationSec);
pTurnOn = min(1, cfg.simulation.timeStepSec / ...
    cfg.resources.meanOffDurationSec);

for timeIdx = 1:timeCount
    if timeIdx > 1
        prior = activeState(timeIdx - 1, :);
        turnOff = prior & rand(1, satelliteCount) < pTurnOff;
        turnOn = ~prior & rand(1, satelliteCount) < pTurnOn;
        activeState(timeIdx, :) = (prior & ~turnOff) | turnOn;
        preferredRank(turnOn) = uint8(randi(candidateCount, 1, nnz(turnOn)));
    end

    allowedNow = squeeze(geometry.zaAllowedBeam( ...
        timeIdx, :, candidateBeams));
    if candidateCount == 1
        allowedNow = allowedNow(:);
    end
    chosen = zeros(1, satelliteCount, "uint8");
    unresolved = activeState(timeIdx, :);
    for offset = 0:(candidateCount - 1)
        rank = mod(double(preferredRank) - 1 + offset, candidateCount) + 1;
        linearIndex = sub2ind([satelliteCount candidateCount], ...
            1:satelliteCount, rank);
        candidateAllowed = reshape(allowedNow(linearIndex), 1, []);
        take = unresolved & candidateAllowed;
        chosen(take) = uint8(candidateBeams(rank(take)));
        unresolved(take) = false;
    end
    selectedBeam(timeIdx, :) = chosen;
end

correlation = exp(-cfg.simulation.timeStepSec / ...
    cfg.resources.meanOnDurationSec);
innovationScale = sqrt(1 - correlation^2) * ...
    cfg.resources.scheduledLoadSigma;
innovations = randn(timeCount, satelliteCount, "single");
innovations(1, :) = 0;
deviation = filter(single(innovationScale), single([1 -correlation]), ...
    innovations, [], 1);
scheduledLoad = min(1, max(cfg.resources.minimumScheduledLoad, ...
    single(cfg.resources.meanScheduledLoad) + deviation));
scheduledLoad(selectedBeam == 0) = 0;
end

function [predictionTimeSec, predictedVisible, protectionRelevant] = ...
    preparePredictionGeometry(cfg, geometry)
% Candidate-independent visibility at ODS command-execution time.

timeCount = numel(geometry.time);
lookAheadSec = cfg.ods.commandLatencySec + cfg.ods.clockErrorSec;
sampleTimeSec = cfg.simulation.timeStepSec;
timeSec = (0:(timeCount - 1)).' * sampleTimeSec;
predictionTimeSec = min(timeSec + lookAheadSec, timeSec(end));
predictedElevationDeg = interp1(timeSec, double(geometry.elevationDeg), ...
    predictionTimeSec, "linear", "extrap");
predictedVisible = predictedElevationDeg >= cfg.site.minimumElevationDeg;
protectionRelevant = geometry.visible | predictedVisible;
end
