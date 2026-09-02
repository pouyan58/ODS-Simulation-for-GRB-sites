function cfg = defaultGrbOdsSpotBeamConfig(site, projectFolder, outputRoot, zaRadiusKm)
%defaultGrbOdsSpotBeamConfig Return the public 10-spot GRB configuration.

arguments
    site (1,1) struct
    projectFolder (1,1) string = string(fileparts(fileparts( ...
        mfilename("fullpath"))))
    outputRoot (1,1) string = projectFolder
    zaRadiusKm (1,1) double = NaN
end

siteName = string(site.Name);
siteId = string(site.SiteId);
if isnan(zaRadiusKm)
    if any(string(site.ReceiverId) == ["fairmont_wv","wallops_va"])
        zaRadiusKm = 30;
    else
        zaRadiusKm = 50;
    end
end
assert(isfinite(zaRadiusKm) && zaRadiusKm >= 0, ...
    "ODS:InvalidZaRadius","ZA radius must be finite and nonnegative.");
diameterM = 4.5;
wavelengthM = physconst("LightSpeed") / 1686.6e6;
diameterWavelengths = diameterM / wavelengthM;

cfg.meta.modelName = "GRB ODS Spot Beam Simulator";
cfg.meta.modelVersion = "1.0.0";
cfg.meta.releaseId = "GRB-ODS-SPOT-BEAM-SIMULATOR-PUBLIC-1.0.0";
cfg.meta.baselineId = "GRB-10-SPOT-60-FOR-30DBW-MHZ";
cfg.meta.experimentId = "GRB-ODS-SPOT-" + upper(siteId);
cfg.meta.caseName = compose("GRB ODS Spot Beam Simulator - %s | %s", ...
    siteName, string(site.Satellite));
cfg.meta.disclaimer = "Public engineering sensitivity simulation, not a coordination " + ...
    "determination or operator commitment. The -27 dBi far-angle receive level, traffic " + ...
    "scheduler, NRAO transfer calibration, and ODS actions require site/operator validation.";
cfg.meta.siteCatalogFile = string(site.SourceCatalog);
cfg.meta.siteCatalogRow = site.SourceRow;

cfg.site.name = siteName;
cfg.site.siteId = siteId;
cfg.site.receiverId = string(site.ReceiverId);
cfg.site.latitudeDeg = site.LatitudeDeg;
cfg.site.longitudeDeg = site.LongitudeDeg;
cfg.site.heightM = site.FeedHeightM;
cfg.site.antennaFeedHeightAglM = site.FeedHeightM;
cfg.site.heightSourceStatus = "Catalog feed height above local ground used as a geodetic-height " + ...
    "proxy; replace with surveyed antenna-reference height when available.";
cfg.site.minimumElevationDeg = 0;
cfg.site.visibilityConvention = "local geometric horizon (0-degree minimum LEO elevation)";
cfg.site.boresight.type = "geoLongitude";
cfg.site.boresight.satelliteAssignment = string(site.Satellite);
cfg.site.boresight.geoLongitudeDeg = site.BoresightLongitudeDeg;
cfg.site.boresight.longitudeSource = string(site.BoresightLongitudeSource);

cfg.receiver.centerFrequencyHz = 1686.6e6;
cfg.receiver.bandwidthHz = 10.9e6;
cfg.receiver.noiseTemperatureK = 120;
cfg.receiver.protectionCriterionDb = -6;
cfg.receiver.allowedExceedancePercent = 0.1;
cfg.receiver.referencePlane = ...
    "GRB antenna terminals, one circular-polarization receiver chain";
cfg.receiver.service = "GOES Rebroadcast (GRB)";
cfg.receiver.modulation = "QPSK DVB-S2 (current NOAA default)";
cfg.receiver.polarization = "Dual circular RHCP and LHCP; simulation " + ...
    "protects one worst-case co-polar chain";
cfg.receiver.minimumRequiredGOverTDbK = 15.2;
cfg.receiver.profileStatus = "NOAA nominal worst-location user receiver " + ...
    "applied at each catalog location";
cfg.receiver.profileSource = "https://www.ospo.noaa.gov/operations/goes/grb/";
cfg.receiver.downlinkSpecificationSource = ...
    "https://goes-r.noaa.gov/users/docs/GRB_downlink.pdf";
cfg.receiver.protectionCriterionStatus = ...
    "Engineering -6 dB I/N criterion at 99.9-percent modeled availability.";

cfg.antenna.model = "ituRS580Aperec015SmallExtension";
cfg.antenna.diameterM = diameterM;
requiredPeakGainDbi = cfg.receiver.minimumRequiredGOverTDbK + ...
    10*log10(cfg.receiver.noiseTemperatureK);
cfg.antenna.efficiency = 10^(requiredPeakGainDbi/10) / ...
    (pi*diameterM/wavelengthM)^2;
cfg.antenna.sidelobeFloorDbi = -10;
cfg.antenna.sidelobeFloorUsage = "unused by the selected S.580/APEREC015 model";
cfg.antenna.pointingBiasDeg = 0;
cfg.antenna.farAngleGainDbiOverride = -27;
cfg.antenna.farAngleOverrideStartDeg = 10^(42/25);
cfg.antenna.farAngleOverrideStatus = "Public sensitivity assumption applied " + ...
    "beyond the S.580 far-angle breakpoint; replace with measured site data.";
cfg.antenna.reference = "ITU-R S.580-6 with attached APEREC015 Appendix 8 extension where applicable";
cfg.antenna.applicabilityNote = compose( ...
    "D/lambda = %.6f; NOAA minimum G/T %.1f dB/K with %.0f K gives " + ...
    "modeled peak gain %.3f dBi and derived aperture efficiency %.4f.", ...
    diameterWavelengths,cfg.receiver.minimumRequiredGOverTDbK, ...
    cfg.receiver.noiseTemperatureK,requiredPeakGainDbi,cfg.antenna.efficiency);

cfg.constellation.shells = [ ...
    struct("name", "DTC-53", "altitudeKm", 340.0, "inclinationDeg", 53.0, ...
        "satellites", 325, "planes", 13, "phasing", 1); ...
    struct("name", "DTC-43", "altitudeKm", 355.0, "inclinationDeg", 43.0, ...
        "satellites", 325, "planes", 13, "phasing", 1)];
cfg.constellation.sourceNote = "Engineering two-shell D2C baseline: 325 satellites at 340 km/53 degrees and 325 at 355 km/43 degrees.";

cfg.beams.numberPerSatellite = 10;
cfg.beams.boresightDirection = "10 snapshot steering positions inside the field of regard";
cfg.beams.transmissionBoundary = "nadirConePackedSpots";
cfg.beams.coneHalfAngleDeg = 60;
cfg.beams.fullConeAngleDeg = 120;
cfg.beams.layoutModel = "optimizedTriangularClosePackedGeodesicCells";
cfg.beams.layoutOrientation = "Earth-fixed local east/north snapshot steering lattice";
cfg.beams.layoutFile = fullfile(projectFolder, "evidence", ...
    "beam_packing_layout.csv");
cfg.beams.nonOverlapToleranceDeg = 1e-7;
cfg.beams.transmitPatternModel = "ituRS1528LeoPerBeam";
cfg.beams.spotHalfPowerRadiusDeg = 2.0;
cfg.beams.spotFullHalfPowerBeamwidthDeg = 2 * cfg.beams.spotHalfPowerRadiusDeg;
cfg.beams.spotWidthInterpretation = "The 60-degree value is steering field of regard, not beamwidth. Each spot has a separately configurable -3 dB angular radius.";
cfg.beams.halfPowerBoundaryDefinition = "Configured satellite-view -3 dB angular radius for each independently steered spot.";
cfg.beams.nearSidelobeCrossPointDb = -6.75;
cfg.beams.farOutSidelobeGainDbi = 0;
cfg.beams.offAxisMaskMarginDb = 0.0;
cfg.beams.offAxisMaskTransitionStartNormalized = 1.0;
cfg.beams.offAxisMaskTransitionEndNormalized = 1.5;
cfg.beams.offAxisMaskMarginStatus = "Optional calibration assumption; zero in the default case so the S.1528 reference envelope is retained without extra satellite-pattern credit.";
cfg.beams.scanLossModel = "projectedApertureCosine";
cfg.beams.scanLossCosineExponent = 1;
cfg.beams.scanLossCompensationEnabled = true;
cfg.beams.scanLossCompensationNote = "Modeled power control compensates projected-aperture scan loss so every spot retains the same peak EIRP-density ceiling.";
cfg.beams.patternReference = "ITU-R S.1528-0 LEO reference pattern plus explicit configurable off-axis mask margin.";
cfg.beams.footprintModel = "10 narrow spots at snapshot steering positions; every visible satellite can couple through spot sidelobes";
cfg.beams.frequencyReuseAssumption = "Four 5-MHz carriers with cyclic " + ...
    "assignment; the GRB passband overlaps three carriers.";
cfg.beams.loadingInterpretation = "Independent two-state scheduler with " + ...
    "at most one active spot on each GRB-overlapping carrier per satellite.";
cfg.beams.cellEdgeInterpretation = "The configured spot half-power radius is independent of the 60-degree field of regard.";
cfg.beams.sourceNote = "The 10-point lattice defines snapshot steering locations only; it does not assert continuous tiling of the field of regard.";

cfg.za.enabled = true;
cfg.za.protectedCoreRadiusKm = 30;
cfg.za.minimumActiveBeamCenterDistanceKm = 30;
cfg.za.definition = "No D2D served user/active beam center inside 30 km of the protected receiver";
cfg.za.sourceStatus = "User-defined 30-km geographic avoidance radius; not a coordination result.";

cfg.resources.carrierBandwidthHz = 5e6;
cfg.resources.carrierEdgesHz = 1675e6:5e6:1695e6;
cfg.resources.carrierCount = 4;
cfg.resources.beamCarrierIndex = mod(0:(cfg.beams.numberPerSatellite - 1), cfg.resources.carrierCount) + 1;
receiverPassbandHz = cfg.receiver.centerFrequencyHz + ...
    0.5*cfg.receiver.bandwidthHz*[-1 1];
carrierLowHz = cfg.resources.carrierEdgesHz(1:end-1);
carrierHighHz = cfg.resources.carrierEdgesHz(2:end);
carrierOverlapHz = max(0,min(carrierHighHz,receiverPassbandHz(2)) - ...
    max(carrierLowHz,receiverPassbandHz(1)));
cfg.resources.receiverPassbandHz = receiverPassbandHz;
cfg.resources.carrierOverlapBandwidthHz = carrierOverlapHz;
cfg.resources.protectedCarrierIndices = find(carrierOverlapHz > 0);
cfg.resources.protectedCarrierIndex = cfg.resources.protectedCarrierIndices(1);
cfg.resources.cochannelBeamIndices = find(ismember( ...
    cfg.resources.beamCarrierIndex,cfg.resources.protectedCarrierIndices));
cfg.resources.maximumSimultaneousCochannelBeamsPerSatellite = ...
    numel(cfg.resources.protectedCarrierIndices);
cfg.resources.schedulerMode = "oneBeamPerProtectedCarrier";
cfg.resources.beamInterferenceBandwidthHz = ...
    carrierOverlapHz(cfg.resources.beamCarrierIndex);
cfg.resources.frequencyPlan = "Four contiguous 5-MHz carriers; the 10.9-MHz " + ...
    "GRB passband overlaps carriers 2, 3, and 4 by 3.85, 5.00, and 2.05 MHz.";
cfg.resources.frequencyPlanStatus = "Exact rectangular spectral-overlap " + ...
    "integration; transmitter and receiver filter skirts are not modeled.";
cfg.resources.meanSatelliteCarrierDutyCycle = 0.20;
cfg.resources.meanOnDurationSec = 30;
cfg.resources.meanOffDurationSec = 120;
cfg.resources.meanScheduledLoad = 0.60;
cfg.resources.scheduledLoadSigma = 0.15;
cfg.resources.minimumScheduledLoad = 0.20;
cfg.resources.trafficModel = "Independent two-state burst scheduler on each " + ...
    "overlapping carrier; one spot per carrier per satellite; ZA-aware reassignment";
cfg.resources.trafficSourceStatus = "Engineering traffic assumption requiring operator validation.";

cfg.rf.studyBandHz = [1675e6 1695e6];
cfg.rf.sourceBandHz = [1990e6 1995e6];
cfg.rf.eirpDensityDbwPerMHz = 30;
cfg.rf.maximumBeamEirpDensityDbwPerMHz = 30;
cfg.rf.eirpBasis = "User-selected 30 dBW/MHz maximum EIRP density for every spot beam.";
cfg.rf.referencePeakAntennaGainDbi = 38;
cfg.rf.perCellEirpConvention = "Every active spot has the same 30 dBW/MHz maximum boresight EIRP density; pattern and load reduce off-axis/instantaneous EIRP only.";
cfg.rf.filedPolarizations = "RHCP and LHCP";
cfg.rf.maximumPfdCapEnabled = true;
cfg.rf.maximumPfdCapDbwM2MHz = -80;
cfg.rf.maximumPfdSource = "Prior FCC engineering-table assumption; applied per beam after range loss.";
cfg.rf.meanActiveProbability = cfg.resources.meanSatelliteCarrierDutyCycle;
cfg.rf.activitySourceStatus = "Compatibility alias for the burst scheduler.";
cfg.rf.polarizationLossDb = 0;
cfg.rf.polarizationStatus = "Worst-case co-polar coupling.";
cfg.rf.atmosphericLossDb = 0;
cfg.rf.atmosphericLossUsage = "zero when the propagation model is enabled";
cfg.rf.implementationLossDb = 1;

cfg.propagation.enabled = true;
cfg.propagation.atmosphericModel = "ituRP618Lookup";
cfg.propagation.atmosphericAnnualExceedancePercent = 5;
cfg.propagation.minimumLookupElevationDeg = 5;
cfg.propagation.belowLookupTreatment = "Use the 5-degree P.618 value for visible 0-to-5-degree links.";
cfg.propagation.atmosphericLookupFile = fullfile(projectFolder, "data", ...
    string(site.ReceiverId) + "_grb_p618_1686p6MHz_4p5m_5pct.csv");
cfg.propagation.atmosphericReference = "ITU-R P.618 lookup generated with MATLAB Satellite Communications Toolbox";
cfg.propagation.clutterModel = "ituRP2108OpenRuralHeightGain";
cfg.propagation.clutterRepresentativeHeightM = 10;
cfg.propagation.receiverAntennaHeightAglM = site.FeedHeightM;
cfg.propagation.clutterReference = "ITU-R P.2108-1 open/rural height-gain model; provisional heights";

cfg.aggregation.selectionMode = "allEligible";
cfg.aggregation.powerCombinationMode = "linearSum";
cfg.aggregation.maxContributingSatellites = sum([cfg.constellation.shells.satellites]);
cfg.aggregation.contributorLimitApplied = false;
cfg.aggregation.selectionMetric = "all visible satellites with one or more " + ...
    "scheduled ZA-permitted spots on the three carriers overlapping GRB";
cfg.aggregation.selectionOrder = "calculate every overlapping-carrier link, " + ...
    "apply per-beam ODS action, integrate exact spectral overlap, and sum powers linearly";
cfg.aggregation.dominantSelectionStage = "notApplicableAllSatelliteLinearSum";
cfg.aggregation.dominantDefinition = "All eligible scheduled satellites contribute.";

cfg.ods.steerMitigationDb = 0;
cfg.ods.muteMitigationDb = 80;
cfg.ods.deepNullMitigationDb = cfg.ods.muteMitigationDb;
cfg.ods.measuredMuteLowerBoundDb = 20;
cfg.ods.outerMinimumBeamCenterDistanceKm = 180;
cfg.ods.minimumOperationalOuterAngleDeg = 0;
cfg.ods.minimumOperationalOuterAngleStatus = "No operating floor; the search selects the unconstrained minimum modeled angle pair.";
cfg.ods.outerRetaskMode = "retask protected-carrier spots to a feasible center at least 180 km from the protected receiver and recompute gain";
cfg.ods.innerProtectedCarrierShutdownEnabled = true;
cfg.ods.actionMode = "physical outer beam retasking and inner protected-carrier shutdown";
cfg.ods.unaffectedBeamOperation = "Outside the outer cone operation is unchanged; in the annulus the protected-carrier spot is retasked; in the inner cone it is shut down.";
cfg.ods.mitigationSourceStatus = "180-km retasking is an NRAO-informed engineering assumption; inner action is exact protected-carrier shutdown; no outer-angle operating floor is applied.";
cfg.ods.demonstratedMuteAngleDeg = 0.5;
cfg.ods.demonstratedSteerExampleAngleDeg = 2;
cfg.ods.commandLatencySec = 15;
cfg.ods.clockErrorSec = 1;
cfg.ods.hysteresisDeg = 0.5;
cfg.ods.failSafeProbability = 0;
cfg.ods.failSafeMitigationDb = 80;

cfg.calibration.enabled = true;
cfg.calibration.evidenceFile = fullfile(projectFolder, "evidence", "nrao_dtc_calibration.csv");
cfg.calibration.sidelobeCorrectionDb = -18;
cfg.calibration.rampStartDistanceKm = cfg.za.minimumActiveBeamCenterDistanceKm;
cfg.calibration.fullCorrectionDistanceKm = 150;
cfg.calibration.sourceFrequencyHz = 1992.5e6;
cfg.calibration.transferMode = "Relative NRAO D2D-test sidelobe correction " + ...
    "transferred from 1990-1995 MHz to 1675-1695 MHz";
cfg.calibration.sourceStatus = "Provisional cross-service engineering calibration; receiver- and band-specific measurements supersede it.";

cfg.uncertainty.powerSigmaDb = 3;
cfg.uncertainty.powerTreatment = "Independent beam uncertainty below the 30 dBW/MHz ceiling; positive excursions clipped at the ceiling.";
cfg.uncertainty.ephemerisAngleSigmaDeg = 0.05;
cfg.uncertainty.pointingSigmaDeg = 0.10;
cfg.uncertainty.angularGeometryTreatment = "Angular errors are clamped to visible sky.";
cfg.uncertainty.loadingSigma = cfg.resources.scheduledLoadSigma;
cfg.uncertainty.temporalCorrelationSec = cfg.resources.meanOnDurationSec;
cfg.uncertainty.robustPercentile = 95;

cfg.reporting.strictAllowedExceedancePercent = 0;
cfg.reporting.practicalAllowedExceedancePercent = 0.1;
cfg.reporting.practicalAvailabilityPercent = 99.9;
cfg.reporting.practicalIOverNPercentile = 99.9;
cfg.reporting.note = "The primary result allows at most 0.1-percent robust " + ...
    "exceedance and requires the robust 99.9th-percentile I/N to be at or " + ...
    "below -6 dB. A zero-exceedance result is retained as a diagnostic.";

cfg.simulation.startTime = datetime(2026, 8, 1, 0, 0, 0, TimeZone="UTC");
cfg.simulation.duration = hours(2);
cfg.simulation.timeStepSec = 2;
cfg.simulation.monteCarloRuns = 4;
cfg.simulation.randomSeed = 16751695;

cfg.search.outerAngleRangeDeg = unique([0 0.5 1 2 3 5:5:145]);
cfg.search.innerAngleRangeDeg = unique([0 0.5 1 2 3 5:5:145]);
cfg.search.adaptiveRefinementEnabled = true;
cfg.search.refinementStepDeg = 1;
cfg.search.refinementHalfWidthDeg = 7;
cfg.search.finalStepDeg = 0.5;
cfg.search.finalHalfWidthDeg = 1.5;
cfg.search.serviceImpactWeight = 0;
cfg.search.serviceImpactNote = "Rank feasible pairs by outer plus inner angle; " + ...
    "report the 99.9-percent-availability pair as primary.";

cfg.output.writeDetailedTimeSeries = true;
cfg.output.folder = fullfile(outputRoot, "results_" + siteId);
cfg = configureGrbOdsZaRadius(cfg,zaRadiusKm);
end
