function plotOdsResults(results, outputFolder)
%plotOdsResults Create baseline ODS result figures.

cfg = results.config;
figureHandle = figure(Visible="off", Color="white", Position=[100 100 1200 650]);
displayFloorDb = -100;
displayIOverNDb = max(results.selected.iOverNDb, displayFloorDb);
plot(results.selected.time, displayIOverNDb, LineWidth=1.1);
yline(cfg.receiver.protectionCriterionDb, "r--", "Protection criterion");
grid on
xlabel("UTC time");
ylabel("Scheduled all-satellite I/N (dB; display floor -100 dB)");
title(sprintf("%s strict I/N, representative MC run: outer %.1f deg, inner %.1f deg", ...
    cfg.site.name, ...
    results.selected.outerAngleDeg, results.selected.innerAngleDeg));
exportgraphics(figureHandle, fullfile(outputFolder, "aggregate_in_over_n.png"), ...
    Resolution=180);
close(figureHandle);

figureHandle = figure(Visible="off", Color="white", Position=[100 100 1200 650]);
displayPracticalIOverNDb = max(results.practical.iOverNDb, displayFloorDb);
plot(results.selected.time, displayPracticalIOverNDb, LineWidth=1.1);
yline(cfg.receiver.protectionCriterionDb, "r--", "Protection criterion");
grid on
xlabel("UTC time");
ylabel("Scheduled all-satellite I/N (dB; display floor -100 dB)");
title(sprintf("%s practical I/N, representative MC run: outer %.1f deg, inner %.1f deg", ...
    cfg.site.name, results.practical.outerAngleDeg, ...
    results.practical.innerAngleDeg));
exportgraphics(figureHandle, fullfile(outputFolder, ...
    "practical_aggregate_in_over_n.png"), Resolution=180);
close(figureHandle);

figureHandle = figure(Visible="off", Color="white", Position=[100 100 800 650]);
valid = results.searchTable.Feasible;
scatter(results.searchTable.OuterAngleDeg(~valid), ...
    results.searchTable.InnerAngleDeg(~valid), 18, ...
    results.searchTable.RobustExceedancePercent(~valid), "filled");
hold on
scatter(results.searchTable.OuterAngleDeg(valid), ...
    results.searchTable.InnerAngleDeg(valid), 32, ...
    results.searchTable.RobustExceedancePercent(valid), "filled", Marker="square");
plot(results.selected.outerAngleDeg, results.selected.innerAngleDeg, ...
    "kp", MarkerSize=15, MarkerFaceColor="yellow");
plot(results.practical.outerAngleDeg, results.practical.innerAngleDeg, ...
    "kd", MarkerSize=11, MarkerFaceColor="cyan");
hold off
grid on
xlabel("Outer angle (deg)");
ylabel("Inner angle (deg)");
title(sprintf("ODS angle search (yellow=strict, cyan=%.1f%% availability)", ...
    cfg.reporting.practicalAvailabilityPercent));
colorbar
exportgraphics(figureHandle, fullfile(outputFolder, "angle_search.png"), ...
    Resolution=180);
close(figureHandle);

figureHandle = figure(Visible="off", Color="white", Position=[100 100 800 650]);
sortedPracticalIOverN = sort(max(results.practical.iOverNDb, displayFloorDb));
practicalProbability = (1:numel(sortedPracticalIOverN)) / ...
    numel(sortedPracticalIOverN);
plot(sortedPracticalIOverN, practicalProbability, LineWidth=1.5);
xline(cfg.receiver.protectionCriterionDb, "r--", "Protection criterion");
grid on
xlabel("Aggregate scheduled I/N (dB)");
ylabel("Empirical cumulative probability");
title("Practical selected-angle I/N CDF, representative MC run");
exportgraphics(figureHandle, fullfile(outputFolder, ...
    "practical_aggregate_in_over_n_cdf.png"), Resolution=180);
close(figureHandle);

figureHandle = figure(Visible="off", Color="white", Position=[100 100 800 650]);
sortedIOverN = sort(max(results.selected.iOverNDb, displayFloorDb));
probability = (1:numel(sortedIOverN)) / numel(sortedIOverN);
plot(sortedIOverN, probability, LineWidth=1.5);
xline(cfg.receiver.protectionCriterionDb, "r--", "Protection criterion");
grid on
xlabel("Aggregate all-satellite I/N (dB)");
ylabel("Empirical cumulative probability");
title("Strict selected-angle I/N CDF, representative MC run");
exportgraphics(figureHandle, fullfile(outputFolder, "aggregate_in_over_n_cdf.png"), ...
    Resolution=180);
close(figureHandle);
end
