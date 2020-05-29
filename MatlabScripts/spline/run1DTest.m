y = [1, 2, 3, 5, 8, 5, 2, 3, 9, 8, 6, 5, 5, 5, 5, 5];

query = 3:0.1:length(y)-2;

[t, centers, controlPts] = interpSP(y, query, 9);
figure
hold on
scatter(1:length(y), y)
scatter(centers, controlPts)
legend("Input points", "Control points")
plot(query, t)
