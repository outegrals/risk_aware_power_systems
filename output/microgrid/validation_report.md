# Synthetic microgrid study

30 paired 48-hour trajectories per setting; 15-minute steps. Cash in synthetic currency units; no real data calibration. Each confidence interval is an unadjusted, approximate paired Student-t interval across independent replicates, not a simultaneous familywise claim.

| n | s | RN cash/operator | fixed ERS | retuned | gain [95% interval] | shortfall kWh/operator RN / ERS / retuned |
|---:|---:|---:|---:|---:|---|---|
| 10 | 0 | 1641.25 | 1733.59 | 1733.62 | 92.36 [79.34, 105.39] | 0.00 / 0.00 / 0.00 |
| 10 | 0.5 | 1642.28 | 1731.50 | 1731.55 | 89.27 [76.82, 101.72] | 31.65 / 5.41 / 5.88 |
| 10 | 1 | 1643.53 | 1724.27 | 1724.38 | 80.85 [69.79, 91.91] | 92.33 / 20.01 / 21.72 |
| 10 | 2 | 1634.55 | 1691.47 | 1691.94 | 57.39 [49.32, 65.47] | 352.86 / 93.71 / 107.98 |
| 50 | 0 | 1575.61 | 1710.39 | 1710.43 | 134.82 [127.26, 142.37] | 0.00 / 0.00 / 0.00 |
| 50 | 0.5 | 1582.38 | 1706.89 | 1706.95 | 124.57 [117.62, 131.52] | 39.50 / 5.35 / 5.84 |
| 50 | 1 | 1587.77 | 1698.68 | 1698.81 | 111.04 [104.86, 117.22] | 98.23 / 18.42 / 20.52 |
| 50 | 2 | 1590.79 | 1665.22 | 1665.79 | 75.00 [70.91, 79.08] | 396.21 / 94.39 / 107.80 |
| 100 | 0 | 1549.35 | 1690.23 | 1690.26 | 140.91 [135.80, 146.03] | 0.00 / 0.00 / 0.00 |
| 100 | 0.5 | 1558.32 | 1687.33 | 1687.39 | 129.06 [124.27, 133.86] | 39.34 / 5.69 / 6.17 |
| 100 | 1 | 1565.08 | 1680.24 | 1680.37 | 115.29 [110.66, 119.92] | 97.26 / 18.52 / 20.84 |
| 100 | 2 | 1571.53 | 1648.52 | 1649.08 | 77.56 [73.44, 81.67] | 404.93 / 98.97 / 113.20 |

All strategies share the same bid-independent physical controller and therefore the same SOC, imports, curtailment, and battery throughput. Shortfall measures contracted export missing at settlement, not unserved local load. Import purchases and battery throughput costs are included in cash. No terminal energy value or recharge requirement is imposed; identical terminal states make the paired comparison valid but absolute cash is not optimized lifetime profit.

Wall time for this validation: 10.40 s; this includes validation and physical simulation and is not a distributed deployment benchmark.
