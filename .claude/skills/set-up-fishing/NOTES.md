# Project notes: set-up-fishing

- The model has two gears, both at effort 0 in the saved states: "Trawl" on
  cod (sigmoid length, l50 = 35 cm) and "Sprat trawl" on sprat (knife edge at
  1 g). Give effort as a named vector or a time × gear array.
- Sprat trawl and cod invasion fitness R0 in the cod-free state:
  - Knife edge 1 g: F = 0.25 → R0 1.18, F = 0.5 → R0 1.55 (fixed points).
    F ≥ 1 sends sprat into limit cycles, and R0 at one moment means nothing.
  - Knife edge 3 g: F = 0.25 → R0 0.2 only. F ≥ 0.5 → sprat cycles.
- Recovery by prey thinning (cod-free state with 1 % cod immigration, sprat
  trawl at F = 0.5 from year 10): 10, 15 or 20 years are not enough, 30 years
  are. F = 0.3 for 50 years also works.
