# Project notes: change-parameters

- Changing cod `w_min` (0.5 g → 1 mg) also changes the default `h`, which
  mizer derives from `age_mat` *and* `w_min` (22.3 → 24.5), and with it `ks`.
  `new_combined_params()` in `R/combined_model.R` therefore passes the original
  `h` and `ks` explicitly.
- A size-dependent alternative food is set with `ext_encounter(params) <-`,
  which freezes the array. `E_ext` in the species parameters is still read by
  `combined_ext_encounter()` to build it, but changing `E_ext` afterwards does
  nothing until the array is rebuilt.
