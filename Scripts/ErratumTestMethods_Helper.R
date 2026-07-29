# TESTING of ARRANGE and desc() vs. - ++++++

# arrange() is effectiveoly ignored.
# desc(), withing group_by(), is needed and is the same as using - on numeric variables


# original - apparently arrange is ignored.
set.seed(128)
genOrd_l2h_O <- 1:reps |> map(~pre_meta_Guild_use |>
                              arrange(desc(generality)) |>
                              group_by(generality)|>
                              sample_frac() |>
                              filter(node != "BASAL NODE") |>
                              pull(node)) |> 
  set_names("genOrd_l2h_O")

set.seed(128)
genOrd_l2h_O1 <- 1:reps |> map(~pre_meta_Guild_use |>
                                arrange(desc(generality)) |>
                                group_by(desc(generality))|>
                                sample_frac() |>
                                filter(node != "BASAL NODE") |>
                                pull(node)) |> 
  set_names("genOrd_l2h_O1")

set.seed(128)
genOrd_l2h_O2 <- 1:reps |> map(~pre_meta_Guild_use |>
                                arrange(generality) |>
                                group_by(generality) |>
                                sample_frac() |>
                                filter(node != "BASAL NODE") |>
                                pull(node)) |> 
  set_names("genOrd_l2h_O2")

set.seed(128)
genOrd_l2h_O3 <- 1:reps |> map(~pre_meta_Guild_use |>
                                 arrange(generality) |>
                                 group_by(desc(generality)) |>
                                 sample_frac() |>
                                 filter(node != "BASAL NODE") |>
                                 pull(node)) |> 
  set_names("genOrd_l2h_O3")

set.seed(128)
genOrd_l2h_O4 <- 1:reps |> map(~pre_meta_Guild_use |>
                                 group_by(desc(generality)) |>
                                 sample_frac() |>
                                 filter(node != "BASAL NODE") |>
                                 pull(node)) |> 
  set_names("genOrd_l2h_O4")

# checking ----
# 0 and 3 are the same (group_by is normal generality)
# 1,3,4 are the same with group_by() either desc() or -.
bind_cols(genOrd_l2h_O, genOrd_l2h_O1, genOrd_l2h_O2, genOrd_l2h_O3, genOrd_l2h_O4)

## - category

set.seed(128)
out1 <- pre_meta_Guild_use |>
  mutate(tiering_rev = forcats::fct_rev(tiering_i2p)) |>
  group_by(tiering_rev) |>
  sample_frac() |>
  filter(node != "BASAL NODE") |>
  pull(node)

set.seed(128)
out2 <- pre_meta_Guild_use |>
  group_by(desc(tiering_i2p)) |>
  sample_frac() |>
  filter(node != "BASAL NODE") |>
  pull(node)

bind_cols(out1, out2)
