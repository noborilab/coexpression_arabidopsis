test_that("example config is valid YAML with required top-level keys", {
  cfg_path <- file.path("..", "..", "config", "example_config.yaml")
  skip_if_not(file.exists(cfg_path), "example_config.yaml not found — run tests from package root")
  cfg <- yaml::read_yaml(cfg_path)
  expect_type(cfg, "list")
  for (key in c("input", "stratum", "estimation", "output")) {
    expect_true(key %in% names(cfg), info = paste("missing required key:", key))
  }
})

test_that("example config estimation mode is a valid value", {
  cfg_path <- file.path("..", "..", "config", "example_config.yaml")
  skip_if_not(file.exists(cfg_path))
  cfg <- yaml::read_yaml(cfg_path)
  expect_true(cfg$estimation$mode %in% c("pseudobulk", "singlecellggm"))
})
