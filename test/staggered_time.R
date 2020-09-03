# test data
y <- data.frame(
  id = 1:3,
  start = c(7, 3, 3),
  end = c(10, 6, 5),
  event = c(0, 0, 1),
  sampling = c(0, 1, 2),
  stringsAsFactors = FALSE
)

# recalculate the end time if start zeroed
y$end_zero <- y$end - y$start

ragged <- SamplingDesignTools::compute_km_weights(
  cohort = y,
  t_name = "end",
  t_start_name = "start",
  y_name = "event",
  id_name = "id",
  sample_stat = y$sampling,
  n_per_case = 1)

zeroed <-SamplingDesignTools::compute_km_weights(
  cohort = y,
  t_name = "end_zero",
  y_name = "event",
  id_name = "id",
  sample_stat = y$sampling,
  n_per_case = 1) 

identical(ragged, zeroed)  # FALSE
