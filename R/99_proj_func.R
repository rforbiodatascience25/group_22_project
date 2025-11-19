join_unique <- function(x) {
  x |>
    unique() |>
    str_trim() |>
    discard(~ is.na(.x) || .x == "") |>
    sort() |>
    (\(v) if (length(v) == 0) NA_character_ else str_c(v, collapse = "; "))()
}
