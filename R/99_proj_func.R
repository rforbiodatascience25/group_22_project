join_unique <- function(x) {
  x |>
    unique() |>
    str_trim() |>
    discard(~ is.na(.x) || .x == "") |>
    sort() |>
    (\(v) if (length(v) == 0) NA_character_ else str_c(v, collapse = "; "))()
}

molecular_weight <- function(formula) {

  atomic_weights <- c(
    H = 1.008, He = 4.003, Li = 6.941, Be = 9.012, B = 10.81,
    C = 12.011, N = 14.007, O = 15.999, F = 19.00, Ne = 20.18,
    Na = 22.99, Mg = 24.305, Al = 26.982, Si = 28.086, P = 30.974,
    S = 32.06, Cl = 35.45, Ar = 39.95, K = 39.098, Ca = 40.078,
    Sc = 44.956, Ti = 47.867, V = 50.942, Cr = 51.996, Mn = 54.938,
    Fe = 55.845, Co = 58.933, Ni = 58.693, Cu = 63.546, Zn = 65.38,
    Ga = 69.723, Ge = 72.630, As = 74.922, Se = 78.971, Br = 79.904,
    Kr = 83.798, Rb = 85.468, Sr = 87.62, Y = 88.906, Zr = 91.224,
    Nb = 92.906, Mo = 95.95, Ru = 101.07, Rh = 102.91, Pd = 106.42,
    Ag = 107.87, Cd = 112.41, In = 114.82, Sn = 118.71, Sb = 121.76,
    Te = 127.60, I = 126.90, Xe = 131.29, Cs = 132.91, Ba = 137.33,
    La = 138.91, Ce = 140.12
  )

  formula |>
    stringr::str_extract_all("[A-Z][a-z]?[0-9]*") |>   # ← FIXED
    purrr::flatten_chr() |>
    purrr::discard(~ is.na(.x) || .x == "") |>
    purrr::map_dbl(function(part) {
      element <- stringr::str_extract(part, "^[A-Z][a-z]?")
      count   <- stringr::str_extract(part, "[0-9]+") |> as.numeric()
      if (is.na(count)) count <- 1
      atomic_weights[[element]] * count
    }) |>
    sum()
}




count_nitrogen <- function(formula) {
  match <- formula |> stringr::str_extract("N[0-9]*")
  if (is.na(match)) return(0)

  count <- match |> stringr::str_remove("N") |> as.numeric()
  if (is.na(count)) 1 else count
}


count_oxygen <- function(formula) {
  match <- formula |> stringr::str_extract("O[0-9]*")
  if (is.na(match)) return(0)

  count <- match |> stringr::str_remove("O") |> as.numeric()
  if (is.na(count)) 1 else count
}


hbond_donors <- function(smiles) {
  nh <- smiles |>
    stringr::str_extract_all("N(H+)") |>
    purrr::map_int(~ stringr::str_count(.x, "H"))

  oh <- smiles |>
    stringr::str_extract_all("O(H+)") |>
    purrr::map_int(~ stringr::str_count(.x, "H"))

  nh + oh
}

hbon_acceptors <- function(formula) {
  count_nitrogen(formula) + count_oxygen(formula)
}



