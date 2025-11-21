join_unique <- function(x) {
  x |>
    unique() |>
    str_trim() |>
    discard(~ is.na(.x) || .x == "") |>
    sort() |>
    (\(v) if (length(v) == 0) NA_character_ else str_c(v, collapse = "; "))()
}

count_hbd <- function(smiles) {

  n_h <- stringr::str_count(smiles, "N(H)")

  # Count O–H donors
  o_h <- stringr::str_count(smiles, "O(H)")

  n_explicit <- stringr::str_count(smiles, "\\[N[Hh][^]]*\\]")
  o_explicit <- stringr::str_count(smiles, "\\[O[Hh][^]]*\\]")


  n_h + o_h + n_explicit + o_explicit
}

count_hba <- function(smiles) {

  n_count <- stringr::str_count(smiles, "N|n")


  o_count <- stringr::str_count(smiles, "O")

  n_count + o_count
}

calc_clogp <- function(smiles) {
  mol <- tryCatch(parse.smiles(smiles)[[1]], error = function(e) NULL)
  if (is.null(mol)) return(NA_real_)

  tryCatch({
    do.aromaticity(mol)
    set.atom.types(mol)
    convert.implicit.to.explicit(mol)
  }, error = function(e) return(NA_real_))

  # Direct descriptor evaluation
  out <- tryCatch(
    eval.desc(
      mol,
      "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor"
    ),
    error = function(e) return(NA_real_)
  )

  as.numeric(out[[1]])
}
