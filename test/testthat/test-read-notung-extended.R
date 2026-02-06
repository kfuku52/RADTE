# Extended read_notung_parsable tests

# --- Single duplication event ---

test_that("read_notung_parsable handles single duplication event", {
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    "1\t0\t0\t0\t3\t5\t3\t0.1,1.0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tn1\troot_sp\tancestor_sp"
  ), tmp)
  on.exit(unlink(tmp))

  result <- read_notung_parsable(tmp, mode = "D")
  expect_equal(nrow(result), 1)
  expect_equal(result$event[1], "D")
  expect_equal(result$gn_node[1], "n1")
  expect_equal(result$lower_sp_node[1], "root_sp")
  expect_equal(result$upper_sp_node[1], "ancestor_sp")
})

# --- Many duplication events ---

test_that("read_notung_parsable handles many duplication events", {
  lines <- c(
    "5\t0\t0\t0\t6\t11\t3\t0.01,2.0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tn1\tsp1\tsp_root",
    "#D\tn2\tsp2\tsp_root",
    "#D\tn3\tsp1\tsp_root",
    "#D\tn4\tsp2\tsp_root",
    "#D\tn5\tsp1\tsp_root"
  )
  tmp <- tempfile(fileext = ".txt")
  writeLines(lines, tmp)
  on.exit(unlink(tmp))

  result <- read_notung_parsable(tmp, mode = "D")
  expect_equal(nrow(result), 5)
  expect_true(all(result$event == "D"))
  expect_equal(result$gn_node, c("n1", "n2", "n3", "n4", "n5"))
})

# --- Duplication with species names containing underscores ---

test_that("read_notung_parsable handles species with underscore names", {
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    "1\t0\t0\t0\t3\t5\t3\t0.1,1.0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tn1\tHomo_sapiens\tn_root"
  ), tmp)
  on.exit(unlink(tmp))

  result <- read_notung_parsable(tmp, mode = "D")
  expect_equal(nrow(result), 1)
  expect_equal(result$lower_sp_node[1], "Homo_sapiens")
  expect_equal(result$upper_sp_node[1], "n_root")
})

# --- Parsable file with extra whitespace lines ---

test_that("read_notung_parsable is robust to trailing blank lines", {
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    "1\t0\t0\t0\t3\t5\t3\t0.1,1.0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tn1\tsp1\tsp_root",
    "",
    "",
    ""
  ), tmp)
  on.exit(unlink(tmp))

  result <- read_notung_parsable(tmp, mode = "D")
  expect_equal(nrow(result), 1)
})

# --- read_generax_nhx: unbalanced parentheses fix ---

test_that("read_generax_nhx handles NHX with extra closing paren", {
  skip_if_not(requireNamespace("treeio", quietly = TRUE), "treeio not installed")

  # Simulate an NHX where there's one extra ) — the code strips ); to ;
  nhx_text <- "((A:0.1[&&NHX:S=sp1:D=N],B:0.2[&&NHX:S=sp2:D=N])n1:0.3[&&NHX:S=root:D=Y]);"
  nhx_file <- tempfile(fileext = ".nhx")
  writeLines(nhx_text, nhx_file)
  on.exit(unlink(nhx_file))

  old_wd <- getwd()
  setwd(tempdir())
  on.exit(setwd(old_wd), add = TRUE)

  result <- read_generax_nhx(nhx_file)
  expect_true(inherits(result, "treedata"))
  expect_equal(length(result@phylo$tip.label), 2)
})
