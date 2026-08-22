test_files <- list.files(file.path("tests", "r"), pattern = "^test_.*[.]R$", full.names = TRUE)
for (test_file in test_files) {
  source(test_file, local = new.env(parent = globalenv()))
}
