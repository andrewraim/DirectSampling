# Mimic fprintf function in C
fprintf = function(file, msg, ...) {
	cat(sprintf(msg, ...), file = file)
}

# Print a vector x as a string "c(x1, ..., xn)"
print_vector = function(x) {
	sprintf("c(%s)", paste(x, collapse = ","))
}

# Print a matrix with latex separators
print_matrix_latex = function(x, fmt = "%f")
{
	m = nrow(x)
	n = ncol(x)

	if (!is.null(colnames(x))) {
		cat(paste(colnames(x), collapse = " & "), "\\\\ \n")
	}

	for (i in 1:m) {
		if (!is.null(rownames(x))) {
			cat(rownames(x)[i], "& ")
		}
		cat(paste(sprintf(fmt, x[i,]), collapse = " & "), "\\\\ \n")
	}
}
