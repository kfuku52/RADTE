read_generax_nhx = function(generax_file) {
    treetext = readLines(generax_file, warn=FALSE)
    count_matches = function(pattern, text_vec) {
        matches = gregexpr(pattern, text_vec, perl=TRUE)
        return(sum(vapply(matches, function(m) {
            if ((length(m)==1) && (m[[1]]==-1)) {
                return(0L)
            }
            return(as.integer(length(m)))
        }, integer(1))))
    }
    num_open = count_matches('\\(', treetext)
    num_close = count_matches('\\)', treetext)
    if ((num_open - num_close)==-1) {
        cat('Number of parentheses in the .nhx is not consistent. Trying to fix.')
        treetext <- gsub("\\)\\s*;", ";", treetext, perl=TRUE)
    }
    tmp_nhx_file = tempfile(pattern='radte_nhx_', fileext='.txt')
    write(treetext, tmp_nhx_file)
    on.exit(file.remove(tmp_nhx_file), add=TRUE)
    nhxtree = treeio::read.nhx(tmp_nhx_file)
    return(nhxtree)
}
