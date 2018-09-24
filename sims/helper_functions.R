QQ_plot = function (p_values, col = rainbow(min(length(p_values), ncol(p_values))), 
    main = "", pch = 20, errors = TRUE, lambda = TRUE, p_thresh = 1e-06, 
    showNames = FALSE, ylim = NULL, xlim = NULL, plot = TRUE, 
    new = TRUE, box.lty = par("lty"), collapse = FALSE, ...) 
{
    if (collapse) {
        p_values = as.vector(unlist(p_values, use.names = FALSE))
    }
    if (!is.list(p_values)) {
        names(p_values) = c()
        keys = colnames(p_values)
        if ("" %in% keys) {
            keys[which("" %in% keys)] = "NULL"
        }
        p_values = as.matrix(p_values)
        p_values_list = list()
        for (i in 1:ncol(p_values)) {
            p_values_list[[i]] = p_values[, i]
        }
        names(p_values_list) = keys
        p_values = p_values_list
        rm(p_values_list)
    }
    rge = range(p_values, na.rm = TRUE)
    if (rge[1] < 0 || rge[2] > 1) {
        stop("p_values outside of range [0,1]")
    }
    if (is.null(names(p_values))) {
        names(p_values) = 1:length(p_values)
    }
    if (is.null(pch)) {
        pch = rep(20, length(p_values))
    }
    if (length(pch) == 1) {
        pch = rep(pch, length(p_values))
    }
    p_values = as.list(p_values)
    ry = 0
    rx = 0
    for (key in names(p_values)) {
        p_values[[key]] = p_values[[key]][which(!is.na(p_values[[key]]))]
        j = which(p_values[[key]] == 0)
        if (length(j) > 0) {
            p_values[[key]][j] = min(p_values[[key]][-j])
        }
        ry = max(ry, -log10(min(p_values[[key]])))
        rx = max(rx, -log10(1/(length(p_values[[key]]) + 1)))
    }
    if (!is.null(ylim)) {
        ry = max(ylim)
    }
    if (!is.null(xlim)) {
        rx = max(xlim)
    }
    r = max(rx, ry)
    xlab = expression(-log[10]("expected p-value"))
    ylab = expression(-log[10]("observed p-value"))
    if (plot && new) {
        plot(1, type = "n", las = 1, pch = 20, xlim = c(0, rx), 
            ylim = c(0, ry), xlab = xlab, ylab = ylab, main = main, 
            ...)
        abline(0, 1)
    }
    lambda_values = c()
    se = c()
    i = 1
    for (key in names(p_values)) {
        observed = sort(p_values[[key]], decreasing = TRUE)
        obs = sort(p_values[[key]][which(p_values[[key]] > p_thresh)], 
            decreasing = TRUE)
        lambda_values[i] = qchisq(median(obs), 1, lower.tail = FALSE)/qchisq(0.5, 
            1)
        if (plot) {
            j = which(observed == 0)
            if (length(j) > 0) {
                observed[j] = min(observed[-j])
            }
            expected = (length(observed):1)/(length(observed) + 
                1)
            p = length(expected)
            if (p < 1e+06) {
                intervals = ceiling(c(1, 0.2 * p, 0.4 * p, 0.7 * 
                  p, 0.9 * p, 0.95 * p, 0.99 * p, p))
                scaling = c(28, 200, 800, 1000, 3000, 5000, p)
            }
            else {
                intervals = ceiling(c(1, 0.2 * p, 0.4 * p, 0.7 * 
                  p, 0.9 * p, 0.95 * p, 0.99 * p, 0.995 * p, 
                  0.999 * p, p))
                scaling = c(28, 200, 800, 1000, 3000, 5000, 10000, 
                  50000, p)
            }
            for (j in 1:length(scaling)) {
                k = seq(intervals[j], intervals[j + 1], intervals[j + 
                  1]/scaling[j])
                points(-log10(expected[k]), -log10(observed[k]), 
                  col = col[i], pch = pch[i])
            }
        }
        i = i + 1
    }
    if (plot && lambda) {
        if (!showNames) {
            namesStrings = rep("", length(names(p_values)))
        }
        else {
            namesStrings = paste(names(p_values), ": ", sep = "")
        }
        legend("topleft", legend = paste(namesStrings, format(lambda_values, 
            digits = 4, nsmall = 3)), col = col, pch = 15, pt.cex = 1.5, 
            title = expression(lambda["GC"]), box.lty = box.lty)
    }
    if (plot && errors) {
        error_quantile = 0.95
        M = length(p_values[[key]])
        alpha = seq(M, M/10 + 1, length.out = 1000)
        if (M/10 > 1) 
            alpha = append(alpha, seq(M/10, M/100 + 1, length.out = 1000))
        if (M/100 > 1) 
            alpha = append(alpha, seq(M/100, M/1000 + 1, length.out = 1000))
        if (M/1000 > 1) 
            alpha = append(alpha, seq(M/1000, M/10000 + 1, length.out = 10000))
        alpha = append(alpha, seq(min(alpha), 1, length.out = 10000))
        alpha = round(alpha)
        beta = M - alpha + 1
        x_top = qbeta(error_quantile, alpha, beta)
        x_bot = qbeta(1 - error_quantile, alpha, beta)
        lines(-log10(alpha/(M + 1)), -log10(x_top))
        lines(-log10(alpha/(M + 1)), -log10(x_bot))
    }
    return(invisible(lambda_values))
}






.makepretty = function(x){
    msg = gsub('\n', ' ', x)
    msg = gsub('    ', '', msg)
    msg
}


.check_error_model = function(extras, paired){

    # make sure it's an available model
    error_model = match.arg(extras$error_model,
        c('uniform', 'illumina4', 'illumina5', 'custom'))

    # check uniform model --> error rate
    if(error_model == 'uniform'){
        if('error_rate' %in% names(extras)){
            error_rate = extras$error_rate
            stopifnot(is.numeric(error_rate))
            stopifnot(error_rate >= 0 & error_rate <= 1)
        }
    }

    # check paths and such for custom model
    if(error_model == 'custom'){
        if(!('model_path' %in% names(extras)) |
            !('model_prefix' %in% names(extras))){
            stop(.makepretty('with custom error models, you must provide both
                the path to the folder that holds your error model
                (model_path) and the prefix of your error model (model_prefix),
                where the prefix is whatever comes before _mate1 and _mate2
                (for paired reads) or _single (for single-end reads). You
                provided prefix when running build_error_models.py.'))
        }
        model_path = extras$model_path
        model_prefix = extras$model_prefix
        if(paired){
            if(!file.exists(paste0(model_path, '/', model_prefix, '_mate1')) |
               !file.exists(paste0(model_path, '/', model_prefix, '_mate2'))){
               stop('could not find error model')
            }
        }else{
            if(!file.exists(paste0(model_path, '/', model_prefix, '_single'))){
                stop('could not find error model')
            }
        }
    }
}

.check_fold_changes = function(fold_changes, num_reps, transcripts){

    # make sure fold change matrix is compatible with experiment size
    if(length(num_reps) == 1 | length(num_reps) == 2){
        stopifnot(is.numeric(fold_changes))
    }else{
        stopifnot(is.matrix(fold_changes))
        if(ncol(fold_changes) != length(num_reps)){
            stop(.makepretty('wrong number of columns in fold change matrix:
                need same number of columns as number of groups.'))
        }
        if(nrow(fold_changes) != length(transcripts)){
            stop(.makepretty('wrong number of rows in fold change matrix: need
                same number of rows as number of simulated transcripts. see
                count_transcripts to find out that number.'))
        }
    }
}

.check_extras = function(extras, paired){

    if(!('distr' %in% names(extras))){
        extras$distr = 'normal'
    }else{
        extras$distr = match.arg(extras$distr, 
            c('normal', 'empirical', 'custom'))
        if(extras$distr == 'custom' & !('custdens' %in% names(extras))){
            stop(.makepretty('to use custom fragment distribution, provide
                "custdens", a logspline object representing the distribution.'))
        }
    }
    if(!('fraglen' %in% names(extras))){
        extras$fraglen = 250
    }
    if(!('fragsd' %in% names(extras))){
        extras$fragsd = 25
    }## I don't love this--these arguments aren't needed unless distr is normal.
    # but we store them anyway. should code better?
    if(!('readlen' %in% names(extras))){
        extras$readlen = 100
    }

    if(!('bias' %in% names(extras))){
        extras$bias = 'none'
    }else{
        extras$bias = match.arg(extras$bias, c('none', 'rnaf', 'cdnaf'))
    }

    if(!('error_model' %in% names(extras))){
        extras$error_model = 'uniform'
    }
    .check_error_model(extras, paired)

    if(!('error_rate' %in% names(extras))){
        extras$error_rate = 0.005
    }
    if(extras$error_model == 'custom'){
        extras$path = paste0(extras$model_path, '/', extras$model_prefix)
    }#this should work beause we already checked stuff.

    if(!('bias' %in% names(extras))){
        extras$bias = 'none'
    }else{
        extras$bias = match.arg(extras$bias, c('none', 'rnaf', 'cdnaf'))
    }

    return(extras)

}


simulate_experiment_justCounts = function (transcripts = NULL, gtf = NULL, seqpath = NULL, outdir = ".", 
    num_reps = c(10, 10), reads_per_transcript = 300, size = NULL, 
    fold_changes, paired = TRUE, reportCoverage = F, simReads = FALSE,...) 
{
    extras = list(...)
    extras = .check_extras(extras, paired)
    if (!("lib_sizes" %in% names(extras))) {
        extras$lib_sizes = rep(1, sum(num_reps))
    }
    else {
        stopifnot(is.numeric(extras$lib_sizes))
        stopifnot(length(extras$lib_sizes) == sum(num_reps))
    }
   
    if (length(num_reps) == 1) {
        fold_changes = rep(1, length(transcripts))
    }
    .check_fold_changes(fold_changes, num_reps, transcripts)
    if ("meanmodel" %in% names(extras)) {
        b0 = -3.0158
        b1 = 0.8688
        sigma = 4.152
        logmus = b0 + b1 * log2(width(transcripts)) + rnorm(length(transcripts), 
            0, sigma)
        reads_per_transcript = 2^logmus - 1
    }
    basemeans = ceiling(reads_per_transcript * fold_changes)
    if (is.null(size)) {
        size = basemeans/3
    }
    else if (class(size) == "numeric") {
        size = matrix(size, nrow = nrow(basemeans), ncol = ncol(basemeans))
    }
    else if (class(size) == "matrix") {
        stopifnot(nrow(size) == nrow(basemeans))
        stopifnot(ncol(size) == ncol(basemeans))
    }
    else {
        stop("size must be a number, numeric vector, or matrix.")
    }
    if ("seed" %in% names(extras)) {
        set.seed(extras$seed)
    }
    group_ids = rep(1:length(num_reps), times = num_reps)
    numreadsList = vector("list", sum(num_reps))
    numreadsList = lapply(1:sum(num_reps), function(i) {
        group_id = group_ids[i]
        NB(as.matrix(basemeans)[, group_id], as.matrix(size)[, 
            group_id])
    })
    readmat = matrix(unlist(numreadsList), ncol = sum(num_reps))
    readmat = ceiling(t(extras$lib_sizes * t(readmat)))
    if ("gcbias" %in% names(extras)) {
        stopifnot(length(extras$gcbias) == sum(num_reps))
        gcclasses = unique(sapply(extras$gcbias, class))
        if (sum(gcclasses %in% c("numeric", "loess")) < length(extras$gcbias)) {
            stop(.makepretty("gc bias parameters must be integers 0 through 7\n                or loess objects"))
        }
        if (any(extras$gcbias != 0)) {
            readmat = add_gc_bias(readmat, extras$gcbias, transcripts)
        }
    }
    sysoutdir = gsub(" ", "\\\\ ", outdir)
    if (.Platform$OS.type == "windows") {
        shell(paste("mkdir", sysoutdir))
    }
    else {
        system(paste("mkdir -p", sysoutdir))
    }
    if(simReads){
        sgseq(readmat, transcripts, paired, outdir, extras, reportCoverage)

        if (!("write_info" %in% names(extras))) {
            write_info = TRUE
        }
        if (write_info) {
            .write_info(extras, transcripts, num_reps, fold_changes, 
                outdir, group_ids)
        }
    }
    readmat
}

