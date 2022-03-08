
# starting message
.onAttach <- function(...) packageStartupMessage("Welcome to Weighted Quantile Sum (WQS) Regression.\nIf you are using a Mac you have to install XQuartz.\nYou can download it from: https://www.xquartz.org/\n")


# function to remove terms from formula
remove_terms <- function(form, term) {
  fterms <- terms(form)
  fac <- attr(fterms, "factors")
  if(term %in% rownames(fac)){
    fac_wqs <- fac[grep(term, rownames(fac)), ]
    if(NCOL(fac_wqs) == 1) idx <- which(as.logical(fac[term, ]))
    else idx <- which(apply(fac_wqs, 2, function(i) any(i==1)))
    new_fterms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
    return(formula(new_fterms))
  }
  else return(form)
}


# function to manage NAs
na_action <- function(data, na.action){
  dtf <- match.call(expand.dots = FALSE)
  m <- match(c("na.action", "data"), names(dtf), 0)
  dtf <- dtf[m]
  names(dtf)[2] <- "object"
  dtf <- eval(dtf, parent.frame())
  return(dtf)
}


# function to create stratified elements in the mixture by levels of factors
stratified_f = function(Q, dtf, stratified, mix_name){
  ls = levels(unlist(dtf[, stratified, drop = FALSE]))
  if(is.null(ls)) stop("'stratified' must be factor\n")
  llsm = lapply(ls, function(ls){
    mat = diag(as.numeric(dtf[, stratified] == ls))
    sub_dt = mat%*%Q[, mix_name]
    colnames(sub_dt) = paste(mix_name, ls, sep = "_")
    return(sub_dt)
  })
  Q = do.call("cbind", llsm)
  mix_name = colnames(Q)
  strtfd_out = list(Q, mix_name)
  names(strtfd_out) = c("Q", "mix_name")

  return(strtfd_out)
}


# function to create variables with quantile of the components
# quantile_f <- function(dtf, mix_name, q){
#
#   if(!is.numeric(q)) stop("'q' must be a number\n")
#   Ql <- lapply(1:length(mix_name), function(i){
#     q_i <- unique(quantile(dtf[[mix_name[i]]], probs = seq(0, 1, by = 1/q), na.rm = TRUE))
#     if(length(q_i) == 1) q_i = c(-Inf, q_i)
#     else{
#       q_i[1] <- -Inf
#       q_i[length(q_i)] <- Inf
#     }
#     q <- cut(dtf[[mix_name[i]]], breaks = q_i, labels = FALSE, include.lowest = TRUE) - 1
#     return(list(q_i, q))
#   })
#   q_i <- lapply(Ql, function(x) x[[1]])
#   Q <- matrix(unlist(lapply(Ql, function(x) x[[2]])), ncol = length(Ql))
#   colnames(Q) <- names(q_i) <- mix_name
#
#   qf_out <- list(Q, q_i)
#   names(qf_out) <- c("Q", "q_i")
#   return(qf_out)
# }


# function to split the dataset
create_rindex <- function(dtf, N, validation, valid_var, m, family){

  if(!m){
    if(!is.numeric(validation) | validation < 0 | validation >= 1) stop("'validation' must be numeric >= 0 and < 1\n")
    # if(!is.numeric(pred) | pred < 0 | pred >= 1) stop("'pred' must be numeric >= 0 and < 1")
    # if(pred + validation >= 1) stop("the sum of pred and validation must be between 0 and 1")
    # if(pred>0 & family$family == "multinomial") stop("The predictive model is not available for multinomial regression")
    groups <- rep(0, N)
    if(validation > 0) groups[sample(1:N, round(N*validation))] <-1
    # if(pred > 0) groups[sample(which(groups!=1), round(N*pred))] <-2
  }
  else{
    groups = unlist(dtf[, valid_var, drop = FALSE])
    if(!any(unique(groups) %in% c(0, 1))) stop("valid_var values must be 0 or 1\n")
    if(!(0 %in% unique(groups))) stop(("0 must identify test dataset\n"))
  }
  it = which(groups == 0)
  iv = which(groups == 1)
  # ip = which(groups == 2)
  if(length(iv) == 0) iv = it
  # if(length(ip) == 0) ip = NULL

  indexl = list(it, iv)
  names(indexl) = c("it", "iv")
  # indexl = list(it, iv, ip)
  # names(indexl) = c("it", "iv", "ip")
  return(indexl)
}


# parameter names in model matrix
parnames <- function(df, formula, form2){
  if(!is.null(form2)){
    mf <- model.frame(form2, df)
    Y <- model.response(mf, "any")
    df$yz <- ifelse(Y == 0, 0, 1)
  }
  mm <- model.matrix(formula, df)
  colnames(mm)
}


# functtion to define the objective function
# objfn <- function(initp, kw, bdtf, Y, offset, Q, kx, Xnames, n_levels, level_names, wqsvars, family, zilink, zero_infl, formula, ff, wghts, stratified, b1_pos, b1_constr){
#
#   if (family$family == "multinomial"){
#     if(b1_constr){
#       par_pos <- which(grepl("wqs", names(initp)))
#       initp[par_pos] <- sapply(1:length(b1_pos), function(i) ifelse(b1_pos[i], abs(initp[par_pos[i]]), -abs(initp[par_pos[i]])))
#     }
#     # w <- matrix(exp(initp[(kx + 1):length(initp)]), kw, n_levels-1)
#     # w <- matrix(1/(1+exp(-(initp[(kx + 1):length(initp)]))), kw, n_levels-1)
#     w <- matrix(initp[(kx + 1):length(initp)]^2, kw, n_levels-1)
#     # w <- matrix(((initp[(kx + 1):(kx + kw)]-min(initp[(kx + 1):(kx + kw)]))/(max(initp[(kx + 1):(kx + kw)])-min(initp[(kx + 1):(kx + kw)]))), kw, n_levels-1)
#     # w <- matrix(abs(initp[(kx + 1):length(initp)]), kw, n_levels-1)
#     # w <- matrix(sqrt(abs(initp[(kx + 1):length(initp)])), kw, n_levels-1)
#     w <- apply(w, MARGIN = 2, FUN = function(i) i/sum(i))
#     bdtf[, wqsvars] <- Q%*%w
#     # fm_l = sapply(wqsvars, function(i) as.formula(gsub("wqs", i, format(formula))))
#     # Xl = lapply(fm_l, function(i) model.matrix(i, data = bdtf))
#     Xl = lapply(wqsvars, function(i){
#       fmi <- as.formula(gsub("wqs", i, format(formula)))
#       model.matrix(fmi, data = bdtf)
#     })
#     X = do.call("cbind", Xl)
#     b_covs = matrix(0, kx, n_levels-1)
#     i = 1:(kx/(n_levels-1))
#     for (j in 0:(n_levels-2)){
#       b_covs[(kx/(n_levels-1))*j+i, j+1] = initp[(kx/(n_levels-1))*j+i]
#     }
#     term = X%*%b_covs + offset
#   }
#   else{
#     if(b1_constr) initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), -abs(initp["wqs"]))
#     # w <- 1/(1+exp(-(initp[(kx + 1):(kx + kw)])))
#     w <- initp[(kx + 1):(kx + kw)]^2
#     # w <- ((initp[(kx + 1):(kx + kw)]-min(initp[(kx + 1):(kx + kw)]))/(max(initp[(kx + 1):(kx + kw)])-min(initp[(kx + 1):(kx + kw)])))
#     # w <- abs(initp[(kx + 1):(kx + kw)])
#     # w <- sqrt(abs(initp[(kx + 1):(kx + kw)]))
#     # w <- exp(initp[(kx + 1):(kx + kw)])
#     w <- w/sum(w)
#     bdtf$wqs <- as.numeric(Q%*%w)
#     X <- model.matrix(formula, bdtf)
#     b_covs <- initp[1:kx]
#     if(family$family == "negbin") theta <- exp(initp[length(initp)])
#     term <- as.numeric(X%*%b_covs) + offset
#   }
#
#   # pen <- l*ifelse(b1_pos, ifelse(initp["wqs"] >= 0, 0, abs(initp["wqs"])), ifelse(initp["wqs"] <= 0, 0, abs(initp["wqs"])))
#   if(family$family == "multinomial") f = -sum((diag(Y%*%t(term)) - log(1 + rowSums(exp(term))))*wghts)
#   else if(family$family == "negbin") f = -sum((suppressWarnings(dnbinom(Y, size = theta, mu = exp(term), log = TRUE)))*wghts)
#   else if(family$family %in% c("poisson", "quasipoisson")) f = -sum(dpois(Y, lambda = exp(term), log = TRUE))
#   else f = sum(family$dev.resids(y = Y, mu = family$linkinv(term), wt = wghts))
#
#   return(f)
# }


# function to define the equality constraint
# linconst = function(initp, kw, bdtf, Q, kx, Xnames, n_levels, level_names, wqsvars, family, zilink, zero_infl, formula, ff, wghts, stratified, b1_pos, b1_constr){
#
#   if (family$family == "multinomial"){
#     w = matrix(initp[(kx + 1):length(initp)], kw, (n_levels-1))
#     wsum = colSums(w)
#   }
#   else{
#     if(is.null(stratified)) wsum = sum(initp[(kx + 1):(kx + kw)]^2)
#     # if(is.null(stratified)) wsum = sum(1/(1+exp(-(initp[(kx + 1):(kx + kw)]))))
#     else wsum = sapply(0:(n_levels-1), function(i) sum(initp[((kx + kw*i/n_levels) + 1):(kx + kw*i/n_levels + kw/n_levels)]))
#   }
#
#   return(wsum)
# }


# function to determine bounded prameters
# bounded_param = function(initp, kw, bdtf, Q, kx, Xnames, n_levels, level_names, wqsvars, family, zilink, zero_infl, formula, ff, wghts, stratified, b1_pos, b1_constr){
#
#   bp = initp[(kx + 1):(kx + kw*ifelse(family$family == "multinomial", (n_levels-1), 1))]
#   # if (b1_constr){
#   #   wqs_site = which(Xnames == "wqs")
#   #   bp = c(initp[wqs_site], bp)
#   # }
#
#   return(bp)
# }


# function to define the lower bounds
# LBound = function(kw, n_levels, family, b1_pos, b1_constr){
#
#   LB = rep(0, kw)
#   if(family$family == "multinomial") LB = rep(LB, n_levels-1)
#   if(b1_constr) LB = c(sapply(b1_pos, function(i) ifelse(i, 0, -Inf)), LB)
#
#   return(LB)
# }


# function to define the upper bounds
# UBound = function(kw, n_levels, family, b1_pos, b1_constr){
#
#   UB = rep(1, kw)
#   if(family$family == "multinomial") UB = rep(UB, n_levels-1)
#   if(b1_constr) UB = c(sapply(b1_pos, function(i) ifelse(i, Inf, 0)), UB)
#
#   return(UB)
# }


# function to define the parameters initial values
values.0 = function(kw, Xnames, kx, n_levels, formula, ff, wghts, bdtf, stratified, stratlev, b1_pos, family, wp, wn, dwqs, zilink, zero_infl){

  # w = rep(1/kw*ifelse(is.null(stratified), 1, n_levels), kw)
  # w = rep(1, kw)
  w = rep(1/kw*stratlev, kw)
  w <- sqrt(w)

  if(family$family == "multinomial"){
    fit = multinom(formula, bdtf, trace = F, weights = wghts)
    bj = c(sapply(1:(n_levels-1), function(i) coef(fit)[i,]))
    names(bj) <- Xnames
    # if(dwqs) val.0 = c(bj, wp, wn)
    if(dwqs) val.0 = c(bj, sqrt(wp), sqrt(wn))
    else{
      w = rep(w, (n_levels-1))
      val.0 = c(bj, w)
    }
  }
  else{
    if(family$family == "negbin") fit = suppressWarnings(glm.nb(formula, bdtf, weights = wghts))
    else fit = glm(formula, bdtf, family = family, weights = wghts)
    bj = coef(fit)
    # if(dwqs) val.0 = c(bj, wp, wn)
    if(dwqs) val.0 = c(bj, sqrt(wp), sqrt(wn))
    else val.0 = c(bj, w)
    if(family$family == "negbin"){
      if(length(attr(terms(formula), "term.labels")) == 1) val.0 = c(val.0, 1)
      else val.0 = c(val.0, log(fit$theta))
    }
  }

  if(dwqs){
    pwqs_site <- which(grepl("pwqs", Xnames))
    nwqs_site <- which(grepl("nwqs", Xnames))
    val.0[pwqs_site] = 0.0001
    val.0[nwqs_site] = -0.0001
  }
  else{
    wqs_site <- which(grepl("wqs", Xnames))
    val.0[wqs_site] = sapply(b1_pos, function(i) ifelse(i, 0.0001, -0.0001))
  }

  return(val.0)
}


# optimization function to estimate the weights
optim.f <- function(i, objfn, Y, Xm, Q, offset, wghts, initp, n_levels, level_names, wqsvars, pwqsvars, nwqsvars,
                    b1_pos, b1_constr, n_vars, dwqs, family, rs, zilink, zero_infl, formula, ff, kx, kw, Xnames,
                    stratified, b, optim.method, control, lambda, stratlev){

  if(rs){
    bindex <- 1:nrow(Xm)
    slctd_vars <- sample(colnames(Q), n_vars, replace=FALSE)
    if(dwqs) initp <- initp[c(rep(T, ncol(Xm)), rep(if(family$family == "multinomial") rep(colnames(Q) %in% slctd_vars, (n_levels-1)) else colnames(Q) %in% slctd_vars, 2))]
    else initp <- initp[c(rep(T, ncol(Xm)), if(family$family == "multinomial") rep(colnames(Q) %in% slctd_vars, (n_levels-1)) else colnames(Q) %in% slctd_vars)]
    kw <- length(slctd_vars)
  }
  else{
    if(b == 1) bindex <- 1:nrow(Xm)
    else bindex <- sample(1:nrow(Xm), nrow(Xm), replace=TRUE)
    slctd_vars <- colnames(Q)
  }
  bXm <- Xm[bindex,]
  bQ <- Q[bindex, slctd_vars]
  bY <- if(family$family == "multinomial")  Y[bindex,] else Y[bindex]
  bwghts <- wghts[bindex]
  boffset <- offset[bindex]
  # Xnames <- parnames(bdtf, formula, NULL)
  # kx <- length(Xnames)
  # if(family$family == "multinomial"){
  #   n_levels <- nlevels(eval(formula[[2]], envir = bdtf))
  #   if(n_levels == 0) stop("y must be of class factor when 'family = \"multinomial\"'\n")
  #   level_names <- levels(eval(formula[[2]], envir = bdtf))
  #   Xnames <- c(sapply(level_names[-1], function(i) paste0(Xnames[1:kx], "_", i, "_vs_", level_names[1])))
  #   kx <- kx*(n_levels-1)
  #   if(dwqs){
  #     pwqsvars = Xnames[grepl("^pwqs_", Xnames)]
  #     nwqsvars = Xnames[grepl("^nwqs_", Xnames)]
  #     bdtf[, c(pwqsvars, nwqsvars)] <- 0
  #     bdtf[, pwqsvars] <- bQ%*%wp
  #     bdtf[, nwqsvars] <- bQ%*%wn
  #     wp <- as.vector(wp)
  #     wn <- as.vector(wn)
  #     wqsvars <- NULL
  #   }
  #   else{
  #     wqsvars = Xnames[grepl("^wqs_", Xnames)]
  #     bdtf[, wqsvars] <- 0
  #     pwqsvars <- nwqsvars <- NULL
  #   }
  # }
  # else {
  #   n_levels <- 2
  #   level_names <- wqsvars <- pwqsvars <- nwqsvars <- NULL
  #   if(dwqs){
  #     bdtf$pwqs <- as.numeric(bQ%*%wp)
  #     bdtf$nwqs <- as.numeric(bQ%*%wn)
  #   }
  # }
  # kw <- dim(bQ)[2]
  # initp <- values.0(kw, Xnames, kx, n_levels, formula, ff, wghts, bdtf, stratified, stratlev, b1_pos, family, wp, wn, dwqs, zilink, zero_infl)
  # mf <- model.frame(formula, bdtf)
  # Y <- model.response(mf, "any")
  # if(family$family == "binomial" & any(class(Y) %in% c("factor", "character"))){
  #   if(class(Y) == "character") Y = factor(Y)
  #   Y <- as.numeric(Y != levels(Y)[1])
  # }
  # if(family$family == "multinomial") Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == level_names[i], 1, 0)))
  # offset <- model.offset(mf)
  # if(is.null(offset)) offset <- 0
  #
  # if(family$family == "multinomial"){
  #   Xl = lapply(wqsvars, function(i){
  #     fmi <- as.formula(gsub("wqs", i, format(formula)))
  #     model.matrix(fmi, data = bdtf)
  #   })
  #   X = do.call("cbind", Xl)
  # }
  # else X <- model.matrix(formula, bdtf)
  opt_res <- optim(par = initp, fn = objfn, method = optim.method, control = control,
                            kw = kw, bXm = bXm, bY = bY, boffset = boffset, bQ = bQ, kx = kx,
                            Xnames = Xnames, n_levels = n_levels, level_names = level_names,
                            wqsvars = wqsvars, pwqsvars = pwqsvars, nwqsvars = nwqsvars, family = family, dwqs = dwqs,
                            zilink = zilink, zero_infl = zero_infl, formula = formula, ff = ff, bwghts = bwghts,
                            stratified = stratified, stratlev = stratlev, b1_pos = b1_pos, b1_constr = b1_constr,
                            lambda = lambda)
  # opt_res <- tryCatch(optim(par = initp, fn = objfn, method = optim.method, control = control,
  #                           kw = kw, X = X, Y = Y, offset = offset, Q = bQ, kx = kx,
  #                           Xnames = Xnames, n_levels = n_levels, level_names = level_names,
  #                           wqsvars = wqsvars, pwqsvars = pwqsvars, nwqsvars = nwqsvars, family = family, dwqs = dwqs,
  #                           zilink = zilink, zero_infl = zero_infl, formula = formula, ff = ff, wghts = wghts,
  #                           stratified = stratified, stratlev = stratlev, b1_pos = b1_pos, b1_constr = b1_constr,
  #                           lambda = lambda), error = function(e) NULL)

  # if(!is.null(opt_res)) {
  #   if(dwqs){
  #     if(family$family == "multinomial") par_opt <- rbind(apply(matrix(abs(opt_res$par[(kx + 1):(kx + kw*(n_levels-1))]), kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)),
  #                                                         apply(matrix(abs(opt_res$par[(kx + kw*(n_levels-1) + 1):(kx + 2*(n_levels-1)*kw)]), kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)))
  #     else par_opt <- c(abs(opt_res$par[(kx + 1):(kx + kw)])/sum(abs(opt_res$par[(kx + 1):(kx + kw)])),
  #                       abs(opt_res$par[(kx + kw + 1):(kx + 2*kw)])/sum(abs(opt_res$par[(kx + kw + 1):(kx + 2*kw)])))
  #   }
  #   else{
  #     if(family$family == "multinomial") par_opt <- apply(matrix(abs(opt_res$par[(kx + 1):length(initp)]), kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i))
  #     else{
  #       if(!is.null(stratified)){
  #         par_opt <- lapply(1:stratlev, function(i){
  #           ncolQ <- ncol(bQ)/stratlev
  #           tmp <- abs(opt_res$par[(kx + 1 + (i-1)*ncolQ):(kx + kw/stratlev*i)])/sum(abs(opt_res$par[(kx + 1 + (i-1)*ncolQ):(kx + kw/stratlev*i)]))
  #           return(tmp)
  #         })
  #         par_opt <- do.call("c", par_opt)
  #       }
  #       else par_opt <- abs(opt_res$par[(kx + 1):(kx + kw)])/sum(abs(opt_res$par[(kx + 1):(kx + kw)]))
  #     }
  #   }
  #   conv <- opt_res$convergence
  #   counts <- opt_res$counts
  #   val <- opt_res$val
  #   mex <- opt_res$message
  # }
  # else{
  #   if(dwqs){
  #     if(family$family == "multinomial") par_opt <- rbind(apply(matrix(abs(initp[(kx + 1):(kx + kw*(n_levels-1))]), kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)),
  #                                                         apply(matrix(abs(initp[(kx + kw*(n_levels-1) + 1):(kx + 2*(n_levels-1)*kw)]), kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)))
  #     else par_opt <- c(abs(initp[(kx + 1):(kx + kw)]^2/sum(initp[(kx + 1):(kx + kw)])),
  #                       abs(initp[(kx + kw + 1):(kx + 2*kw)]^2/sum(initp[(kx + kw + 1):(kx + 2*kw)])))
  #   }
  #   else{
  #     if(family$family == "multinomial") par_opt <- apply(matrix(abs(initp[(kx + 1):length(initp)]), kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i))
  #     else par_opt <- abs(initp[(kx + 1):(kx + kw)])/sum(abs(initp[(kx + 1):(kx + kw)]))
  #   }
  #   conv <- 1
  #   counts <- val <- mex <- NULL
  # }
  if(!is.null(opt_res)) {
    if(dwqs){
      if(family$family == "multinomial") par_opt <- rbind(apply(matrix(opt_res$par[(kx + 1):(kx + kw*(n_levels-1))]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)),
                                                          apply(matrix(opt_res$par[(kx + kw*(n_levels-1) + 1):(kx + 2*(n_levels-1)*kw)]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)))
      else par_opt <- c(opt_res$par[(kx + 1):(kx + kw)]^2/sum(opt_res$par[(kx + 1):(kx + kw)]^2),
                        opt_res$par[(kx + kw + 1):(kx + 2*kw)]^2/sum(opt_res$par[(kx + kw + 1):(kx + 2*kw)]^2))
    }
    else{
      if(family$family == "multinomial") par_opt <- apply(matrix(opt_res$par[(kx + 1):length(initp)]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i))
      else{
        if(!is.null(stratified)){
          par_opt <- lapply(1:stratlev, function(i){
            ncolQ <- ncol(bQ)/stratlev
            tmp <- (opt_res$par[(kx + 1 + (i-1)*ncolQ):(kx + kw/stratlev*i)]^2)/sum(opt_res$par[(kx + 1 + (i-1)*ncolQ):(kx + kw/stratlev*i)]^2)
            return(tmp)
          })
          par_opt <- do.call("c", par_opt)
        }
        else par_opt <- (opt_res$par[(kx + 1):(kx + kw)]^2)/sum(opt_res$par[(kx + 1):(kx + kw)]^2)
      }
    }
    conv <- opt_res$convergence
    counts <- opt_res$counts
    val <- opt_res$val
    mex <- opt_res$message
  }
  else{
    if(dwqs){
      if(family$family == "multinomial") par_opt <- rbind(apply(matrix(initp[(kx + 1):(kx + kw*(n_levels-1))]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)),
                                                          apply(matrix(initp[(kx + kw*(n_levels-1) + 1):(kx + 2*(n_levels-1)*kw)]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)))
      else par_opt <- c(initp[(kx + 1):(kx + kw)]^2/sum(initp[(kx + 1):(kx + kw)]^2),
                        initp[(kx + kw + 1):(kx + 2*kw)]^2/sum(initp[(kx + kw + 1):(kx + 2*kw)]^2))
    }
    else{
      if(family$family == "multinomial") par_opt <- apply(matrix(initp[(kx + 1):length(initp)]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i))
      else par_opt <- (initp[(kx + 1):(kx + kw)]^2)/sum(initp[(kx + 1):(kx + kw)]^2)
    }
    conv <- 1
    counts <- val <- mex <- NULL
  }
  if(any(is.infinite(par_opt))){
    if(dwqs){
      if(family$family == "multinomial"){
        par_opt[which(is.infinite(par_opt[1:kw,]), arr.ind = T)] <- 1/colSums(is.infinite(par_opt[1:kw,])[,which(is.infinite(par_opt[1:kw,]), arr.ind = T)[,2]])
        par_opt[which(!is.infinite(par_opt[1:kw,]), arr.ind = T)] <- 0
        par_opt[(kw+1):2*kw,][which(is.infinite(par_opt[(kw+1):2*kw,]), arr.ind = T)] <- 1/colSums(is.infinite(par_opt[(kw+1):2*kw,])[,which(is.infinite(par_opt[(kw+1):2*kw,]), arr.ind = T)[,2]])
        par_opt[(kw+1):2*kw,][which(!is.infinite(par_opt[(kw+1):2*kw,]), arr.ind = T)] <- 0
      }
      else{
        par_opt[is.infinite(par_opt[1:kw])] <- 1/sum(is.infinite(par_opt[1:kw]))
        par_opt[!(is.infinite(par_opt[1:kw]))] <- 0
        par_opt[(kw+1):2*kw][is.infinite(par_opt[(kw+1):2*kw])] <- 1/sum(is.infinite(par_opt[(kw+1):2*kw]))
        par_opt[(kw+1):2*kw][!(is.infinite(par_opt[(kw+1):2*kw]))] <- 0
      }
    }
    else{
      if(family$family == "multinomial"){
        par_opt[which(is.infinite(par_opt), arr.ind = T)] <- 1/colSums(is.infinite(par_opt)[,which(is.infinite(par_opt), arr.ind = T)[,2]])
        par_opt[which(!is.infinite(par_opt), arr.ind = T)] <- 0
      }
      else{
        par_opt[is.infinite(par_opt)] <- 1/sum(is.infinite(par_opt))
        par_opt[!(is.infinite(par_opt))] <- 0
      }
    }
  }
  # mfit <- model.fit(w = par_opt, bdtf = bdtf, bQ = bQ, family = family, dwqs = dwqs, zilink = zilink, formula = formula,
  #                   ff = ff, wghts = wghts, stratified = stratified, b1_pos = b1_pos, zero_infl = zero_infl)
  mfit <- model.fit(w = par_opt, dt = bXm, bQ = bQ, Y = bY, family = family, dwqs = dwqs, zilink = zilink,
                    formula = formula, ff = ff, wghts = bwghts, offset = boffset, initp = initp, Xnames = Xnames,
                    n_levels = n_levels, level_names = level_names, wqsvars = wqsvars, pwqsvars = pwqsvars,
                    nwqsvars = nwqsvars, stratified = stratified, b1_pos = b1_pos, zero_infl = zero_infl, kx = kx, kw = kw)

  out <- list(par_opt = par_opt, conv = conv, counts = counts, val = val, mex = mex, mfit = mfit, bindex = bindex, slctd_vars = slctd_vars)
  return(out)
}


# function that fit the wqs model
# model.fit <- function(w, bdtf, bQ, family, dwqs, zilink, formula, ff, wghts, stratified, b1_pos, zero_infl){
model.fit <- function(w, dt, bQ, Y, family, dwqs, zilink, formula, ff, wghts, offset, initp, Xnames, n_levels, level_names,
                      wqsvars, pwqsvars, nwqsvars, stratified, b1_pos, zero_infl, kx, kw){

  if(is.matrix(dt)){
    dtf <- as.data.frame(dt)
    formula <- ff <- Y ~ 0 + .
    zero_infl <- FALSE
  }
  else{
    if(family$family == "multinomial"){
      if(dwqs){
        tmp <- sapply(1:length(pwqsvars), function(i) gsub("pwqs", pwqsvars[i], format(formula)))
        fm_l <- sapply(1:length(nwqsvars), function(i) as.formula(gsub("nwqs", nwqsvars[i], tmp[i])))
      }
      else fm_l <- sapply(wqsvars, function(i) as.formula(gsub("wqs", i, format(formula))))
      dtfl <- lapply(fm_l, function(i) model.matrix(i, data = dt))
      dtf <- do.call("cbind", dtfl)
    }
    else dtf <- dt
  }

  # kw <- dim(bQ)[2]
  if (family$family == "multinomial"){
    # mf <- model.frame(formula, bdtf)
    # Y <- model.response(mf, "any")
    # n_levels <- nlevels(Y)
    # level_names <- levels(Y)
    # Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == level_names[i], 1, 0)))
    # offset <- model.offset(mf)
    # if(is.null(offset)) offset <- 0
    # Xnames <- parnames(bdtf, formula, NULL)
    # Xnames <- c(sapply(level_names[-1], function(i) paste0(Xnames[1:length(Xnames)], "_", i, "_vs_", level_names[1])))
    # kx <- length(Xnames)
    if(dwqs){
      wp <- w[1:kw, 1:(n_levels-1)]
      wn <- w[(kw+1):(2*kw), 1:(n_levels-1)]
      # wp <- matrix(w[1:kw*(n_levels-1)], kw, n_levels-1)
      # wn <- matrix(w[(kw*(n_levels-1)+1):2*kw*(n_levels-1)], kw, n_levels-1)
      # pwqsvars <- Xnames[grepl("^pwqs_", Xnames)]
      # nwqsvars <- Xnames[grepl("^nwqs_", Xnames)]
      pwqs <- dtf[, pwqsvars] <- bQ%*%wp
      nwqs <- dtf[, nwqsvars] <- bQ%*%wn
      colnames(pwqs) <- paste0("pwqs_", level_names[-1], "_vs_", level_names[1])
      colnames(nwqs) <- paste0("nwqs_", level_names[-1], "_vs_", level_names[1])
      # fm_l <- sapply(1:(n_levels-1), function(i){
      #   tmp <- gsub("pwqs", pwqsvars[i], format(formula))     # CONTROLLARE !!!!!!!!!!!!!!!!!!!!!!!
      #   tmp <- gsub("nwqs", nwqsvars[i], format(formula))
      #   as.formula(tmp)
      # })
      wqs <- NULL
    }
    else{
      w <- matrix(w, kw, n_levels-1)
      # wqsvars <- Xnames[grepl("^wqs_", Xnames)]
      wqs <- dtf[, wqsvars] <- bQ%*%w
      colnames(wqs) <- paste0("wqs_", level_names[-1], "_vs_", level_names[1])
      # fm_l <- sapply(wqsvars, function(i) as.formula(gsub("wqs", i, format(formula))))
      pwqs <- nwqs <- NULL
    }
    # Xl <- lapply(fm_l, function(i) model.matrix(i, data = bdtf))
    # X <- do.call("cbind", Xl)
    # initp <- values.0(kw, Xnames, kx, n_levels, formula, ff, wghts, bdtf, stratified, b1_pos, family, zilink, zero_infl)
    initp <- initp[1:kx]
  }
  else{
    if(dwqs){
      pwqs <- dtf[,"pwqs"] <- as.numeric(bQ%*%w[1:kw])
      nwqs <- dtf[,"nwqs"] <- as.numeric(bQ%*%w[(kw+1):(2*kw)])
      wqs <- NULL
    }
    else{
      wqs <- dtf[,"wqs"] <- as.numeric(bQ%*%w)
      form_terms <- attr(terms.formula(formula, data = dtf), "term.labels")
      form_terms <- gsub("^`|`$", "", form_terms)
      wqsint <- any(grepl("wqs:", form_terms))
      intwqs <- any(grepl(":wqs", form_terms))
      intpos <- ifelse(wqsint, which(grepl("wqs:", form_terms)),
                       ifelse(intwqs, which(grepl(":wqs", form_terms)), NA))
      if(wqsint | intwqs){
        int_vars <- unlist(strsplit(form_terms[intpos], ":"))
        intvar <- int_vars[int_vars != "wqs"]
        # if(is.numeric(dtf[, intvar])) dtf[,"wqs"] <- as.numeric(scale(dtf[,"wqs"]))
      }
      pwqs <- nwqs <- NULL
    }
  }

  if(zero_infl) m_f <- zeroinfl(ff, bdtf, dist = family$family, link = zilink$name)
  else{
    if(family$family == "multinomial") {
      LL <- function(p, dtfm){
        # b_covs = matrix(0, kx, n_levels-1)
        # i = 1:(kx/(n_levels-1))
        # for (j in 0:(n_levels-2)){
        #   b_covs[(kx/(n_levels-1))*j+i, j+1] = p[(kx/(n_levels-1))*j+i]
        # }
        # term = X%*%b_covs + offset
        tmp <- lapply(1:(n_levels-1), function(i) as.matrix(p[((i-1)*kx/(n_levels-1)+1):(i*kx/(n_levels-1))]))
        b_covs <- as.matrix(bdiag(tmp))
        term = dtfm%*%b_covs + offset
        -sum((diag(Y%*%t(term)) - log(1 + rowSums(exp(term))))*wghts)
      }
      nlm_out = nlm(LL, initp, hessian = TRUE, dtfm = as.matrix(dtf))

      Estimate = nlm_out$estimate
      Standard_Error = sqrt(diag(solve(nlm_out$hessian)))
      stat = Estimate/Standard_Error
      p_value = 2*pnorm(-abs(stat))
      coefficients <- data.frame(Estimate = Estimate, Standard_Error = Standard_Error, stat = stat, p_value = p_value)
      rownames(coefficients) = Xnames
      m_f = list(nlm_out, coefficients)
      names(m_f) = c("nlm_out", "coefficients")
    }
    else if(family$family == "negbin") m_f = glm.nb(formula, data = dtf, weights = wghts)
    else m_f = glm(formula, data = dtf, family = family, weights = wghts)
  }

  mf_out = list(wqs = wqs, pwqs = pwqs, nwqs = nwqs, m_f = m_f)
  return(mf_out)
}

# function to sample from data rownames
# sample_f = function(i, dtf){
#
#   bindex = sample(1:nrow(dtf), nrow(dtf), replace=TRUE)
#
#   return(bindex)
# }

# function that call optim.f to estimate parameters for each bootstrap sample
# estimate_param = function(bindex, dtf, Q, b1_pos, b1_constr, family, zilink, zero_infl, formula, ff, wghts, stratified, control){
#
#   param = optim.f(dtf[bindex,], Q[bindex,], b1_pos, b1_constr, family, zilink, zero_infl, formula, ff, wghts[bindex], stratified, control)
#
#   return(param)
# }


# function to be passed to future_lapply to fit the models
# model.fit_f = function(i, param, bindex, dtf, Q, b1_pos, family, zilink, formula, ff, wghts, stratified, zero_infl){
#
#   b_fit <- model.fit(param[[i]]$par_opt, dtf[bindex[[i]],], Q[bindex[[i]],], family, zilink, formula, ff, wghts, stratified, b1_pos, zero_infl)
#
#   return(b_fit)
# }

# function that calls the optimization function and the function to fit the model for each bootstrap sample
# par.modl.est <- function(dtf, Q, formula, ff, wghts, b, b1_pos, b1_constr, family, zilink, zero_infl, plan_strategy, stratified, control){
#
#   if (family$family %in% c("gaussian", "quasipoisson")) ts = "t"
#   else if (family$family %in% c("binomial", "poisson", "multinomial", "negbin")) ts = "z"
#
#   if(!is.numeric(b)) stop("'b' must be a number\n")
#
#   bindex = lapply(X=1:b, FUN = sample_f, dtf = dtf)
#
#   plan(plan_strategy)
#   if(control$trace) cat("start opt\n")
#   param <- lapply(X = bindex, FUN = optim.f, bdtf = dtf, bQ = Q, b1_pos = b1_pos, b1_constr = b1_constr,
#                          family = family, zilink = zilink, zero_infl = zero_infl, formula = formula, ff = ff,
#                          wghts = wghts, stratified = stratified, control = control)
#   # param <- future_lapply(X = bindex, FUN = estimate_param, dtf = dtf, Q = Q, b1_pos = b1_pos, b1_constr = b1_constr,
#   #                        family = family, zilink = zilink, zero_infl = zero_infl, formula = formula, ff = ff,
#   #                        wghts = wghts, stratified=stratified, control = control, lp = lp, ln = ln, future.seed = FALSE)
#
#   # b_fit <- lapply(X = 1:b, FUN = model.fit_f, param = param, bindex = bindex, dtf = dtf, Q = Q, b1_pos = b1_pos,
#   #                        family = family, zilink = zilink, formula = formula, ff = ff, wghts = wghts,
#   #                        stratified = stratified, zero_infl = zero_infl)
#
#   conv <- c(sapply(param, function(i) i$conv))
#   counts <- c(sapply(param, function(i) i$counts))
#   val <- c(sapply(param, function(i) i$val))
#   mex <- lapply(param, function(i) i$mex)
#   if(family$family == "multinomial"){
#     n_levels <- dim(param[[1]]$par_opt)[2]+1
#     wqs_site <- which(grepl("^wqs_", rownames(param[[1]]$mfit$m_f$sum_stat)))
#     wght_matrix <- lapply(1:(n_levels-1), function(j) do.call("rbind", lapply(param, function(i) i$par_opt[,j])))
#     b1 <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$nlm_out$estimate[j]))
#     se <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$sum_stat$Standard_Error[j]))
#     stat <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$sum_stat$stat[j]))
#     p_val <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$sum_stat$p_value[j]))
#   }
#   else{
#     wght_matrix <- do.call("rbind", lapply(param, function(i) i$par_opt))
#     if(zero_infl){
#       b1_count <- sapply(param, function(i) i$mfit$m_f$coefficients$count["wqs"])
#       se_count <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$count["wqs", "Std. Error"])
#       stat_count <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$count["wqs", paste0(ts, " value")])
#       p_val_count <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$count["wqs", gsub("x", ts, "Pr(>|x|)")])
#       if("wqs" %in% names(param[[1]]$mfit$m_f$coefficients$zero)){
#         b1_zero <- sapply(param, function(i) i$mfit$m_f$coefficients$zero["wqs"])
#         se_zero <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$zero["wqs", "Std. Error"])
#         stat_zero <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$zero["wqs", paste0(ts, " value")])
#         p_val_zero <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$zero["wqs", gsub("x", ts, "Pr(>|x|)")])
#       }
#       else b1_zero <- se_zero <- stat_zero <- p_val_zero <- NULL
#     }
#     else{
#       b1 <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", "Estimate"])
#       se <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", "Std. Error"])
#       stat <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", paste0(ts, " value")])
#       p_val <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", gsub("x", ts, "Pr(>|x|)")])
#     }
#     n_levels <- 1
#   }
#
#   n_non_conv = sum(conv == 1)
#   if(n_non_conv == 0 & control$trace) cat(paste0("The optimization function always converged\n"))
#   else if(n_non_conv == b) stop("The optimization function never converged\n")
#   else if(control$trace) cat(paste0("The optimization function did not converge ", n_non_conv, " time/times\n"))
#   if(control$trace) cat(paste0("There are ", ifelse(b1_pos, sum(b1 >= 0, na.rm = T), sum(b1 <= 0, na.rm = T)),
#                                ifelse(b1_pos, " positive", " negative"), " bootstrapped b1 out of ", b, "\n"))
#
#   # estimate mean weight for each component (exclude weights from iterations with failed convergence)
#   if (family$family == "multinomial"){
#     bres <- Map(cbind, wght_matrix, b1, se, stat, p_val)
#     bres <- lapply(bres, as.data.frame)
#     bres <- lapply(bres, setNames, c(colnames(Q), "b1", "Std_Error", "stat", "p_val"))
#     strata_names <- gsub("wqs_", "", rownames(param[[1]]$mfit$m_f$sum_stat)[wqs_site])
#     names(bres) <- strata_names
#   }
#   else {
#     if(zero_infl){
#       if(is.null(b1_zero)){
#         bres <- as.data.frame(cbind(wght_matrix, b1_count, se_count, stat_count, p_val_count))
#         names(bres) <- c(colnames(Q), "b1_count", "Std_Error_count", "stat_count", "p_val_count")
#       }
#       else{
#         bres <- as.data.frame(cbind(wght_matrix, b1_count, se_count, stat_count, p_val_count, b1_zero, se_zero, stat_zero, p_val_zero))
#         names(bres) <- c(colnames(Q), "b1_count", "Std_Error_count", "stat_count", "p_val_count", "b1_zero", "Std_Error_zero", "stat_zero", "p_val_zero")
#       }
#     }
#     else{
#       bres <- as.data.frame(cbind(wght_matrix, b1, se, stat, p_val))
#       names(bres) <- c(colnames(Q), "b1", "Std_Error", "stat", "p_val")
#     }
#     strata_names <- NULL
#   }
#
#   par_model_out <- list(bres = bres, conv = conv, val = val, bindex = bindex, counts = counts, n_levels = n_levels, strata_names = strata_names, mex = mex)
#   return(par_model_out)
# }
#
#
# # function to estimate mean weights for each component
# mean_weight_f = function(mix_name, bres, conv, b1_pos, family, n_levels, strata_names, zero_infl){
#
#   if (family$family == "multinomial"){
#     mean_weight <- lapply(1:(n_levels-1), function(i){
#       if(b1_pos[i]) w_t = apply(bres[[i]][bres[[i]]$b1 > 0 & conv == 0, mix_name], 2, weighted.mean, abs(bres[[i]][bres[[i]]$b1 > 0 & conv == 0, "stat"]))
#       else if(!b1_pos[i]) w_t = apply(bres[[i]][bres[[i]]$b1 < 0 & conv == 0, mix_name], 2, weighted.mean, abs(bres[[i]][bres[[i]]$b1 < 0 & conv == 0, "stat"]))
#       if (all(is.nan(w_t)))
#         stop(paste0("There are no ", ifelse(b1_pos[i], "positive", "negative"), " b1 in the bootstrapped models for ", strata_names[i], "\n"))
#       return(w_t)
#     })
#     mean_weight <- list.cbind(mean_weight)
#   }
#   else{
#     if(zero_infl){
#       if(b1_pos){
#         if(is.null(bres$b1_zero)) mean_weight = apply(bres[bres$b1_count > 0 & conv == 0, mix_name], 2, weighted.mean, abs(bres[bres$b1_count > 0 & conv == 0, "stat_count"]))
#         else mean_weight = apply(bres[bres$b1_count > 0 & bres$b1_zero > 0 & conv == 0, mix_name], 2, weighted.mean, abs(bres[bres$b1_count > 0 & bres$b1_zero > 0 & conv == 0, "stat_count"]))
#       }
#       else{
#         if(is.null(bres$b1_zero)) mean_weight = apply(bres[bres$b1_count < 0 & conv == 0, mix_name], 2, weighted.mean, abs(bres[bres$b1_count < 0 & conv == 0, "stat_count"]))
#         else mean_weight = apply(bres[bres$b1_count < 0 & bres$b1_zero < 0 & conv == 0, mix_name], 2, weighted.mean, abs(bres[bres$b1_count < 0 & bres$b1_zero < 0 & conv == 0, "stat_count"]))
#       }
#     }
#     else{
#       if(b1_pos) mean_weight = apply(bres[bres$b1 > 0 & conv == 0, mix_name], 2, weighted.mean, abs(bres[bres$b1 > 0 & conv == 0, "stat"]))
#       else mean_weight = apply(bres[bres$b1 < 0 & conv == 0, mix_name], 2, weighted.mean, abs(bres[bres$b1 < 0 & conv == 0, "stat"]))
#     }
#     if(all(is.nan(mean_weight)))
#       stop("There are no ", ifelse(b1_pos, "positive", "negative"), " b1 in the bootstrapped models\n")
#   }
#
#   return(mean_weight)
# }


# function to build predictive model in case of binomial dist
# predict_f = function(Q, dtf, w, m_f, formula){
#
#   dtf$wqs = as.numeric(Q%*%w)
#   pred = predict(m_f, newdata = dtf, type = "response")
#   y = model.response(model.frame(formula, dtf), "any")
#   df_pred = data.frame(y = y, p_y = pred)
#
#   return(df_pred)
# }
#
# plots <- function(data_plot, y_adj_wqs_df, q, mix_name, mean_weight, fit, family, n_levels, strata_names, df_roc, zero_infl){
#
#   if(family$family == "multinomial"){
#     data_plot = data_plot[order(data_plot[, strata_names[1]]),]
#     pos = match(data_plot$mix_name, sort(mix_name))
#     data_plot$mix_name = factor(data_plot$mix_name, levels(data_plot$mix_name)[pos])
#     data_plot_l = melt(data_plot, id.vars = "mix_name")
#     bar_plot_h = ggplot(data_plot_l, aes(x = mix_name, y = value, fill = mix_name)) +
#       facet_wrap(~ variable)
#   }
#   else {
#     data_plot <- data_plot[order(data_plot$mean_weight),]
#     data_plot$mix_name <- factor(data_plot$mix_name, levels = data_plot$mix_name)
#     bar_plot_h <- ggplot(data_plot, aes(x = mix_name, y = mean_weight, fill = mix_name))
#   }
#
#   bar_plot_h = bar_plot_h + geom_bar(stat = "identity", color = "black") + theme_bw() +
#     theme(axis.ticks = element_blank(),
#           axis.title = element_blank(),
#           axis.text.x = element_text(color='black'),
#           legend.position = "none") + coord_flip()
#
#   print(bar_plot_h)
#
#   y_labs = ifelse(family$family %in% c("multinomial", "binomial"), "y", "y_adj")
#
#   yadj_vs_wqs = ggplot(y_adj_wqs_df, aes_string("wqs", y_labs)) +
#     geom_point() + stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw()
#
#   if(family$family == "multinomial") yadj_vs_wqs = yadj_vs_wqs + facet_wrap(~ level)
#
#   print(yadj_vs_wqs)
#
#   if(!(family$family %in% c("binomial", "multinomial"))){
#     if(zero_infl) fit_df = data.frame(.fitted = fit$fitted.values, .resid = fit$residuals)
#     else fit_df = augment(fit)
#     res_vs_fitted = ggplot(fit_df, aes_string(x = ".fitted", y = ".resid")) + geom_point() + theme_bw() +
#       xlab("Fitted values") + ylab("Residuals")
#     print(res_vs_fitted)
#   }
#
#   if(n_levels == 3){
#     w1_vs_w2 = ggplot(data_plot, aes_string(names(data_plot)[2], names(data_plot)[3])) + geom_point() +
#       theme_bw() + xlab(names(data_plot)[2]) + ylab(names(data_plot)[3]) + geom_abline(linetype = 2) +
#       geom_text_repel(aes(label=mix_name))
#     print(w1_vs_w2)
#   }
#
#   if(!is.null(df_roc) & family$family == "binomial"){
#     if(class(df_roc$y) == "character") df_roc$y = factor(df_roc$y)
#     if(class(df_roc$y) == "factor") df_roc$y <- as.numeric(df_roc$y != levels(df_roc$y)[1])
#     gg_roc = suppressWarnings(ggplot(df_roc, aes_string(d="y", m="p_y")) + geom_roc(n.cuts = 0) +
#                                 style_roc(xlab = "1 - Specificity", ylab = "Sensitivity"))
#     auc_est = calc_auc(gg_roc)
#     gg_roc = gg_roc + annotate("text", x=0.75, y=0.25, label=paste0("AUC = ", round(auc_est[, "AUC"], 3)))
#
#     print(gg_roc)
#   }
# }

# Function for creating the file containing the tables
# tables <- function(final_weight, mf, family, n_levels, zero_infl){
#
#   if(family$family == "multinomial"){
#     final_weight[, c(2:n_levels)] = signif(final_weight[, c(2:n_levels)], 3)
#     mf_df = signif(mf$sum_stat, 3)
#   }
#   else{
#     final_weight <- data.frame(Mix_name = final_weight$mix_name, Final_weight = signif(final_weight$mean_weight, 3))
#     if(zero_infl){
#       mf_df_count = as.data.frame(signif(coef(summary(mf))$count, 3))
#       mf_df_zero = as.data.frame(signif(coef(summary(mf))$zero, 3))
#       rownames(mf_df_zero) = paste0("z_", rownames(mf_df_zero))
#       mf_df = rbind(mf_df_count, mf_df_zero)
#     }
#     else mf_df = as.data.frame(signif(coef(summary(mf)), 3))
#   }
#
#   if(zero_infl) print(kable_styling(kable(mf_df, row.names = TRUE)) %>%
#     group_rows("Count model", 1, dim(mf_df_count)[1]) %>%
#     group_rows("Zero-inflation model", dim(mf_df_count)[1]+1, dim(mf_df_count)[1]+dim(mf_df_zero)[1]))
#   else print(kable_styling(kable(mf_df, row.names = TRUE)))
#   print(kable_styling(kable(final_weight, row.names = FALSE)))
# }


predictmultinom <- function(object, data, sumtype, type){
  fm_l <- sapply(colnames(object$wqs), function(i) as.formula(gsub("wqs", i, format(object$formula))))
  Xl <- lapply(fm_l, function(i) model.matrix(i, data = data))
  predl <- lapply(1:length(Xl), function(i) exp(Xl[[i]]%*%coef(object, sumtype)[i,]))
  predl <- lapply(predl, function(i) i/(1+Reduce("+", predl)))
  pred <- do.call("cbind", predl)
  pred <- cbind(1-rowSums(pred), pred)
  y <- model.response(model.frame(fm_l[[1]], data), "any")
  if(type == "response"){
    Ylevels <- levels(y)
    pred <- factor(Ylevels[apply(pred, 1, which.max)], levels = Ylevels)
  }
  else if(!(type %in% c("response", "prob"))) stop("If family is \"multinomial\" then predict type must be \"response\" or \"prob\"\n")
  return(list(pred = pred, y = y))
}

set_par_names <- function(i, slctd_vars, param, q_name, family, dwqs){

  temp <- param[[i]]$par_opt
  if(family$family == "multinomial"){
    param[[i]]$par_opt <- matrix(NA, length(q_name), dim(temp)[2])
    param[[i]]$par_opt[which(q_name %in% slctd_vars[[i]]),] <- temp
    rownames(param[[i]]$par_opt) <- q_name
  }
  else{
    param[[i]]$par_opt <- rep(NA, ifelse(dwqs, 2, 1)*length(q_name))
    names(param[[i]]$par_opt) <- rep(q_name, ifelse(dwqs, 2, 1))
    if(dwqs){
      param[[i]]$par_opt[1:length(q_name)][slctd_vars[[i]]] <- temp[1:(length(temp)/2)]
      param[[i]]$par_opt[(length(q_name)+1):(2*length(q_name))][slctd_vars[[i]]] <- temp[(length(temp)/2+1):length(temp)]
    }
    else param[[i]]$par_opt[slctd_vars[[i]]] <- temp
  }

  return(param[[i]])
}

