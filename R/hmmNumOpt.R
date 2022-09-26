hmmNumOpt <- function(Dat,par0,stationary,verbose,itmax,bicm,
                      rhovals,npar,optimiser,optimMethod,hessian=FALSE,...) {

# Check whether matrices of auxilliary predictors are being used.
npred <- length(attr(Dat,"prednames"))
Xused <- npred > 1

# Check on the "stationary" argument.
if(Xused & !stationary) {
    stop("When auxiliary predictors are use the model must be stationary.\n")
}

# Pick out the useAnalGrad argument, if it's in "...".
useAnalGrad <- list(...)$useAnalGrad
if(is.null(useAnalGrad)) useAnalGrad <- FALSE
if(useAnalGrad & Xused) {
    if(verbose) {
        cat("Cannot currently use an analytic gradient when\n")
        cat("auxilliary predictors (\"X\") are used.\n")
        cat("Re-setting \"useAnalGrad\" to FALSE.\n")
    }
    useAnalGrad <- FALSE
}

pars0 <- reparam(par0,stationary=stationary)
K     <- ncol(par0$tpm)

# Check on the length of pars0.  In particular note that the value
# of pars0 could have come from fitting a model that did not use
# auxiliary predictors.
nval  <- length(attr(Dat,"lvls"))
npp   <- (npred - 1) * K * (nval-1)
np0   <- K*(K-1) + K*(nval-1)
if(length(pars0) == np0) {
   pars0 <- c(pars0,rep(0,npp))
} else if(length(pars0) != np0 + npp) {
   stop("Initial parameter vector is of the wrong length.\n")
}

.cache4nlmPars     <- new.env()
assign("kount",0,envir=.cache4nlmPars)

if(useAnalGrad) {
    objFun <- function(pars,Dat,K) {
        prev.pars <- .cache4nlmPars$last.pars
        .cache4nlmPars$prev.pars <- prev.pars
        .cache4nlmPars$last.pars <- pars
        xxx <- try(get.gl(pars,K,Dat),silent=TRUE)
        if(inherits(xxx,"try-error")) {
            stop("Log likelihood problem.\n")
        }
        vvv <- -xxx$ll
        attr(vvv,"gradient") <- -xxx$grad
        vvv
    }
    gradFun <- function(pars,Dat,K) {
        -get.gl(pars,K,Dat)$grad
    }
} else {
    if(Xused) {
        objFun <- function(pars,Dat,K) {
            prev.pars <- .cache4nlmPars$last.pars
            .cache4nlmPars$prev.pars <- prev.pars
            .cache4nlmPars$last.pars <- pars
            tpm <- getTpm(pars,K,stationary=TRUE)
            Rho <- getRho(pars,K,rhovals=attr(Dat,"lvls"),
                          stationary=TRUE,prednames=attr(Dat,"prednames"))
            xxx <- try(logLikHmm(y=Dat,tpm=tpm,Rho=Rho))
            if(inherits(xxx,"try-error")) {
                stop("Log likelihood problem.\n")
            }
            -xxx
        }
    } else {
        objFun <- function(pars,Dat,K) {
        kount <- .cache4nlmPars$kount+1
cat("count=",kount,"pars=",pars,"\n")
.cache4nlmPars$kount <- kount
            prev.pars <- .cache4nlmPars$last.pars
            .cache4nlmPars$prev.pars <- prev.pars
            .cache4nlmPars$last.pars <- pars
            xxx <- try(get.l(pars,K,Dat),silent=TRUE)
            #xxx <- try(get.gl(pars,K,Dat)$ll,silent=TRUE)
            if(inherits(xxx,"try-error")) {
                stop("Log likelihood problem.\n")
            }
            -xxx
        }
    }
    gradFun <- NULL
}

if(optimiser=="nlm") {
    useOptim <- FALSE
    PL <- if(verbose) 2 else 0
    argh <- list(f=objFun,p=pars0,print.level=PL,iterlim=itmax,
                 Dat=Dat,K=K,hessian=hessian)
    dotz <- list(...)
    dotz$useAnalGrad <- NULL
    dotz$keep.x <- NULL
    dotz$keep.y <- NULL
    argh <- c(argh,dotz)
    mrslt <- try(do.call(nlm,argh),silent=TRUE)
    if(inherits(mrslt,"try-error")) {
        if(verbose) {
            cat("The nlm() optimiser has crashed.\n")
            cat("Trying the optim() optimiser instead.\n")
        }
        pars0 <- .cache4nlmPars$prev.pars
        if(is.null(pars0)) stop("Starting values not working!\n")
        useOptim <- TRUE
    } else {
        ll    <- -mrslt$minimum
        pars  <- mrslt$estimate
        converged <- mrslt$code==1
        attr(converged,"code") <- mrslt$code
        nstep <- mrslt$iterations
    }
} else useOptim <- TRUE

if(useOptim) {
    trace <- if(verbose) 6 else 0
    method <- optimMethod
    if(is.null(method)) method <- "BFGS"
    mrslt <- optim(pars0,objFun,gr=gradFun,method=method,hessian=hessian,
                   control=list(trace=trace,maxit=itmax),Dat=Dat,K=K)
    ll <- -mrslt$value
    pars <- mrslt$par
    converged <- mrslt$convergence==0
    attr(converged,"convergence") <- mrslt$convergence
    nstep <- sum(mrslt$counts)
}
Hess  <- mrslt$hessian
grad  <- mrslt$gradient
tpm   <- getTpm(pars,K,stationary)
Rho   <- getRho(pars,K,rhovals,stationary,prednames=attr(Dat,"prednames"))
ispd  <- if(stationary) revise.ispd(tpm) else getIspd(pars,K)
AIC   <- -2*ll + 2*npar
BIC   <- -2*ll + bicm*npar
rslt  <- list(Rho=Rho,tpm=tpm,ispd=ispd,log.like=ll,par0=par0,npar=npar,
              bicm=bicm,converged=converged,nstep=nstep,prior.emsteps=0,
              AIC=AIC,BIC=BIC)
if(!is.null(grad)) {
    naft <- grep("log.like",names(rslt))
    names(grad) <- names(pars0)
    rslt <- append(rslt,list(grad=grad),after=naft)
    naft <- naft+1
} else naft <- NULL

if(!is.null(Hess)) {
    if(is.null(naft)) naft <- grep("log.like",names(rslt))
    dimnames(Hess) <- list(names(pars0),names(pars0))
    rslt <- append(rslt,list(hessian=Hess),after=naft)
}
rslt
}
