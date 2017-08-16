#' @title Print assay summary
#' @method print assay
#' @description \code{print.assay} simply prints a summary of the HTS assay
#' @param x The assay you want to print
#' @param plate The plate number (Default:1)
#' @param ... Ellipsis to be passed to the default \code{print()} function
#' @return None
#' @export
print.assay<-function(x,...,plate=1){
  assay=x
  if(class(assay)!="assay")
    stop("Error: x is not an assay.")
  cat("HTS assay (",dim(assay$m)[1],"rows x ",dim(assay$m)[2],"columns x ",dim(assay$m)[3]," plates):\nPlate ",plate,"\n")
  print(assay$m[,,plate],...) # Print first plate raw measurements
  cat("Minimum value: ",min(assay$m)," Maximum value: ",max(assay$m))
}
#' @title Plot assay plate-wise
#' @method plot assay
#' @description  \code{plot.assay} plots a hit map of the assay (only one plate at a time)
#' @param x The assay you want to plot
#' @param plate The plate number (Default:1)
#' @param type Either "R" - raw assay, "C" - corrected assay (if it exists) or "P" - spatial bias position
#' @param ... Ellipsis to be passed to the \code{levelplotplot()} function from the \code{lattice} package
#' @return None
#' @export
plot.assay<-function(x,...,plate=1,type="R"){
  assay=x
  if(class(assay)!="assay")
    stop("Error: x is not an assay.")
  rotate <- function(x) t(apply(x, 2, rev)) # Helper method to put rows back horizontally and columns vertically
  theme <- latticeExtra::custom.theme(region=rev(brewer.pal(n=11, 'RdBu'))) # Custom red/blue theme
  if(type=="R")
    lattice::levelplot(rotate(assay$m[,,plate]),xlab="Columns",ylab="Rows",main=paste("Plate",plate,"raw measurements"),par.settings=theme,margin=F,pretty=T,scales=list(x=list(alternating=2,tck=c(0,0)),y=list(tck=c(0,0))),...)
  else if(type=="C" & 'mCorrected'%in%names(assay))
    lattice::levelplot(rotate(assay$mCorrected[,,plate]),xlab="Columns",ylab="Rows",main=paste("Plate",plate,"raw measurements"),par.settings=theme,margin=F,pretty=T,scales=list(x=list(alternating=2,tck=c(0,0)),y=list(tck=c(0,0))),...)
  else  if(type=="P")
    lattice::levelplot(rotate(assay$biasPositions[,,plate]),xlab="Columns",ylab="Rows",main=paste("Plate",plate,"raw measurements"),par.settings=theme,margin=F,pretty=T,scales=list(x=list(alternating=2,tck=c(0,0)),y=list(tck=c(0,0))),...)
  else
    stop("Error. Check your assay")
}
#' @title Create a new \code{assay}
#' @description \code{create_assay} makes a new object of class assay. You should pass this object to \code{detect_bias()} and \code{correct_bias()} methods
#' @param m The assay you want to be corrected
#' @param ctrl An optional boolean array of the same dimensions as \code{m}. Each entry is 1 if the well is a control well, 0 otherwise. All control wells are excluded from all computations
#' @return assay The created assay object.
#' It containts the following fields:
#'
#' \code{n} The HTS matrix of raw measurements
#'
#' \code{ctrl} The binary matrix of control wells
#'
#' \code{biasPositions} The binary matrix where 1:well is biased, 0:well is unbiased, as suggested by Mann-Whitney test
#'
#' \code{mCorrected} The HTS matrix of corrected measurements, initilized to a zero array, and subsequently storing the corrected version of \code{m} via \code{correct_bias()}
#'
#' \code{biasType} Vector of length p, where p is the number of plates. It tells, for each plate of the assay, A:Additive trend, M:Multiplicative trend, U:Undetermined trend and C:Error-free plate.
#'
#' \code{biasModel} Vector of length p, where p is the number of plates. It tells, for each plate of the assay, the most likely spatial bias model (1 through 6)
#'
#' \code{biasConf} Vector of length p, where p is the number of plates. It tells, for each plate of the assay, the confidence in the model, (0 - lowest to 3-  highest). It is computed by counting the number of bias models (additive or mutliplicative) which agree together.
#' @examples
#' # Fictive 8x12x5 assay
#' assay<-create_assay(m)
#' # Plate 7 taken from Carralot et al. 2012
#' assay<-create_assay(plate7)
#' @export
create_assay<-function(m,ctrl=NA){
  if(!(class(m) %in% c("array","matrix")))
    stop("Error: m is not an array.")
  if(is.na(ctrl)){ # If no control pattern supplied, assuming none
    ctrl=m
    ctrl[]=0
  }
  if(any(!(ctrl%in%0:1))) # If non zero/one values detected
    stop("Control array must be binary")
  if(any(is.na(m)) || any(is.infinite(m))) # If values are missing or are infinite
    stop("NA in input assay")
  assay=NULL
  Rows=dim(m)[1]
  Columns=dim(m)[2]
  rownames(m)=if(Rows<=26)LETTERS[1:Rows]else{warning("Too many rows, naming them by numbers instead of letters.");sapply(1:Rows,toString)}
  colnames(m)=sapply(1:Columns,toString)
  rownames(ctrl)=rownames(m)
  colnames(ctrl)=colnames(m)
  assay$m=m
  assay$ctrl=ctrl
  assay$biasType=rep(NA,dim(m)[3])
  assay$biasModel=rep(NA,dim(m)[3])
  assay$biasConf=rep(NA,dim(m)[3])
  biasPositions=ctrl
  biasPositions[]=0
  assay$biasPositions=biasPositions
  assay$mCorrected=biasPositions
  class(assay)="assay"
  return(assay)
}
#' @title Detect the type of bias present in the assay
#' @description \code{detect}  (1) identifies rows and columns of all plates of the assay affected by spatial bias (following the results of the Mann-Whitney U test); (2) identifies well locations (i.e., well positions scanned across all plates of a given assay) affected by spatial bias (also following the results of the Mann-Whitney U test).
#' @param assay The assay to be corrected. Has to be an \code{assay} object.
#' @param alpha Significance level threshold (defaults to 0.01)
#' @param type \code{P}:plate-specific, \code{A}:assay-specific, \code{PA}:plate then assay-specific, \code{AP}:assay then plate-specific
#' @param test \code{KS}:Kolmogorov-Smirnov (1933), \code{AD}:Anderson-Darling (1952), \code{CVM}:Cramer-von-Mises (1928)
#' @examples
#' assay<-create_assay(m)
#' detected<-detect_bias(assay)
#' @return The corrected assay (\code{assay} object)
#' @export
detect_bias<-function(assay,alpha=0.01,type="P",test="AD"){
  if(class(assay)!="assay")
    stop("Error: This is not an assay.")
  ac=as.character
  m=assay$m
  ctrl=assay$ctrl
  biasType=assay$biasType
  biasModel=assay$biasModel
  biasConf=assay$biasConf
  PMPmapping=c(1,4,5,3,2,6) # Due to a change in the order of methods, need to convert publication notation to code one
  m.E<-new.env()
  for (model in 1:6){
    m.E[[ac(model)]]=ctrl # Initialize empty corrected assay
  }
  dimensions=dim(m)
  Depth=dimensions[3]
  .test.f=NULL
  if(test=="KS")
    .test.f=function(x,y)ks.test(x,y)$p.value
  else if(test=="AD")
    .test.f=function(x,y)as.numeric(tail(strsplit(capture.output(kSamples::ad.test(x,y))[17]," ")[[1]],1))
  else if(test=="CVM")
    .test.f=function(x,y)RVAideMemoire::CvM.test(x,y)$p.value
  else{
    stop("Error: This is not a valid test. Please use KS, AD or CVM")
  }
  if(type=="AP") # If we want assay->plate correction, first apply assay-wise correction
    m=.assay(m,ctrl,alpha)
  for (k in 1:Depth){
    for (model in 1:6){
      m.E[[ac(model)]][,,k]=try(.PMP(m[,,k],ctrl[,,k],PMPmapping[model],alpha)) # Correct the plate k using PMP
      if(class(m.E[[ac(model)]][,,k])=="try-error")
        stop("PMP encountered a problem") # Problem in PMP - check your data
    }
    mww=(m.E[[ac(1)]][,,k]!=m[,,k])*1 # Determine which rows and columns the Mann-Whitney test detected as biased (both technologies have same bias locations, hence using additive)
    assay$biasPositions[,,k]=mww # Save the bias positions suggested by the Mann-Whitney test
    biased.E=new.env()
    unbiased=list()
    for(i in 1:dimensions[1]){
      for(j in 1:dimensions[2]){
        if(mww[i,j]&!ctrl[i,j,k]){ # If the MW test flaged this row/column AND this cell is not a control well
          for(model in 1:6){
            biased.E[[ac(model)]]=c(biased.E[[ac(model)]],m.E[[ac(model)]][i,j,k]) # Add well (i,j) corrected by PMP to set of corrected wells
          }
        }
        else if(!mww[i,j]&!ctrl[i,j,k]) # If the MW did not flag this row/column AND this cell is not a control well
          unbiased=c(unbiased,m.E[[ac(model)]][i,j,k]) # Add this well to set of unbiased wells
      }
    }
    for(model in 1:6){
      biased.E[[ac(model)]]=unlist(biased.E[[ac(model)]])
    }
    unbiased=unlist(unbiased)
    if(Reduce(function(x,y)length(biased.E[[ac(y)]])*x,1:6,1)==0){ # 100% unbiased (computed via fold left)
        biasType[k]='C'
        next
    }
    pvalue.E=new.env()
    for(model in 1:6){
      pvalue.E[[ac(model)]]=.test.f(biased.E[[ac(model)]],unbiased)
    }
    aMethods=1:3
    mMethods=4:6
    p=function(model)pvalue.E[[ac(model)]]
    if(all(sapply(aMethods,p) < alpha) & any(sapply(mMethods,p) > alpha)){ # mPMP did better
      biasType[k]='M' # Bias is multiplicative
      biasModel[k]=3+which.max(sapply(mMethods,p)) # Model number (between 4 and 6)
      biasConf[k]=sum(sapply(mMethods,p) > alpha) # Number of good corrections
    }
    else if(all(sapply(mMethods,p) < alpha) & any(sapply(aMethods,p) > alpha)){ # aPMP did better
      biasType[k]='A' # Bias is additive
      biasModel[k]=which.max(sapply(aMethods,p)) # Model number (between 1 and 3)
      biasConf[k]=sum(sapply(aMethods,p) > alpha) # Number of good corrections
    }
    else if(all(sapply(mMethods,p) < alpha) & all(sapply(aMethods,p) < alpha)){ # Undetermined, both are bad
      biasType[k]='U'
    }
    else if(all(sapply(mMethods,p) > alpha) & all(sapply(aMethods,p) > alpha)){ # Undetermined, both are good
      biasModel[k]=which.max(sapply(c(aMethods,mMethods),p)) # Model number (between 1 and 6)
      biasConf[k]=sum(sapply(c(aMethods,mMethods),p) > alpha) # Number of good corrections
      biasType[k]=ifelse(biasModel[k]%in%aMethods,'A', # Bias is additive
                         'M' # Bias is multiplicative
                         )
    }
    else{
      biasType[k]='U' # Undetermined type of bias, not enough agreement
    }
  }
  assay$biasType=biasType # Write the resulting vectors back
  assay$biasModel=biasModel
  assay$biasConf=biasConf
  return(assay)
}
#' @title Correct the bias present in the assay, previously detected by the \code{detect_bias()} method
#' @description \code{correct_bias()} (1) uses either of the three additive or either of the three multiplicative PMP (Partial Mean Polish) methods (the most appropriate spatial bias model can be either specified or determined by the program following the results of the Kolmogorov-Smirnov, Anderson-Darling or Cramer-von-Mises two-sample test) to correct the assay measurements if the plate-specific correction is specified; (2) carries out the assay-specific correction if specified.
#' @param assay The assay to be corrected. Has to be an \code{assay} object.
#' @param method \code{NULL}:autodetect (default), \code{1}:additive, \code{2}:multiplicative
#' @param alpha Significance level threshold (defaults to 0.05)
#' @param type \code{P}:plate-specific, \code{A}:assay-specific, \code{PA}:plate then assay-specific, \code{AP}:assay then plate-specific
#' @return The corrected assay (\code{assay} object)
#' @examples
#' assay<-create_assay(m)
#' detected<-detect_bias(assay)
#' corrected<-correct_bias(detected,method=2)
#' @export
correct_bias<-function(assay,method=NULL,alpha=0.05,type="PA"){
  if(class(assay)!="assay")
    stop("Error: This is not an assay.")
  m=assay$m
  ctrl=assay$ctrl
  biasType=assay$biasType
  biasModel=assay$biasModel
  dimensions=dim(m)
  Depth=dimensions[3]
  if(is.null(biasType))
    stop("Run detect() first.")
  mCorrected=m # Initialize the corrected assay to the raw one
  if(type=="P"){
    for (k in 1:Depth){
      if(!is.null(method)) # Correct using given method
        mCorrected[,,k]=try(.PMP(m[,,k],ctrl[,,k],method,alpha))
      else{ # Autodetect method for each plate
        if(biasType[k]!='U' & !is.na(biasModel[k]))
          mCorrected[,,k]=try(.PMP(m[,,k],ctrl[,,k],biasModel[k],alpha))
        # else, if the bias is undefined, we cannot apply the correction algorithm, so we skip it
      }
    }
  }
  else if(type=="A"){
    mCorrected=.assay(m,ctrl,alpha)
  }
  else if(type=="PA"){
    # Part 1: Plate-wise
    for (k in 1:Depth){
      if(!is.null(method)) # Correct using given method
        mCorrected[,,k]=try(.PMP(m[,,k],ctrl[,,k],method,alpha))
      else{ # Autodetect method for each plate
        if(biasType[k]!='U' & !is.na(biasModel[k]))
          mCorrected[,,k]=try(.PMP(m[,,k],ctrl[,,k],biasModel[k],alpha))
        # else, if the bias is undefined, we cannot apply the correction algorithm, so we skip it
      }
    }
    # Part 2: Assay-wise
    mCorrected=.assay(m,ctrl,alpha)
  }
  else if(type=="AP"){
    # Part 1: Assay-wise
    mCorrected=.assay(m,ctrl,alpha)
    # Part 2: Plate-wise
    for (k in 1:Depth){
      if(!is.null(method)) # Correct using given method
        mCorrected[,,k]=try(.PMP(m[,,k],ctrl[,,k],method,alpha))
      else{ # Autodetect method for each plate
        if(biasType[k]!='U' & !is.na(biasModel[k]))
          mCorrected[,,k]=try(.PMP(m[,,k],ctrl[,,k],biasModel[k],alpha))
        # else, if the bias is undefined, we cannot apply the correction algorithm, so we skip it
      }
    }
  }
  else{
    stop("Unknown correction type.")
  }
  assay$mCorrected=mCorrected
  return(assay)
}
