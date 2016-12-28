PMP_EPSILON = 0.0001
PMP_MAX_ITERATIONS = 1000
THRESHOLD=2.5
.zify=function(m)(m-mean(m))/ifelse(!sd(m),PMP_EPSILON,sd(m))
.rify=function(m)(m-median(m))/ifelse(!mad(m),PMP_EPSILON,mad(m))
.PMP=function(m,ctrl,method,alpha=0.05){
  bak=m
  ctrl[is.na(ctrl)] <- 1 # All missing controls and reals are excluded (flagged as control wells)
  ctrl[is.na(m)] <- 1
  Rows=dim(m)[1]
  Columns=dim(m)[2]
  Loop=1
  tmp.c=colSums(ctrl)==Rows
  tmp.r=rowSums(ctrl)==Columns
  if(any(tmp.c)){m=m[,-which(tmp.c)];ctrl=ctrl[,-which(tmp.c)]}
  if(any(tmp.r)){m=m[-which(tmp.r),];ctrl=ctrl[-which(tmp.r),]}
  Rows.b=Rows
  Columns.b=Columns
  Rows=dim(m)[1]
  Columns=dim(m)[2]
  bak2=m
  CFlag=rep(FALSE,Columns)
  RFlag=rep(FALSE,Rows)
  n=0
  ctrlS=sum(ctrl)
  while (Loop< min(Rows-2,Columns-2) && (n < .5*(length(RFlag)*length(CFlag)-ctrlS))){
    p.values=data.frame(r=numeric(),c=numeric(),p=numeric())
    Loop=Loop+1
    Cind=which(CFlag)
    Rind=which(RFlag)
    for (i in 1:Columns){
      if (i %in% Cind)next
      tmp1=list()
      tmp2=list()
      for(j in 1:Rows){
        for(k in 1:Columns){
          if(!(j %in% Rind) & (k == i) & !ctrl[j,k]) tmp1=c(tmp1,m[j,k])
          if(!(j %in% Rind) & (k != i) & !ctrl[j,k] & !(k %in% Cind)) tmp2=c(tmp2,m[j,k])
        }
      }
      tmp1=unlist(tmp1)
      tmp2=unlist(tmp2)
      if(!length(tmp1)||!length(tmp2))next
      p.v=try(wilcox.test(tmp1,tmp2,correct = FALSE)$p.value,TRUE)
      if(class(p.v)!="try-error")p.values=rbind(p.values,data.frame(c=i,r=-1,p=ifelse(is.nan(p.v),1,p.v)))
    }
    for (j in 1:Rows){
      if (j %in% Rind)next
      tmp1=list()
      tmp2=list()
      ct=0
      for(i in 1:Rows){
        for(k in 1:Columns){
          if(!(k %in% Cind) & (j == i) & !ctrl[i,k]) tmp1=c(tmp1,m[i,k])
          if(!(k %in% Cind) & (j != i) & !ctrl[i,k] & !(i %in% Rind)) tmp2=c(tmp2,m[i,k])
          if(ctrl[i,k]) ct=ct+1
        }
      }
      tmp1=unlist(tmp1)
      tmp2=unlist(tmp2)
      if(!length(tmp1)||!length(tmp2))next
      p.v=try(wilcox.test(tmp1,tmp2,correct = FALSE)$p.value,TRUE)
      if(class(p.v)!="try-error")p.values=rbind(p.values,data.frame(c=-1,r=j,p=ifelse(is.nan(p.v),1,p.v)))
    }
    # Sorted p-values in increasing order (most evidence to reject -> least evidence to reject)
    p.values=p.values[order(p.values$p),]
    # If all p-values are greater than alpha, break
    if (all(p.values$p >=alpha))break
    if(p.values[1,]$c==-1){
      i=p.values[1,]$r
      RFlag[i]=TRUE
    }else{
      j=p.values[1,]$c
      CFlag[j]=TRUE
    }
    n=(sum(RFlag)*length(CFlag)+(length(RFlag)-sum(RFlag))*sum(CFlag))
  }
  # Number of row/column/total errors
  NR=max(0,sum(RFlag,na.rm=TRUE))
  NC=max(0,sum(CFlag,na.rm=TRUE))
  N=NR+NC
  if(N==0){
    i.m=1
    for (i in 1:Rows.b){
      j.m=1
      for (j in 1:Columns.b){
        if(!tmp.r[i]&!tmp.c[j])
          bak[i,j]=m[i.m,j.m]
        if(!tmp.c[j])
          j.m=j.m+1
      }
      if(!tmp.r[i])
        i.m=i.m+1
    }
    return(bak)
  }
  Mu=list()
  for(i in 1:Rows){
    if (RFlag[i])next
    for(j in 1:Columns)
    {
      if (CFlag[j] || ctrl[i,j])next
      Mu=c(Mu,m[i,j]) # median
    }
  }
  # Mu=median(unlist(Mu))
  Mu=mean(unlist(Mu))
  Rmu=rep(0,Rows)
  Cmu=rep(0,Columns)
  Loop=1
  Converge=0
  repeat{
    Rmu=rep(0,Rows)
    Cmu=rep(0,Columns)
    Diff=0
    Converge=0
    for (i in 1:Rows){
      intersection=0
      q=list()
      for (j in 1:Columns){
        if(ctrl[i,j])next
        if(RFlag[i]&&CFlag[j]){
          intersection=intersection+1
          next
        }
        q=c(q,m[i,j]) # median
      }
      # Rmu[i]=median(unlist(q))
      Rmu[i]=mean(unlist(q))
    }
    for (j in 1:Columns){
      intersection=0
      q=list()
      for (i in 1:Rows){
        if(ctrl[i,j])next
        if(RFlag[i]&&CFlag[j]){
          intersection=intersection+1
          next
        }
        q=c(q,m[i,j]) # median
      }
      # Cmu[j]=median(unlist(q))
      Cmu[j]=mean(unlist(q))
    }
    if (method %in% c(1,2)){
      for (i in 1:Rows){
        if(!RFlag[i]) next
        Diff=Mu-Rmu[i]
        Converge=Converge+abs(Diff)
        for (j in 1:Columns){
          if(ctrl[i,j])next
          m[i,j]=switch(method,
                        m[i,j]+Diff, # Method 1
                        m[i,j]*abs(Mu/Rmu[i]) # Method 2
          )
        }
      }

      for (j in 1:Columns){
        if(!CFlag[j]) next
        Diff=Mu-Cmu[j]
        Converge=Converge+abs(Diff)
        for (i in 1:Rows){
          if(ctrl[i,j])next
          m[i,j]=switch(method,
                        m[i,j]+Diff, # Method 1
                        m[i,j]*abs(Mu/Cmu[j]) # Method 2
          )
        }
      }}
    if (method %in% c(3,4,5,6)){
      for (i in 1:Rows){
        for (j in 1:Columns){
          if(ctrl[i,j])next
          if (CFlag[j] && !RFlag[i]){
            Diff=Mu-Cmu[j]
            Converge=Converge+abs(Diff)
            m[i,j]=switch(method-2,abs(Mu*m[i,j]/(Cmu[j])), # Method 3
                          m[i,j]-(Cmu[j]-Mu), # Method 4
                          m[i,j]-(Cmu[j]-Mu), # Method 5
                          abs(Mu*m[i,j]/Cmu[j])) # Method 6
          }
          else if (!CFlag[j] && RFlag[i]){
            Diff=Mu-Rmu[i]
            Converge=Converge+abs(Diff)
            m[i,j]=switch(method-2,abs(Mu*m[i,j]/(Rmu[i])), # Method 3
                          m[i,j]-(Rmu[i]-Mu), # Method 4
                          m[i,j]-(Rmu[i]-Mu), # Method 5
                          abs(Mu*m[i,j]/Rmu[i])) # Method 6
          }
        }
      }
    }
    # Convergence condition
    if(class(Converge)!="numeric")Converge=.Machine$double.xmax
    Loop=Loop+1
    if(! (Converge >PMP_EPSILON && Loop<PMP_MAX_ITERATIONS) ){
      break
    }
  }
  # Intersection polishing
  if(method %in% c(3,4,5,6)){
    Rmu=rep(0,Rows)
    Cmu=rep(0,Columns)
    for (i in 1:Rows){
      intersection=0
      q=list()
      for (j in 1:Columns){
        if(ctrl[i,j])next
        if(RFlag[i]&&CFlag[j]){
          intersection=intersection+1
          next
        }
        z=switch(method-2,
                 bak2[i,j]/m[i,j], # Method 3
                 bak2[i,j]-m[i,j], # Method 4
                 bak2[i,j]-m[i,j], # Method 5
                 bak2[i,j]/m[i,j]) # Method 6
        q=c(q,z) # avg
      }
      # Rmu[i]=median(unlist(q)) # median
      Rmu[i]=mean(unlist(q)) # median
    }
    for (j in 1:Columns){
      intersection=0
      q=list()
      for (i in 1:Rows){
        if(ctrl[i,j])next
        if(RFlag[i]&&CFlag[j]){
          intersection=intersection+1
          next
        }
        z=switch(method-2,
                 bak2[i,j]/m[i,j], # Method 3
                 bak2[i,j]-m[i,j], # Method 4
                 bak2[i,j]-m[i,j], # Method 5
                 bak2[i,j]/m[i,j]) # Method 6
        q=c(q,z) # median
      }
      # Cmu[j]=median(unlist(q))
      Cmu[j]=mean(unlist(q))
    }
    for (i in 1:Rows){
      for (j in 1:Columns){
        if(ctrl[i,j])next
        if (CFlag[j] && RFlag[i]){
          Diff=2*Mu-Cmu[j]-Rmu[i]
          Converge=Converge+abs(Diff)
          m[i,j]=switch(method-2,
                        m[i,j]/abs(Rmu[i]+Cmu[j]),          # Method 3
                        m[i,j]-(Rmu[i]*Cmu[j]),          # Method 4
                        m[i,j]-(Rmu[i]+Cmu[j])/2,        # Method 5
                        m[i,j]/sqrt(abs(Rmu[i]*Cmu[j]))) # Method 6
        }
      }
    }
  }
  i.m=1
  for (i in 1:Rows.b){
    j.m=1
    for (j in 1:Columns.b){
      if(!tmp.r[i]&!tmp.c[j]){
        bak[i,j]=m[i.m,j.m]
      }
      if(!tmp.c[j])
        j.m=j.m+1
    }
    if(!tmp.r[i])
      i.m=i.m+1
  }
  return(bak)
}
.assay=function(m,ctrl,alpha){
  ctrl[is.na(ctrl)] <- 1
  ctrl[is.na(m)] <- 1
  dims=dim(m)
  Depth=dims[3]
  raw.bak=m
  for(k in 1:Depth){
    m[,,k]=.rify(m[,,k])
  }
  bak=m
  bak2=ctrl
  m=array(dim=c(dims[3],dims[1]*dims[2]))
  ctrl=array(dim=c(dims[3],dims[1]*dims[2]))
  for(i in 1:Depth){
    m[i,]=as.vector(bak[,,i])
    ctrl[i,]=as.vector(bak2[,,i])
  }
  rm(bak,bak2)
  Columns=dim(m)[2]
  CFlag=rep(FALSE,Columns)
  Loop=0
  n=0
  ctrlS=sum(ctrl)
  CFlag=rep(TRUE,Columns) # AD HOC
  if(!sum(CFlag))return(raw.bak) # If no errors, return non normalized matrix
  if(Depth>1){
    for(i in 1:length(CFlag)){ # Normalize biased wells to correct error
      if(!CFlag[i]) next
      mu=mean(m[,i],na.rm = TRUE)
      sd=sd(m[,i],na.rm = TRUE)
      tmp=array(dim=c(length(m[,i])))
      for(j in 1:length(m[,i])){
        if (!is.na(m[j,i])&(!(m[j,i]>=mu+THRESHOLD*sd || m[j,i]<=mu-THRESHOLD*sd)))tmp[j]=m[j,i]
      }
      tmp=(tmp-mean(tmp,na.rm=TRUE))/ifelse(!sd(tmp,na.rm = TRUE),PMP_EPSILON,sd(tmp,na.rm = TRUE))
      for(j in 1:length(m[,i])){
        if (!is.na(tmp[j]))m[j,i]=tmp[j]
      }
    }
  }
  m=array(as.vector(t(m)),dim=dims)
  return(m)
}
