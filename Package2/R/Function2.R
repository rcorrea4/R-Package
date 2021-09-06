rm(list = ls())

library(rdrobust)
#library(np) # To perform LLR
library(np)

largo_lista<-function(list_a,list_b){
  if (length(list_a)==length(list_b))
    return(TRUE)
  else
    return(FALSE)
}

#Function for verifying
verify_is_vector=function(x){
  if(is.vector(x)==FALSE){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

verify_vector_numeric=function(x){
  for (i in x){
    if (is.numeric(i)==FALSE){
      return(FALSE)
    }
  }
  return(TRUE)
}


same_length=function(a,b){
  if(length(a)!=length(b)){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

verify_c=function(x){
  if(is.numeric(x)==FALSE || length(x)!=1){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

a=c(1,2,3)
b=c(1,2,3)
d=1

verificator_list=list(verify_is_vector(a),verify_is_vector(b),verify_vector_numeric(a),verify_vector_numeric(b),same_length(a,b),verify_c(d))

error_list=list("R no es vector", "X no es vector", "R no es num?rico", "X no es num?rico", "R y X no son del mismo largo", "C debe ser solo un n?mero")
cont=1
for (i in verificator_list){
  if (i==FALSE){
    print(error_list[cont])
  }
  cont=cont+1
}

myweights<-function(r, x, c) {
  library(rdrobust)
  #library(np) # To perform LLR
  library(np)
  verificator_list=list(verify_is_vector(r),verify_is_vector(x),verify_vector_numeric(r),verify_vector_numeric(x),same_length(r,x),verify_c(c))
  error_list=list("R no es vector", "X no es vector", "R no es num?rico", "X no es num?rico", "R y X no son del mismo largo", "C debe ser solo un n?mero")
  errors=list()
  for (i in verificator_list){
    if (i==FALSE){
      print(error_list[i])
    }
  }







  # estimacion de bandwidths
  bw_xr=npcdensbw(x ~ r, cykertype="epanechnikov", cxkertype="epanechnikov") # bw comun a todos

  #ptm=proc.time()

  ######################
  # Preparing data
  ######################
  h_x=bw_xr$ybw
  ## 4 casos para bw de r ##
  # 1: 0,75 de h_r de f(x|r)
  # 2: h_r de f(x|r)
  # 3: 1,25 de h_r de f(x|r)
  # 4: h_r de cat_RD
  h_r=(bw_xr$xbw)
  n=length(r)
  w=array(NaN,dim=c(n,1))
  #generar n, largo de r
  N_ef_w=array(NaN,dim=c(1,1))


  #l=2
  ########################################
  # weight 1 para quienes estan lejos de c
  ########################################
  condition1=(abs((r-c)/h_r)<sqrt(5)) # verdadero quienes estan cerca
  w[condition1==FALSE]=1
  N_ef_w=sum(condition1) # numero efectivo de observaciones que se usan en la estimacion
  #################################################
  # para los weights no necesitamos todos los datos
  #################################################
  condition2=(abs((r-c)/h_r)<=2*sqrt(5)) # el igual es por si acaso
  id=seq(1:n)
  nece=data.frame(id=id[condition2], r=r[condition2], x=x[condition2], estos=condition1[condition2])

  ################################################
  # estimacion del weight solo para los necesarios
  ################################################
  w_si=rep(NaN, sum(condition2)) # para for loop
  for (ii in (which.max(nece$estos)):(sum(nece$estos)+which.max(nece$estos)-1)) {
    # 1. prepare data
    condition3=(abs((nece$r-nece$r[ii])/h_r)<sqrt(5)) # solo los cercanos a ii
    ker_x=(0.75*(1-0.2*((nece$x[condition3]-nece$x[ii])/h_x)^2)/sqrt(5))*(abs((nece$x[condition3]-nece$x[ii])/h_x)<sqrt(5))*(1/h_x)
    R=cbind(rep(1,sum(condition3)),(nece$r[condition3]-nece$r[ii]))
    ker_r=(0.75*(1-0.2*((nece$r[condition3]-nece$r[ii])/h_r)^2)/sqrt(5))
    # 2. contrafactual conditional density through LLR:
    GinvC=chol2inv(chol(crossprod(sqrt(ker_r)*R)))
    num=c(1, 0)%*%GinvC%*%crossprod(R*ker_r,ker_x) # sin limite inferior
    # 3. empirical conditional density through LLR:
    if (nece$r[ii]<c) {
      ker_rL=(0.75*(1-0.2*((nece$r[condition3]-nece$r[ii])/h_r)^2)/sqrt(5))*(nece$r[condition3]<c)
      GinvL=chol2inv(chol(crossprod(sqrt(ker_rL)*R)))
      den2=c(1, 0)%*%GinvL%*%crossprod(R*ker_rL,ker_x)
      den=ifelse(abs(den2)<0.0001,0.0001,den2) # si es muy chico
    } else {
      ker_rR=(0.75*(1-0.2*((nece$r[condition3]-nece$r[ii])/h_r)^2)/sqrt(5))*(nece$r[condition3]>=c)
      GinvR=chol2inv(chol(crossprod(sqrt(ker_rR)*R)))
      den2=c(1, 0)%*%GinvR%*%crossprod(R*ker_rR,ker_x)
      den=ifelse(abs(den2)<0.0001,0.0001,den2) # si es muy chico
    }
    # 4. weight estimation
    w_si[ii]=num/den

  }
  send_w=data.frame(id=nece$id, w=w_si) # link de weights con el id

  # reemplazar los valores estimados en el vector de weights
  w[(nece$id[which.max(nece$estos)]):(nece$id[sum(nece$estos)+which.max(nece$estos)-1])]=send_w$w[!is.na(send_w$w)]


  #proc.time()-ptm

  # output
  out=list(w=w, h_x=h_x, h_r=h_r, N_ef_w=N_ef_w)
  return(out)
}

myweights(a,b,d)

verificator_list=list(verify_is_vector(a),verify_is_vector(b),verify_vector_numeric(a),verify_vector_numeric(b),same_length(a,b),verify_c(d))
