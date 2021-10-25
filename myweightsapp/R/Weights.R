rm(list = ls())

#libraries
#install.packages('rdrobust')
#install.packages('np')
library(rdrobust)
library(np)

#Function to check if the input is a vector
verify_is_vector=function(x){
  if(is.matrix(x)==FALSE || nrow(x)<2){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}
verify_is_matrix=function(x){
  if(is.vector(x)==FALSE || length(x)<2){
    if(is.matrix(x)==FALSE || nrow(x)<2){
      return(FALSE)
    }
    else{
      return(TRUE)
    }
  }
  else{
    return(TRUE)
  }

}

#Function to check if the input is numeric
verify_vector_numeric=function(x){
  for (i in x){

    if (is.numeric(i)==FALSE){
      return(FALSE)
    }
  }
  return(TRUE)
}

#Function to check if the input length of one vector is equal to the other
same_length=function(a,b){
  if(nrow(a)!=nrow(b)){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

#Function to check that the input c is numeric
verify_c=function(x){
  if(is.numeric(x)==FALSE || length(x)!=1){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

verify_c_range=function(r,c){
  if (c<min(r)||c>max(r)){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

verify_is_factor=function(r){
  if(is.factor(r)){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}






my_weights<-function(r, x, c) {

  ri=cbind(r)
  xi=cbind(x)
  #Libraries for the function
  #library(rdrobust)
  #library(np) # To perform LLR
  #library(np)

  #Verify input
  count5=1
  for (i in ri){
    if (is.nan(i)||is.na(i)){
      ri=cbind(ri[-count5])
      xi=cbind(xi[-count5,])
    }
    count5=count5+1
  }
  count6=1
  while(count6<ncol(xi)){
    count7=1
    for (i in xi[,count6]){
      if (is.nan(i)||is.na(i)){
        ri=cbind(ri[-count7])
        xi=cbind(xi[-count7,])
      }
      count7=count7+1
    }
    count6=count6+1
  }

  verificator_list=list(verify_is_vector(ri),
                        verify_is_matrix(xi),
                        verify_vector_numeric(ri),
                        verify_vector_numeric(xi),
                        same_length(ri,xi),
                        verify_c(c),
                        verify_c_range(ri,c),
                        verify_is_factor(r))


  #verificar que r este en el rango

  error_list=list("R must be a vector with length greater than one",
                  "X must be a matrix with length greater than one",
                  "R must be a numeric vector",
                  "X must be a numeric matrix",
                  "R and X must have the same number of rows",
                  "C must be numeric of length one",
                  "C must be in the range of R",
                  "R must be a continous variable")
  r=ri
  x=xi


  problem=FALSE
  counter=1
  for (i in verificator_list){
    if (i==FALSE){
      message(error_list[counter])
      problem=TRUE
    }
    counter=counter+1
  }
  if (problem==TRUE){
    message('There is an input problem')
    return('')
  }

  if (problem==TRUE){
    message('There is an input problem')
    return('')
  }


  # estimation of bandwidths
  bw_xr=npcdensbw(ydat=x,xdat= r, cykertype="epanechnikov", cxkertype="epanechnikov") # bw comun a todos

  #ptm=proc.time()



  ######################
  # Preparing data
  ######################
  #hx es vector hr escalar
  h_x=bw_xr$ybw

  h_r=min(bw_xr$xbw,1000)#agregar en recomendaciones





  ## 4 casos para bw de r ##
  # : 0,75 de h_r de f(x|r)
  # 2: h_r de f(x|r)
  # 3: 1,25 de h_r de f(x|r)
  # 4: h_r de cat_RD

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
  count4=1
  condition2=(abs((r-c)/h_r)<=2*sqrt(5)) # el igual es por si acaso
  for (i in condition2){
    if(is.na(i)){
      condition2[count4]=TRUE
    }
    count4=count4+1
  }
  id=seq(1:n)

  #
  if(ncol(x)==1){
    nece=data.frame(id=id[condition2], r=r[condition2], x=x[condition2], estos=condition1[condition2])

  }
  nece=data.frame(id=id[condition2], r=r[condition2], x=x[condition2,], estos=condition1[condition2])

  ################################################
  # estimacion del weight solo para los necesarios
  ################################################



  w_si=rep(NaN, sum(condition2)) # para for loop

  #Lista para guardar todos los ker_x y ker_x_result que es luego de la multiplicacion
  ker_xs=c()
  ker_xs_results=c()


  for (ii in (which.max(nece$estos)):(sum(nece$estos)+which.max(nece$estos)-1)) {

    #Ir cambiando ii hasta 3
    #  prepare data

    condition3=(abs((nece$r-nece$r[ii])/h_r)<sqrt(5)) # solo los cercanos a ii

    #Crear while y multiplicar kernels
    count2=1
    count3=1
    #ker_xs=c()
    #ker_xs_results=c()
    ker_x_result=1
    for (i in nece){
      if (count2>2 && count2<ncol(nece)){
        ker_x=(0.75*(1-0.2*((i[condition3]-i[ii])/h_x[count3])^2)/sqrt(5))*(
          abs((i[condition3]-i[ii])/h_x[count3])<sqrt(5))*(1/h_x[count3])
        ker_xs=c(ker_xs,ker_x)
        ker_x_result=ker_x_result*ker_x
        count3=count3+1
      }
      count2=count2+1
    }
    ker_xs_results=c(ker_xs_results,ker_x_result)
    #ker_x=1
    #for(i in ker_xs){
    #  if (!is.na(i)){
    #    print(i)
    #    ker_x=ker_x*i
    #  }

    #}

    #Ker_x = multiplicaci?n de los Ker_x

    R=cbind(rep(1,sum(condition3)),(nece$r[condition3]-nece$r[ii]))
    ker_r=(0.75*(1-0.2*((nece$r[condition3]-nece$r[ii])/h_r)^2)/sqrt(5))
    # 2. contrafactual conditional density through LLR:
    GinvC=chol2inv(chol(crossprod(sqrt(ker_r)*R)))
    #Se ocup? el ker_x_result para la formula de abajo
    num=c(1, 0)%*%GinvC%*%crossprod(R*ker_r,ker_x_result) # sin limite inferior
    # 3. empirical conditional density through LLR:

    if (nece$r[ii]<c) {
      ker_rL=(0.75*(1-0.2*((nece$r[condition3]-nece$r[ii])/h_r)^2)/sqrt(5))*(nece$r[condition3]<c)
      GinvL=chol2inv(chol(crossprod(sqrt(ker_rL)*R)))
      den2=c(1, 0)%*%GinvL%*%crossprod(R*ker_rL,ker_x_result)
      den=ifelse(abs(den2)<0.0001,0.0001,den2) # si es muy chico
    } else {
      ker_rR=(0.75*(1-0.2*((nece$r[condition3]-nece$r[ii])/h_r)^2)/sqrt(5))*(nece$r[condition3]>=c)
      GinvR=chol2inv(chol(crossprod(sqrt(ker_rR)*R)))
      #La linea tira error para ii=3
      den2=c(1, 0)%*%GinvR%*%crossprod(R*ker_rR,ker_x_result)
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
  plot(r,w)
  out=list(w=w, h_x=h_x, h_r=h_r)
  return(out)
}


