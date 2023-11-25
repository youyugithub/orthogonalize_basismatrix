# orthogonalize_basismatrix
```
orthogonalize_basismatrix_stable<-function(basismatrix){
  
  normalize<-function(x)x/sqrt(sum(x^2)/ntimepoints)
  
  nbasis<-ncol(basismatrix)
  ntimepoints<-nrow(basismatrix)
  if(nbasis==1){
    Y0<-rep(1,ntimepoints)
    return(matrix(Y1))
  }
  if(nbasis==2){
    Y0<-rep(1,ntimepoints)
    Y1<-seq(-1,1,length.out=ntimepoints)
    Y1<-Y1/sqrt(sum(Y1^2/ntimepoints))
    return(cbind(Y0,Y1))
  }
  
  pattern<-rbind(
    c( 1,-2, 1),
    c(-2, 4,-2),
    c( 1,-2, 1))
  Q<-matrix(0,ntimepoints,ntimepoints)
  for(jj in 2:(ntimepoints-1)){
    Q[jj+(-1:1),jj+(-1:1)]<-Q[jj+(-1:1),jj+(-1:1)]+pattern
  }
  
  basismatrix_ortho<-cbind(1,seq(-1,1,length.out=ntimepoints),basismatrix[,-(1:2)])
  for(jj in 1:nbasis)basismatrix_ortho[,jj]<-normalize(basismatrix_ortho[,jj])
  basismatrix_ortho[,1]<-normalize(basismatrix_ortho[,1])
  for(gg in 2:nbasis){
    projection<-sum(basismatrix_ortho[,gg]*basismatrix_ortho[,1]/ntimepoints)
    basismatrix_ortho[,gg]<-basismatrix_ortho[,gg]-basismatrix_ortho[,1]*projection
    basismatrix_ortho[,gg]<-normalize(basismatrix_ortho[,gg])
  }
  basismatrix_ortho[,2]<-normalize(basismatrix_ortho[,2])
  for(gg in 3:nbasis){
    projection<-sum(basismatrix_ortho[,gg]*basismatrix_ortho[,2]/ntimepoints)
    basismatrix_ortho[,gg]<-basismatrix_ortho[,gg]-basismatrix_ortho[,2]*projection
    basismatrix_ortho[,gg]<-normalize(basismatrix_ortho[,gg])
  }
  
  for(jj in 3:nbasis){
    mat_denominator<-t(basismatrix_ortho[,jj:nbasis])%*%basismatrix_ortho[,jj:nbasis]
    mat_denominator<-0.5*(mat_denominator+t(mat_denominator))
    mat_numerator<-t(basismatrix_ortho[,jj:nbasis])%*%Q%*%basismatrix_ortho[,jj:nbasis]
    mat_numerator<-0.5*(mat_numerator+t(mat_numerator))
    sqrt_denominator<-expm::sqrtm(mat_denominator)
    sqrt_denominator<-0.5*(sqrt_denominator+t(sqrt_denominator))
    mat_temp<-sqrt_denominator%*%solve(mat_numerator)%*%sqrt_denominator
    mat_temp<-0.5*(mat_temp+t(mat_temp))
    eigen_temp<-eigen(mat_temp)
    
    sqrt_denominator_bb<-eigen_temp$vectors[,1]
    bb<-c(solve(sqrt_denominator,sqrt_denominator_bb))
    basismatrix_ortho[,jj]<-c(basismatrix_ortho[,jj:nbasis,drop=F]%*%bb)
    basismatrix_ortho[,jj]<-sign(basismatrix_ortho[ntimepoints,jj])*basismatrix_ortho[,jj]/sqrt(sum(basismatrix_ortho[,jj]^2)/ntimepoints)
    if(jj<nbasis){
      for(gg in (jj+1):nbasis){
        projection<-sum(basismatrix_ortho[,gg]*basismatrix_ortho[,jj]/ntimepoints)
        basismatrix_ortho[,gg]<-basismatrix_ortho[,gg]-basismatrix_ortho[,jj]*projection
        basismatrix_ortho[,gg]<-normalize(basismatrix_ortho[,gg])
      }
    }
  }
  return(basismatrix_ortho)
}

orthogonalize_basismatrix_twice<-function(basismatrix){
  basismatrix_transform<-orthogonalize_basismatrix_stable(basismatrix)
  transform<-MASS::ginv(basismatrix)%*%basismatrix_transform
  basismatrix_transform2<-orthogonalize_basismatrix_stable(basismatrix%*%transform)
  return(basismatrix_transform2)
}

orthogonalize_basismatrix<-function(basismatrix){
  
  nbasis<-ncol(basismatrix)
  ntimepoints<-nrow(basismatrix)
  if(nbasis==1){
    Y0<-rep(1,ntimepoints)
    return(matrix(Y1))
  }
  if(nbasis==2){
    Y0<-rep(1,ntimepoints)
    Y1<-seq(-1,1,length.out=ntimepoints)
    Y1<-Y1/sqrt(sum(Y1^2/ntimepoints))
    return(cbind(Y0,Y1))
  }
  
  pattern<-rbind(
    c( 1,-2, 1),
    c(-2, 4,-2),
    c( 1,-2, 1))
  Q<-matrix(0,ntimepoints,ntimepoints)
  for(jj in 2:(ntimepoints-1)){
    Q[jj+(-1:1),jj+(-1:1)]<-Q[jj+(-1:1),jj+(-1:1)]+pattern
  }
  
  basismatrix_ortho<-cbind(1,seq(-1,1,length.out=ntimepoints),basismatrix[,-(1:2)])
  basismatrix_ortho[,1]<-basismatrix_ortho[,1]/sqrt(sum(basismatrix_ortho[,1]^2)/ntimepoints)
  for(gg in 2:nbasis){
    projection<-sum(basismatrix_ortho[,gg]*basismatrix_ortho[,1]/ntimepoints)
    basismatrix_ortho[,gg]<-basismatrix_ortho[,gg]-basismatrix_ortho[,1]*projection
  }
  basismatrix_ortho[,2]<-basismatrix_ortho[,2]/sqrt(sum(basismatrix_ortho[,2]^2)/ntimepoints)
  for(gg in 3:nbasis){
    projection<-sum(basismatrix_ortho[,gg]*basismatrix_ortho[,2]/ntimepoints)
    basismatrix_ortho[,gg]<-basismatrix_ortho[,gg]-basismatrix_ortho[,2]*projection
  }
  
  for(jj in 3:nbasis){
    mat_denominator<-t(basismatrix_ortho[,jj:nbasis])%*%basismatrix_ortho[,jj:nbasis]
    mat_denominator<-0.5*(mat_denominator+t(mat_denominator))
    mat_numerator<-t(basismatrix_ortho[,jj:nbasis])%*%Q%*%basismatrix_ortho[,jj:nbasis]
    mat_numerator<-0.5*(mat_numerator+t(mat_numerator))
    sqrt_denominator<-expm::sqrtm(mat_denominator)
    sqrt_denominator<-0.5*(sqrt_denominator+t(sqrt_denominator))
    sqrtinv_denominator<-solve(mat_denominator)
    sqrtinv_denominator<-0.5*(sqrtinv_denominator+t(sqrtinv_denominator))
    mat_temp<-sqrtinv_denominator%*%mat_numerator%*%sqrtinv_denominator
    mat_temp<-0.5*(mat_temp+t(mat_temp))
    eigen_temp<-eigen(mat_temp)
    
    sqrt_denominator_bb<-eigen_temp$vectors[,ncol(eigen_temp$vectors)]
    bb<-c(sqrtinv_denominator%*%sqrt_denominator_bb)
    basismatrix_ortho[,jj]<-c(basismatrix_ortho[,jj:nbasis,drop=F]%*%bb)
    basismatrix_ortho[,jj]<-basismatrix_ortho[,jj]/sqrt(sum(basismatrix_ortho[,jj]^2)/ntimepoints)
    if(jj<nbasis){
      for(gg in (jj+1):nbasis){
        projection<-sum(basismatrix_ortho[,gg]*basismatrix_ortho[,jj]/ntimepoints)
        basismatrix_ortho[,gg]<-basismatrix_ortho[,gg]-basismatrix_ortho[,jj]*projection
      }
    }
  }
  return(basismatrix_ortho)
}


##### test
# a_bs<-splines::bs(1:ntimepoints,df=8,intercept=T)
# matplot(a_bs,type="l")
# matplot(orthogonalize_basismatrix_stable(a_bs),type="l")
# matplot(orthogonalize_basismatrix_twice(a_bs),type="l")
# basismatrix_ortho<-orthogonalize_basismatrix_stable(a_bs)
# t(basismatrix_ortho)%*%basismatrix_ortho
# temp<-orthogonalize_basismatrix_stable(orthogonalize_basismatrix_stable(a_bs))
# t(temp)%*%temp
# temp
# matplot(orthogonalize_basismatrix_stable(basismatrix),type="l")
```
