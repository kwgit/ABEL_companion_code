# Functions----
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

Blocks <- function(n.use, M, L){
  Q <- floor((n.use-M)/L)+1
  starting.points <- seq(1,1+(Q-1)*L,by = L )
  end.points <- seq(M,M+(Q-1)*L, by =L)
  blocks.ind <- apply(matrix( c(starting.points,end.points), nrow = 2, byrow = T), 2, function(p){
    seq(p[1], p[2])
  })
  return( list(blocks.ind, Q)) 
}

Matrix.power <- function (mat, pow) {

  mat.eigVal <- eigen(mat)$values
  mat.eigVec <- eigen(mat)$vectors
  
  mat.pow <- mat - mat  # 0 matrix
  for (i in 1:length(mat.eigVal)) {
    mat.pow = mat.pow + mat.eigVal[i]^pow * mat.eigVec[,i] %*% t(mat.eigVec[,i])
  }
  
  return (mat.pow)
}


Ti.fun = function(est_fun, blocks.ind, Q, M){
  
  if (is.null(dim(est_fun))) {
    Ti <- est_fun[1:Q]
    for (i in 1:Q) {
      Ti[i] <- sum(est_fun[blocks.ind[,i]])/M
    }
  }else{
    Ti <- est_fun[1:Q,]
    for (i in 1:Q) {
      Ti[i,] <- apply(est_fun[blocks.ind[,i],], 2, sum)/M
    }
  }
  
  return (Ti)
}

Ti.fun.pro <- function (est_fun, block.list, kw.n) {
  Q <- length(block.list)
  if (is.null(dim(est_fun))) {
    Ti <- est_fun[1:Q]
    for (i in 1:Q) {
      Ti[i] <- mean(est_fun[block.list[[i]][which(block.list[[i]]<=kw.n)]])
    }
  }else{
    Ti <- est_fun[1:Q,]
    for (i in 1:Q) {
      Ti[i,] <- apply(est_fun[block.list[[i]][which(block.list[[i]]<=kw.n)],], 2, mean)
    }
  }
  return (Ti)
}


Prod.sum <- function(ijk, lag.par, kw.Q, kw.y) {
  if (lag.par >= 0) {
    result <- sum(sapply(seq(1, kw.Q-lag.par), function(comp.ind){
      yyy <- 1
      counter <- 0
      for (i in ijk) {
        counter <- counter + 1
        if (counter == length(ijk)){
          yyy <- yyy * kw.y[comp.ind+lag.par, i]
        }else{
          yyy <- yyy * kw.y[comp.ind, i]
        }
      }
      return (yyy)  
    }))
  }else{
    result <- sum(sapply(seq(1-lag.par, kw.Q), function(comp.ind){
      yyy <- 1
      counter <- 0
      for (i in ijk) {
        counter <- counter + 1
        if (counter == length(ijk)){
          yyy <- yyy * kw.y[comp.ind+lag.par, i]
        }else{
          yyy <- yyy * kw.y[comp.ind, i]
        }
      }
      return (yyy) 
    }))
  }
  return (result)
}

Prod.sum2 <- function(ijk, lag.par, kw.Q, kw.y) {
 
  if (lag.par >= 0) {
    result <- sum(sapply(seq(1, kw.Q-lag.par), function(comp.ind){
      yyy <- 1
      counter <- 0
      for (i in ijk) {
        counter <- counter + 1
        if (counter >= length(ijk)-1){
          yyy <- yyy * kw.y[comp.ind+lag.par, i]
        }else{
          yyy <- yyy * kw.y[comp.ind, i]
        }
      }
      return (yyy)  
    }))
  }else{
    result <- sum(sapply(seq(1-lag.par, kw.Q), function(comp.ind){
      yyy <- 1
      counter <- 0
      for (i in ijk) {
        counter <- counter + 1
        if (counter >= length(ijk)-1){
          yyy <- yyy * kw.y[comp.ind+lag.par, i]
        }else{
          yyy <- yyy * kw.y[comp.ind, i]
        }
      }
      return (yyy) 
    }))
  }
  return (result)
}

Prod.sum3 <- function(ijk, lag.par, kw.Q, kw.y) {
 
  if (lag.par >= 0) {
    result <- sum(sapply(seq(1, kw.Q-lag.par-1), function(comp.ind){
      yyy <- 1
      counter <- 0
      for (i in ijk) {
        counter <- counter + 1
        if (counter == length(ijk)){
          yyy <- yyy * kw.y[comp.ind+lag.par+1, i]
        }else if (counter == length(ijk)-1) {
          yyy <- yyy * kw.y[comp.ind+lag.par, i]
        }else{
          yyy <- yyy * kw.y[comp.ind, i]
        }
      }
      return (yyy)  
    }))
  }else{
    result <- sum(sapply(seq(2-lag.par, kw.Q), function(comp.ind){
      yyy <- 1
      counter <- 0
      for (i in ijk) {
        counter <- counter + 1
        if (counter == length(ijk)){
          yyy <- yyy * kw.y[comp.ind+lag.par-1, i]
        }else if (counter == length(ijk)-1) {
          yyy <- yyy * kw.y[comp.ind+lag.par, i]
        }else{
          yyy <- yyy * kw.y[comp.ind, i]
        }
      }
      return (yyy) 
    }))
  }
  return (result)
}

Find.aii <- function(ijk.mat, kw.M, kw.Q, kw.n, kw.q, kw.y){
  
  t1a <- 1/3 * sum(apply(ijk.mat, 1, function(ind) {
    ikj <- c(ind[1], ind[3], ind[2])
    kji <- c(ind[3], ind[2], ind[1])
    kw.M^2/kw.Q * Prod.sum(ind, lag.par = 0, kw.Q, kw.y) * kw.M^2/kw.Q * ( Prod.sum(ind, lag.par = 0, kw.Q, kw.y) + Prod.sum(ind, lag.par=1, kw.Q, kw.y) + Prod.sum(ikj, lag.par=1, kw.Q, kw.y) + Prod.sum(kji, lag.par=1, kw.Q, kw.y) + Prod.sum(ind, -1, kw.Q, kw.y) + Prod.sum(ikj, -1, kw.Q, kw.y) + Prod.sum(kji, -1, kw.Q, kw.y)) 
  }))
  
  t1b <- 3/8 * sum(apply(ijk.mat, 1, function(ind) {
    ijk = ind
    kji = c(ind[3], ind[2], ind[1])
    kw.M^2/kw.Q * ( Prod.sum(ijk, lag.par=0, kw.Q, kw.y) + Prod.sum(ijk, lag.par = 1, kw.Q, kw.y) + Prod.sum(ijk, lag.par = -1, kw.Q, kw.y) ) * kw.M^2/kw.Q * ( Prod.sum(kji, lag.par=0, kw.Q, kw.y) + Prod.sum(kji, lag.par = 1, kw.Q, kw.y) + Prod.sum(kji, lag.par = -1, kw.Q, kw.y) ) 
  })) - 5/6 * sum(apply(ijk.mat, 1, function(ind) {
    kw.M^2/kw.Q * Prod.sum(ind, lag.par = 0, kw.Q, kw.y) * kw.M^2/kw.Q * ( Prod.sum(ind, lag.par=0, kw.Q, kw.y) + Prod.sum(ind, lag.par = 1, kw.Q, kw.y) + Prod.sum(ind, lag.par = -1, kw.Q, kw.y) ) 
  })) - 5/6 * sum(apply(ijk.mat, 1, function(ind) {
    kji = c(ind[3], ind[2], ind[1])
    kw.M^2/kw.Q * Prod.sum(ind, lag.par = 0, kw.Q, kw.y) * kw.M^2/kw.Q * ( Prod.sum(kji, lag.par=0, kw.Q, kw.y) + Prod.sum(kji, lag.par = 1, kw.Q, kw.y) + Prod.sum(kji, lag.par = -1, kw.Q, kw.y) ) 
  })) + 8/9 * sum(apply(ijk.mat, 1, function(ind) {
    kw.M^2/kw.Q * Prod.sum(ind, lag.par = 0, kw.Q, kw.y) * kw.M^2/kw.Q * Prod.sum(ind, lag.par = 0, kw.Q, kw.y)
  }))
  
  t1c <- -5/12 * sum(apply(ijk.mat, 1, function(ind) {
    kw.M^2/kw.Q * Prod.sum(ind, lag.par = 0, kw.Q, kw.y) * kw.M^2/kw.Q * (Prod.sum(ind, lag.par=0, kw.Q, kw.y) + Prod.sum(ind, lag.par = 1, kw.Q, kw.y) +   Prod.sum(ind, lag.par = -1, kw.Q, kw.y) ) 
  })) + 2/9 * sum(apply(ijk.mat, 1, function(ind) {
    kw.M^2/kw.Q * Prod.sum(ind, 0, kw.Q, kw.y) * kw.M^2/kw.Q * Prod.sum(ind, 0, kw.Q, kw.y)
  }))
  
  
  t2a <- 3/8 * sum(apply(ijk.mat, 1, function(ind) {
    ikk = c(ind[1], ind[3], ind[3])
    ill = c(ind[1], ind[2], ind[2])
    kw.M^2/kw.Q * (Prod.sum(ikk, 0, kw.Q, kw.y) + Prod.sum(ikk, 1, kw.Q, kw.y) + Prod.sum(ikk, -1, kw.Q, kw.y)) * kw.M^2/kw.Q * (Prod.sum(ill, 0, kw.Q, kw.y) + Prod.sum(ill, 1, kw.Q, kw.y) +   Prod.sum(ill, -1, kw.Q, kw.y)) 
  })) + 4/9 * sum(apply(ijk.mat, 1, function(ind) {
    iil <- c(ind[1], ind[1], ind[2])
    lkk <- c(ind[2], ind[3], ind[3])
    kw.M^2/kw.Q * Prod.sum(iil, 0, kw.Q, kw.y) * kw.M^2/kw.Q * Prod.sum(lkk, 0, kw.Q, kw.y)
  })) - 5/6 * sum(apply(ijk.mat, 1, function(ind) {
    iik <- c(ind[1], ind[3], ind[3])
    kll <- c(ind[3], ind[2], ind[2])
    kw.M^2/kw.Q * Prod.sum(iik, 0, kw.Q, kw.y) * kw.M^2/kw.Q * ( Prod.sum(kll,0, kw.Q, kw.y) +   Prod.sum(kll, 1, kw.Q, kw.y) + (kw.Q-1)/kw.Q* Prod.sum(kll, -1, kw.Q, kw.y))
  }))
  
  t2b <- 1/4 * sum(apply(ijk.mat, 1, function(ind) {
    ikk = c(ind[1], ind[2], ind[2])
    ill = c(ind[1], ind[3], ind[3])
    kw.M^2/kw.Q * (Prod.sum(ikk, 0, kw.Q, kw.y) +   Prod.sum(ikk, 1, kw.Q, kw.y) +   Prod.sum(ikk, -1, kw.Q, kw.y)) * kw.M^2/kw.Q * (Prod.sum(ill, 0, kw.Q, kw.y) +   Prod.sum(ill, 1, kw.Q, kw.y) +   Prod.sum(ill, -1, kw.Q, kw.y))
  })) - 1/3 * sum(apply(ijk.mat, 1, function(ind) {
    ikk <- c(ind[1], ind[3], ind[3])
    ill <- c(ind[1], ind[2], ind[2])
    kw.M^2/kw.Q* Prod.sum(ikk, 0, kw.Q, kw.y) * kw.M^2/kw.Q * ( Prod.sum(ill,0, kw.Q, kw.y) +   Prod.sum(ill, 1, kw.Q, kw.y) + (kw.Q-1)/kw.Q* Prod.sum(ill, -1, kw.Q, kw.y))
  })) + 1/9 * sum(apply(ijk.mat, 1, function(ind) {
    ikk <- c(ind[1], ind[3], ind[3])
    ill <- c(ind[1], ind[2], ind[2])
    kw.M^2/kw.Q * Prod.sum(ikk, 0, kw.Q, kw.y) * kw.M^2/kw.Q * Prod.sum(ill, 0, kw.Q, kw.y)
  }))
  
  t3a <- -1/2 * sum(apply(ijk.mat, 1, function(ind) {
    ikki <- c(ind[1], ind[3], ind[3], ind[1])
    ikik <- c(ind[1], ind[3], ind[1], ind[3])
    kw.M^3/kw.Q * ( Prod.sum2(ikki, 0, kw.Q, kw.y) +   Prod.sum(ikki, 1, kw.Q, kw.y) +   Prod.sum(ikik, 1, kw.Q, kw.y) +   Prod.sum2(ikki, 1, kw.Q, kw.y) +   Prod.sum(ikki, -1, kw.Q, kw.y) +   Prod.sum(ikik, -1, kw.Q, kw.y) +   Prod.sum2(ikki, -1, kw.Q, kw.y) + (kw.Q-2)/kw.Q * Prod.sum2(ikki, 2, kw.Q, kw.y) + (kw.Q-2)/kw.Q * Prod.sum2(ikki, -2, kw.Q, kw.y) )
  }))
  
  t3b <- 3/8 * sum(apply(ijk.mat, 1, function(ind) {
    ikki <- c(ind[1], ind[3], ind[3], ind[1])
    kw.M^3/kw.Q * (Prod.sum2(ikki, 0, kw.Q, kw.y) + 2* (kw.Q-1)/kw.Q*Prod.sum2(ikki, 1, kw.Q, kw.y) )
  })) + sum(apply(ijk.mat, 1, function(ind) {
    iill <- c(ind[1], ind[1], ind[2], ind[2])
    kw.M^2/kw.Q * (Prod.sum(iill, 0, kw.Q, kw.y) +   Prod.sum(iill, 1, kw.Q, kw.y) +   Prod.sum(iill, -1, kw.Q, kw.y)) })) - 3/4 * sum(apply(ijk.mat, 1, function(ind) {
      iikk <- c(ind[1], ind[1], ind[3], ind[3])
      kw.M^3/kw.Q * Prod.sum(iikk, 0, kw.Q, kw.y)
    }))
  
  t3c <- 1/4 * sum(apply(ijk.mat, 1, function(ind) {
    ikki <- c(ind[1], ind[3], ind[3], ind[1])
    kw.M^3/kw.Q * (Prod.sum2(ikki, 0, kw.Q, kw.y) + 2* (kw.Q-1)/kw.Q*Prod.sum2(ikki, 1, kw.Q, kw.y) 
    )}))
  
  
  aii <- 1/kw.q * (2*t1a + 2*t1b + t1c + 2*t2a + t2b + 2*t3a + 2*t2b + t3c)
  
  return (aii/2 * kw.Q/kw.n/kw.q)
}

Alpha.hat <- function (yi, rst) {
  return (mean(yi^sum(rst)))
}

Alpha.rr <- function (s.size, yi, rr) {
  return (s.size/(s.size-1) * Alpha.hat(yi, rr))
}

Alpha.rrrr <- function (s.size, yi, rrrr) {
  return (1/(s.size-4) * (s.size * Alpha.hat(yi, rrrr) - 6 * Alpha.rr(s.size, yi, rrrr[1:2])^2))
}


Alpha.rr2 <- function (s.size, yi, rrrr) {
  return (Alpha.rr(s.size, yi, rrrr[1:2])^2 - 1/s.size * Alpha.rrrr(s.size, yi, rrrr))
}

Alpha.rrr <- function (s.size, yi, rrr) {
  return (s.size/(s.size-3) * Alpha.hat(yi, rrr))
}

Alpha.rrr2 <- function (s.size, yi, rrrrrr) {
  (1+1/s.size) * Alpha.rrr(s.size, yi, rrrrrr[1:3])^2 - 1/s.size * Alpha.hat(yi, rrrrrr)
}

Alpha.rr3 <- function (s.size, yi, rr) {
  Alpha.rr(s.size, yi, rr)^3
}


Alpha.rrss <- function (s.size, yi, rrss) {
  1/(s.size-4) * (s.size * Alpha.hat(yi, rrss) - 2 * Alpha.rr(s.size, yi, rrss[1:2]) * Alpha.rr(s.size, yi, rrss[3:4]))
}

Alpha.rr.Alpha.ss <- function (s.size, yi, rrss) {
  Alpha.rr(s.size, yi, rrss[1:2]) * Alpha.rr(s.size, yi, rrss[3:4]) - 1/s.size * Alpha.rrss(s.size, yi, rrss)
}

Alpha.rss2 <- function (s.size, yi, rrssss) {
  (1+1/s.size)*Alpha.rrr(s.size, yi, rrssss[2:4])^2 - 1/s.size * (Alpha.hat(yi, rrssss))
}

Alpha.rr.Alpha.ss2 <- function (s.size, yi, rrsstt) {
  Alpha.rr(s.size, yi, rrsstt[1:2]) * Alpha.rr(s.size, yi, rrsstt[3:4]) * Alpha.rr(s.size, yi, rrsstt[5:6])
}

b11.sumand <- function (r, yi, s.size) {
  rrrrrr <- rep(r, 6)
  return (Alpha.rrrr(s.size, yi, rrrrrr[1:4])/(2*Alpha.rr2(s.size, yi, rrrrrr[1:4])) - Alpha.rrr2(s.size, yi, rrrrrr)/(3*Alpha.rr3(s.size, yi, rrrrrr[1:2])))
}

b11 <- function (d, yi, s.size) {
  mean( sapply( seq(1,d), function(ind){
    b11.sumand(ind, yi, s.size)
  } ) )
}

b12.sumand <- function (rs, yi, s.size) {
  rrssss <- c(rep(rs[1],2), rep(rs[2], 4))
  Alpha.rrss(s.size, yi, rrssss[1:4])/Alpha.rr.Alpha.ss(s.size, yi, rrssss[1:4]) - Alpha.rss2(s.size, yi, rrssss)/Alpha.rr.Alpha.ss2(s.size, yi, rrssss)
}

b12 <- function (rs.ind, yi, s.size) {
  mean(apply(rs.ind, 2, function(ind){
    b12.sumand(ind, yi, s.size)
  }))
}

b21.sumand <- function (rs, yi, s.size) {
  rrssss <- c(rep(rs[1],2), rep(rs[2], 4))
  Alpha.rss2(s.size, yi, rrssss)/Alpha.rr.Alpha.ss2(s.size, yi, rrssss)
}

b21 <- function (rs.ind, yi, s.size) {
  mean(apply(rs.ind, 2, function(ind) {
    b21.sumand(ind, yi, s.size)
  }))
}

b22.sumand <- function (rst, yi, s.size) {
  rrsstt <- c(rep(rst[1],2), rep(rst[2],2), rep(rst[3],2))
  Alpha.rss2(s.size, yi, rrsstt)/Alpha.rr.Alpha.ss2(s.size, yi, rrsstt)
}

b22 <- function (rst.ind, yi, s.size) {
  mean(apply(rst.ind, 2, function(ind) {
    b22.sumand(ind, yi, s.size)
  })) * 2
}

b1 <- function (d, yi, s.size) {
  b11(d, yi, s.size) + b12(combn(d, 2), yi, s.size)
}

b2 <- function (d, yi, s.size) {
  if (d==2) {
    return ( b21(combn(d,2), yi, s.size) )
  }else{
    return ( b21(combn(d, 2), yi, s.size) + b22(combn(d, 3), yi, s.size))
  }
}

Boot.replication <- function(ijk.mat, kw.M, kw.Q, kw.n, kw.q, kw.y, bt.blocks, bt.num.blocks){
    bpf.sample.group.index <- sample(seq(1, bt.num.blocks), bt.num.blocks, replace = TRUE)
  bpf.sample.index <- c(bt.blocks[,bpf.sample.group.index])  # bootstrap sample
  bpf.y <- kw.y[bpf.sample.index,]
  
  Find.aii(ijk.mat, kw.M, bt.num.blocks, kw.n, kw.q, bpf.y)
  
}

Boot.plug.formula <- function (ijk.mat, kw.M, kw.Q, kw.n, kw.q, kw.y, bt.blocks, bt.num.blocks, bt.B, kw.a) {
  bpf.a <- replicate(bt.B, Boot.replication(ijk.mat, kw.M, kw.Q, kw.n, kw.q, kw.y, bt.blocks, bt.num.blocks))
  
  bpf.bias <- mean(bpf.a) - kw.a
  
  bpf.tune <- kw.a - bpf.bias
  
  bpf.sd <- sd(bpf.a)
  
  return (c(bpf.tune, bpf.sd, bpf.bias))
}


Extra.Point <- function (kw.a, new.point) {
  if (kw.a > 0){
    return( -kw.a*new.point )
  }else{
    return( c(-kw.a*new.point, (kw.a+2)*new.point) ) 
  }
}

ABELregression <- function (Y, X, null.coef, 
                            kw.q = nrow(null.coef),
                            block.length=floor(sqrt(nrow(Y))), 
                            bt.block.length=2, 
                            bt.block.overlap=2,
                            bt.B=200,
                            pro.block=TRUE) {
  library(emplik)
  library(gtools)
  est.fun <- t(X) * matrix(rep((Y-X%*%null.coef)[,1], ncol(X)), nrow = ncol(X), byrow = TRUE)

  if (pro.block == TRUE) {
    kw.n <- nrow(Y)-1
    
    l1 <- floor((sqrt(4*kw.n+1)-1)/2)
    l2 <- ceiling((sqrt(4*kw.n+1)-1)/2)
    
    if (kw.n-l1*(l1+1) < l2*(l2+1) - kw.n) {
      kw.Q <- l1
    }else{
      kw.Q <- l2
    }
    
    block.list <- list()
    for (i in 1:kw.Q){
      block.list[[i]] <- seq((i-1)*i+1, i*(i+1))
    }
    
    ti <- Ti.fun.pro(t(est.fun), block.list, kw.n)
    kw.scale <- floor(sqrt(kw.n))
    kw.M <- floor(sqrt(kw.n))
    
  }else{
    blocks.info <- Blocks(ncol(est.fun), block.length, block.length)
    kw.Q <- blocks.info[[2]]
    kw.M <- block.length
    kw.n <- nrow(Y)-1
    
    ti <- Ti.fun(t(est.fun), blocks.ind = blocks.info[[1]], Q = kw.Q,
                 M = block.length)
    kw.scale <- block.length
  }
  
  V.hat <- 1/kw.Q * (kw.Q -1) * var(ti)
  V.hat.half <- Matrix.power(kw.scale * V.hat, -1/2)
  kw.y <- t(apply(ti, 1, function(ind){ V.hat.half %*% ind}))
  
  kw.y.bar <- apply(kw.y, 2, mean)
  
  maha.y.bar <- sqrt(t(kw.y.bar) %*% solve((kw.Q-1)/kw.Q*var(kw.y)) %*% kw.y.bar)[1,1]
  
  bart.new.point <- kw.y.bar / (1+0.1*maha.y.bar^2)
  
  ijk.mat <- permutations(n=kw.q, r=3, v=seq(1,kw.q), repeats.allowed = T)
  kw.a <- Find.aii(ijk.mat, kw.M, kw.Q, kw.n, kw.q, kw.y)
  
  bt.blocks <- Blocks(kw.Q, bt.block.length, bt.block.overlap)[[1]]
  bt.num.blocks <- Blocks(kw.Q, bt.block.length, bt.block.overlap)[[2]]
  
  bt.bias.sd <- Boot.plug.formula(ijk.mat, kw.M, kw.Q, kw.n, kw.q, kw.y, 
                                  bt.blocks, bt.num.blocks, bt.B, kw.a)
  
  bt.a <- bt.bias.sd[1]

  if (bt.a >0){add.num <- 1}else{add.num <- 2}
  kw.Ti.adj.bart <- rbind( kw.y, matrix(Extra.Point(bt.a, bart.new.point),
                                        nrow = add.num, byrow = T) )
  
  abel.ratio <- el.test(kw.Ti.adj.bart, rep(0, kw.q))$"-2LLR"
  
  bel.ratio <- el.test(ti, rep(0, kw.q))$"-2LLR"
  
  
  return (list(abel=round(abel.ratio,3), bel=round(bel.ratio,3)))
  
}


# Example
library(MASS)
b <- c(1,2,3,4,0)

set.seed(100)
x <- mvrnorm(40, rep(0,5), diag(rep(1,5)))
y <- x %*% b + rchisq(40,2)

# Test H_0: b_5 = 0 v.s. H_1: b_5 != 0
# \hat{b} under H_0
b.hat <- matrix(c(lm(y~x[,1:4])$coefficients[-1], 0))
# Find ABEL and BEL test statistic and compare to qchisq(0.95, 1)
ABELregression(y, x, b.hat)
