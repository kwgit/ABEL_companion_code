# ABEL simulation

# Functions----

Ar1.fun <- function (rho, inova) {
  kw.z <- inova
  for (i in 2:dim(inova)[1]) {
    kw.z[i,] <- rho %*% kw.z[i-1, ] + inova[i, ]
  }
  return (kw.z)
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


Ti.fun.pro <- function (est_fun, block.list) {
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

Boot.belr.md <- function(ori.data, bt.blocks, kw.blocks, kw.Q, kw.M, tune, kw.mean) {

  bt.block.ind <- sample(seq(1, bt.num.blocks), bt.sample.num.blocks)

  bt.data.ind <- c(bt.blocks[,bt.block.ind])

  bt.data <- ori.data[bt.data.ind, ]
  

  bt.Ti <- Ti.fun(bt.data, kw.blocks, kw.Q, kw.M) 
  
  bt.mu.hat <- apply(bt.data, 2, ave)
  
  bt.Ti.hat <- Ti.fun(bt.data - bt.mu.hat, kw.blocks, kw.Q, kw.M)
  
  V.hat <- 1/kw.Q * (kw.Q-1) * var(bt.Ti.hat)
  
  
  V.hat.half <- Matrix.power(kw.M * V.hat, -1/2)
  bt.y <- t(apply(bt.Ti, 1, function(ind){ V.hat.half %*% ind}))
  bt.y.mu.hat <- apply(bt.y, 2, ave)
  
  maha.y.bar <- try(sqrt(t(bt.y.mu.hat[1,]) %*% solve(var(bt.y)) %*% bt.y.mu.hat[1,])[1])
  
  pre.new.point <- apply(bt.y, 2, mean) * 1/(1+0.1* maha.y.bar^2)
  
  new.point <- -tune * pre.new.point
  
  bt.y.adj.sudo <- rbind(bt.y, new.point)
  
  est.elr <- el.test(bt.y.adj.sudo, kw.mean)$"-2LLR"
  
  return (est.elr)
}

Boot.tune <- function(tune.data, bt.blocks, kw.blocks, kw.Q, kw.M, tune.par, kw.mean,boot.num){
  ini.elrs <- rep(0, boot.num)
  for (i in 1:boot.num){
    ini.elrs[i] <- Boot.belr.md(tune.data, bt.blocks, kw.blocks, kw.Q, kw.M, tune.par, kw.mean)
  }
  
  ini.bart <- abs(mean(ini.elrs)-kw.q)
  
  ini.tune <- ini.bart * kw.Q /2 /kw.q
  
  return (ini.tune)
}


Boot.belr.1d <- function(ori.data, bt.blocks, kw.blocks, kw.Q, kw.M) {
  bt.block.ind <- sample(seq(1, bt.num.blocks), bt.sample.num.blocks)
  bt.data.ind <- c(bt.blocks[,bt.block.ind])
  bt.data <- ori.data[bt.data.ind]
  
  bt.Ti <- Ti.fun(bt.data, kw.blocks, kw.Q, kw.M)
  
  bt.mu.hat <- ave(bt.data)
  
  bt.Ti.hat <- Ti.fun(bt.data - bt.mu.hat, kw.blocks, kw.Q, kw.M)
  
  V.hat <- 1/kw.Q * (kw.Q-1) * var(bt.Ti.hat)
  
  maha.Ti.bar <- try(sqrt(t(bt.mu.hat[1]) %*% solve(V.hat) %*% bt.mu.hat[1])[1])
  
  V.hat.half <- Matrix.power(kw.M * V.hat, -1/2)
  
  bt.y <- sapply(bt.Ti, function(ind){ V.hat.half %*% ind})
  
  new.point <- mean(bt.y) * 1/(1+0.1* maha.Ti.bar^2)
  
  bt.y.adj_logn <- c(bt.y, Extra.Point(log(kw.n)/2, new.point))
  
  est.elr <- el.test(bt.y.adj_logn, kw.mu)$"-2LLR"
  
  return (est.elr)
}


Boot.replication <- function(ijk.mat, kw.M, kw.Q, kw.n, kw.q, kw.y, bt.blocks, bt.num.blocks){
  
  bpf.sample.group.index <- sample(seq(1, bt.num.blocks), bt.num.blocks, replace = TRUE)
  bpf.sample.index <- c(bt.blocks[,bpf.sample.group.index])
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

Matrix.power <- function (mat, pow) {
  mat.eigVal <- eigen(mat)$values
  mat.eigVec <- eigen(mat)$vectors
  
  mat.pow <- mat - mat
  for (i in 1:length(mat.eigVal)) {
    mat.pow = mat.pow + mat.eigVal[i]^pow * mat.eigVec[,i] %*% t(mat.eigVec[,i])
  }
  
  return (mat.pow)
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

# Simulation----
library(MASS)
library(emplik)
library(parallel)
library(gtools)

args <- commandArgs(TRUE)

args = as.integer(args)
print(args)

cores <- 64
num.runs <- 1000

sample.size <- c(100, 400)
proce.coef <- c(-0.8, -0.5, -0.2, 0.2, 0.5, 0.8) 
para.dim <- c(2, 3, 4, 5, 10)
job.dim <- length(sample.size) * length(proce.coef) * length(para.dim)

job.matrix <- matrix( rep(0, job.dim*3), ncol = 3)
row.count <- 0
for (i in sample.size){
  for (j in proce.coef) {
    for (k in para.dim) {
      row.count = row.count + 1
      job.matrix[row.count,] = c(i, j, k)
    }
  }
}

(job.id <- job.matrix[args,])

kw.q <- job.id[3] 
kw.m <- rep(0, kw.q)
kw.var <- diag(kw.q)
kw.rho <- diag(kw.q) * job.id[2]
ijk.mat <- permutations(n=kw.q, r=3, v=seq(1,kw.q), repeats.allowed = T)
kw.n = job.id[1]

bt.block.length <- 2
bt.block.overlap <- 2
bt.B <- 200

Sim.fun <- function (dumy.arg) {
  
  kw.inova <- mvrnorm(kw.n + 200, kw.m, kw.var)
  ar1 <- Ar1.fun(kw.rho, kw.inova)[201:(kw.n + 200), ]
  elr.matrix <- matrix(rep(0, 10*16), ncol = 16)
  
  block.lengths <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
  
  for (kw.M in block.lengths) {
    
    if (kw.M == 17) {
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
      
      kw.Ti <- Ti.fun.pro(ar1, block.list)
      
    }else{
      kw.L = kw.M 
      kw.blocks <- Blocks(kw.n, kw.M, kw.L)[[1]]
      kw.Q <- dim(kw.blocks)[2]
      
      kw.Ti <- Ti.fun(ar1, kw.blocks, kw.Q, kw.M)
    }
    
    if (is.null(dim(ar1))) {
      kw.mu.hat <- ave(ar1)
      if (kw.M == 17) {
        kw.Ti.hat <- Ti.fun.pro(ar1-kw.mu.hat, block.list)
        kw.scale <- floor(sqrt(kw.n))
      }else{
        kw.Ti.hat <- Ti.fun(ar1-kw.mu.hat, kw.blocks, kw.Q, kw.M)
        kw.scale <- kw.M
      }
      
      V.hat <- 1/kw.Q * (kw.Q-1) * var(kw.Ti.hat)
      
      maha.Ti.bar <- (t(kw.mu.hat[1]) %*% (1/V.hat) %*% kw.mu.hat[1])[1]
      
      new.point <- kw.mu.hat[1] * 1/(1+0.1* maha.Ti.bar)
      
      kw.Ti.adj.log <- c( kw.Ti, Extra.Point(log(kw.n)/2, new.point) ) 
      
      kw.Ti.adj.a05 <- c( kw.Ti, Extra.Point(0.5, new.point) ) 
      
      kw.Ti.adj.a08 <- c( kw.Ti, Extra.Point(0.8, new.point) )
      
      kw.Ti.adj.a1 <- c( kw.Ti, Extra.Point(1, new.point) )
    
      V.hat.half <- Matrix.power(kw.scale * V.hat, -1/2)
      
      kw.y <- matrix(kw.Ti * V.hat.half[1,1])
      
      kw.y.bar <- apply(kw.y, 2, mean)
      
      bt.blocks <- Blocks(kw.Q, bt.block.length, bt.block.overlap)[[1]]
      
      bt.num.blocks <- Blocks(kw.Q, bt.block.length, bt.block.overlap)[[2]]
      
      bart.new.point <- kw.y.bar / (1+0.1*maha.Ti.bar)
      
      kw.a <- kw.n^(-2/3)*kw.Q*kw.q
      
      kw.Ti.adj.bart <- c( kw.y, Extra.Point(kw.a, bart.new.point) )
      
    }else{
      kw.mu.hat <- apply(ar1, 2, ave)
      
      if (kw.M == 17) {
        kw.Ti.hat <- Ti.fun.pro(ar1-kw.mu.hat, block.list)
        kw.scale <- floor(sqrt(kw.n))
      }else{
        kw.Ti.hat <- Ti.fun(ar1-kw.mu.hat, kw.blocks, kw.Q, kw.M)
        kw.scale <- kw.M
      }

      V.hat <- 1/kw.Q * (kw.Q-1) * var(kw.Ti.hat)
      
      maha.Ti.bar <- try(sqrt(t(kw.mu.hat[1,]) %*% solve(V.hat) %*% kw.mu.hat[1,])[1])
      
      if (class(maha.Ti.bar) != "try-error") {
        V.hat.half <- Matrix.power(kw.scale * V.hat, -1/2)
        
        kw.y <- t(apply(kw.Ti, 1, function(ind){ V.hat.half %*% ind}))
        
        kw.y.bar <- apply(kw.y, 2, mean)
        
        maha.y.bar <- sqrt(t(kw.y.bar) %*% solve((kw.Q-1)/kw.Q*var(kw.y)) %*% kw.y.bar)[1,1]
        
        bart.new.point <- kw.y.bar / (1+0.1*maha.y.bar^2)
        
        kw.a <- Find.aii(ijk.mat, kw.M, kw.Q, kw.n, kw.q, kw.y)
        
        bt.blocks <- Blocks(kw.Q, bt.block.length, bt.block.overlap)[[1]]
      
        bt.num.blocks <- Blocks(kw.Q, bt.block.length, bt.block.overlap)[[2]]
        
        bt.bias.sd <- Boot.plug.formula(ijk.mat, kw.M, kw.Q, kw.n, kw.q, kw.y, bt.blocks, bt.num.blocks, bt.B, kw.a)
        
        bt.a <- bt.bias.sd[1]
        
        new.point <- kw.mu.hat[1,] * 1/(1+0.1* maha.Ti.bar^2) 
        
        kw.Ti.adj.log <- rbind( kw.Ti, matrix(Extra.Point(log(kw.n)/2, new.point), nrow = 1, byrow = T) )
        
        kw.Ti.adj.a05 <- rbind( kw.Ti, matrix(Extra.Point(0.5, new.point), nrow = 1, byrow = T) )
        
        kw.Ti.adj.a08 <- rbind( kw.Ti, matrix(Extra.Point(0.8, new.point), nrow = 1, byrow = T) )
        
        kw.Ti.adj.a1 <- rbind( kw.Ti, matrix(Extra.Point(1, new.point), nrow = 1, byrow = T) )
        
        if (bt.a >0){add.num <- 1}else{add.num <- 2}
        kw.Ti.adj.bart <- rbind( kw.y, matrix(Extra.Point(bt.a, bart.new.point), nrow = add.num, byrow = T) )       
        
      }else{
        kw.Ti.adj.log <- 0
        
        kw.Ti.adj.a05 <- 0
        
        kw.Ti.adj.a08 <- 0
        
        kw.Ti.adj.a1 <- 0
        
        kw.Ti.adj.bart <- 0
      }
    } 
    
    belr <- try( el.test(kw.Ti, kw.m)$'-2LLR')
    if (class(belr) == "try-error") {belr = 0}
    abelr.log <- try(el.test(kw.Ti.adj.log, kw.m)$"-2LLR")
    if (class(abelr.log) == "try-error") {abelr.log = 0}
    abelr.a05 <- try(el.test(kw.Ti.adj.a05, kw.m)$"-2LLR")
    if (class(abelr.a05) == "try-error") {abelr.a05 = 0}
    abelr.a08 <- try(el.test(kw.Ti.adj.a08, kw.m)$"-2LLR")
    if (class(abelr.a08) == "try-error") {abelr.a08 = 0}
    abelr.a1 <- try(el.test(kw.Ti.adj.a1, kw.m)$"-2LLR")
    if (class(abelr.a1) == "try-error") {abelr.a1 = 0}
    abelr.bart <- try(el.test(kw.Ti.adj.bart, kw.m)$"-2LLR")
    if (class(abelr.bart) == "try-error") {abelr.bart = 0}
    try(elr.matrix[,(kw.M-1)] <- c(belr, abelr.log, abelr.a05, abelr.a08, abelr.a1, abelr.bart, kw.a, bt.a, bt.bias.sd[2],bt.bias.sd[3]))
  }
  
  return (elr.matrix)
}

llr.list.blocks <- mclapply(seq(1,num.runs), Sim.fun, mc.cores = cores)
llr.blocks <- do.call(rbind, llr.list.blocks) 

for (i in 1:dim(llr.blocks)[2]) {
  
  name.block <- i + 1
  
  llr <- matrix(llr.blocks[,i], ncol = 10, byrow = T)
  
  cov.prob <- sapply(c(0.75, 0.85, 0.9, 0.95, 0.99), function(ind){ 
    chi.quantile <- qchisq(ind, kw.q)
    return ( c(sum(llr[,1] < chi.quantile), sum(llr[,2] < chi.quantile), sum(llr[,3] < chi.quantile),sum(llr[,4] < chi.quantile), sum(llr[,5] < chi.quantile), sum(llr[,6] < chi.quantile) )/num.runs )
  })
  
  cov.prob <- data.frame(cov.prob)
  names(cov.prob) <- c('75%', '85%', '90%', '95%', '99%')
  bel_method <- c('BEL', 'ABEL_log', 'ABEl_0.5', 'ABEL_0.8', 'ABEL_1', 'ABEL_bart')
  cov.prob <- cbind(bel_method, cov.prob)
  
  write.csv(cov.prob, file = paste("Coverage_Prob_AR1", job.id[1], job.id[2], job.id[3], name.block, '.csv', sep = '_') )
  
  llr <- data.frame(llr)
  names(llr) <- c('BEL', 'ABERL_log', 'ABEL_0.5', 'ABEL_0.8', 'ABEL_1', 'ABEL_bart', 'Plugin_tuning', 'Bias_corrected_tuning', 'sd_plugin', 'bias_plugin')
  
  write.csv(llr, file = paste("Data_AR1", job.id[1], job.id[2], job.id[3], name.block, ".csv", sep = '_'))
  
}


