align_quantiles = function(P1,P2,quantiles=100)
{
  dims = nrow(P1)
  for (i in 1:dims)
  {
    warp_func = approxfun(quantile(P2[i,],seq(0,1,by=1/quantiles)),quantile(P1[i,],seq(0,1,by=1/quantiles)))
    P2[i,] = warp_func(P2[i,])
  }
  return(list(P1,P2))
}