INLA_run = function(virus, wmel_proportion, in.project, prediction_type){

  grid_res = 0.5 #km


  if(virus == "dengue"){
    case.grid = dengue.grid
  }
  
  if(virus == "chik"){
    case.grid = chik.grid
    t = (1 + ncol(dengue.grid) - ncol(chik.grid)):ncol(dengue.grid)
    aegypti.captured.grid = aegypti.captured.grid[,t]
    in.project.grid = in.project.grid[,t]
    release.grid = release.grid[,t]
    wMel.grid = wMel.grid[,t]
    }
  
  #Remove places where Population=0 from the analysis
  case.grid = case.grid[covariates.df$Population>0,]
  in.project.grid = in.project.grid[covariates.df$Population>0,]
  wMel.grid = wMel.grid[covariates.df$Population>0,]
  covariates.df = covariates.df[covariates.df$Population>0,]
  

  #Create the input data frame-------------------------
  t = 1:ncol(case.grid)
  month = ((t-1)%%12)+1
  
  df.st = data.frame(obs = as.vector(case.grid[,t]),
                     x = rep(covariates.df$x, length(t)),
                     y = rep(covariates.df$y, length(t)),
                     t = rep(t,each=(nrow(covariates.df))),
                     month = rep(month,each=(nrow(covariates.df))),
                     In.project = as.double(as.vector(in.project.grid[,t])),
                     Wolbachia = as.vector(wMel.grid[,t]),
                     Wolbachia.bin = as.vector(wMel.grid.bin[,t]),
                     Population = rep(covariates.df$Population, length(t)))
  
  #Held_out data----------------
  #Full set - All cases areused for data prediciton
  if(prediction_type==1){ 
  } 
  #Half set - 50% of cases are removed at random
  if(prediction_type==2){ 
    out = sample(nrow(df.st), size=round(0.5*nrow(df.st)), replace=F)
    df.st$obs[out] = NA
  }

  #Spatial set - whole years removed(20%)
  if(prediction_type==3){ 
    radius_chunk=3/grid_res #Approx.
    out = c()
    nstep.in.year  = 12/(1+1*tm)
    for(i in 1:(ncol(case.grid)/nstep.in.year)){
      total_chunks = round(0.2*nrow(covariates.df)/(pi*radius_chunk^2))
      ratio.check = 0
      while(ratio.check<0.195){ #We want to approx. hold out 20% of the data per time step
        kern = sample(nrow(covariates.df), size=total_chunks)
        kern.coords = covariates.df[kern,c("x","y")]
        kern.neighbours = get.knnx(covariates.df[,c("x","y")], kern.coords, (pi*radius_chunk^2)+1)
        ratio.check = sum((1:nrow(covariates.df) %in% kern.neighbours$nn.index))/nrow(covariates.df)
        total_chunks = total_chunks+1
      }
      out = c(out, rep((1:nrow(covariates.df) %in% kern.neighbours$nn.index), nstep.in.year))
    }
    
    df.st$obs[out] = NA
  } 
  
  
  #Creation of the mesh and the SPDE----------
  spatial.range = 22000

  mesh = inla.mesh.2d(loc=covariates.df[,c("x","y")], max.edge=c(spatial.range/5, spatial.range),
                           cutoff = 1500,
                           offset=c(spatial.range/5, diff(range(covariates.df$x))/3))
  
  spde = inla.spde2.matern(mesh=mesh,alpha=2)
  
  #Formula-------------
  if(wmel_proportion){
    diag.param = 1e-1 #increasing "diagonal" parameter to avoid numerical singularity in wMel.grid - check https://groups.google.com/forum/#!topic/r-inla-discussion-group/8zFb6oi9O-U for more detail
    length.wb = length(levels(as.factor(inla.group(wMel.grid.bin, method='cut', n=100)))) - 1
    
    formula = obs ~ -1 + Intercept + offset(log(Population)) +  f(t, model = 'ar1') +   f(spatial.field, model = spde) +
      f(inla.group(Wolbachia.bin, method='cut', n=100), model=rand_model, diagonal=diag.param, constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length.wb)), nrow=1), e=rep(0,1)))
  }
  
  if(in.project){
    
    formula = obs ~ -1 + Intercept + offset(log(Population)) +  f(t, model = 'ar1') + f(spatial.field, model = spde) + In.project
    
  }


  #Stack creation and inla analysis-------------------------------
  s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde) 
  
  A.est = inla.spde.make.A(mesh=mesh,
                           loc=as.matrix(df.st[,c("x","y")]))

  stack.est = inla.stack(data=list(obs=df.st$obs),
                         A=list(A.est,
                                1),
                         effects=list(c(s.index, list(Intercept=1)),
                                      list(df.st[,-1])),
                         tag='stdata')
  
  output = inla(formula,
                data=inla.stack.data(stack.est),
                family="poisson",
                control.predictor=list(A=inla.stack.A(stack.est), compute = F),
                control.fixed=list(expand.factor.strategy='inla'),
                verbose=T)

  return(output)
  
}
