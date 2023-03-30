INLA_run = function(data=df.st ,virus, wmel_proportion, in.project, prediction_type){
  
  grid_res = 0.5 #km
  rand_model = "rw1"
  
  #Create the input data frame-------------------------
  if(virus == "dengue"){
    data = data |> 
      select(dengue_cases,x,y,t,month,In.project,wmel_pctg,Population,cell.id) |>
      rename(obs = dengue_cases, Wolbachia = wmel_pctg) |>
      mutate(Wolbachia.bin = cut(Wolbachia, right=T, include.lowest=T, breaks=c(0,seq(1e-5,0.7,by=0.1),1)))
  }
  
  if(virus == "chik"){
    data = data |> mutate(month=as.Date(month)) |>
      mutate(t=1+as.integer(month)-min(as.integer(month))) |>
      select(chik_cases,x,y,t,month,In.project,wmel_pctg,Population,cell.id) |>
      rename(obs = chik_cases, Wolbachia = wmel_pctg) |>
      mutate(Wolbachia.bin = cut(Wolbachia, right=T, include.lowest=T, breaks=c(0,seq(1e-5,0.7,by=0.1),1)))
  }
  
  
  
  #Held_out data----------------
  #Full set - All cases areused for data prediciton
  if(prediction_type==1){ 
  } 
  #Half set - 50% of cases are removed at random
  if(prediction_type==2){ 
    out = sample(nrow(data), size=round(0.5*nrow(data)), replace=F)
    data$obs[out] = NA
  }
  
  #Spatial set - whole years removed(20%)
  if(prediction_type==3){ 
    radius_chunk=3/grid_res #Approx.
    out = c()
    nstep.in.year  = 12
    for(i in 1:(max(data$t)/nstep.in.year)){
      total_chunks = round(0.2*length(unique(data$cell.id))/(pi*radius_chunk^2))
      ratio.check = 0
      while(ratio.check<0.195){ #We want to approx. hold out 20% of the data per time step
        kern = sample(length(unique(data$cell.id)), size=total_chunks)
        kern.coords = data[!duplicated(data$cell.id),c("x","y")][kern,]
        kern.neighbours = get.knnx(data[!duplicated(data$cell.id),c("x","y")], kern.coords, (pi*radius_chunk^2)+1)
        ratio.check = sum((1:length(unique(data$cell.id)) %in% kern.neighbours$nn.index))/length(unique(data$cell.id))
        total_chunks = total_chunks+1
      }
      out = c(out, rep((1:length(unique(data$cell.id))  %in% kern.neighbours$nn.index), nstep.in.year))
    }
    
    data$obs[out] = NA
  } 
  
  
  #Creation of the mesh and the SPDE----------
  spatial.range = 22000
  
  mesh = inla.mesh.2d(loc=data[!duplicated(data$cell.id),c("x","y")], max.edge=c(spatial.range/5, spatial.range),
                      cutoff = 1500,
                      offset=c(spatial.range/5, diff(range(data$x))/3))
  
  spde = inla.spde2.matern(mesh=mesh,alpha=2)
  
  #Formula-------------
  if(wmel_proportion){
    diag.param = 1e-1 #increasing "diagonal" parameter to avoid numerical singularity in wMel.grid - check https://groups.google.com/forum/#!topic/r-inla-discussion-group/8zFb6oi9O-U for more detail
    length.wb = length(levels(data$Wolbachia.bin)) - 1
    
    formula = obs ~ -1 + Intercept + offset(log(Population)) +  f(t, model = 'ar1') +   f(spatial.field, model = spde) +
      f(Wolbachia.bin, model=rand_model, diagonal=diag.param, constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length.wb)), nrow=1), e=rep(0,1)))
  }
  
  if(in.project){
    
    formula = obs ~ -1 + Intercept + offset(log(Population)) +  f(t, model = 'ar1') + f(spatial.field, model = spde) + In.project
    
  }
  
  
  #Stack creation and inla analysis-------------------------------
  s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde) 
  
  A.est = inla.spde.make.A(mesh=mesh,
                           loc=as.matrix(data[,c("x","y")]))
  
  stack.est = inla.stack(data=list(obs=data$obs),
                         A=list(A.est,
                                1),
                         effects=list(c(s.index, list(Intercept=1)),
                                      list(data[,-1])),
                         tag='stdata')
  
  output = inla(formula,
                data=inla.stack.data(stack.est),
                family="poisson",
                control.predictor=list(A=inla.stack.A(stack.est), compute = F),
                control.fixed=list(expand.factor.strategy='inla'),
                verbose=T)
  
  return(output)
  
}
