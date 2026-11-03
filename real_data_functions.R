# Data helpers for the two applications. Design recipes use covariates only.
# in validation, every recipe is learned on the training fold and applied to all sites.

# Ames: field checks and the basic main-effects design.
# Amenity areas a: H = I(a > 0) and L = log(a / g) if a > 0, else 0 (g: geometric mean of positive a).
check_ames_fields <- function(d, response=TRUE) {
  numeric_fields <- c("Lot_Area","Gr_Liv_Area","Total_Bsmt_SF","Garage_Cars",
    "Fireplaces","Year_Built","Year_Sold","Full_Bath","Half_Bath","Bedroom_AbvGr",
    "TotRms_AbvGrd","Wood_Deck_SF","Open_Porch_SF","Longitude","Latitude")
  if (response) numeric_fields <- c(numeric_fields,"Sale_Price")
  needed <- c(numeric_fields,"Central_Air","source_row")
  if (!all(needed %in% names(d))) stop("Required Ames fields are missing.")
  if (anyNA(d[,needed])) stop("Missing analysis fields: review them before fitting.")
  if (!all(vapply(d[,numeric_fields],is.numeric,logical(1))))
    stop("An analysis field has the wrong numeric type.")
  if (!all(is.finite(as.matrix(d[,numeric_fields])))) stop("Non-finite analysis data.")
  if (anyDuplicated(d$source_row)) stop("Duplicate source identifiers.")
  stopifnot(all(d$Lot_Area>0),all(d$Gr_Liv_Area>0),
    all(d$Year_Built<=d$Year_Sold),all(d$Year_Built>0),
    all(abs(d$Latitude)<=90),all(abs(d$Longitude)<=180))
  nonnegative <- c("Total_Bsmt_SF","Garage_Cars","Fireplaces","Full_Bath","Half_Bath",
    "Bedroom_AbvGr","TotRms_AbvGrd","Wood_Deck_SF","Open_Porch_SF")
  stopifnot(all(as.matrix(d[,nonnegative])>=0))
  counts <- c("Garage_Cars","Fireplaces","Full_Bath","Half_Bath","Bedroom_AbvGr",
    "TotRms_AbvGrd","Year_Built","Year_Sold")
  stopifnot(all(abs(as.matrix(d[,counts])-round(as.matrix(d[,counts])))<1e-8))
  if (response) stopifnot(all(d$Sale_Price>0))
  air <- as.character(d$Central_Air)
  if (!all(air %in% c("N","Y"))) stop("Unexpected Central_Air coding.")
  stopifnot(all(d$source_row>0),all(d$source_row==as.integer(d$source_row)))
  invisible(TRUE)
}

ames_features <- function(d, area_encoding, area_reference) {
  check_ames_fields(d,response=FALSE)
  Z <- cbind(log_Lot_Area=log(d$Lot_Area),log_Gr_Liv_Area=log(d$Gr_Liv_Area))
  if (area_encoding=="log1p")
    Z <- cbind(Z,log1p_Total_Bsmt_SF=log1p(d$Total_Bsmt_SF))
  Z <- cbind(Z,Garage_Cars=d$Garage_Cars,Fireplaces=d$Fireplaces,
    Age=d$Year_Sold-d$Year_Built,Full_Bath=d$Full_Bath,Half_Bath=d$Half_Bath,
    Bedroom_AbvGr=d$Bedroom_AbvGr,TotRms_AbvGrd=d$TotRms_AbvGrd,
    Central_Air=as.numeric(as.character(d$Central_Air)=="Y"))
  if (area_encoding=="log1p") {
    Z <- cbind(Z,log1p_Wood_Deck_SF=log1p(d$Wood_Deck_SF),
      log1p_Open_Porch_SF=log1p(d$Open_Porch_SF))
  } else {
    for (v in names(area_reference)) {
      positive <- d[[v]]>0
      L <- numeric(nrow(d))
      L[positive] <- log(d[[v]][positive])-area_reference[[v]]
      added <- cbind(as.numeric(positive),L)
      colnames(added) <- paste0(c("Has_","LogSize_"),v)
      Z <- cbind(Z,added)
    }
  }
  rownames(Z) <- as.character(d$source_row)
  Z
}

fit_ames_recipe <- function(d, area_encoding="presence_size") {
  stopifnot(area_encoding %in% c("presence_size","log1p"))
  check_ames_fields(d,response=FALSE)
  areas <- c("Total_Bsmt_SF","Wood_Deck_SF","Open_Porch_SF")
  ref <- vapply(areas,function(v) {
    positive <- d[[v]]>0
    if (!any(positive)) stop(paste("No positive training areas:",v))
    mean(log(d[[v]][positive]))
  },numeric(1))
  Z <- ames_features(d,area_encoding,ref)
  centers <- colMeans(Z)
  scales <- sqrt(colMeans(sweep(Z,2,centers,"-")^2))
  if (any(scales<=0)) stop("Constant training predictor: specification needs review.")
  list(area_encoding=area_encoding,area_reference=ref,centers=centers,
       scales=scales,columns=colnames(Z),training_source_rows=d$source_row)
}

# Ames weights: symmetrized k-nearest neighbours from great-circle distances, row-standardized.
make_observed_weights <- function(xy, k) {
  n <- nrow(xy)
  lon <- xy[, 1] * pi / 180
  lat <- xy[, 2] * pi / 180
  h <- sin(outer(lat, lat, "-") / 2)^2 +
       outer(cos(lat), cos(lat)) * sin(outer(lon, lon, "-") / 2)^2
  distance <- 2 * 6371008.8 * asin(sqrt(pmin(pmax(h, 0), 1)))
  diag(distance) <- Inf
  adjacency <- matrix(0, n, n)
  ties <- integer(n)
  for (i in seq_len(n)) {
    nearest <- order(distance[i, ], seq_len(n))
    adjacency[i, nearest[seq_len(k)]] <- 1
    ties[i] <- sum(distance[i, ] == distance[i, nearest[k]]) > 1L
  }
  adjacency <- 1 * ((adjacency + t(adjacency)) > 0)
  diag(adjacency) <- 0
  degree <- rowSums(adjacency)
  W <- adjacency / degree
  component <- integer(n)
  ncomp <- 0L
  while (any(component == 0L)) {
    ncomp <- ncomp + 1L
    todo <- which(component == 0L)[1]
    component[todo] <- ncomp
    while (length(todo)) {
      current <- todo[1]
      todo <- todo[-1]
      new <- which(adjacency[current, ] > 0 & component == 0L)
      component[new] <- ncomp
      todo <- c(todo, new)
    }
  }
  es <- eigen(adjacency / sqrt(outer(degree, degree)), symmetric = TRUE)
  vectors <- es$vectors
  left <- vectors / sqrt(degree)
  right <- t(vectors) * rep(sqrt(degree), each = n)
  K <- crossprod(vectors, degree * vectors) *
       crossprod(vectors, vectors / degree)
  stopifnot(max(abs(rowSums(W) - 1)) < 1e-12, all(diag(W) == 0),
            max(abs(W - left %*% (es$values * right))) < 1e-10)
  list(W = W, values = es$values, left = left, right = right, K = K,
       degree = degree, xy = xy, n = n, kind = paste0("observed_symmetric_", k, "NN"),
       components = ncomp, boundary_ties = sum(ties), adjacency = adjacency)
}


ames_broad_features <- function(d, recipe=NULL, min_count=10L) {
  fitting <- is.null(recipe)
  if(fitting) {
    base_recipe <- fit_ames_recipe(d)
    areas <- c("Mas_Vnr_Area","Garage_Area","Enclosed_Porch","Screen_Porch")
    area_reference <- vapply(areas,function(v) {
      stopifnot(is.numeric(d[[v]]),all(is.finite(d[[v]])),all(d[[v]]>=0))
      if(!any(d[[v]]>0)) return(0)
      mean(log(d[[v]][d[[v]]>0]))
    },numeric(1))
    recipe <- list(base=base_recipe,area_reference=area_reference,
      categorical=list(),min_count=min_count)
  }
  B <- ames_features(d,recipe$base$area_encoding,recipe$base$area_reference)
  additional <- c("Year_Remod_Add","Bsmt_Full_Bath","Bsmt_Half_Bath","Kitchen_AbvGr")
  stopifnot(all(vapply(d[,additional],is.numeric,logical(1))),
    all(is.finite(as.matrix(d[,additional]))),all(as.matrix(d[,additional])>=0),
    all(d$Year_Remod_Add<=d$Year_Sold))
  B <- cbind(B,Remodel_Age=d$Year_Sold-d$Year_Remod_Add,
    Bsmt_Full_Bath=d$Bsmt_Full_Bath,Bsmt_Half_Bath=d$Bsmt_Half_Bath,
    Kitchen_AbvGr=d$Kitchen_AbvGr)
  for(v in names(recipe$area_reference)) {
    x <- d[[v]]; stopifnot(all(is.finite(x)),all(x>=0))
    positive <- x>0; z <- numeric(length(x))
    z[positive] <- log(x[positive])-recipe$area_reference[[v]]
    extra <- cbind(as.numeric(positive),z)
    colnames(extra) <- paste0(c("Has_","LogSize_"),v)
    B <- cbind(B,extra)
  }
  factors <- c("MS_Zoning","Alley","Lot_Shape","Land_Contour","Lot_Config",
    "Land_Slope","Neighborhood","Bldg_Type","House_Style","Overall_Cond",
    "Roof_Style","Exterior_1st","Exter_Cond","Foundation","Bsmt_Exposure",
    "BsmtFin_Type_1","Heating_QC","Electrical","Functional","Garage_Type",
    "Garage_Finish","Paved_Drive","Fence")
  for(v in factors) {
    x <- as.character(d[[v]]); stopifnot(!anyNA(x))
    if(fitting) {
      counts <- table(x); rare <- names(counts)[counts<min_count]
      mapped <- x; mapped[x %in% rare] <- "Other"
      counts2 <- table(mapped)
      omitted <- length(counts2)<2L || any(counts2<min_count)
      lev <- names(counts2)[order(-as.integer(counts2),names(counts2))]
      spec <- list(rare=rare,observed=names(counts),levels=lev,reference=lev[1],
        omit=omitted,reason=if(omitted) "Insufficient category support after pooling" else "",counts=counts2)
      recipe$categorical[[v]] <- spec
    } else {
      spec <- recipe$categorical[[v]]
      mapped <- x; mapped[x %in% spec$rare] <- "Other"
      unseen <- !x %in% spec$observed
      if(!spec$omit && any(unseen) && !"Other" %in% spec$levels)
        stop(paste("Unseen category without learned Other group:",v))
      mapped[unseen] <- "Other"
    }
    if(spec$omit) next
    for(lev in spec$levels[-1]) {
      z <- matrix(as.numeric(mapped==lev),ncol=1)
      colnames(z) <- paste(v,make.names(lev),sep="__")
      B <- cbind(B,z)
    }
  }
  stopifnot(all(is.finite(B)),!anyDuplicated(colnames(B)))
  if(fitting) {
    spread <- apply(B,2,sd)
    constant <- names(spread)[spread<1e-12]
    candidate <- B[,spread>=1e-12,drop=FALSE]
    candidate <- scale(candidate)
    # Keep the first independent column; no y is used in this step.
    keep <- character(); Q <- matrix(1,nrow(B),1)
    for(v in colnames(candidate)) {
      trial <- cbind(Q,candidate[,v])
      if(qr(trial,tol=1e-9)$rank==ncol(trial)) {keep <- c(keep,v); Q <- trial}
    }
    recipe$columns <- keep
    recipe$dropped <- data.frame(term=c(constant,setdiff(colnames(candidate),keep)),
      reason=c(rep("Constant",length(constant)),rep("Exact linear dependency",ncol(candidate)-length(keep))))
    recipe$centers <- colMeans(B[,keep,drop=FALSE])
    recipe$scales <- sqrt(colMeans(sweep(B[,keep,drop=FALSE],2,recipe$centers,"-")^2))
    recipe$training_source_rows <- d$source_row
  }
  B <- B[,recipe$columns,drop=FALSE]
  X0 <- cbind(Intercept=1,B)
  X <- cbind(Intercept=1,sweep(sweep(B,2,recipe$centers,"-"),2,recipe$scales,"/"))
  rownames(X) <- rownames(X0) <- as.character(d$source_row)
  list(X=X,X0=X0,recipe=recipe,centers=recipe$centers,scales=recipe$scales)
}

# Apply a training . Unseen levels go to a learned Other group if available,
apply_ames_validation_recipe <- function(d,recipe) {
  modified <- d
  audit <- data.frame(source_row=integer(),variable=character(),observed=character(),mapped_to=character())
  for(v in names(recipe$categorical)) {
    spec <- recipe$categorical[[v]]
    if(spec$omit) next
    x <- as.character(d[[v]]); unknown <- !x %in% spec$observed
    if(!any(unknown)) next
    target <- if("Other" %in% spec$levels) "Other" else spec$reference
    audit <- rbind(audit,data.frame(source_row=d$source_row[unknown],variable=v,
      observed=x[unknown],mapped_to=target))
    if(target!="Other") {x[unknown] <- target; modified[[v]] <- x}
  }
  design <- ames_broad_features(modified,recipe)
  list(design=design,unseen=audit)
}

# Chicago: 15 predictors, centred and scaled with training rows only.
chicago_design <- function(d,expanded=FALSE,train=seq_len(nrow(d)),recipe=NULL) {
  base <- c("EP_NOHSDP","Unemp14","EP_AGE17","Ovr6514P","EP_DISABL",
    "EP_SNGPNT","EP_CROWD","EP_LIMENG","EP_NOVEH","EP_MUNIT","EP_GROUPQ",
    "CTA14P","CmTm14","log_density","MEANMI_14")
  d$log_density <- log(d$PDENS14)
  v <- as.matrix(d[,base])
  stopifnot(all(is.finite(v)),length(train)>0,!anyDuplicated(train),
    all(train %in% seq_len(nrow(d))))
  if(is.null(recipe)) {
    parent_centers <- colMeans(v[train,,drop=FALSE])
  } else {
    expanded <- recipe$expanded
    stopifnot(identical(base,recipe$base))
    parent_centers <- recipe$parent_centers
  }
  vc <- sweep(v,2,parent_centers,"-")
  parents <- list(education_sq="EP_NOHSDP",unemployment_sq="Unemp14",density_sq="log_density",
    education_unemployment=c("EP_NOHSDP","Unemp14"),
    crowding_singleparent=c("EP_CROWD","EP_SNGPNT"),
    noveh_transit=c("EP_NOVEH","CTA14P"),
    multiunit_density=c("EP_MUNIT","log_density"))
  if(expanded) {
    extra <- vapply(parents,function(pp) {
      if(length(pp)==1) vc[,pp]^2 else vc[,pp[1]]*vc[,pp[2]]
    },numeric(nrow(v)))
    v <- cbind(v,extra)
  }
  if(is.null(recipe)) {
    centers <- colMeans(v[train,,drop=FALSE])
    scales <- sqrt(colMeans(sweep(v[train,,drop=FALSE],2,centers,"-")^2))
    recipe <- list(expanded=expanded,base=base,parent_centers=parent_centers,
      centers=centers,scales=scales,columns=colnames(v),parents=if(expanded) parents else list())
  }
  stopifnot(identical(colnames(v),recipe$columns),all(recipe$scales>0))
  X0 <- cbind(Intercept=1,v)
  X <- cbind(Intercept=1,sweep(sweep(v,2,recipe$centers,"-"),2,recipe$scales,"/"))
  rownames(X0) <- rownames(X) <- sprintf("%.0f",d$geoid10)
  list(X=X,X0=X0,recipe=recipe,centers=recipe$centers,scales=recipe$scales)
}

# Chicago weights: queen contiguity, row-standardized.
queen_weights <- function(polygons,ids) {
  # Needs the sf package.
  # A_ij=1 when tract boundaries touch (edge or vertex); W=D^(-1)A.
  nb <- st_touches(polygons)
  n <- length(nb)
  A <- matrix(0,n,n)
  for(i in seq_len(n)) A[i,nb[[i]]] <- 1
  stopifnot(all(diag(A)==0),max(abs(A-t(A)))==0)
  degree <- rowSums(A)
  if(any(degree==0)) stop(paste("Isolated tracts require review:",paste(ids[degree==0],collapse=", ")))
  component <- integer(n); count <- 0L
  for(i in seq_len(n)) if(component[i]==0L) {
    count <- count+1L; queue <- i; component[i] <- count; head <- 1L
    while(head<=length(queue)) {
      j <- queue[head]; head <- head+1L
      new <- nb[[j]][component[nb[[j]]]==0L]
      component[new] <- count; queue <- c(queue,new)
    }
  }
  W <- A/degree
  dimnames(W) <- list(ids,ids)
  B <- A/sqrt(outer(degree,degree))
  es <- eigen(B,symmetric=TRUE)
  V <- es$vectors

  left <- V/sqrt(degree)
  right <- t(V)*rep(sqrt(degree),each=n)
  K <- crossprod(V,degree*V)*crossprod(V,V/degree)
  stopifnot(max(abs(rowSums(W)-1))<1e-12,max(abs(es$values))<1+1e-10,
    max(abs(W-left%*%(es$values*right)))<1e-10)
  list(W=W,B=B,values=es$values,left=left,right=right,K=K,
    degree=degree,components=count,component=component,adjacency=A,n=n,kind="queen")
}
