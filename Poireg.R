Poireg <- function (formula,data,estimate_se=T,init=NULL) {
vars <- all.vars(formula)
S1 <- data[,paste0(vars[1],1)]
S2 <- data[,paste0(vars[1],2)]
X1 <- data.frame(data[,paste0(vars[-1],1)])
X2 <- data.frame(data[,paste0(vars[-1],2)])
cat <- sapply(X1,is.factor) #which variable is a categorical variable
cont <- 1-cat #which variable is a continuous variable
cat <- which(cat)
cont <- which(cont==1)
nd <- length(cat) #number of categorical variable
nc <- length(cont) #number of continuous variable
n <- which(S1*S2!=0)
S1 <- S1[n]
S2 <- S2[n]
X1 <- data.frame(X1[n,])
X2 <- data.frame(X2[n,])
n <- length(S1)
if (nd>0) {
	c <- sapply(cat,function (i) length(levels(X1[,i]))) #number of categories in each categorical variable, the first category of each variable is the reference category
}
X1 <- data.matrix(X1)
X2 <- data.matrix(X2)
if (nd>0) {
	X1[,cat] <- X1[,cat]+t(matrix(rep(c(0,c[-nd])),nd,n))
	X2[,cat] <- X2[,cat]+t(matrix(rep(c(0,c[-nd])),nd,n))
}
nll <- function (p) {
	if (nd>0) {
		b <- numeric(sum(c))
		ref <- cumsum(c(1,c))
		ref <- ref[-length(ref)]
		b[ref] <- 1
		b[-ref] <- exp(p[nc+(1:sum(c-1))])
	}
	r1 <- exp(p[length(p)+c(-(n-1):0)])
	r2 <- r1
	if (nc>0) {
		X <- X2[,cont]-X1[,cont]
		if (nc==1) {
			tmp <- X*p[1:nc]
		} else {
			tmp <- X%*%p[1:nc]
		}
		r1[tmp!=0] <- r1[tmp!=0]*(1-exp(-tmp[tmp!=0]/2))/tmp[tmp!=0]
		r2[tmp!=0] <- -r2[tmp!=0]*(1-exp(tmp[tmp!=0]/2))/tmp[tmp!=0]
	}
	if (nd>0) {
		r1 <- r1*sapply(1:n,function (j) prod(sapply(cat,function (i) b[X1[j,i]])))
		r2 <- r2*sapply(1:n,function (j) prod(sapply(cat,function (i) b[X2[j,i]])))
	}
out <- numeric(n)
for (i in which(r1!=0)) {
out[i] <- S1[i]*log(r1[i])-lfactorial(S1[i])
}
for (i in which(r2!=0)) {
out[i] <- out[i]+S2[i]*log(r2[i])-lfactorial(S2[i])
}
sum(r1+r2)-sum(out)
}
if (is.null(init)) {
	out <- nlm(nll,p=c(numeric(nc)+0.05,numeric(ifelse(nd>0,sum(c-1),0)+n)),hessian=T)
} else {
	out <- nlm(nll,p=c(init,numeric(ifelse(nd>0,sum(c-1),0)+n)),hessian=T)
}
nll <- function (p) {
n <- length(S1)
r <- exp(p)
r1 <- r
r2 <- r
out <- numeric(n)
for (i in 1:n) {
out[i] <- S1[i]*log(r1[i])+S2[i]*log(r2[i])-lfactorial(S1[i])-lfactorial(S2[i])
}
sum(r1+r2)-sum(out)
}
null <- nlm(nll,p=numeric(n),hessian=T)
mean<-out$estimate[1:(nc+ifelse(nd>0,sum(c-1),0))]
if (estimate_se) {
	se<-sqrt(diag(solve(out$hessian)))[1:(nc+ifelse(nd>0,sum(c-1),0))]
	t.se<-se*qt(0.025,2*length(S1)-length(out$estimate))
	CI<-rbind(mean-t.se,mean+t.se)
} else {
	se <- NA
	t.se <- NA
	CI <- NA
}
likratio <- 2*(null$minimum-out$minimum)
pval <- 1-pchisq(likratio,df=nc+ifelse(nd>0,sum(c-1),0))
Rs.full<--sum(S1[S1!=0])-
sum(S2[S2!=0])+sum(S1[S1!=0]*log(S1[S1!=0]))+sum(S2[S2!=0]*log(S2[S2!=0]))-
sum(lfactorial(S1[S1!=0]))-sum(lfactorial(S2[S2!=0]))
Rs.hat<--out$minimum
sbar<-exp(null$estimate)
Rs.bar<--sum(sbar[sbar!=0])-
sum(sbar[sbar!=0])+sum(S1[sbar!=0]*log(sbar[sbar!=0]))+sum(S2[sbar!=0]*log(sbar[sbar!=0
]))-sum(lfactorial(S1[sbar!=0]))-sum(lfactorial(S2[sbar!=0]))
Rsquare<-1-(Rs.full-Rs.hat)/(Rs.full-Rs.bar)
list(mean=mean,se=se,CI=CI,likratio=likratio,pval=pval,Rsquare=Rsquare)
}
