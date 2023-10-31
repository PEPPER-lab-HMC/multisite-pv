#How to use:
	#1) run/source this code
	#2) set your working directory
	#3) call the pv() function. Minimally just supply a sample name e.g. pv("sample 1") 

#PV curve function - this function asks for pressure and mass and stores them in a new .csv file or adds them if the file already exists in the working directory. Data is then plotted.
#A time stamp was initially added to each entry, but I found that data aren't always entered as soon as they're created, so the timestamp becomes meaningless. 
#2016-06-07: Addition: getting regressions from last point through every other point and showing the r2 next to each point. 
#2016-06-20: Make and use a column in saved data for whether to use a point when calculating r2. The point will still be shown but in grey. 
#2023-05-16: added time back in and made it flexible for units.

## VERSION "2"
#2023-05-24: There are two goals in this reversion. 
	#1) is to ask for offset mass where you'd add the mass of the gasket and make any other adjustements as you go. Then there would be columns for observed mass, offset mass and leaf mass. Then leaf mass would be plotted. 
	#2) is to (try to) fit a model to the data and estimate the location of the TLP. This would be done once there are at least three points (if not more). It would also be good to record estimated TLP as you go and then show how it is changing with each measurement.
#2023-10-18: Forced pressure to be positive. Changed plots to read -Psi. Removed "showr2" argument

pv=function(nom,plotonly=F,unit="MPa",logfit=T){
	tm=Sys.time()
	#setwd("~/desktop")
	par(mgp=c(2,0.5,0),mar=c(4,4,1,1))
	fls=list.files()
	fnom=paste(nom,"pvdata.csv")
	fls=fls[fls==fnom]
	#fls=fls[grep(nom,fls)]
	pnom=paste0("P.",unit)
	if(length(fls)>0){
		dat=read.csv(paste(nom,"pvdata.csv"))
		dat$time=as.character(dat$time)
		names(dat)[names(dat)==pnom]="p.bar"
		nob=dim(dat)[1]+1
		lst=dat[nob-1,"p.bar"]
		lst.m=dat[nob-1,"mass.g"]
		lst.om=dat[nob-1,"offset.mass.g"]
	}else{
		cat("creating NEW file:\n\t",fnom,"\n\tin working directory\n",sep="")
		lst="NA"
		lst.m="NA"
		lst.om=0
		nob=1
	}
	if(!plotonly){
		p.bar=as.numeric(readline(paste0("pressure (",unit,")? - last = ",lst,": ")))
		while(is.na(p.bar)){
			p.bar=as.numeric(readline(paste0("invalid answer: pressure (",unit,")? - last = ",lst,": ")))
		}
		p.bar=abs(p.bar)
		m.g=as.numeric(readline(paste0("total mass (g)? - last = ",lst.m,": ")))
		while(is.na(m.g)){
			m.g=as.numeric(readline("invalid answer: total mass (g)? "))
		}
		om.g=as.numeric(readline(paste0("new offset mass (g)? ")))
		if(is.na(om.g)) om.g=0
		om.g=om.g+ lst.om
		lm.g=m.g-om.g #leaf mass
		
		if(length(fls)>0){	
			dat[nob,"p.bar"]=p.bar
			dat[nob,"mass.g"]=m.g
			dat[nob,"leaf.mass.g"]=lm.g
			dat[nob,"offset.mass.g"]=om.g
			dat[nob,"keep"]=T
			dat[nob,"time"]=as.character(tm)
		}else{
			dat=data.frame(p.bar=p.bar,mass.g=m.g,leaf.mass.g= lm.g,offset.mass.g=om.g,keep=T,time=tm,note=" ",eTLP=NA,xTLP=NA,a=NA,b=NA,c=NA, d=NA)
		}
	}
	if(length(fls)==0 & plotonly){
		print(paste("no data to plot in current directory:",getwd()))
	}else{
		p=dat$p.bar
		ip=1/dat$p.bar
		ml=dat$leaf.mass.g[1]-dat$leaf.mass.g	
		#Fit the model if possible
		if(sum(dat$keep)>2){
			if(logfit){
				fit=try(nls(log(ip)~a*exp(-ml*b)+cc+d*ml,start=list(a=5,b=150,cc=0.5,d=-1.5),subset=dat$keep),silent=T)
			}else{
				fit=try(nls(ip~a*exp(-ml*b)+cc+d*ml,start=list(a=5,b=150,cc=0.5,d=-1.5),subset=dat$keep),silent=T)
			}
			if(!inherits(fit,"try-error")){
				fit=summary(fit)$coef
				if(logfit){
					tpf=function(x,a,b,cc,d) exp(a*exp(-x*b)+cc+d*x) #total component 
					otpf=function(x,cc,d) exp(cc+d*x) #osmotic component
					ptpf=function(x,a,b,cc) exp(a*exp(-b*x)+cc) #pressure component
					xTLP=-log(-fit["d",1]/(fit["a",1]*fit["b",1]))/fit["b",1]
				
				}else{
					tpf=function(x,a,b,cc,d) (a*exp(-x*b)+cc+d*x) #total component 
					otpf=function(x,cc,d) (cc+d*x) #osmotic component
					ptpf=function(x,a,b,cc) (a*exp(-b*x)+cc) #pressure component
					xTLP=-log(-fit["d",1]/(fit["a",1]*fit["b",1]))/fit["b",1]
				}
				eTLP=1/tpf(xTLP,fit["a",1],fit["b",1],fit["cc",1], fit["d",1])
				dat[nob-plotonly,"xTLP"]=xTLP
				dat[nob-plotonly,"eTLP"]=eTLP
				dat[nob-plotonly,c("a","b","c","d")]=fit[c("a","b","cc", "d"),1]
				is.fit=T
			}else {is.fit=F;cat("Unable to fit data\n")}
		}else is.fit=F
		if(!is.fit) {xTLP=NA;dat[nob,c("a","b","c","d","eTLP","xTLP")]=NA}
		## Inverse curve
		par(mfrow=c(1,2))
		plot(ip~ml,pch=19,ylab=as.expression(substitute(1/-psi~(u^-1),list(u=unit))),xlab="mass lost (g)",type="n")
		if(is.fit){ 
			abline(v=dat$xTLP,lwd=2,col="grey")
			abline(v=xTLP,lwd=2)
			corner(1,c("est TLP:",round(eTLP,3)),cex=0.6)
		}
		points(ip[dat$keep]~ml[dat$keep],pch=19)
		points(ip[!dat$keep]~ml[!dat$keep],pch=19,col="grey")
		points(ip[nob]~ml[nob],pch=19,cex=0.8,col="cyan")
		if(is.fit){ 
			curve(tpf(x,fit["a",1],fit["b",1],fit["cc",1],fit["d",1]),add=T, col="red")
			curve(otpf(x,fit["cc",1],fit["d",1]),add=T,lty=2)
			curve(ptpf(x,fit["a",1],fit["b",1],fit["cc",1]),add=T,lty=2)
		}
		## Non inverse curve
		plot(p~ml,pch=19,ylab=as.expression(substitute(-psi~(u^-1), list(u=unit))),xlab="mass lost (g)",type="n")
		if(is.fit){ 
			abline(v=dat$xTLP,lwd=2,col="grey")
			abline(v=xTLP,lwd=2)
		}
		points(p[dat$keep]~ml[dat$keep],pch=19)
		points(p[!dat$keep]~ml[!dat$keep],pch=19,col="grey")
		points(p[nob]~ml[nob],pch=19,cex=0.8,col="cyan")
		if(is.fit){ 
			curve(1/tpf(x,fit["a",1],fit["b",1],fit["cc",1],fit["d",1]),add=T, col="red")
			curve(1/otpf(x,fit["cc",1],fit["d",1]),add=T,lty=2)
			curve(1/ptpf(x,fit["a",1],fit["b",1],fit["cc",1]),add=T,lty=2)	
		}
		mtext(nom,3,-1,outer=T,adj=0.01)
	}	
	if(plotonly) tst=readline("save data? (y) ") else tst=""
	if(!plotonly | tst=="y"){
		#points(ip[nob]~ml[nob],pch=19,cex=0.8,col="cyan")
		names(dat)[names(dat)=="p.bar"]=pnom
		note=readline("note? ")
		if(note=="") note=" "
		dat$note[nob-plotonly]=note
		write.csv(dat,paste(nom,"pvdata.csv"),row.names=F)
		print(paste("data saved to",getwd()))
		print("open .csv file and change keep=TRUE to keep=FALSE to exclude points")
	}
}