# 2021-07-07
# return the xy location of the specified corner and the x and y justification. Code taken from corner() and this function now replaces that code. Made a separate function to easily make use of it in functions like euclegend and gradleg.
# changed xyadj so its a fraction of axis rather than absolute units.
# NOTE: if the axes are logged, you still need to convert. I don't convert in this function because corner() makes an adjustment for multiline text before converting from log.
# 2021-12-07 I was looking at par(mfg) to decide how to scale. But it's easier to just check par(cex) and doing so allows me to change cex if needed (which prompted this all)

c.loc=function(xy,ofst=0.01,wrt=c(NULL,"x","y"),xyadj=c(0,0),outer=F) {
	#xy:	corner number (think quadrants) or description e.g. "topleft"
	#ofst:	offset from chosen corner, as fraction of short axis length
	#cex:	size of text  
	#wrt:	with respect to which axis (short axis by default)
	#xyadj:	xy adjustment
	#outer:	whether to go in the margins or not
	
	usr=par("usr") #(xlo,xhi,ylo,yhi)
	dimen=par("pin") #width, height
	shorty=which(dimen==min(dimen))[1] #which is the short axis?
	if(wrt[1]=="x" & shorty==2) shorty=1
	if(wrt[1]=="y" & shorty==1) shorty=2
	xconv=(usr[2]-usr[1])/dimen[1] #units per inch (conversion)
	yconv=(usr[4]-usr[3])/dimen[2] #units per inch (conversion)
	if(outer){
		din=par("din")
		mar=par("mar")
		oma=par("oma")
		mfg=par("mfg") #panels
		fac=0.2
		#if(mfg[3]==2 & mfg[4]==2) fac=fac*0.83
		#if(any(mfg[3:4]>2)) fac=fac*0.66	
		fac=fac*par("cex")
		#for each usr, add or substract inner and outer margin(s) and panels (if any) and marginS of those panels
		usr[1]=usr[1]-(mar[2]+oma[2])*fac*xconv-dimen[1]*xconv*(mfg[2]-1)-(mar[2]+mar[4])*fac*xconv*(mfg[2]-1) 
		usr[2]=usr[2]+(mar[4]+oma[4])*fac*xconv+dimen[1]*xconv*(mfg[4]-mfg[2])+(mar[2]+mar[4])*fac*xconv*(mfg[4]-mfg[2])
		usr[3]=usr[3]-(mar[1]+oma[1])*fac*yconv-dimen[2]*yconv*(mfg[3]-mfg[1])-(mar[1]+mar[3])*fac*yconv*(mfg[3]-mfg[1]) 
		usr[4]=usr[4]+(mar[3]+oma[3])*fac*yconv-dimen[2]*yconv*(mfg[1]-1)-(mar[1]+mar[3])*fac*yconv*(mfg[1]-1)
	}
	if(shorty==1){ #x is short
		xofst=ofst
		xinch=ofst*(usr[2]-usr[1])/xconv #how many inches from the axis?
		yunits=xinch*yconv #how many y units is xinch?
		yofst=yunits/(usr[4]-usr[3]) #what proportion of y-axis is xinch?	
	}else{	#y is short
		yofst=ofst
		yinch=ofst*(usr[4]-usr[3])/yconv #how many inches from the axis?
		xunits=yinch*xconv #how many x units is yinch?
		xofst=xunits/(usr[2]-usr[1]) #what proportion of x-axis is yinch?
	}
	
	if(xy=="topleft" | xy=="topright" | xy==1 | xy==2) {
		y=(1-yofst)*usr[4]+yofst*usr[3]	#top
		ud=1 #text below x,y (updown)
	}else{
		y=(yofst)*usr[4]+(1-yofst)*usr[3] #bottom
		ud=0 #text above x,y
	}
	if(xy=="topleft" | xy=="bottomleft" | xy==2 | xy==3) {
		x=(xofst)*usr[2]+(1-xofst)*usr[1]	#left
		lr=0 #text to the right of x,y (leftright)
	}else{
		x=(1-xofst)*usr[2]+xofst*usr[1]	#right
		lr=1 #text to the left of x,y
	}
	
	x=x+xyadj[1]*(usr[2]-usr[1])
	y=y+xyadj[2]*(usr[4]-usr[3])
	
	return(c(x,y,lr,ud))
}
