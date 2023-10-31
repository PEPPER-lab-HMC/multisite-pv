#this function adds text to the corner of a plot by specifying which corner and an offset
#this is much easier and more precise than manually choosing x and y values
#the offset is a fraction of the total spans in the x and y directions.
#only one offset is specified, so the text must be in the corner. 

#2012-10-01 update: having ofst as a proportion of each axis leads to funny results in
#plots that are far from square. Make it so ofst is only based on short axis
#now resizing quartz requires that plot be remade for correct placement

#2012-11-20 update: the offset from the corner wasn't right when text was
#multiple characters long (text could run off of plot). Within text() 'pos' was used at first.
#This specifies the text as to the left, right, top, or bottom side of the x,y point
#Using 'adj' is better because it can specify corners instead of sides.
#Also, compatability with log scales was left out initially. Now works fine.
#Finally, multiline capability was added such that "text" can be of length > 1

#2013-05-05 update: added ability to specify quadrant number for which corner  

#2013-12-06 update: added 'wrt' which allows you to specify which axis to offset from. 
														#Shorty is default though 
#2018-10-22 update: added 'xyadj' so you can shift where the reference corner is. In plot units. An alternate to this would be adding option that ofst can be length 2

#2021-07-06 update: added 'outer' option by redefining usr as the limits of the plot. This should only work for single panel figures but shouldn't be difficult to add that option. Would just need to know which panel you're in and the size of each panel.

#2021-07-06 an hour later: Now it works for multi panel figures. Wasn't as 'easy' as expected. Unexpected need to account for different scale in mult panel figures (i.e. not always 0.2 inches per margin unit) and account for oma and mar of other panels. Don't know about log scale

#2021-07-07: move the guts to the c.loc function. Changed xyadj so that it's a fraction rather than absolute units

#2022-03-23: at the end, set xpd to the original value - rather than F by default

corner=function(xy,text,ofst=0.01,cex=1,test=F,wrt=c(NULL,"x","y"),spcadj=1,xyadj=c(0,0),outer=F,...) {
	#xy:	corner number (think quadrants) or description e.g. "topleft"
	#text:	character string(s)
	#ofst:	offset from chosen corner, as fraction of short axis length
	#cex:	size of text  
	#test:	should ablines be drawn to show the selected point relative to the text?
	
	loc=c.loc(xy,ofst,wrt,xyadj,outer)
	x=loc[1]
	y=loc[2]
	lr=loc[3]
	ud=loc[4]
	oxpd=par("xpd")
	
	if(length(text)>1){
		usr=par("usr") #(xlo,xhi,ylo,yhi)
		dimen=par("pin") #width, height
		xconv=(usr[2]-usr[1])/dimen[1] #units per inch (conversion)
		yconv=(usr[4]-usr[3])/dimen[2] #units per inch (conversion)
		# ^^ c.loc calculated xconv too. But need it again here if mult lines
		cin=par("cin")[2]*spcadj #character height
		if(ud==1){ #top, values must go down FROM y
			y=y-cex*cin*yconv*1*(1:length(text)-1)
		}else{	#bottom, values go down TO y
			y=sort(y+cex*cin*yconv*1*(1:length(text)-1),T)
		}
	}
	
	if(par("xlog")) x=10^x
	if(par("ylog")) y=10^y
	#x=ifelse(par()$xlog==T,10^x,x) #adjust to log scale if needed
	#y=ifelse(par()$ylog==T,10^y,y) #ditto
	par(xpd=NA)
	text(x,y,text,cex=cex,adj=c(lr,ud),...)
	if(test==T) {abline(v=x,h=y,lwd=0.5)} #add lines to check placement
	par(xpd=oxpd)
	###test using multiple lines
	#cin=par("cin")[2]
	#text(x,y-cex*cin*yconv*0.9,"test",cex=cex,adj=c(lr,ud),col=col)
	#text(x,y-2*cex*cin*yconv*0.9,"grrrr!",cex=cex,adj=c(lr,ud),col=col)
	########
	
	#points(x,y,col="red",pch=19,cex=0.2)
	
}
