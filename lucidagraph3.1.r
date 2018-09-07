
lucidaxy<-choose.files(caption="Seleziona il file LucidaXY",multi=FALSE)
setwd(dirname(lucidaxy))
dati<-read.csv(lucidaxy)

borderWidth=150


# hist(dati$PositionX)
library(polyclip)
library (sp)
library(rgeos)
#Il sistema prevede un file CSV da Mark&Find dove siano state rimosse le prime righe
#il file deve iniziare con le intestazioni delle colonne (Comments	PositionX	PositionY...)
#il contorno del cervello deve essere marcato "contorno", i foci tumorali con stringhe contententi "tumor"
#eventuale canale di iniezione deve essere etichettato con una stringa contenente "iniez"
#le strutture vuote (strappi, ventricoli, cavità varie possono essere nominate a piacere,
#ognuna con un nome differente
#le cellule devono essere marcate con una o più stringhe da indicare nella variabile labCell 
#qui sotto
#########################################################
#inserire qui le etichette (classification) delle cellule###
labCell<-c("CD4")
########################################################


area<-function(pol){ ####pol deve essere di tipo polygons 
	np<-length(pol)
	aparz<-0
	for(i in 1:np){
		tx<-pol[[i]]$x
		ty<-pol[[i]]$y
		nlt<-length(tx)
		tx[nlt+1]<-tx[1]
		ty[nlt+1]<-ty[1]
		s1<-sum(tx[1:nlt]*ty[2:(nlt+1)])
		s2<-sum(ty[1:nlt]*tx[2:(nlt+1)])
		#aparz<-aparz+(abs(s1-s2)/2)
		aparz<-aparz+(s1-s2)/2
	}
return(aparz)}
#######################################
#Utilità:#
#i<-i+1;points(da$PositionX[i],da$PositionY[i],pch=19,col=2,cex=.8)
##########################################


if(winDialog(type = c("yesno"), "Ci sono cellule localizzate su foto?")=="YES"){
##############################################################
#################Sitema di lettura cellule da file############
##############################################################

foto<-read.csv(choose.files(default=paste(dirname(lucidaxy),"*.*",sep="/"),caption="Seleziona il file CellXY",multi=FALSE))
centrifoto<-read.csv(choose.files(default=paste(dirname(lucidaxy),"*.*",sep="/"),caption="Seleziona il file MetadatiFotoXY",multi=FALSE))

PositionX<-(foto$Index)*0
PositionY<-(foto$Index)*0
PositionZ<-(foto$Index)*0
for(ind in foto$Index){
	
	Z<-foto$Z[ind]
	CentroFotoX<-centrifoto$centrox[Z]
	CentroFotoY<-centrifoto$centroy[Z]
	RicentrCellX<-(foto$X[ind]-694)*centrifoto$Scale[Z]
	RicentrCellY<-(foto$Y[ind]-520)*centrifoto$Scale[Z]
	PositionX[ind]<-CentroFotoX+RicentrCellX
	PositionY[ind]<-CentroFotoY-RicentrCellY
	PositionZ[ind]<-Z
	}
	
Comments<-(rep("Cells da foto",length(PositionX)))
Color<-(rep("Blue",length(PositionX)))
Classification<-(rep(labCell,length(PositionX)))
CellDaFoto<-data.frame(Comments,PositionX,PositionY,PositionZ,Color,Classification)
##############################################################
##############################################################


da<-rbind(dati,CellDaFoto)
}else
{da<-dati}

#centratura mediante traslazione del contorno del cervello
w<-max(da$PositionX)-min(da$PositionX)
h<-max(da$PositionY)-min(da$PositionY)
cx<-min(da$PositionX)+w/2
cy<-min(da$PositionY)+h/2
da$PositionX<-da$PositionX-cx
da$PositionY<-da$PositionY-cy

# cont<-da$Classification=="contorno"
cont<-unique((da$Classification[grep("contorn",da$Classification,ignore.case = T)]))
tumori<-unique((da$Classification[grep("tum|iniez",da$Classification,ignore.case = T)]))
buchi<-unique((da$Classification[grep(paste("tum|iniez|contorn",paste(labCell,collapse="|"),sep="|"),da$Classification,invert=TRUE,ignore.case = T)]))


tipiOggetti<-unique (da$Classification)
contornoGlob<-da$Classification%in%cont

colori<-9+(as.integer(tipiOggetti)*0)
colori[tipiOggetti%in%cont]<-1
colori[tipiOggetti%in%tumori]<-3
colori[tipiOggetti%in%buchi]<-2
colori[tipiOggetti%in%labCell]<-match(tipiOggetti[tipiOggetti%in%labCell],labCell)+3


plot(da$PositionX[contornoGlob],da$PositionY[contornoGlob],asp=1,type="n",xlab="micron",ylab="")


rims<-unique(da$Classification[grep(paste(labCell,collapse="|"),da$Classification,invert=TRUE,ignore.case = T)])
coloRims<-colori[tipiOggetti%in%rims]

for(i in (1:length(rims))){
	
	lin<-da$Classification==rims[i]
	polygon(da$PositionX[lin],da$PositionY[lin],border= coloRims[i] )
}

cell<- !(da$Classification %in% rims)
cellpiede<-cell&!(da$Comments=="Cells da foto")
cellfoto<-da$Comments=="Cells da foto"

points(da$PositionX[cellpiede],da$PositionY[cellpiede],col=colori[match(da$Classification[cellpiede],tipiOggetti)],cex=.5,pch=19)
points(da$PositionX[cellfoto],da$PositionY[cellfoto],col=colori[match(da$Classification[cellfoto],tipiOggetti)],pch=".")



tumlist<-list()
if(length(tumori)>0)
for(i in 1:length(tumori)){
	tumlist[[i]]<-list(x=da$PositionX[da$Classification==tumori[i]],y=da$PositionY[da$Classification==tumori[i]])
	
	}else{
	tumlist[[1]]<-list(x=c(0,0),y=c(0,0))
}

buchilist<-list()
if(length(buchi)>0)
for(i in 1:length(buchi)){
	buchilist[[i]]<-list(x=da$PositionX[da$Classification==buchi[i]],y=da$PositionY[da$Classification==buchi[i]])
	}else{
	buchilist[[1]]<-list(x=c(0,0),y=c(0,0))
}

brain<-list()
for(i in 1:length(cont)){
	brain[[i]]<-list(x=da$PositionX[da$Classification==cont[i]],y=da$PositionY[da$Classification==cont[i]])
}

# brain<-list(list(x=da$PositionX[cont],y=da$PositionY[cont]))
brainmass<-polyclip(brain,buchilist,op="minus")
tumormass<-polyclip(tumlist,buchilist,op="minus")
#area(polyclip(list(list(x=da$PositionX[da$Classification=="contorno1"],y=da$PositionY[da$Classification=="contorno1"])),so1,op="union"))
#tumormass2<-polyclip(polyclip(brain,tumlist,op="intersection"),buchilist,op="minus")
##########
tum<-da$Classification%in%tumori

##Calcolo cellule nel tumore
inTum<- rep(F,sum(cell))
for(i in (1:length(tumlist))){
	inTum<- inTum|(point.in.polygon(da$PositionX[cell],da$PositionY[cell],tumlist[[i]]$x,tumlist[[i]]$y)>0)
}
##Calcolo cellule nei buchi
inBuchi<-rep(F,sum(cell))
for(i in (1:length(buchilist))){
	inBuchi<- inBuchi|(point.in.polygon(da$PositionX[cell],da$PositionY[cell],buchilist[[i]]$x,buchilist[[i]]$y)>0)
}


##Calcolo cellule nei bordi del tumore (con l'operazione di intersezione viene eliminato dal
##calcolo della superficie del bordo, tutto quello che sta all'esterno del cervello. Questo però
##non elimina dal bordo la porzione di area della fascia interna di un tumore il quale termini su un confine del cervello)
onRim<-rep(F,sum(cell))
areaRim<-0
# borderWidth=150

for(i in (1:length(tumlist))){
	
	Otum<-SpatialPolygons(list(Polygons(list(Polygon(cbind(tumlist[[i]]$x,tumlist[[i]]$y))),i) ))
	
	Gtum<-gBuffer(Otum,width=borderWidth)
	
	Gx<-Gtum@polygons[[1]]@Polygons[[1]]@coords[,1]
	Gy<-Gtum@polygons[[1]]@Polygons[[1]]@coords[,2]
	inG<-point.in.polygon(da$PositionX[cell],da$PositionY[cell],Gx,Gy)>0
	NetinG<-polyclip(list(list(x=Gx,y=Gy)),brainmass,op="intersect") #questo elimina dal calcolo del bordo ciò che è fuori dal cervello
    areaEXT<-area(NetinG)
	# for (k in 1:length(NetinG)){
		# polygon(NetinG[[k]],border=7)
	# }
	
	Stum<-gBuffer(Otum,width=-borderWidth)
	if(is.null(Stum)){
		inS<-rep(FALSE,length(da$PositionX[cell]))
		areaINT<-0
		}else{
		
		Sx<-Stum@polygons[[1]]@Polygons[[1]]@coords[,1]
		Sy<-Stum@polygons[[1]]@Polygons[[1]]@coords[,2]
		inS<-point.in.polygon(da$PositionX[cell],da$PositionY[cell],Sx,Sy)>0
		NetinS<-polyclip(list(list(x=Sx,y=Sy)),brainmass,op="intersect")
		areaINT<-area(NetinS)
		# for (k in 1:length(NetinS)){
			# polygon(NetinS[[k]],border=7)
		# }
	}
	
	areaRim<-areaRim+(areaEXT-areaINT)
	onRim<- onRim|(inG&!inS)
}





ris<-0*(1:11)
ris[1]<-sum(cell) #totale cellule 
ris[2]<-sum(inTum&!inBuchi) #totale cellule in tumori 
ris[3]<-sum(inBuchi) #totale cellule in buchi
ris[4]<-ris[1]-(ris[2]+ris[3]) #totale cellule nel parenchima sano
ris[5]<-sum(onRim) #totale Cellule sui bordi del tumore
ris[6]<-round(area(brainmass)/10^6,2) #area (mm2) cervello senza buchi 
ris[7]<-round(area(tumormass)/10^6,2) #area (mm2)tumori senza buchi 
ris[8]<-round(ris[7]/ris[6]*100,0) # %area (mm2) tumori 
ris[9]<-ris[2]/ris[7] #densità cellule in mm^2 ditumore 
ris[10]<-ris[4]/(ris[6]-ris[7]) #densità cellule in mm^2 di parechima 
ris[11]<-sum(onRim)/(areaRim/10^6)
ris[12]<-ris[9]/ris[10]


reportlab<-c(paste("Cellule",labCell),paste("Cellule",labCell,"nel tumore"),"Cellule nei buchi","Cellule nel parenchima sano","Cellule sui bordi del tumore","Area netta cervello  mm2","Area netta Tumore  mm2","Percentuale area tumore",paste("Densità",labCell,"in tumore"),paste("Densità",labCell,"nel parenchima sano"),paste("Densità",labCell,"nei bordi"),"Rapporto densità tum/par")
a<-data.frame(signif(ris,2))
rownames(a)<-reportlab
colnames(a)<-"valori"
write.table(a, "clipboard", sep="\t") 





