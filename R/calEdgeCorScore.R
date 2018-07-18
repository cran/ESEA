calEdgeCorScore <- function(dataset, class.labels,controlcharactor, edgesbackgrand, k=3, k_iter_max=10) {
#This function calculate differential Mutual information
#Inputs:
#     dataset: Marix of gene expression values (rownames are genes, columnnames are samples) 
#     class.labels: Vector of binary labels.
#     controlcharactor: Charactor of control in the class labels.
#     edgesbackgrand: Marix of the edges' backgrand   
#Outputs:
#     EdgeCorScore:Vector of the aberrant correlation in phenotype P based on mutual information (MI) for each edge
#        
       location<-matrix(0,length(edgesbackgrand[,1]),2)
	   location[,1]<-match(edgesbackgrand[,1],rownames(dataset))
	   location[,2]<-match(edgesbackgrand[,2],rownames(dataset))
	   location<-na.omit(location)
	   controlloca<-which(class.labels==controlcharactor)
	   dataset.1<-dataset[location[,1],]
	   dataset.2<-dataset[location[,2],]
	   Cexpress.1<-dataset[location[,1],controlloca]
	   Cexpress.2<-dataset[location[,2],controlloca]
	   EdgeCorScore<-c()
	   for(i in 1:length(location[,1])){
	   		dataset_1_dis = arules::discretize(x=as.numeric(dataset.1[i,]), method="cluster", centers=k,iter.max=k_iter_max)
			dataset_2_dis = arules::discretize(x=as.numeric(dataset.2[i,]), method="cluster", centers=k,iter.max=k_iter_max)
			dataset_MI = infotheo::mutinformation(dataset_1_dis,dataset_2_dis)

			Cexpress_1_dis = arules::discretize(x=as.numeric(Cexpress.1[i,]), method="cluster", centers=k,iter.max=k_iter_max)
			Cexpress_2_dis = arules::discretize(x=as.numeric(Cexpress.2[i,]), method="cluster", centers=k,iter.max=k_iter_max)
			Cexpress_MI = infotheo::mutinformation(Cexpress_1_dis,Cexpress_2_dis)

		   	EdgeCorScore[i]<-dataset_MI-Cexpress_MI
		  }   
		EdgeID<-matrix(0,length(location[,1]),2)
		EdgeID[,1]<-row.names(dataset.1)
		EdgeID[,2]<-row.names(dataset.2)
		colnames(EdgeID)<-c("Edge1","Edge2")
		Paste<-function(x,c1,c2) paste(x[c1],x[c2],sep="|")
		names(EdgeCorScore)<-apply(EdgeID,1,Paste,c1="Edge1",c2="Edge2")
		return(EdgeCorScore)
	   }