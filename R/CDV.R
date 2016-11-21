####################################################################

# CDV package

#' Run the analysis
#'
#' This function gets a fasta or multi-fasta file containing full-length
#' hemmaglutinin nucleotide sequences from CDV and performs the assignment
#' to the currently described lineages and sublineages, according to their
#' phylogenetic position.
#' The function outputs various results including: i) the new alignment, ii) 
#' the new phylogenetic tree, iii) a table describing the assignments for each
#' input sequence and iii) a PDF with a schematic representation of the tree.
#'
#' @param test.sequences is the input fasta or multi-fasta file
#' @param output.prefix is a prefix for naming the output files
#' @export

CDV<-function(test.sequences,output.prefix){
	
	options(warn=-1)
	options(getClass.msg=FALSE)
			   	
	# Dependencies
	
	print('Loading dependencies...')
	
	suppressMessages(require(msa))
	suppressMessages(require(seqinr))
	suppressMessages(require(ape))
	suppressMessages(require(phangorn))
	suppressMessages(require(phylobase))
	suppressMessages(require(ggtree))
			   				   	
	# Generate new dataset 

	print('Generating new dataset...')

	wdir<-getwd()
			   
	reference<-system.file('data','cdv_reference.fasta',
						   package='CDV')
						             
	rfasta<-read.fasta(reference)
	rsequs<-lapply(getSequence(rfasta),toupper)
	rnames<-getName(rfasta)
			   
	tfasta<-read.fasta(test.sequences)
	tsequs<-lapply(getSequence(tfasta),toupper)
	tnames<-getName(tfasta)
			   
	nfasta<-c(rsequs,tsequs)
	nnames<-c(rnames,tnames)
			   
	seqname<-paste(output.prefix,'.sequences.fasta',sep='')
			   
	write.fasta(nfasta,names=nnames,file=seqname)

	# Align new dataset
	
	print('Aligning sequences...')
			   
	sequences<-readDNAStringSet(seqname)
	alignment<-msa(sequences,method='Muscle')

	ali<-msaConvert(alignment,'seqinr::alignment')
	    	   
	aliname<-paste(output.prefix,'.alignment.fasta',sep='')
	
	write.fasta(as.list(ali$seq),names=ali$nam,file=aliname)
	
	# Extract lineage markers
	
	print('Extracting markers...')
	
	markers<-system.file('data','aa_markers.csv',package='CDV')
	amimark<-read.table(markers,sep='\t',fill=F,header=T)
	posmark<-as.numeric(as.vector(amimark$Site))
	
	wtest<-which(ali$nam%in%tnames)
	stest<-ali$seq[wtest]
	ntest<-ali$nam[wtest]
	
	tlate<-lapply(lapply(stest,s2c),function(x){translate(x,numcode=1)})
	
	nmarkers<-NULL
		
	for (e in 1:length(tlate)){
		
		tmarkers<-tlate[[e]][posmark]
		uniques<-apply(amimark[,-1],1,unique)
		
		for (u in 1:length(uniques)){
			
			ift<-tmarkers[u]%in%uniques[[u]]
		
			if (ift==FALSE){
				
				tmarkers[u]<-''
			}
		}
		
		nmarkers<-cbind(nmarkers,tmarkers)
	}
	
	colnames(nmarkers)<-ntest
	
	nmarkers<-cbind(amimark,nmarkers)
	
	markname<-paste(output.prefix,'.markers.csv',sep='')
	
	write.table(nmarkers,sep='\t',quote=F,row.names=F,na='',file=markname)
	
	# Make the tree
	
	print('Making phylogeny...')
	
	phydat<-msaConvert(alignment,'phangorn::phyDat')
	distam<-dist.ml(phydat,model='F81')
	phytre<-NJ(distam)

	edges<-phytre$edge.length
	edges[which(edges<=0)]<-1e-10
	phytre$edge.length<-edges

	# Define reference lineages
	
	print('Identifying lineages...')

	reflin<-system.file('data',
						'reference_lineages.Rdata',
						package='CDV')
	load(reflin)
	
	lineages<-paste('L',1:7,sep='')
	sublineages<-c('L1.1','L1.2','L2.1','L2.2','L3.1','L3.2',
				   'L3.3','L4.1','L4.2','L5.1','L5.2','L6.1',
				   'L6.2','L6.4','L7.1','L7.2')

	rl<-reference_lineages

	for (l in lineages){

		aux<-as.vector(rl[which(rl$lineage==l),'strain'])
		assign(l,aux)		
	}
	
	for (s in sublineages){
		
		aux<-as.vector(rl[which(rl$sublineage==s),'strain'])
		assign(s,aux)	
	}

	# Refine the tree
	
	n7<-MRCA(phytre,L7)
	phytre2<-reroot(phytre,node=n7)
	
	treename<-paste(output.prefix,'.tree.nwk',sep='')

	write.tree(phytre2,treename)
	
	# Define ancestors

	for (l in lineages){
		
		aux<-MRCA(phytre2,get(l))
		assign(gsub('L','n',l),aux)
	}

	for (s in sublineages){
		
		aux<-MRCA(phytre2,get(s))
		assign(gsub('L','n',s),aux)	
	}
	
	# Look for descendants
	
	phytre4<-as(phytre2,'phylo4')

	lindes<-list()
	subdes<-list()
	
	for (x in 1:length(lineages)){

		n<-gsub('L','n',lineages[x])
		lindes[[x]]<-names(descendants(phytre4,node=get(n)))
	}
	
	for (y in 1:length(sublineages)){
	
		s<-gsub('L','n',sublineages[y])
		subdes[[y]]<-names(descendants(phytre4,node=get(s)))
	}
	
	names(lindes)<-lineages
	names(subdes)<-sublineages	

	# Classify sequences
	
	print('Classifying sequences..')
	
	result<-NULL
	
	for (e in tnames){
		
		classification<-c(e,NA,NA,NA)
		
		gln<-grep(e,lindes)
		gsb<-grep(e,subdes)
		lle<-length(gln)
		sle<-length(gsb)
		
		if ((lle>0 & sle>0)){
			
			lin<-toupper(paste(names(lindes)[gln],collapse='/'))
			sli<-toupper(paste(names(subdes)[gsb],collapse='/'))
			
			classification[2]<-lin
			classification[3]<-sli
			classification[4]<-'A'
		
		} else if ((lle>0 & sle==0)){
			
			lin<-toupper(paste(names(lindes)[gln],collapse='/'))
			classification[2]<-lin
			classification[4]<-'NS'
			
		} else if ((lle==0)){
			
			classification[4]<-'NL'
		}
		
		result<-rbind(result,classification)
		rownames(result)<-NULL
	}
	
	result<-as.data.frame(result)
	colnames(result)<-c('strain','lineage',
						'sublineage','assignment')
	
	tabname<-paste(output.prefix,'.table.csv',sep='')
	
	write.table(result,quote=F,row.names=F,sep='\t',file=tabname)
	
	# Draw phylogeny
	
	print('Generating output files...')
	
	a<-as.vector(result[which(result$assignment=='A'),'strain'])
	ns<-as.vector(result[which(result$assignment=='NS'),'strain'])
	nl<-as.vector(result[which(result$assignment=='NL'),'strain'])
	
	nodes<-c(n1,n2,n3,n4,n5,n6,n7)
	  
	tiplist<-list(a,ns,nl)
	names(tiplist)<-c('A','NS','NL')
	  
	p<-ggtree(phytre2,col='grey90') + xlim(NA,.1)

	ncolors<-c('cyan','orange','salmon','grey50',
			   'green','yellow','purple')

	phyloname<-paste(output.prefix,'.tree.pdf',sep='')
		
	pdf(phyloname)

	p + 

	geom_cladelabel(node=n1,label='Lineage 1',align=T,col=ncolors[1]) +
	geom_cladelabel(node=n2,label='Lineage 2',align=T,col=ncolors[2]) +
	geom_cladelabel(node=n3,label='Lineage 3',align=T,col=ncolors[3]) +
	geom_cladelabel(node=n4,label='Lineage 4',align=T,col=ncolors[4]) +
	geom_cladelabel(node=n5,label='Lineage 5',align=T,col=ncolors[5]) +
	geom_cladelabel(node=n6,label='Lineage 6',align=T,col=ncolors[6]) +
	geom_cladelabel(node=n7,label='Lineage 7',align=T,col=ncolors[7]) +

	geom_point2(aes(label=node,subset=(node%in%nodes)),
		alpha=.5,color=ncolors[c(1,3,2,5,4,6,7)], size=6) +

	geom_point2(aes(label=label,subset=(label%in%a)),
		size=1.5,shape=23,color='blue') +

	geom_point2(aes(label=label,subset=(label%in%nl)),
		size=1.5,shape=24,color='magenta') +

	geom_point2(aes(label=label,subset=(label%in%ns)),
		size=1.5,shape=21,color='black') + scale_shape_identity()

	dev.off() 
	
	print('Job completed!')
}

####################################################################