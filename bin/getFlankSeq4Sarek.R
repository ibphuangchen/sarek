#!/usr/bin/env Rscript

library(seqinr)
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-m", "--maf", type='character',default=TRUE,
                    help="MAF file")
parser$add_argument("-d", "--databse", type="character",
                    help="Protein and CDS database")
parser$add_argument("-n", "--num", type="integer", default=13,
                    help="the length of flanking sequence")
parser$add_argument("-o", "--out", type="character",
                    help="output")
args <- parser$parse_args()

getFlankSeq=function(totalMaf, ensemblSeq, flankAACount=13){
  flankList=list()
  totalMaf = totalMaf[Variant_Classification %in% c('Missense_Mutation','Nonsense_Mutation',
                                                    'Nonstop_Mutation','In_Frame_Del','Frame_Shift_Del',
                                                    'In_Frame_Ins','Frame_Shift_Ins')]
  for(i in 1:nrow(totalMaf)){
    ##For frameshift mutations, HGVSp_Short and Protein_position can be different; Protein_position and CDS_position are always the same
    ##For duplication from inframe insertion, HGVSp_Short always annotates the last AA;
    ##However, there are still inconsistency that I can't resolve, like COSV57128135 records are different from: p.A390_A391insC       RBM23
    aaPos = gsub(totalMaf$Protein_position[i],pattern = '^([0-9]+).*',replacement = '\\1')
    aaPos = as.numeric(aaPos)

		if(totalMaf$HGVSp_Short[i]==''|is.na(aaPos)) {
			warning(paste0('Cannot resolve at line ',i))
			next
		}

    startPos = ifelse(aaPos<(flankAACount+1), 1, aaPos-flankAACount) #if aaPos is less or equal flankAACount (13)
    protein = s2c(ensemblSeq[totalMaf$Transcript_ID[i],proteinSeq,on='TxID'])
    cdna = s2c(ensemblSeq[totalMaf$Transcript_ID[i],coding,on='TxID']) #this is only coding
    if(length(cdna)%%3 !=0) #ensembl seq bug
      next
    varType = totalMaf$Variant_Classification[i]
    if(varType=="Missense_Mutation"){
      refAA=gsub(totalMaf$Amino_acids[i],pattern = '/.*',replacement = '')
      newAA=gsub(totalMaf$Amino_acids[i],pattern = '.*/',replacement = '')
      
			if(aaPos>length(protein) | 
				 protein[aaPos]!=refAA |
				 aaPos==as.numeric(gsub(totalMaf$HGVSp_Short[i],
																pattern = '.*?([0-9]+).*',
																replacement = '\\1'))) {
        warning(paste0("reference sequence does not match maf file, nrow=",i,"type Missense_Mutation"))
        next
      }
      
			endPos = ifelse(aaPos> length(protein)-flankAACount, length(protein), aaPos+flankAACount)
      refSeq= c2s(protein[startPos:endPos])
      newSeq= c2s(c(protein[startPos:(aaPos-1)], newAA, protein[(aaPos+1):endPos]))
      
      aaChangeStartPepNew = aaPos - startPos +1
      aaChangeEndPepNew = aaChangeStartPepNew
      aaChangeStartPepRef = aaChangeStartPepNew
      aaChangeEndPepRef = aaChangeEndPepNew
      
    }
    #check ENST00000493964 ENST00000634670
    else if(varType %in% c('In_Frame_Ins', 'Frame_Shift_Ins')){
      insPos = as.numeric(gsub(totalMaf$CDS_position[i],pattern = '^([0-9]+).*',replacement = '\\1'))
      if (totalMaf$Reference_Allele[i]!='-'){
        warning(paste0('cannot resolve insert position at: '),i)
        next
      }
      insDNAseq=s2c(totalMaf$Tumor_Seq_Allele2[i])
      if(totalMaf$STRAND_VEP[i]==-1) insDNAseq=rev(comp(insDNAseq))
      aaChangeStartPepRef = aaPos - startPos +1
      aaChangeStartPepNew = aaChangeStartPepRef
      
      if(varType=='In_Frame_Ins'){
        newDNA=c(cdna[1:insPos],insDNAseq,cdna[(insPos+1):length(cdna)])
        newProtein=translate(newDNA)

        if(newProtein[length(newProtein)]!='*')
					warning(paste0('Possible error in line ',i,", type In_Frame_Ins"))

        endPos = ifelse(aaPos> length(protein)-flankAACount, length(protein), aaPos+flankAACount)
        refSeq=protein[startPos:endPos]
        newSeq=newProtein[startPos:(endPos+floor(length(insDNAseq)/3))]
        
        aaChangeEndPepNew = aaChangeStartPepNew+length(insDNAseq)/3
        aaChangeEndPepRef = aaChangeStartPepRef
      }else if(varType=='Frame_Shift_Ins'){
        cdna=s2c(ensemblSeq[totalMaf$Transcript_ID[i],cds3utr,on='TxID'])
        newDNA=c(cdna[1:insPos],insDNAseq,cdna[(insPos+1):length(cdna)])
        newProtein=seqinr::translate(newDNA)
        if(!grepl(totalMaf$HGVSp_Short[i],pattern = ".*fs\\*[0-9]+")) 
          next
        ##e.g. p.V156Sfs*6 will be 156+6
        ##for the frameshift, using the aaPos in HGVsp
        aaPos=as.numeric(gsub(totalMaf$HGVSp_Short[i],
                              pattern = '.*(\\*|[A-Z])([0-9]+)[A-Z]fs.*',replacement = '\\2'))
        endPos = as.numeric(gsub(totalMaf$HGVSp_Short[i],
                                 pattern = '.*fs\\*',replacement = ''))+aaPos-1
        startPos = ifelse(aaPos <= flankAACount, 1, aaPos-flankAACount)
        newSeq=newProtein[startPos:endPos]
        if(!grepl(c2s(newSeq),pattern = '\\*$')|grepl(c2s(newSeq),pattern = '\\*.+')){
          warning(paste0("error in resolving line ",i,", type Frame_Shift_Ins"))
          next
        }
        refSeq=protein[startPos:ifelse(endPos>length(protein),length(protein),endPos)]
        aaChangeEndPepNew = length(newSeq)
        aaChangeEndPepRef = length(refSeq)
      }
    }
    else if(varType %in% c('In_Frame_Del','Frame_Shift_Del')){ 
      ##ENST00000394143 #ENST00000390311 #ENST00000366522 #ENST00000371923
      delBaseStart = as.numeric(gsub(totalMaf$CDS_position[i], 
                                     pattern = '^([0-9]+).*',
                                     replacement = '\\1'))
      if(grepl(totalMaf$CDS_position[i],pattern = '-'))
        delBaseEnd = as.numeric(gsub(totalMaf$CDS_position[i], 
                                     pattern = '.*-([0-9]+)/.*',
                                     replacement = '\\1'))
      else
        delBaseEnd=delBaseStart
      delSeq=s2c(totalMaf$Reference_Allele[i])
      if(totalMaf$STRAND_VEP[i]==-1) delSeq=rev(comp(delSeq))
      
			if(all(toupper(delSeq)==cdna[delBaseStart:delBaseEnd])){
				warning(paste0('error in line ',i,', type ',varType))
				next
			}
      
      aaChangeStartPepNew = aaPos - startPos +1
      aaChangeStartPepRef = aaChangeStartPepNew
      
      if(varType=='In_Frame_Del'){
        newDNA=c(cdna[1:(delBaseStart-1)],cdna[(delBaseEnd+1):length(cdna)])
        newProtein=translate(newDNA)
        if (newProtein[length(newProtein)]!='*')
					warning(paste0('Possible error in line ',i,', type In_Frame_Del'))
        
        endPos = ifelse(aaPos> length(protein)-flankAACount, length(protein), aaPos+flankAACount)
        refSeq=protein[startPos:endPos]
        newSeq=newProtein[startPos:(endPos-floor(length(delSeq)/3))]
        
        aaChangeEndPepNew = aaChangeStartPepNew
        aaChangeEndPepRef = aaChangeStartPepRef+length(delSeq)/3
      }else if(varType=='Frame_Shift_Del'){
        cdna=s2c(ensemblSeq[totalMaf$Transcript_ID[i],cds3utr,on='TxID'])
        newDNA=c(cdna[1:(delBaseStart-1)],cdna[(delBaseEnd+1):length(cdna)])
        newProtein=translate(newDNA)
        
        if(!grepl(totalMaf$HGVSp_Short[i],pattern = ".*fs\\*[0-9]+")) 
          next #will skip for p.N751* or p.M1?
        ##e.g. p.V156Sfs*6 will be 156+6
        aaPos=as.numeric(gsub(totalMaf$HGVSp_Short[i],
                              pattern = '.*(\\*|[A-Z])([0-9]+)[A-Z]fs.*',replacement = '\\2'))
        endPos = as.numeric(gsub(totalMaf$HGVSp_Short[i],
                                 pattern = '.*fs\\*',replacement = ''))+aaPos-1
        startPos = ifelse(aaPos <= flankAACount, 1, aaPos-flankAACount)
        newSeq=newProtein[startPos:endPos]
        
        if(!grepl(c2s(newSeq),pattern = '\\*$')|grepl(c2s(newSeq),pattern = '\\*.+')){
          warning(paste0("error in resolving line ",i,", type Frame_Shift_Del"))
          next
        }
        refSeq=protein[startPos:ifelse(endPos>length(protein),length(protein),endPos)]
        aaChangeEndPepNew = length(newSeq)
        aaChangeEndPepRef = length(refSeq)
      }
    }
    else if(varType == "Nonsense_Mutation"){
      
      endPos = ifelse(aaPos> length(protein)-flankAACount, length(protein), aaPos+flankAACount)
      refSeq=protein[startPos:endPos]
      newSeq=protein[startPos:aaPos]
      aaChangeStartPepRef = aaPos - startPos +1
      aaChangeEndPepRef = length(refSeq)
      aaChangeStartPepNew = aaChangeStartPepRef
      aaChangeEndPepNew = aaChangeStartPepRef
    }
    else if(varType == "Nonstop_Mutation"){
      cdna=s2c(ensemblSeq[totalMaf$Transcript_ID[i],cds3utr,on='TxID'])
      changePos=gsub(totalMaf$HGVSc[i],pattern = 'c\\.([0-9]+).*',replacement = '\\1')
      changePos=as.numeric(changePos)
      oriBase=gsub(totalMaf$HGVSc[i],pattern = '.*[0-9]+([A-Z])>.*',replacement = '\\1')
      alterBase=gsub(totalMaf$HGVSc[i],pattern = '.*>([A-Z])$',replacement = '\\1')
      if(cdna[changePos]!=oriBase){
        warning(paste0("cannot resolve at line: ",changePos,", type NonStop_Mutation"))
        next
      }
      newDNA=c(cdna[1:(changePos-1)],alterBase,cdna[(changePos+1):length(cdna)])
      newProtein=translate(newDNA)
      startPos = ifelse(length(protein)>flankAACount, length(protein)-flankAACount, 1)
      refSeq=protein[startPos:length(protein)]
      if(length(which(newProtein=='*'))>0) 
        newEnds = which(newProtein=='*')[1]
      else 
        newEnds=length(newProtein)
      
      newSeq=newProtein[startPos:newEnds]
      aaChangeStartPepRef = length(refSeq)
      aaChangeEndPepRef = length(refSeq)
      aaChangeStartPepNew = length(refSeq)
      aaChangeEndPepNew = length(newSeq)
    }else next
    flankList[[i]]=c(c2s(refSeq), aaChangeStartPepRef,aaChangeEndPepRef,
                     c2s(newSeq), aaChangeStartPepNew,aaChangeEndPepNew,
                     totalMaf$Variant_Classification[i],
                     totalMaf$Hugo_Symbol[i],
                     totalMaf$Transcript_ID[i],
                     totalMaf$Gene[i],
                     totalMaf$dbSNP_RS[i],
                     totalMaf$Existing_variation[i],
                     totalMaf$Amino_acids[i],
                     totalMaf$HGVSp_Short[i],
                     totalMaf$CDS_position[i],
                     totalMaf$Protein_position[i],
                     totalMaf$gnomAD_AF[i],
                     totalMaf$n_ref_count[i],
                     totalMaf$n_alt_count[i],
                     totalMaf$n_depth[i],
                     totalMaf$t_ref_count[i],
                     totalMaf$t_alt_count[i],
                     totalMaf$t_depth[i])
  }
  flankList = do.call(rbind,flankList)
  flankList=data.table(flankList)
  colnames(flankList) = c('RefPep', 'RefPepStart','RefPepEnd',
                          'AltPep', 'AltPepStart', 'AltPepEnd',
                          'Variant_Class','Gene_Symbol','Transcript_ID','Gene_ID',
                          'dbSNP_ID','Existing','Amino_acids','HGVSp_Short','CDS_pos','Protein_pos',
                          'GnomeAD_AF',
                          'RefCountN','AltCountN','DepthN',
                          'RefCountT','AltCountT','DepthT')
  flankList$RefCountN=as.numeric(flankList$RefCountN)
  flankList$AltCountN=as.numeric(flankList$AltCountN)
  flankList$DepthN=as.numeric(flankList$DepthN)
  flankList$RefCountT=as.numeric(flankList$RefCountT)
  flankList$AltCountT=as.numeric(flankList$AltCountT)
  flankList$DepthT=as.numeric(flankList$DepthT)
  flankList$GnomeAD_AF=as.numeric(flankList$GnomeAD_AF)
  return(flankList)
}

flankSeqDt = getFlankSeq(fread(args$m),
                         fread(args$d),
                         args$n)

flankSeqDt[,dedup:=paste0(dbSNP_ID,'_',Transcript_ID,'_',HGVSp_Short)]
flankSeqDt_Dedup = flankSeqDt[!duplicated(dedup)]
flankSeqDt_Dedup$dedup=NULL

fwrite(flankSeqDt_Dedup, sep = '\t', quote=F, file = args$o)
