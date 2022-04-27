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

getFlankSeq=function(totalMaf, cdnaProteinsEnsemble, flankAACount=13){
  flankList=list()
  for(i in 1:nrow(totalMaf)){
    if(totalMaf$HGVSp_Short[i]=='') next
    ##For fs mutations, HGVSp_Short and Protein_position can be different; Protein_position and CDS_position are always the same
    ##For duplication from inframe insertion, HGVSp_Short always annotates the last AA;
    ##However, there are still inconsistency that I can't resolve, like COSV57128135 records are different from: p.A390_A391insC       RBM23
    aaPos = gsub(totalMaf$Protein_position[i],pattern = '^([0-9]+).*',replacement = '\\1')
    aaPos = as.numeric(aaPos)
    startPos = ifelse(aaPos<(flankAACount+1), 1, aaPos-flankAACount) #if aaPos is less or equal flankAACount (13)
    protein = seqinr::s2c(cdnaProteinsEnsemble[totalMaf$Transcript_ID[i],proSeq,on='txID_s'])
    varType = totalMaf$Variant_Classification[i]
    if(varType=="Missense_Mutation"){
      if(totalMaf$Variant_Type[i]!='SNP'){
        warning(paste0("cannot resolve the mutation type at #row:", i))
        next
      }
      refAA=gsub(totalMaf$Amino_acids[i],pattern = '/.*',replacement = '')
      newAA=gsub(totalMaf$Amino_acids[i],pattern = '.*/',replacement = '')
      #stopifnot(protein[aaPos]==refAA)
      if(aaPos>length(protein) | protein[aaPos]!=refAA) {
        warning(paste0("reference sequence does not match maf file, nrow=",i))
        next
      }
      stopifnot(aaPos==as.numeric(gsub(totalMaf$HGVSp_Short[i],
                                       pattern = '.*?([0-9]+).*',
                                       replacement = '\\1')))
      endPos = ifelse(aaPos> length(protein)-13, length(protein), aaPos+13)
      refSeq= c2s(protein[startPos:endPos])
      newSeq= c2s(c(protein[startPos:(aaPos-1)], newAA, protein[(aaPos+1):endPos]))
      
      aaChangeStartPepNew = aaPos - startPos +1
      aaChangeEndPepNew = aaChangeStartPepNew
      aaChangeStartPepRef = aaChangeStartPepNew
      aaChangeEndPepRef = aaChangeEndPepNew
      
    }
    #check ENST00000493964 ENST00000634670
    else if(varType %in% c('In_Frame_Ins', 'Frame_Shift_Ins')){
      cds = seqinr::s2c(cdnaProteinsEnsemble[totalMaf$Transcript_ID[i],cds,on='txID_s']) 
      cdsUTR = seqinr::s2c(cdnaProteinsEnsemble[totalMaf$Transcript_ID[i],cds_3utr,on='txID_s'])
      if(length(cds)%%3 !=0) #ensembl seq bug
        next
      insPos = as.numeric(gsub(totalMaf$CDS_position[i],pattern = '^([0-9]+).*',replacement = '\\1'))
      stopifnot(totalMaf$Reference_Allele[i]=='-')
      insDNAseq=s2c(totalMaf$Tumor_Seq_Allele2[i])
      if(totalMaf$STRAND_VEP[i]==-1) insDNAseq=rev(comp(insDNAseq))
      aaChangeStartPepRef = aaPos - startPos +1
      aaChangeStartPepNew = aaChangeStartPepRef
      
      if(varType=='In_Frame_Ins'){
        newDNA=c(cds[1:insPos],insDNAseq,cds[(insPos+1):length(cds)])
        newProtein=seqinr::translate(newDNA)
        #stopifnot(newProtein[length(newProtein)]=='*')
        
        endPos = ifelse(aaPos> length(protein)-13, length(protein), aaPos+13)
        refSeq=protein[startPos:endPos]
        newSeq=newProtein[startPos:(endPos+floor(length(insDNAseq)/3))]
        
        aaChangeEndPepNew = aaChangeStartPepNew+length(insDNAseq)/3
        aaChangeEndPepRef = aaChangeStartPepRef
      }
      else{#frameshift case
        newDNA=c(cdsUTR[1:insPos],insDNAseq,cdsUTR[(insPos+1):length(cdsUTR)])
        newProtein=seqinr::translate(newDNA)
        
        if(!grepl(totalMaf$HGVSp_Short[i],pattern = ".*fs\\*[0-9]+")) next #will skip for p.N751* or p.M1?
        ##e.g. p.V156Sfs*6 will be 156+6
        ##for the frameshift, using the aaPos in HGVsp
        aaPos=as.numeric(gsub(totalMaf$HGVSp_Short[i],pattern = '.*(\\*|[A-Z])([0-9]+)[A-Z]fs.*',replacement = '\\2'))
        endPos = as.numeric(gsub(totalMaf$HGVSp_Short[i],pattern = '.*fs\\*',replacement = ''))+aaPos-1
        startPos = ifelse(aaPos<(flankAACount+1), 1, aaPos-flankAACount)
        newSeq=newProtein[startPos:endPos]
        refSeq=protein[startPos:ifelse(endPos>length(protein),length(protein),endPos)]
        aaChangeEndPepNew = length(newSeq)
        aaChangeEndPepRef = length(refSeq)
      }
    }
    else if(varType %in% c('In_Frame_Del','Frame_Shift_Del')){ 
      cds = seqinr::s2c(cdnaProteinsEnsemble[totalMaf$Transcript_ID[i],cds,on='txID_s']) 
      cdsUTR = seqinr::s2c(cdnaProteinsEnsemble[totalMaf$Transcript_ID[i],cds_3utr,on='txID_s'])
      if(length(cds)%%3 !=0) #ensembl seq bug
        next
      ##based on CDS_position column
      delBaseStart1 = as.numeric(gsub(totalMaf$CDS_position[i], 
                                     pattern = '^([0-9]+).*',
                                     replacement = '\\1'))
      if(grepl(totalMaf$CDS_position[i],pattern = '-'))
        delBaseEnd1 = as.numeric(gsub(totalMaf$CDS_position[i], 
                                     pattern = '.*-([0-9]+)/.*',
                                     replacement = '\\1'))
      else
        delBaseEnd1=delBaseStart1
      
      ##based on HGVSc
      delBaseStart2 = as.numeric(gsub(totalMaf$HGVSc[i], 
                                      pattern = 'c\\.([0-9]+).*',
                                      replacement = '\\1'))
      if(grepl(totalMaf$HGVSc[i],pattern = '_'))
        delBaseEnd2 = as.numeric(gsub(totalMaf$HGVSc[i], 
                                      pattern = '.*_([0-9]+)del',
                                      replacement = '\\1'))
      else
        delBaseEnd2=delBaseStart2
      
      delSeq=s2c(totalMaf$Reference_Allele[i])
      if(totalMaf$STRAND_VEP[i]==-1) delSeq=rev(comp(delSeq))
      if(all(toupper(delSeq)==cds[delBaseStart1:delBaseEnd1])){
        delBaseStart=delBaseStart1
        delBaseEnd=delBaseEnd1
      }else if(all(toupper(delSeq)==cds[delBaseStart2:delBaseEnd2])){
        delBaseStart=delBaseStart2
        delBaseEnd=delBaseEnd2
      }else{
        warning(paste0('cannot resolve delete seq at #row: ',i))
        next
      }
      
      aaChangeStartPepRef = aaPos - startPos +1
      aaChangeStartPepNew = aaChangeStartPepRef
      
      if(varType=='In_Frame_Del'){
        newDNA=c(cds[1:(delBaseStart-1)],cds[(delBaseEnd+1):length(cds)])
        newProtein=seqinr::translate(newDNA)
        #stopifnot(newProtein[length(newProtein)]=='*')
        
        endPos = ifelse(aaPos> length(protein)-13, length(protein), aaPos+13) #for the ref, AA endPos is predictable
        refSeq = protein[startPos:endPos]
        newSeq=newProtein[startPos:(endPos-floor(length(delSeq)/3))]
        
        aaChangeEndPepNew = aaChangeStartPepNew
        aaChangeEndPepRef = aaChangeStartPepRef+length(delSeq)/3
      }else{ #frame shift case
        newDNA=c(cdsUTR[1:(delBaseStart-1)],cdsUTR[(delBaseEnd+1):length(cdsUTR)])
        newProtein=seqinr::translate(newDNA)
        
        if(!grepl(totalMaf$HGVSp_Short[i],pattern = ".*fs\\*[0-9]+")) next #will skip for p.N751* or p.M1?
        ##e.g. p.V156Sfs*6 will be 156+6
        aaPos=as.numeric(gsub(totalMaf$HGVSp_Short[i],pattern = '.*(\\*|[A-Z])([0-9]+)[A-Z]fs.*',replacement = '\\2'))
        endPos = as.numeric(gsub(totalMaf$HGVSp_Short[i],pattern = '.*fs\\*',replacement = ''))+aaPos-1
        startPos = ifelse(aaPos<(flankAACount+1), 1, aaPos-flankAACount)
        newSeq=newProtein[startPos:endPos]
        refSeq=protein[startPos:ifelse(endPos>length(protein),length(protein),endPos)]
        aaChangeEndPepNew = length(newSeq)
        aaChangeEndPepRef = length(refSeq)
      }
    }
    else if(varType == "Nonsense_Mutation"){
      endPos = ifelse(aaPos> length(protein)-13, length(protein), aaPos+13)
      refSeq=protein[startPos:endPos]
      newSeq=""
      aaChangeStartPepRef = aaPos
      aaChangeEndPepRef = length(refSeq)
      aaChangeStartPepNew = 0
      aaChangeEndPepNew = 0
    }else
			next
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

flankSeqDt = getFlankSeq(totalMaf = fread(args$m),
                         cdnaProteinsEnsemble = fread(args$d),
                         flankAACount = args$n)

fwrite(flankSeqDt, sep = '\t', file = args$o)
