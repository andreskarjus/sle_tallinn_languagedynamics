# Measuring frequency change bias in semantic change measures

# Code to go with Karjus et al., "Two problems and solutions in 
# evolutionary corpus-based language dynamics research", 
# presented at SLE 2018 in Tallinn, the semantic change simulation
# in particular.
# See here for the code to run the topical-cultural advection model:
# https://github.com/andreskarjus/topical_cultural_advection_model

# Running these models takes about 40h on a modern laptop -
# retraining GLoVe for all 10x replicates is a massive bottleneck; 
# disabling that (just comment out or reduce replicates) 
# brings it down to 4h or so. The code is not particularly optimized.
# The lists sims_all and ranks_all (both of length 5, for each of the 5 methods) 
# that will record the (mean) self-similarities and similarity ranks 
# respectively in matrices (latter: how the original w ranks in terms of 
# similarity to the downsampled w', if it is 1st (closest synonym), or 
# lower than that (w' has moved closer to something else)).

# Packages
library(fastmatch)
library(text2vec)
library(Matrix)
library(grr)
library(compiler)
library(plotly)
library(tidyr)

#### Load data ####
# The data should be a character vector, tokenized into words, but not into documents
# A few million words at least, preferably more
perioddata = ""   # e.g. a period, unlisted, as a character vector
#perioddata = .Internal(unlist(perioddata, FALSE, FALSE)) # use if data in list format (eg multiple docs)
minc = 100   # min freq to include a word in model
nounregex = "^S:" # regex to use to filter nouns (or other POS of interest) from the lexicon for the target sample, if words are tagged; set to "." to bypass
thereps = 1:10  # how many times to replicate (resample+rerun) each word*downsample combo



#### Prepares original corpus ####
# The following code creates preparatory term matrices and does other preparations
# Feel free to just run everything from here onwards (but expect it to take time)


# functions
fullppmi = function(pmat, N,colp,rowp, positive=T){
  library(Matrix, quietly = T)
  pmat = Matrix::t(pmat)
  #pmat = (tcmfornews)
  #set.seed(1)
  #pmat = matrix(sample(c(0,0,0,0,0,0,1,10),5*10,T), 5,10, byrow=T)
  #pmat = Matrix::t(Matrix(pmat, sparse=T))
  
  # tcmrs = Matrix::rowSums(pmat)
  # tcmcs = Matrix::colSums(pmat)
  # N = sum(tcmrs)
  # colp = tcmcs/N
  # rowp = tcmrs/N
  
  pp = pmat@p+1
  ip = pmat@i+1
  tmpx = rep(0,length(pmat@x))
  for(i in 1:(length(pmat@p)-1) ){ 
    #for(i in 1:100 ){ 
    #not0 = which(pmat[, i] > 0)
    ind = pp[i]:(pp[i+1]-1)
    not0 = ip[ind]
    icol = pmat@x[ind]
    #print(icol)
    #tmp = log( (pmat[not0,i]/N) / (rowp[not0] * colp[i] ))
    tmp = log2( (icol/N) / (rowp[not0] * colp[i] ))
    tmpx[ind] = tmp
    #print(tmp)
    # tmp = ifelse(tmp < 0, 0, tmp)
    #pmat2[not0, i] = tmp
  }
  if(positive){
    tmpx[tmpx<0] = 0 
  }
  pmat2 = pmat
  pmat2@x = tmpx
  pmat2 = Matrix::t(pmat2) # tagasi sõnad ridadeks
  pmat2 = drop0(pmat2,is.Csparse=T) # nullid: võtab oma 100 mega maha, väga hea
  return(pmat2)
}

someppmi = function(NN, tcm1=tcm[words,],
                    normvoc_all=normvoc, 
                    normvoc_sub=normvoc[words]){
  ptcm = tcm1/NN   # P(w,c)
  for(w in 1:nrow(ptcm)){
    tmp = pmax(0, log2( ptcm[w, ] / 
                          ((normvoc_sub[w]*normvoc_all))  #P(w)*P(c)
    ) )
    tmp[is.nan(tmp)]=0
    ptcm[w,] = tmp
  }
  return(ptcm)
}
someppmi = cmpfun(someppmi, options=list(optimize=3))

dorelevants = function(ppmimat,relevancy_threshold){
  relevants = list()
  relevants = list(); length(relevants) = nrow(ppmimat)
  ppmimat = Matrix::as.matrix(ppmimat)
  for(x in 1:nrow(ppmimat)){
    #y=sort(ppmimat[x,-(union(which(ppmimat[x,]==0), excludeself[x]))], decreasing=T)
    # self-ppmi is zero, fixed in ppmimat (old: thought text2vec handles itself)
    tmp = ppmimat[x,-which(ppmimat[x,]==0)]
    y=tmp[rev(order2(tmp))] # sort2 is decreasing=F BUT MESSES UP NAMES, USE order2
    y=y[1:min(length(y),relevancy_threshold)] # but this need top ones, so rev
    relevants[[x]] = y; names(relevants[[x]]) = names(y)
  }
  return(relevants)
}
dorelevants = cmpfun(dorelevants, options=list(optimize=3))
dosim =  function(relr, relc){
  #if(length(relr)>0 & length(relc)>0){   # doesn't seem to occur, always some context
  # requires fastmatch
  r1 = fmatch(relc, relr, nomatch = 0) # the intersection part; remove nomatch 0s below
  r2 = fmatch(relr[r1], relc, nomatch = 0)
  return(sum( 1/(colMeans(rbind(r1[r1>0],r2))))) #
  #   } else return(NA)
}


# do vocab
maxn = length(perioddata) # COHA last decade is 11 023 006
it = itoken(list(perioddata), progressbar = FALSE)
voc = prune_vocabulary(create_vocabulary(it), term_count_min = minc)
rownames(voc) = voc$term
vectorizer = vocab_vectorizer(voc) # needed later too
tcm = create_tcm(it, vectorizer, skip_grams_window = 5)
tcm = tcm + Matrix::t(triu(tcm))
Matrix::diag(tcm) = 0
#nrow(voc)

#### Chooses words ####
# Randomly samples from the lexicon, but from equally spaced log freq bands
vtmp = voc[voc$term_count>499 & grepl(nounregex, voc$term),]
vtmp1 = vtmp$term_count; names(vtmp1) = rownames(vtmp); rm(vtmp)
vc = cut(log(vtmp1), 116) # adjust to get different number of words (some bands may be empty)
words=c()
for(v in seq_along(levels(vc))){
  try({words = c(words, sample(names(vtmp1[vc%in%levels(vc)[v]]), 1))})
}
#plot(vtmp1[words])
print(paste("Sampled", length(words), "words" ))

words_prime = paste0(words, "_prime")
subsetsizes = c(0,  seq(0.1,7,0.3) ) # starts with 0-change sanity test
actualns = voc[sapply(words, match, voc$term), "term_count"]
randomns = vector("list", length(words))
for(j in seq_along(words)){
  for(r in seq_along(subsetsizes)){
    tmp = round(exp( log(actualns[j]) - subsetsizes[r])) # calculates subsample size
    if(tmp >= 10) { 
      randomns[[j]][r] = tmp
      names(randomns[[j]])[r] = paste0(words_prime[j], "_", subsetsizes[r])
    }
  }
}
#length(unlist(randomns))

#### Prepares LSA, PPMI matrices ####

# LSA full:
lsa_1 = LatentSemanticAnalysis$new(300)
d1 = lsa_1$fit_transform(tcm) # embeds (symmetric tcm, so no difference which way)

# PPMI for chosen vectors:
# 10085666 after n<100 pruning
N = sum(voc$term_count)
normvoc = voc$term_count/N; names(normvoc) = voc$term
ptcm = someppmi(N,tcm1=tcm[words,], normvoc_all=normvoc, normvoc_sub=normvoc[words])
relevants = dorelevants(ptcm[words,], 100); names(relevants) = words
maxapsyn = dosim(1:100, 1:100) # divide by that to get [0,1] measure in apsyn

# for ppmi ranks
ptcm_full = fullppmi(pmat=tcm, N=N, colp=normvoc,rowp=normvoc)
allrelevants = dorelevants(ptcm_full, 100); names(allrelevants) = rownames(ptcm_full)

# since full tcm is created, just fetch +/5 snippets with prime in the middle
# don't need to replace, can just create temp tcm and rename row, and add

# prep dataframes:
allwords = voc$term  # full vocab
sims_all = vector("list", 5)
ranks_all = vector("list", 5)
for(l in 1:5){
  sims_all[[l]] = matrix(NA, ncol=length(subsetsizes), nrow=length(words), dimnames = list(words,subsetsizes))
  ranks_all[[l]] = matrix(NA, ncol=length(subsetsizes), nrow=length(words), dimnames = list(words,subsetsizes))
}
names(sims_all) = c("plaincos_all" ,
                    "ppmicos_all" ,
                    "ppmiapsyn_all", 
                    "lsacos_all" ,
                    "glovecos_all" )
names(ranks_all) = names(sims_all)
# template:
emptymat = as(matrix(nrow=0,ncol=ncol(tcm)), "dgCMatrix")


##### Run models #####
# (one ugly for-loop)
thewords = seq_along(words) 
gc()
startt = Sys.time()
print("")
for(i in thewords ){  # for each word
  cat(paste(" ", i, words[i]) )
  
  # actual corpus positions:
  actualmatches = .Internal(which(perioddata %fin% words[i]))
  # matrices for reps:
  repsims = replicate(5,matrix(NA, ncol=length(subsetsizes), nrow=length(thereps)), simplify = F)
  names(repsims) = 
    c("plaincos" ,
      "ppmicos" ,
      "ppmiapsyn", 
      "lsacos",
      "glovecos") 
  repranks = replicate(5,matrix(NA, ncol=length(subsetsizes), nrow=length(thereps)), simplify = F)
  names(repranks) = 
    c("plaincos" ,
      "ppmicos" ,
      "ppmiapsyn", 
      "lsacos",
      "glovecos") 
  
  for(reps in thereps){   # 10 replications each
    cat(".")
    tcmrows = emptymat
    for(r in seq_along(randomns[[i]])){  # not all words have the same n of subsample values
      try({
        samp = actualmatches[ sample.int(actualns[i], randomns[[i]][r] ) ] # from 1:actualns[i]
        # sample from 1:actualns[i] randomly to replace
        sampfrom = pmax(samp-5, 1); sampto = pmin(samp+5, maxn) 
        # for manual windows; but fixes out-of-bounds indexes
        samplist = vector("list", randomns[[i]][r] )
        for(s in seq_along(samplist)){
          samplist[[s]] = sampfrom[s] : sampto[s] #perioddata[sampfrom[s] : sampto[s]]
        }
        samplist = sort(unique(unlist(samplist))) # indices
        # indices are concatenated but unique'd and sorted since close-by occurrences can create doubling
        # when their windows overlap; removing overlaps messes up the co-occurrence
        # counts of other words in the temp tcm, but it does not matter
        # since only the target _prime row is relevant anyway
        # tested: the vector from the 0-removed sanity sample comes out 100% like the original
        
        it_tmp = itoken(list(perioddata[samplist]), progressbar = FALSE)
        # "create dtm using same vectorizer, but new iterator"
        #  won't include low freq words (but based on global freqs not from this sample!)
        tcm_tmp = create_tcm(it_tmp, vectorizer, skip_grams_window = 5)
        # need to make matrix symmetric/full for this to work:
        tcm_tmp = tcm_tmp + Matrix::t(triu(tcm_tmp))
        Matrix::diag(tcm_tmp) = 0
        tcm_tmp = tcm_tmp[words[i], ,drop=F] # extract relevant row
        rownames(tcm_tmp) = names(randomns[[i]][r]) # just one now
        tcmrows = rbind(tcmrows, tcm_tmp)
      })
    } # end subset s
    
    # test integrity 
    #all(tcmrows[1,] == tcm[words[i],,drop=T])
    #which(tcmrows[1,] != tcm[words[i],,drop=T])
    
    try({
      # plain + cosine
      repsims$plaincos[reps,1:length(randomns[[i]]) ] = sim2(tcm[words[i],,drop=F], tcmrows)[1,]
      simstmp = sim2(tcm, tcmrows) # rows=allwords, cols=subsamplewords; takes a moment
      repranks$plaincos[reps,1:length(randomns[[i]]) ] = 
        sapply(seq_along(randomns[[i]]), function(x) 
          fmatch(words[i], allwords[rev(order2(simstmp[,x])) ] ) 
          # finds sim rank of word_i to word_i_prime_x
        )
      
      #ppmi + cosine
      ppmitmp = someppmi(N, tcm1=tcmrows,normvoc_all=normvoc, normvoc_sub=randomns[[i]]/N)
      repsims$ppmicos[reps, 1:length(randomns[[i]])] = sim2(ptcm[words[i],,drop=F], ppmitmp)[1,]
      simstmp = sim2(ptcm_full, ppmitmp) # rows=allwords, cols=subsamplewords; takes a moment
      repranks$ppmicos[reps,1:length(randomns[[i]]) ] = 
        sapply(seq_along(randomns[[i]]), function(x) {
          #print(allwords[rev(order2(simstmp[,x])) ][1]) # debug
          fmatch(words[i], allwords[rev(order2(simstmp[,x])) ] ) 
          # finds sim rank of word_i to word_i_prime_x
        }
        )
      
      # ppmi + apsyn
      apsyntmp = dorelevants(ppmitmp, 100)
      for(w in 1:length(apsyntmp)){
        repsims$ppmiapsyn[reps, w] = dosim(names(apsyntmp[[w]]), names(relevants[[i]]))/maxapsyn # [0,1] range
        # for ranks:
        tmpalls = rep(NA, length(allrelevants)); names(tmpalls)=allwords
        for(w2 in 1:length(allrelevants)){
          tmpalls[w2] = dosim(names(apsyntmp[[w]]), names(allrelevants[[w2]]))/maxapsyn
        }
        repranks$ppmiapsyn[reps, w] = fmatch(words[i], allwords[rev(order2(tmpalls)) ] ) 
      }
      
      # lsa + cosine
      # non-self mean=0.05; 1st closest non-self, mean=0.61 - not very dense.
      lsatmp = lsa_1
      x = transform(tcmrows, lsatmp)
      repsims$lsacos[reps, 1:length(randomns[[i]]) ] = sim2(d1[words[i],,drop=F], x)[1,]
      
      simstmp = sim2(d1, x) # rows=allwords, cols=subsamplewords; takes a moment
      repranks$lsacos[reps,1:length(randomns[[i]]) ] = 
        sapply(seq_along(randomns[[i]]), function(x) 
          fmatch(words[i], allwords[rev(order2(simstmp[,x])) ] ) 
          # finds sim rank of word_i to word_i_prime_x
        )
      
      # run only a few replicates if very time-consuming
      #if(reps %in% 1){
        # glove+cosine
        emptycols = as(matrix(0, nrow=nrow(tcmrows),ncol=nrow(tcmrows)), "dgCMatrix")
        colnames(emptycols) = rownames(tcmrows)
        gtcm = cbind(rbind(tcm,tcmrows), rbind(Matrix::t(tcmrows),emptycols ) )
        # this seems the right concatenation, but the sanity sample does not come
        # out with 100% identity, cosine is slightly (-0.0001) below 1. maybe a glove problem?
        # also apparently high freq words with 6-7 log reduction get negative cosines
        # (that should be normal)
        glove_model = GloVe$new(word_vectors_size = 100,
                                vocabulary = c(voc$term,rownames(tcmrows)), 
                                x_max = 10)
        word_vectors_main = glove_model$fit_transform(gtcm, n_iter = 20L, 
                                                      n_check_convergence=5L, convergence_tol= 0.006) # takes time
        word_vectors_context = glove_model$components
        word_vectors = word_vectors_main + t(word_vectors_context) # the recommended addition step
        # prime sims:
        repsims$glovecos[reps, 1:length(randomns[[i]]) ] = 
          sim2(word_vectors[words[i],,drop=F], word_vectors[rownames(tcmrows),])[1,]
        # prime ranks:
        simstmp = sim2(word_vectors[1:(nrow(word_vectors)-length(randomns[[i]])),  ], # all w
                       word_vectors[(nrow(word_vectors)-length(randomns[[i]])+1):nrow(word_vectors) ,  ] ) # vs prime w's
        repranks$glovecos[reps,1:length(randomns[[i]]) ] = 
          sapply(seq_along(randomns[[i]]), function(x) 
            fmatch(words[i], allwords[rev(order2(simstmp[,x])) ] ) 
            # finds sim rank of word_i to word_i_prime_x
          )
        
        # fasttext+cosine
          # expects single line of words separated by spaces (so unlisted is fine)
          # also need a re-do for EACH subsample, cannot use tcm, need to replace in
          # actual every time, and then rerun fasttext, meaning 100*25 models (*10 replicates), so maybe later
      #} # end glove
    })
  } # end rep
  
  # record scores:
  sims_all$plaincos_all[i,] = colMeans(repsims$plaincos, na.rm = T)
  sims_all$ppmicos_all[i,] = colMeans(repsims$ppmicos, na.rm = T)
  sims_all$ppmiapsyn_all[i,] = colMeans(repsims$ppmiapsyn, na.rm = T)
  sims_all$lsacos_all[i,] = colMeans(repsims$lsacos, na.rm = T)
  sims_all$glovecos_all[i,] = colMeans(repsims$glovecos, na.rm = T)
  
  ranks_all$plaincos_all[i,] = colMeans(repranks$plaincos, na.rm = T)
  ranks_all$ppmicos_all[i,] = colMeans(repranks$ppmicos, na.rm = T)
  ranks_all$ppmiapsyn_all[i,] = colMeans(repranks$ppmiapsyn, na.rm = T)
  ranks_all$lsacos_all[i,] = colMeans(repranks$lsacos, na.rm = T)
  ranks_all$glovecos_all[i,] = colMeans(repranks$glovecos, na.rm = T)
  
} # end words loop

print(paste("Started at", startt))
print(paste("Finished at", Sys.time()))



