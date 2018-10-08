#### Measuring frequency change bias in semantic change measures ####

# Code to go with Karjus et al., "Two problems and solutions in 
# evolutionary corpus-based language dynamics research", 
# presented at SLE 2018 in Tallinn, the semantic change simulation
# in particular.
# See here for the code to run the topical-cultural advection model:
# https://github.com/andreskarjus/topical_cultural_advection_model

# Running these models takes about 40h on a modern laptop -
# retraining GLoVe for all 10x replicates is a massive bottleneck; 
# reducing that the replicates for that brings it down to 4h or so. 
# The code is not particularly well optimized.
# Running/sourcing this script will produce two lists (see below for a parameter
# combo that runs fast for 1 word, to see an example).
# The lists sims_all and ranks_all (both of length 5, for each of the 5 methods) 
# will record the (mean) self-similarities and similarity ranks 
# respectively in matrices (latter: how the original w ranks in terms of 
# similarity to the downsampled w', if it is 1st (closest synonym), or 
# lower than that (w' has moved closer to something else)).
# If using RStudio to run this, turn on Soft wrap in Global Options->Code.

# Loads and if needed installs required packages (prints package version but does not check for updates of present packages; update if necessary):
invisible(sapply( c( "fastmatch","text2vec","Matrix","grr","compiler" ), function(x) if(require(x, character.only = T)){print(paste(x, packageVersion(x)),quote=F)}else{install.packages(x);library(x)} ))


#### Load data, define params ####
#
# The data should be a corpus as a single character vector (not list, ie no doc or sentence tokenization), tokenized into words, of any length (but enough to do semantics, preferrably millions of tokens; stopwords, punctuation and such should be removed beforehand)
# The corpus object should be called 'perioddata'
# Use readLines & strsplit or such to load text data from file
# or use load("coha2000.RData") to load the pre-parsed coha, if provided
# or define here:
# perioddata = c("")   
#
minc = 100   # min freq threshold to include any word in the vocab -> term co-occurrence -> any models
# (for target words, the mindownsample parameter overrides this)
nounregex = "^S:" # regex to use to filter nouns (or other POS of interest) from the lexicon for the target sample, if words are tagged, e.g. "N:noun" -> "^N:"
# "." bypasses (matches everything), target sample will contain any POS.
nreps = 10  # how many times to replicate (resample+rerun) each word*downsample combo
targetmin = 500 # lowest frequency bound for the target words (should be large enough to allow at least some downsampling levels)
targetmax = Inf # highest frequency bound for the target words
nsample = 116   # number of target (log) frequency bands to use for sampling from the lexicon; adjust to get different number of words (some bands may be empty!)
subsetsizes = c(0,  seq(0.1,7,0.3) ) # downsample sizes in log change units (see slides); starts with 0-change for sanity test; this sequence gives 25 values, choose less values to save time
mindownsample = 10 # minimum sample size (frequency) for downsampling (i.e. will not downsample a word, given a subsetsizes[] value, if the resulting sample would have lower frequency than that)
targetwords = NULL # alternatively, specify a vector of target words (make sure tags match if words affixed with tag in corpus); if using random sampling from frequency bands, set this as NULL
winsize = 2   # +- context window size for tcm construction (linearly weighted by default)
ngloves = 10 # how many times to replicate GloVe (should be =< nreps), use 2-3 if in a hurry

# Parameters for a quick debug/example run:
# uncomment the following line and select-all+run/source entire script
#  nreps=2; subsetsizes=c(0,3);targetwords=c("smart");ngloves=1
# have a look at the sims_all & ranks_all objects afterwards 


#### Define functions ####

# PPMI weighting, using actual term counts and total
fullppmi = function(pmat, N, colp,rowp, positive=T){
  library(Matrix, quietly = T)
  pmat = Matrix::t(pmat)

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
    tmp = log2( (icol/N) / (rowp[not0] * colp[i] ) )
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
  pmat2 = Matrix::t(pmat2)
  pmat2 = drop0(pmat2,is.Csparse=T)
  return(pmat2)
}; fullppmi = cmpfun(fullppmi, options=list(optimize=3))

# # Another PPMI function for smaller bits
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
}; someppmi = cmpfun(someppmi, options=list(optimize=3))

# ppmi as matrix weighting
# makeppmi = function(pmat, positive=T){
#   library(Matrix, quietly = T)
#   pmat = Matrix::t(pmat)
#   
#   tcmrs = Matrix::rowSums(pmat)
#   tcmcs = Matrix::colSums(pmat)
#   N = sum(tcmrs)
#   colp = tcmcs/N
#   rowp = tcmrs/N
#   
#   pp = pmat@p+1
#   ip = pmat@i+1
#   tmpx = rep(0,length(pmat@x))
#   for(i in 1:(length(pmat@p)-1) ){ 
#     #for(i in 1:100 ){ 
#     #not0 = which(pmat[, i] > 0)
#     ind = pp[i]:(pp[i+1]-1)
#     not0 = ip[ind]
#     icol = pmat@x[ind]
#     #print(icol)
#     #tmp = log( (pmat[not0,i]/N) / (rowp[not0] * colp[i] ))
#     tmp = log2( (icol/N) / (rowp[not0] * colp[i] ))
#     tmpx[ind] = tmp
#     #print(tmp)
#     # tmp = ifelse(tmp < 0, 0, tmp)
#     #pmat2[not0, i] = tmp
#   }
#   if(positive){
#     tmpx[tmpx<0] = 0 
#   }
#   pmat2 = pmat
#   pmat2@x = tmpx
#   pmat2 = Matrix::t(pmat2) # words back to rows
#   pmat2 = drop0(pmat2,is.Csparse=T) 
#   return(pmat2)
# }; makeppmi = cmpfun(makeppmi, options=list(optimize=3))

# Functions for the APSyn similarity metric
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

# target word finder function, also calculates sizes of the downsamples for every word x downsample
sampletargets = function(voc,targetmin,targetmax,nounregex,nsample,subsetsizes,mindownsample,targetwords=NULL){
  if(is.null(targetwords)){
    vtmp = voc[voc$term_count>=targetmin & voc$term_count<=targetmax & grepl(nounregex, voc$term),] 
    vtmp1 = vtmp$term_count; names(vtmp1) = rownames(vtmp); rm(vtmp)
    vc = cut(log(vtmp1), nsample) 
    words=c()
    for(v in seq_along(levels(vc))){
      try({words = c(words, sample(names(vtmp1[vc%in%levels(vc)[v]]), 1))},
          silent=T) # will throw error if band is empty, ignore
    }
    print(paste(paste("Sampled", length(words), "target words:   " ),
                paste0(words, " (",vtmp1[words], ")",  collapse=", ") ),quote=F)
  } else {
    words=targetwords
  }
  words_prime = paste0(words, "_prime")
  actualns = voc[sapply(words, match, voc$term), "term_count"]
  randomns = vector("list", length(words))
  for(j in seq_along(words)){
    for(r in seq_along(subsetsizes)){
      tmp = round(exp( log(actualns[j]) - subsetsizes[r])) # calculates subsample size
      if(tmp >= mindownsample) { 
        randomns[[j]][r] = tmp
        names(randomns[[j]])[r] = paste0(words_prime[j], "_", subsetsizes[r])
      }
    }
  }
  return(list(words=words,words_prime=words_prime,actualns=actualns,randomns=randomns))
}

# simlex test function
dosimlex = function(simlexpath, tcm, ptcm_full,d1,allrelevants, nountag="S:"){
  simlex = read.table(file=simlexpath, header = 1, sep = "\t", stringsAsFactors = F)
  # make tags match (edit to work for more tags) 
  simlex$word1 = ifelse(simlex$POS=="N", paste0(nountag, simlex$word1), simlex$word1)
  simlex$word2 = ifelse(simlex$POS=="N", paste0(nountag, simlex$word2), simlex$word2)
  simlex = simlex[which(simlex[,1] %in% rownames(tcm) & simlex[,2] %in% rownames(tcm)   ),] # can only test words that are in the corpus lexicon
  # nrow(simlex) # 974 # 950
  
  # test simlex for the 5 methods
  # get the similarities:
  tmp = sim2(tcm[simlex[,1],,drop=F], tcm[simlex[,2],,drop=F]  )
  simlex$plaincos = Matrix::diag(tmp)
  tmp = sim2(ptcm_full[simlex[,1],,drop=F], ptcm_full[simlex[,2],,drop=F]  )
  simlex$ppmicos = Matrix::diag(tmp)
  tmp = sim2(d1[simlex[,1],,drop=F], d1[simlex[,2],,drop=F]  ) # lsa
  simlex$lsacos = Matrix::diag(tmp)
  tmp=c()
  for(i in 1:nrow(simlex)){
    tmp[i] = dosim(names(allrelevants[[simlex[i,1]]]), names(allrelevants[[simlex[i,2]]]))
  }
  simlex$apsyn = tmp
  
  # train full glove model, takes some time
  glove_model = GloVe$new(word_vectors_size = 100, vocabulary = rownames(tcm), x_max = 10)
  invisible(capture.output( # this disables the glove iteration messages
    word_vectors_main <- glove_model$fit_transform(tcm, n_iter = 20L,
                          n_check_convergence=5L, convergence_tol= 0.006) # takes time
  ))
  word_vectors_main = glove_model$fit_transform(tcm, n_iter = 50L)
  word_vectors_context = glove_model$components
  word_vectors = word_vectors_main + t(word_vectors_context) # the recommended addition step
  tmp = sim2(word_vectors[simlex[,1],,drop=F], word_vectors[simlex[,2],,drop=F]  )
  simlex$glovecos = Matrix::diag(tmp)
  
  # do the correlations:
  corrs = c(
    plaincos = cor(simlex$plaincos, simlex$SimLex999, method = "spearm"), # 0.1039354 # 0.12 window=2
    ppmicos = cor(simlex$ppmicos, simlex$SimLex999, method = "spearm"), #  -0.01462576 # 0.0984 # 0.218 nouns only | 0.16879 newppmi,window=2
    apsyn = cor(simlex$apsyn, simlex$SimLex999, method = "spearm"),   # 0.3372357
    lsacos = cor(simlex$lsacos, simlex$SimLex999, method = "spearm") #,  # 0.38 (300dim), 0.33 (100dim) #  0.46 nouns | 0.341 100dim+newppmi | 0.379 300dim+newppmi | 0.23 oldppmi+300! | 0.37 window=2,counts
    #glovecos = cor(simlex$glovecos, simlex$SimLex999, method = "spearm") # 0.2481524 for 20x, 0.2513886 with 50x iterations
  )
  return(corrs)
}


#### Choose words, calculate downsample sizes ####
#
# iterate over corpus, create vocabulary object
maxn = length(perioddata) # COHA last decade, content words only (excluding proper nouns) is 11 023 006
# requires text2vec
it = itoken(list(perioddata), progressbar = FALSE)
voc = prune_vocabulary(create_vocabulary(it), term_count_min = minc); rownames(voc) = voc$term # includes only >=minc words
#
# Sample target words from the lexicon, but from equally spaced log freq bands:
# (some high-freq bands are likely empty)
# (alternatively, may choose words manually, set the targetwords parameter)
targetdata = sampletargets(voc,targetmin,targetmax,nounregex,nsample,subsetsizes,mindownsample, targetwords)
words = targetdata$words
words_prime =targetdata$words_prime
actualns = targetdata$actualns
randomns = targetdata$randomns


#### Prepare reference corpus ####
#
# This is the reference term-term co-occurrence matrix:
vectorizer = vocab_vectorizer(voc) # needed later too!
# tcm from windows of size 5, linearly weighted counts within window; should maybe experiment with other window sizes (e.g. Levy's recommendation of 2 - but maybe too small for low freq)
tcm = create_tcm(it, vectorizer, skip_grams_window = winsize); tcm = tcm + Matrix::t(triu(tcm)) # makes symmetric
Matrix::diag(tcm) = 0 # removes self-co-occurrences


#### Prepares LSA, PPMI matrices ####
#
# LSA full:
lsa_1 = LatentSemanticAnalysis$new(300) # used 300 dimensions
# NB: this uses count matrix as input to LSA; apparently usually PPMI matrix is used; but seems to makes no difference with Simlex though
d1 = lsa_1$fit_transform(tcm) # embeds (symmetric tcm, so no difference which way)
# should incorporate Levy's "don't use SVD "correctly""?, i.e. without eigenvector weighting (p=0.5).


# PPMI
N = sum(tcm) # denominator for P(w,c)
Nwords = sum(voc$term_count)
normvoc = voc$term_count/Nwords; names(normvoc) = voc$term
ptcm_full = fullppmi(pmat=tcm, N=N, colp=normvoc,rowp=normvoc) # for ppmi and apsyn
ptcm = ptcm_full[words,,drop=F]
relevants = dorelevants(ptcm, 100); names(relevants) = words
maxapsyn = dosim(1:100, 1:100) # value to divide by to get [0,1] similarity value in apsyn
allrelevants = dorelevants(ptcm_full, 100); names(allrelevants) = rownames(ptcm_full) # for apsyn
# should maybe incorporate context distribution smoothing ("raise unigram distribution to the power of alpha=0.75 for context words"), apparently would improve


### optional: simlex test ####
# see how the reference models perform with Simlex999 (takes some time, need to train glove)
# need to specify path to where the SimLex-999.txt file with the similarity scores is:
# dosimlex(simlexpath = "SimLex-999.txt",tcm, ptcm_full,d1,allrelevants, nountag="S:")




#### Prepare structures to hold the similarity values from downsampling ####
#
allwords = voc$term  # full vocab
sims_all = vector("list", 5)
ranks_all = vector("list", 5)
for(l in 1:5){
  sims_all[[l]] = matrix(NA, ncol=length(subsetsizes), nrow=length(words), dimnames = list(words,subsetsizes))
  ranks_all[[l]] = matrix(NA, ncol=length(subsetsizes), nrow=length(words), dimnames = list(words,subsetsizes))
}
names(sims_all) = c("plaincos_all" ,"ppmicos_all" ,"ppmiapsyn_all", "lsacos_all" ,"glovecos_all" )
names(ranks_all) = names(sims_all)
emptymat = as(matrix(nrow=0,ncol=ncol(tcm)), "dgCMatrix") #  template




##### Run models #####
#
# (one simple for-loop)
gc(); startt = Sys.time(); print("")
for(i in seq_along(words)  ){  # for each word_prime
  cat(paste(" ", i, words[i]) )
  
  # find corpus positions for target word i;
  # since full tcm is created, just fetch +/5 snippets with word_prime in the middle
  # don't need to replace, can just create temp tcm and rename row, and add
  actualmatches = .Internal(which(perioddata %fin% words[i]))
  # temp matrices for replications:
  repsims = replicate(5,matrix(NA, ncol=length(subsetsizes), nrow=nreps), simplify = F)
  names(repsims) = c("plaincos" ,"ppmicos" ,"ppmiapsyn", "lsacos","glovecos") 
  repranks = replicate(5,matrix(NA, ncol=length(subsetsizes), nrow=nreps), simplify = F)
  names(repranks) = c("plaincos" ,"ppmicos" ,"ppmiapsyn", "lsacos","glovecos") 
  
  for(reps in 1:nreps){   # do 10 replications each
    cat(".")
    tcmrows = emptymat # initialize tcm for word_prime_downsample
    # this will contain a row per each word_prime * downsample size (ie 25 rows)
    for(r in seq_along(randomns[[i]])){  # for each downsample size of word i
      # (not all words have the same n of subsample values)
      try({
        ## Downsampling ##
        #
        # Sample replacement positions randomly
        # the number of times indicated in randomns[[i]][r] 
        # from among the pool of possible target word occurrences (ie from 1...actualns[i])
        # this currently works only if the input corpus is a single vector 
        # (not list, ie no doc/sentence borders)
        samplelocations = sample.int(n=actualns[i], size=randomns[[i]][r] )
        samp = actualmatches[ samplelocations ] 
        # gets the context window indices; but makes sure indices are not out of corpus bounds:
        sampfrom = pmax(samp-winsize, 1); sampto = pmin(samp+winsize, maxn) 
        samplist = vector("list", randomns[[i]][r] )
        for(s in seq_along(samplist)){
          samplist[[s]] = sampfrom[s] : sampto[s]
        }
        samplist = sort(unique(unlist(samplist))) # indices
        # indices are concatenated but unique'd and sorted since close-by occurrences can create doubling,
        # when their windows overlap; removing overlaps messes up the co-occurrence
        # counts of other words in the temporary tcm, but it does not matter
        # since only the target _prime row is relevant anyway
        # tested: the tcm vector of word_prime_0 from the nothing-removed sanity sample comes out 100% like the original
        
        ## Create the temp tcm from the downsamples ##
        #
        it_tmp = itoken(list(perioddata[samplist]), progressbar = FALSE)
        # "create dtm using same vectorizer, but new iterator"
        #  won't include low freq words (but based on global freqs, not this sample!)
        tcm_tmp = create_tcm(it_tmp, vectorizer, skip_grams_window = winsize) # uses vectorizer from above
        # need to make matrix symmetric/full for this to work:
        tcm_tmp = tcm_tmp + Matrix::t(triu(tcm_tmp))
        Matrix::diag(tcm_tmp) = 0 # make sure no self-co-occurrences
        tcm_tmp = tcm_tmp[words[i], ,drop=F] # extract relevant row
        rownames(tcm_tmp) = names(randomns[[i]][r]) # just one now
        tcmrows = rbind(tcmrows, tcm_tmp)    # append to temporary tcm
      }) # end try
    } # end downsample r
    
    # test integrity 
    #all(tcmrows[1,] == tcm[words[i],,drop=T])
    #which(tcmrows[1,] != tcm[words[i],,drop=T])
    
    
    ## Calculate similarity between reference word and downsampled word_prime ##
    #
    try({
      #
      ## Plain counts + cosine
      #
      repsims$plaincos[reps,1:length(randomns[[i]]) ] = sim2(tcm[words[i],,drop=F], tcmrows)[1,] # cosine similarity function
      simstmp = sim2(tcm, tcmrows) # rows=allwords, cols=subsamplewords; takes a moment
      # will also find how does word_prime rank in terms of similarity to reference word, 
      # compared to all other words; the vector spaces have very different densities in terms of
      # closest neighbors (mean closest may be 0.99 in one and 0.3 in another), making them incomparable.
      repranks$plaincos[reps,1:length(randomns[[i]]) ] = 
        sapply(seq_along(randomns[[i]]), function(x) 
          fmatch(words[i], allwords[rev(order2(simstmp[,x])) ] ) 
          # finds sim rank of word_i to word_i_prime_x
        ) # i.e. "1" would mean reference word is the closest synonym/neigbour to word_prime
      
      ## PPMI + cosine
      #
      ppmitmp = someppmi(NN=N, tcm1=tcmrows,normvoc_all=normvoc, normvoc_sub=randomns[[i]]/Nwords)
      repsims$ppmicos[reps, 1:length(randomns[[i]])] = sim2(ptcm[words[i],,drop=F], ppmitmp)[1,]
      simstmp = sim2(ptcm_full, ppmitmp) # rows=allwords, cols=subsamplewords; takes a moment
      repranks$ppmicos[reps,1:length(randomns[[i]]) ] = 
        sapply(seq_along(randomns[[i]]), function(x) {
          #print(allwords[rev(order2(simstmp[,x])) ][1]) # debug
          fmatch(words[i], allwords[rev(order2(simstmp[,x])) ] ) 
          # finds sim rank of word_i to word_i_prime_x
        }
        )
      
      ## PPMI + apsyn
      #
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
      
      ## LSA + cosine
      #
      lsatmp = lsa_1
      x = transform(tcmrows, lsatmp)
      repsims$lsacos[reps, 1:length(randomns[[i]]) ] = sim2(d1[words[i],,drop=F], x)[1,]
      
      simstmp = sim2(d1, x) # rows=allwords, cols=subsamplewords; takes a moment
      repranks$lsacos[reps,1:length(randomns[[i]]) ] = 
        sapply(seq_along(randomns[[i]]), function(x) 
          fmatch(words[i], allwords[rev(order2(simstmp[,x])) ] ) 
          # finds sim rank of word_i to word_i_prime_x
        )
      
      ## GloVe+cosine
      #
      if(reps %in% 1:ngloves){ # controls how many times to replicate GloVe (bottleneck)
        emptycols = as(matrix(0, nrow=nrow(tcmrows),ncol=nrow(tcmrows)), "dgCMatrix")
        colnames(emptycols) = rownames(tcmrows)
        gtcm = cbind(rbind(tcm,tcmrows), rbind(Matrix::t(tcmrows),emptycols ) )
        # this seems the right concatenation, but the sanity sample does not come
        # out with 100% identity, cosine is slightly (-0.0001) below 1. maybe a glove problem?
        # also apparently high freq words with 6-7 log reduction get negative cosines
        # (ok in principle, but maybe model not trailing long enough?)
        glove_model = GloVe$new(word_vectors_size = 100,
                                vocabulary = c(voc$term,rownames(tcmrows)), 
                                x_max = 10)
        invisible(capture.output( # this disables the glove iteration messages
          word_vectors_main <- glove_model$fit_transform(gtcm, n_iter = 20L, 
                                  n_check_convergence=5L, convergence_tol= 0.006) # takes time
          ))
        word_vectors_context = glove_model$components
        word_vectors = word_vectors_main + t(word_vectors_context) # the recommended addition step
        # word_prime sims:
        repsims$glovecos[reps, 1:length(randomns[[i]]) ] = 
          sim2(word_vectors[words[i],,drop=F], word_vectors[rownames(tcmrows),,drop=F])[1,]
        # word_prime ranks:
        simstmp = sim2(word_vectors[1:(nrow(word_vectors)-length(randomns[[i]])),,drop=F  ], # all w
           word_vectors[(nrow(word_vectors)-length(randomns[[i]])+1):nrow(word_vectors) ,,drop=F  ] ) # vs prime w's
        repranks$glovecos[reps,1:length(randomns[[i]]) ] = 
          sapply(seq_along(randomns[[i]]), function(x) {
              fmatch(words[i], allwords[rev(order2(simstmp[,x])) ] ) 
              # finds sim rank of word_i to word_i_prime_x
            }
          )
      }  # end ngloves
        
      # word2vec+cosine
          # probably not fasttext, which uses subword information (_prime suffix a problem?)
          # if fasttext + r wrapper:
          # expects single vector of words separated by spaces (so unlisted is fine)
          # also need a re-do for EACH subsample, cannot use tcm, need to replace in
          # actual corpus every time, and then rerun fasttext, meaning 100*25 models (*10 replicates)
      
    }) # end try()
  } # end replication loop
  
  # record mean values for the word_prime across all its downsamples
  # similarities to reference word
  sims_all$plaincos_all[i,] = colMeans(repsims$plaincos, na.rm = T)
  sims_all$ppmicos_all[i,] = colMeans(repsims$ppmicos, na.rm = T)
  sims_all$ppmiapsyn_all[i,] = colMeans(repsims$ppmiapsyn, na.rm = T)
  sims_all$lsacos_all[i,] = colMeans(repsims$lsacos, na.rm = T)
  sims_all$glovecos_all[i,] = colMeans(repsims$glovecos, na.rm = T)
  
  # ranks to reference word, compared to all other words
  ranks_all$plaincos_all[i,] = colMeans(repranks$plaincos, na.rm = T)
  ranks_all$ppmicos_all[i,] = colMeans(repranks$ppmicos, na.rm = T)
  ranks_all$ppmiapsyn_all[i,] = colMeans(repranks$ppmiapsyn, na.rm = T)
  ranks_all$lsacos_all[i,] = colMeans(repranks$lsacos, na.rm = T)
  ranks_all$glovecos_all[i,] = colMeans(repranks$glovecos, na.rm = T)
  
} # end words loop

print(paste("Started at", startt))
print(paste("Finished at", Sys.time()))

# Output:
# sims_all: list, each element is a matrix for each of the similarity methods, containing the similarity scores, for all words across all the downsamples (less frequent words are of course NA for higher downsamples, depend on the mindownsample value)
# ranks_all: list, same, but ranks
# in both lists: each of the 5 elements is a matrix (one for each similarity method), each matrix contains the similarity/rank values of the downsampled words (rows), with the log difference/downsampling values (from parameter subsetsizes) as columns.


