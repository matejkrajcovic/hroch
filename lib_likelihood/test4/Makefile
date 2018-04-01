BIN = /projects/dups-cornell/dups/dups/mcmc/scripts/
BIN2 = /projects/dups-cornell/dups/dups/mcmc/bin/
BIN3 = /projects/dups-atoms/scripts/

ALL:
	echo "Hello"

.SECONDARY:
# prevents marker files (and all targets) from being deleted

.DELETE_ON_ERROR:
#delete targets if a rule fails

#generate fasta, atom and history files for each sample
#according to a config file
%.done-sample : %.config
	$(BIN)/generate.pl $<
	touch $@

#for each history we can generate gene trees
%.trees : %.history %.atoms *.tree
	$(BIN)/compute_gene_trees.pl $^ > $@

#generate gene trees for all histories in one config file
%.done-trees : %.done-sample
	ls $**.history > $@.tmp
	$(BIN2)/run-make.pl '*.history' '*.trees' < $@.tmp
	rm $@.tmp
	touch $@

#filter out all atoms shorter than 500bp
f-%.atoms f-%.fasta : test-%.atoms test-%.fasta
	$(BIN)/filter_short_atoms.pl test-$* 500 f-$*
	cp test-$*.fasta f-$*.fasta

#filter out short and compteley deleted atoms from trees, ratoms and histories
f-%.ratoms f-%.trees f-%.history : test-%.ratoms test-%.history test-%.history f-%.atoms
	$(BIN)/filter_short_atoms.pl -o test-$* 500 f-$*

# filtering of short atoms for all samples in one config file
f.done-sample : test.done-sample
	ls test-*.fasta > $@.tmp
	$(BIN2)/run-make.pl 'test-*' 'f-*' < $@.tmp
	$(BIN2)/run-make.pl 'test-*.fasta' 'f-*.ratoms' < $@.tmp
	rm $@.tmp
	touch $@

#for each atom and fasta file, we can sample gene trees by mrbayes
%.sampled-trees : %.fasta %.atoms 
	$(BIN3)/sample_trees.pl $^ > $@

#generate sampled gene trees for all samples in one config file
%.done-sampled-trees : %.done-sample
	ls $**.atoms > $@.tmp
	$(BIN2)/run-make.pl '*.atoms' '*.sampled-trees' < $@.tmp
	rm $@.tmp
	touch $@

# hisotry in a new format
%.nhistory : %.atoms %.history test.tree
	../reformat.pl $^ $@.map > $@
	perl -lane 'die $$_ unless $$F[1] eq $$F[0]' $@.map

#new history for all samples
%.done-nhistory : %.done-sample
	ls $**.history > $@.tmp
	$(BIN2)/run-make.pl '*.history' '*.nhistory' < $@.tmp
	rm $@.tmp
	touch $@

#prune (collapse) short branches from mr bayes as well as from real trees
%.sampled-trees-pruned : %.sampled-trees
	$(BIN)/collapse_short_branches.pl $*.atoms 5 < $< > $@ 2>$@.log

%.trees-pruned : %.trees
	$(BIN)/collapse_short_branches.pl $*.atoms 5 < $< > $@ 2>$@.log

#pruning for all samples
%.done-sampled-trees-pruned : %.done-sampled-trees
	ls $**.sampled-trees > $@.tmp
	$(BIN2)/run-make.pl '*' '*-pruned' < $@.tmp
	rm $@.tmp
	touch $@

#see how good are sampled trees
%.sampled-trees-pruned-quality : %.sampled-trees-pruned %.trees
	$(BIN)/unify_tree_topologies.pl $^ >$@

%.sampled-trees-quality : %.sampled-trees %.trees
	$(BIN)/unify_tree_topologies.pl $^ >$@

#quality for all samples
%.done-sampled-trees-quality : %.done-sampled-trees %.done-trees
	ls $**.sampled-trees > $@.tmp
	$(BIN2)/run-make.pl '*' '*-quality' < $@.tmp
	rm $@.tmp
	touch $@


#get summary of event number per edge for simulated data
%.event-summary : %.history
	$(BIN)/summarize_history.pl -1 param/hky-primate.tree $< > $@

%.done-event-summary : %.done-sample
	ls $**.history > $@.tmp
	$(BIN2)/run-make.pl '*.history' '*.event-summary' < $@.tmp
	rm $@.tmp
	touch $@

###################################################
# Stats
results/%.event-summary : %.history results/%.*.samples
	$(BIN)/summarize_history.pl 2500 param/hky-primate.tree results/$*.*.samples > $@

results/%.done-event-summary : %.done-sample results/%*.samples
	ls $**.history > $@.tmp
	$(BIN2)/run-make.pl '*.history' 'results/*.event-summary' < $@.tmp
	rm $@.tmp
	touch $@

