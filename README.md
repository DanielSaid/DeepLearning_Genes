# DeepLearning_Genes
MaCSBio project on predicting gene expressions across different domains

Making predictions:
Existing sets (i.e. random, nested or orthologs):
-	Call file PredictEncoding.py from job script
-	Scroll to bottom of PredictEncoding (“main”)
-	Change output domain as desired
-	Change parameters of loops as desired
-	Change method as desired: untrained_ae, autoencoder (-> pre-trained AE), rf, knn, elastic or cnn
-	In order to run things in parallel, create separate job scripts AND separate copies of PredictEncoding (at least I never tried editing a running file)

Gene list provided (i.e. NAFLD, Steatosis etc.):
-	Call file PredictEncoding_preselected.py from job script
-	Scroll to bottom of PredictEncoding_preselected (“main”)
-	Change gene list as desired
-	Change domain and method as desired

Creating new gene sets:
-	Call file RandomGenes.py from job script
-	Change output domain as desired

Other remarks:
-	In order to investigate a data file, use OpenData.sh and OpenData.py
-	data files (e.g. data_X20_1.p) from RandomGenes.py will have format
[og_X, data_compounds, gene_list_x, gene_variance],
where og_X are the actual gene expression values
-	compounds used can be found in file CompoundLists.py (default is UNGENERAL_45)
-	Keras needs to be installed inside cluster; follow: https://doc.itc.rwth-aachen.de/display/CC/tensorflow (and replace “tensorflow” by “keras”)
