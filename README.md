General Idea:

I have taken the full hitting data (large time series matrices of markered motion capture data) from Driveline's OpenBiomechanics Dataset, and am performing dimensionality reduction on it to get to a more easily modellable, understandable place

A lot of the ideas come from Liu et al, 2006 ("Human motion estimation from a reduced marker set") and Barbic et al, 2004 ("Segmenting Motion Capture Data into Distinct Behaviors")

The core idea of the work is to first perform principal feature analysis (via clustering), then probabilistic principal component analysis, to distinguish signal and noise, and build covariance matrices of the principal features/markers. Then, learn a section of the swing data, and check the rows immediately following the learning window, to see if they are similar in cov structure or different beyond a certain threshold. If different, create a segment break at that point. Following this, runs a kmeans clustering algorithm to link together similar behaviours within swings throughout the dataset. This informs me of who swings similarly to who else, who has a very repeatable swing, and most importantly, which covariance structures of body parts create the best hitters outcomes (highest exit velos and bat speeds). Also, I can then reconstruct the full swing using only the feature vectors of principal markers, which worked incredibly well

The majority of the work is done and explained in the DimensionalityReductionBiomechData.R file, allowing a user to take one individual swing through the entire pipeline to inspect every element. There are separate files to run the pipeline on multiple swings, or the entire dataset altogether
There is a 2d plotting script which allows a user to quickly sanity check the motions happening at a particular part of the swing
There is also a c3d to csv conversion script (in python) since the original Driveline provided files were in c3d format

If you'd like to get in touch with me or ask about my work, you can do so here:

ryangunther1 on bluesky

ryan_gunther1 on X

https://www.linkedin.com/in/ryan-gunther-68a156195/

ryangunther98@gmail.com
