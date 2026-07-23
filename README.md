# merge_lanes branch

Custom branch to merge across lanes, which isn't needed that often. 

In the example where a sample was sequenced, and then the same sample re-sequenced to a higher depth, this branch can be used to run the full processing pipeline whilst also merging between the sequencing runs. This is done after running reachtools as sometimes R2 length varies between sequencing lanes. 

In our example we actually need to do 2 things:
- merge across lanes
- merge between lanes on the same run

We have the naming convention {SAMPLE}_SP{1-2}_S{1-8}_S1{1-6}:
- SP{1-2} represents separate plates for splitpool barcoding
- The first S (S1-S8) represents sublibraries - distinct cell populations
- The second S (S1-6) represents the same cells sequenced over different 

The aim is to collapse S1-6 together, as they are the same cells, so we are left with SP{1-2}_S{1-8}. 

Note: if we want to merge across lanes, we just need to include the other folders/files into the pipeline input, and it will collapse all SP1-S1, SP1-S2 etc, across different sequencing lanes together. 

This only works if files are slightly differently named, for example the file name may include sequencing ID information, but in the case where the files are identically named, but in different folders (e.g. Lane1/Sample1_SP1_S1_S1.fq.gz, Lane2/Sample1_SP1_S1_S2.fq.gz), it will not work by default, as parts of the pipeline will just see "Sample1_SP1_S1_S1.fq.gz". **To assist here I have added the add_lane_to_filename.sh script**. To run: simply go to the dir you want to drop the files, run "bash add_lane_to_filename.sh <PATH/TO/FILES/DIR>" and it will create a folder called input_data and softlink to your files (the new files will be input_data/L1_Sample1_SP1_S1_S1.fq.gz, input_data/L2_Sample1_SP1_S1_S1.fq.gz).
