# beta_sheets
These are python and xml scripts developed for designing all-beta protein structures.

The Blueprint class eases the generation and manipulation of protein topologies, and is used for loop analysis and setting up design calculations.
The loop_analysis folder includes the PyRosetta script to analyze beta-arch loop connections of a non-redundant database of natural protein structures ("pisces_pdblist_seqid30_Res2.0.txt")
The design folder includes the python script to automatically generate combinations of all-beta protein topologies. This automatically generates blueprint files to be used with the Blueprint builder mover specified in the accompanying xml file.
