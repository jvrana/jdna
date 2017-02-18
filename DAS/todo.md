# TODO
1. Use different parser, store all files as json files for communication
1. (ALMOST) Update features on plasmid map viewer
2. Update j5 Interface (correct the XMLRPCLIB errors
4. Update BLAST so that temp databases are used instead
5. Merge j5 parameter file with cost function file
5. Multithreading j5 submission
5. Add suggested sequencing primers
5. Deploy Django user interface
6. Deploy JBEI integration
7. Deploy JBEI on EC2 instance (or on Biofab AWS)
8. !!!!!!! Fix subject start and end for products
8. !!!!!!! Implement RC query searches as well if query is on opposite strand
9. Add true circular subject search
9. Fix remove redundant contigs ('except contig_type')
10. Primer ends fix
11. Separate out Fragments (used directly in a reaction) from plasmids
(that must be pcr amplified)
12. Find sequencing primers for junctions, account for that in the cost
13. Use node.js for live update on assembly
14. Heirarchical assemblies (primer layers, prc'ing off of gibsons)
15. Fix disagreements between query start in blast and start index 0
16. Speed improvments for high complexity (high Query, low length) assemblies
17. Implement costs and part quality for backtracing
18. Implement 'robustness' or assemblies that are robust against failings