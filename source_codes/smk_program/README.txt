Generates the SMK file of each individual node (node of type “user”) in the input IIN graph.

Commands:
1. make
3. execute

	./alg <IIN Graph> <Influence Graph> <distance> <output directory>

	Ex: ./compute_SMK ../../data/iin.txt ../../data/iin_influence.txt 3 ../../data/smk_data/

Parameters:
	(1) IIN input graph File

	(2) IIN influence File
                In this file, each line represents the directed influence between a pair of individual in the IIN that share a social relationship. Each line has the format <user_id1>,<user_id2>,<influence prob.> where <influence prob.> is the probability of influence from <user_id1> to <user_id2>

	(3) distance 
                If the distance to source_node is larger than this value, the node will not be computed SMK

	(4) Output directory, 
                The output SMK file of each individual node will be put in this directory
The output file name follows the format => SMK_n<node_id>.txt
        EX: SMK_n15396.txt
