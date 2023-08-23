Generates the answer to the input queries.

1. make
2. execute:

	./alg <IIN Graph> <Influence Graph> <meta-path/structure file name> <input user> <input item file> <source node file> <delta> <output result File> <number of threads> <pruning option> <SMK files directory> <CET files directory> <distance> <flag> <spatial distance> <SIM files directory>

Parameters:
	Ex: ./alg ../../data/iin.txt ../../data/iin_influence.txt ../../data/meta_path_file/meta_path.txt 15396 ../../data/input_item_file/7204408.txt ../../data/source_node.txt 0.0  ../../data/output_file/15396_7204408_result.txt 10 0 ../../data/smk_data/ ../../data/cet_data/ 3 0 5 ../../data/sim_data/

	(1) IIN input graph File

	(2) IIN Influence file

	(3) meta-paragraph/structure file name

	(4) input user id: An integer

	(5) input items lists: This File lists the input items, each line is an item ID

	(6) source node file: This File lists the source node, each line is an Source_Node ID

	(7) delta: threshold of omega (float number of 0.0~1.0). 
        Ex: if delta = 0.2, the omega value less than 0.2 will be taken as 0.

	(8) result File: List all the possible source node with their distance to the input user, and SIMER value
					
	(9) number of threads, to execute the program using multiple thread

	(10) prune switch: 0 doesn't prune, 1 prune

	(11) SMK file directory: the folder which contains SMK files

	(12) CET file directory: the folder containing the CET files of users with respected to specific meta_structures and    all users

	(13) The distance from the input user node. If the distance is larger than this value, it will not compute SIMER.

	(14) flag switch: if larger than 0, the program only consider the instance that end with the input item
				  if 0, consider all the instances.
	(15) Spatail distance: If the distance between CET's leaf node and CET's root larger than this value, this CET's leaf node will not be added in.

	(16) SIM file directory: the folder which contains SIM files