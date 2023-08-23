Generates the answer to the input queries.

1. make
2. execute:

	./alg <IIN Graph> <Influence Graph> <meta-path/structure file name> <input user> <input item file> <k> <delta> <output top K user File> <output result File> <number of threads> <pruning option> <SMK files directory> <CET files directory> <distance> <flag> <spatial distance>

Parameters:
	Ex: ./alg ../../data/iin.txt ../../data/iin_influence.txt ../../data/meta_path_file/meta_path.txt 15396 ../../data/input_item_file/7204408.txt 3 0.0 ../../data/output_file/15396_7204408_topK.txt ../../data/output_file/15396_7204408_result.txt 10 0 ../../data/smk_data/ ../../data/cet_data/ 3 0 5

	(1) IIN input graph File

	(2) IIN Influence file

	(3) meta-paragraph/structure file name

	(4) input user id: An integer

	(5) input items lists: This File lists the input items, each line is an item ID

	(6) K: An integer, there will be at most k node in output result

	(7) delta: threshold of omega (float number of 0.0~1.0). 
        Ex: if delta = 0.2, the omega value less than 0.2 will be taken as 0.

	(8) output topK user file contains output the top-k node with highest SIMER

	(9) result File: List all the possible source node with their distance to the input user, and SIMER value
					
	(10) number of threads, to execute the program using multiple thread

	(11) prune switch: 0 doesn't prune, 1 prune

	(12) SMK file directory: the folder which contains SMK files

	(13) CET file directory: the folder containing the CET files of users with respected to specific meta_structures and    all users

	(14) The distance from the input user node. If the distance is larger than this value, it will not compute SIMER.

	(15) flag switch: if larger than 0, the program only consider the instance that end with the input item
				  if 0, consider all the instances.
	(16) Spatail distance: If the distance between CET's leaf node and CET's root larger than this value, this CET's leaf node will not be added in.