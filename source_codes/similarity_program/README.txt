Generates Similarity Files for each type of Nodes in IIN_map.

Commands:
1. make
2. execute
	./alg <IIN graph> <NodeNum> <output directory>

	EX: ./alg ../../data/iin.txt 1 ../../data/sim_data/

Parameters:
	(1) IIN input graph file

	(2) Node Num
		The num of Node you want to calculate the similarity

	(3) output directory
		The output Similarity file for each type of Node will be put into this directory.

The file name follows the format   =>   SIM_<IIN file name>_t<node Type>.txt
	EX: CET_iin.txt_n15396.txt