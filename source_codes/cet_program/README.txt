Generates CET Files for each meta-path/structure from each node matching the source of the meta-path/structure or for each node in the node list.

Commands:
1. make
2. execute
	./alg <IIN graph> <meta-path/sturcture list file> <output directory> <node location file> <node list file>
	[The last argument is optional]
	EX: ./alg ../../data/iin.txt ../../data/meta_path_list_file/meta_path_list.txt ../../data/cet_data/ ../../data/node_location_file/node_location.txt
	EX: ./alg ../../data/iin.txt ../../data/meta_path_list_file/meta_path_list.txt ../../data/cet_data/ ../../data/node_location_file/node_location.txt ../../data/example_of_node_list.txt

Parameters:
	(1) IIN input graph file

	(2) meta-path/structure list file 
		In this file, each line lists a meta-path/structure file name.

	(3) output directory
		The output CET file for each meta-path/structure and the relevant nodes will be put into this directory.

	(4) node location file
The file name follows the format   =>   CET_<meta-path/structure file name>_n<node ID>.txt
	EX: CET_metapath.txt_n15396.txt
