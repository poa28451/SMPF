 package ca.pfv.spmf.algorithms.frequentpatterns.fpgrowth;
 
 /* This file is copyright (c) 2008-2015 Philippe Fournier-Viger
 * 
 * This file is part of the SPMF DATA MINING SOFTWARE
 * (http://www.philippe-fournier-viger.com/spmf).
 * 
 * SPMF is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * SPMF is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with
 * SPMF. If not, see <http://www.gnu.org/licenses/>.
 */


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ca.pfv.spmf.patterns.itemset_array_integers_with_count.Itemset;
import ca.pfv.spmf.patterns.itemset_array_integers_with_count.Itemsets;
import ca.pfv.spmf.tools.MemoryLogger;

/** 
 * This is an implementation of the FPGROWTH algorithm (Han et al., 2004).
 * FPGrowth is described here:
 * <br/><br/>
 * 
 * Han, J., Pei, J., & Yin, Y. (2000, May). Mining frequent patterns without candidate generation. In ACM SIGMOD Record (Vol. 29, No. 2, pp. 1-12). ACM
 * <br/><br/>
 * 
 * This is an optimized version that saves the result to a file
 * or keep it into memory if no output path is provided
 * by the user to the runAlgorithm method().
 *
 * @see FPTree
 * @see Itemset
 * @see Itemsets
 * @author Philippe Fournier-Viger
 */
public class AlgoFPGrowth {

	// for statistics
	private long startTimestamp; // start time of the latest execution
	private long endTime; // end time of the latest execution
	private int transactionCount = 0; // transaction count in the database
	private int itemsetCount; // number of freq. itemsets found
	
	// parameter
	public int minSupportRelative;// the relative minimum support
	
	BufferedWriter writer = null; // object to write the output file
	
	// The  patterns that are found 
	// (if the user want to keep them into memory)
	protected Itemsets patterns = null;
		
	// This variable is used to determine the size of buffers to store itemsets.
	// A value of 50 is enough because it allows up to 2^50 patterns!
	final int BUFFERS_SIZE = 2000;
	
	// buffer for storing the current itemset that is mined when performing mining
	// the idea is to always reuse the same buffer to reduce memory usage.
	private String[] itemsetBuffer = null;
	// another buffer for storing fpnodes in a single path of the tree
	private FPNode[] fpNodeTempBuffer = null;
	
	// This buffer is used to store an itemset that will be written to file
	// so that the algorithm can sort the itemset before it is output to file
	// (when the user choose to output result to file).
	private String[] itemsetOutputBuffer = null;
	
	/**
	 * Editted ver.
	 */
	//private final String SEPARATOR = " ";
	private final String SEPARATOR = ",";
	private double pruneConf;
	
	/**
	 * Constructor
	 */
	public AlgoFPGrowth() {
		
	}

	/**
	 * Method to run the FPGRowth algorithm.
	 * @param input the path to an input file containing a transaction database.
	 * @param output the output file path for saving the result (if null, the result 
	 *        will be returned by the method instead of being saved).
	 * @param minsupp the minimum support threshold.
	 * @return the result if no output file path is provided.
	 * @throws IOException exception if error reading or writing files
	 */
	public Itemsets runAlgorithm(String input, String output, double minsupp, double pruneConf) throws FileNotFoundException, IOException {
		// record start time
		startTimestamp = System.currentTimeMillis();
		// number of itemsets found
		itemsetCount = 0;
		
		//initialize tool to record memory usage
		MemoryLogger.getInstance().reset();
		MemoryLogger.getInstance().checkMemory();
		
		// if the user want to keep the result into memory
		if(output == null){
			writer = null;
			patterns =  new Itemsets("FREQUENT ITEMSETS");
	    }else{ // if the user want to save the result to a file
			patterns = null;
			writer = new BufferedWriter(new FileWriter(output)); 
			itemsetOutputBuffer = new String[BUFFERS_SIZE];
		}
		
		// (1) PREPROCESSING: Initial database scan to determine the frequency of each item
		// The frequency is stored in a map:
		//    key: item   value: support
		final Map<String, Integer> mapSupport = scanDatabaseToDetermineFrequencyOfSingleItems(input); 
		
		/*for(Map.Entry<String, Integer> entry: mapSupport.entrySet()){
			System.out.println(entry.getKey() + " " + entry.getValue());
		}*/
		
		// convert the minimum support as percentage to a
		// relative minimum support
		this.minSupportRelative = (int) Math.ceil(minsupp * transactionCount);
		
		/**
		 * Print support counts of every items
		 */
		System.out.println("132:Support Count: (minsup = " + minSupportRelative + " (" + minsupp + "%))\n" + mapSupport);
		
		// (2) Scan the database again to build the initial FP-Tree
		// Before inserting a transaction in the FPTree, we sort the items
		// by descending order of support.  We ignore items that
		// do not have the minimum support.
		FPTree tree = new FPTree();
		
		// read the file
		BufferedReader reader = new BufferedReader(new FileReader(input));
		String line;
		// for each line (transaction) until the end of the file
		while( ((line = reader.readLine())!= null)){ 
			// if the line is  a comment, is  empty or is a
			// kind of metadata
			if (line.isEmpty() == true ||	line.charAt(0) == '#' || line.charAt(0) == '%'
				|| line.charAt(0) == '@') {
				continue;
			}
			
			//line = line.trim();
			
			String[] lineSplited = line.split(SEPARATOR);
//			Set<Integer> alreadySeen = new HashSet<Integer>();
			List<String> transaction = new ArrayList<String>();
			
			// for each item in the transaction
			for(String itemString : lineSplited){
				//Integer item = Integer.parseInt(itemString);
				// only add items that have the minimum support
				if(mapSupport.get(itemString) >= minSupportRelative){
					transaction.add(itemString);	
				}
			}
			
			/**
			 * By this point, this one transaction has been converted into Unordered Frequent Items
			 * The next thing to do is to sort it, make it into Ordered Frequent Items
			 */
			
			// sort item in the transaction by descending order of support
			Collections.sort(transaction, new Comparator<String>(){
				public int compare(String item1, String item2){
					// compare the frequency
					int compare = mapSupport.get(item2) - mapSupport.get(item1);
					// if the same frequency, we check the lexical ordering!
					if(compare == 0){ 
						return item1.compareTo(item2);
					}
					// otherwise, just use the frequency
					return compare;
				}
			});
			// add the sorted transaction to the fptree.
			/**
			 * Add this Ordered Frequent Items transaction into the FP-tree
			 */
			tree.addTransaction(transaction);
		}
		// close the input file
		reader.close();		
		
		// We create the header table for the tree using the calculated support of single items
		/**
		 * i.e. this is the list of every items that Support > minsup, sorted in order from highest sup to lowest sup. 
		 */
		tree.createHeaderList(mapSupport);
		System.out.println("204:Ordered Items: (" + tree.headerList.size() + " items)\n " + tree.headerList + "\n");
		

		/**
		 * At this point, FP-tree has been fully constructed. What's left is to find Conditional Pattern Bases
		 * 	and further find Conditional FP-trees which are used to create Frequent patterns
		 */
		
		/**
		 * Prune the FP-tree
		 */
		if(tree.headerList.size() > 0) {
			this.pruneConf = pruneConf;
			fptreePruning(tree, mapSupport);
		}
		
		// (5) We start to mine the FP-Tree by calling the recursive method.
		// Initially, the prefix alpha is empty.
		// if at least an item is frequent
		if(tree.headerList.size() > 0) {
			// initialize the buffer for storing the current itemset
			itemsetBuffer = new String[BUFFERS_SIZE];
			// and another buffer
			fpNodeTempBuffer = new FPNode[BUFFERS_SIZE];
			// recursively generate frequent itemsets using the fp-tree
			// Note: we assume that the initial FP-Tree has more than one path
			// which should generally be the case.
			fpgrowth(tree, itemsetBuffer, 0, transactionCount, mapSupport);
		}
		
		// close the output file if the result was saved to a file
		if(writer != null){
			writer.close();
		}
		// record the execution end time
		endTime= System.currentTimeMillis();
		
		// check the memory usage
		MemoryLogger.getInstance().checkMemory();
		
		// return the result (if saved to memory)
		return patterns;
	}

	private void fptreePruning(FPTree tree, Map<String, Integer> mapSupport){
		traverseTree(tree.root, tree, mapSupport);
		
		/**
		 * Try to cut tree branches here
		 */
		cutBranches(tree.root, tree);
		System.out.println("-------------------------");
		testTraverse(tree.root);
	}
	
	private void cutBranches(FPNode curNode, FPTree tree){
		for(int i=0; i<curNode.childs.size(); i++){
			cutBranches(curNode.childs.get(i), tree);
		}
		
		//Cut a branch
		System.out.println(curNode.itemID + " " + curNode.maxConf);
		if(curNode.itemID == "-1"){//Ignore the node if it's the root
			return;
		}
		
		//If it's not the upper most node (node that's one level next to root) and its maxConf is lower than the specified Conf from user)
		if(curNode.maxConf != -1 && curNode.maxConf < this.pruneConf){
			curNode.parent.childs.remove(curNode);//Remove the node from the tree
			FPNode tempNode = tree.mapItemNodes.get(curNode.itemID);
			while(tempNode.nodeLink != null){
				if(tempNode.nodeLink == curNode){//Also remove the node from the linked list of nodes with the same item ID
					tempNode.nodeLink = curNode.nodeLink;
					break;
				}
				tempNode = tempNode.nodeLink;
			}
		}
	}
	
	private void testTraverse(FPNode curNode){
		for(int i=0; i<curNode.childs.size(); i++){
			testTraverse(curNode.childs.get(i));
		}
		System.out.println(curNode.itemID + " " + curNode.maxConf);
	}
	
	private void traverseTree(FPNode curNode, FPTree tree, Map<String, Integer> mapSupport){		
		//Traverse recursively to the left-most child node first.
		//This is Pre-order traversal (read the data in the parent first, then traverse from the left-most child to the right-most child)
		
		//Do whatever you want
		if(curNode.itemID != "-1"){
			double nodeConf = calculateMaxConf(curNode, tree, mapSupport);
			curNode.maxConf = nodeConf;
		}
		
		for(int i=0; i<curNode.childs.size(); i++){
			traverseTree(curNode.childs.get(i), tree, mapSupport);
		}
	}
	
	private double calculateMaxConf(FPNode curNode, FPTree tree, Map<String, Integer> mapSupport){
		if(curNode.parent.itemID == "-1"){
			//Cannot find a max confidence of a single item
			return -1;
		}

		// Store all the nodes in this path up to the root, make it into current subtree
		List<FPNode> subtree = new ArrayList<>();
		FPNode temp = curNode;		
		while (temp.itemID != "-1") {
			subtree.add(temp);
			temp = temp.parent;
		}
		
		//Min Conf = Sup(parents+currentNode)/Sup(x), where x has 1 item to give max sup
		//Max Conf = Sup(parents+currentNod)/Sup(x), where x has k-1 items to give min sup, k = number of items in (parents+currentNod)
		//We want to find max conf
		
		//Find the whole subtree's Sup
		int subtreeSupport = calculateSubtreeSupport(subtree, tree, mapSupport, true);
		
		//Find the minimum possible Sup of the subtree
		int subtreeMinSup = transactionCount;
		List<List<FPNode>> subtreeCombination = findK_1Permutation(subtree); //K-1 subset of this subtree will always have the minimum Sup
		for(List<FPNode> combi : subtreeCombination){
			int combiSup = calculateSubtreeSupport(combi, tree, mapSupport, false);
			subtreeMinSup = (subtreeMinSup < combiSup) ? subtreeMinSup : combiSup;
			
			if(subtreeMinSup == 1){//There's no need to calculate Sup anymore if it reaches 1, because that's the minimum value possible
				break;
			}
		}
		System.out.print(" // minsup = " + subtreeMinSup);
		System.out.println(" // Max conf = " + subtreeSupport + "/" + subtreeMinSup);
		
		return (double) subtreeSupport / subtreeMinSup;
	}
	
	private int calculateSubtreeSupport(List<FPNode> subtree, FPTree tree, Map<String, Integer> mapSupport, boolean isPrint){
		int subtreeSupport = 0;
		/*List<FPNode> subtree = new ArrayList<>();
		FPNode temp = curNode;
		
		// Store all the nodes in this subtree, make it into current subtree
		while(temp.itemID != "-1"){
			subtree.add(temp);
			temp = temp.parent;
		}*/

		// We are trying to count the support by comparing this subtree with a subtree from FP-tree
		//	by back-tracking from the lowest node upto the root
		
		if(subtree.size() == 1){
			subtreeSupport = mapSupport.get(subtree.get(0).itemID);
			//System.out.println(subtree.get(0).itemID + " sup=" + subtreeSupport);
			return subtreeSupport;
		}
		
		FPNode temp = tree.mapItemNodes.get(subtree.get(0).itemID);
		//System.out.println(temp.itemID + ": ");
		int matchedParents = 0;
		while(temp != null){// Do a loop until we have already considered every nodes with this itemID
			int i = 1;
			int nodeSupport = temp.counter; //Remember this FP-tree node's support
			FPNode parent = temp.parent;
			//System.out.print(parent.itemID + " ");
			
			while(parent.itemID != "-1" && i < subtree.size()){ //If it isn't a root and the whole of current subtree hasn't been found yet
				//System.out.println(parent.itemID + " + " + subtree.get(i).itemID + "(" + (parent.itemID == subtree.get(i).itemID) + ")");
				
				if(parent.itemID.equals(subtree.get(i).itemID)){ // If the node from FP-tree has the same item as the current node from this subtree
					// Count the matching node and move the pointer (i) to the next node in current subtree
					matchedParents++;
					i++;
				}
				
				if(matchedParents == subtree.size()-1){ // If the current subtree is completely matched with the branch from FP-tree
					// Aggregate this subtree's support
					subtreeSupport += nodeSupport;
					break;
				}
				parent = parent.parent; // Move the pointer to the next node in FP-tree's subtree
			}
			matchedParents = 0;
			temp = temp.nodeLink; // Move onto the next node in FP-tree with the same item
			//System.out.print(" " + (temp == null) + " / ");
		}
		//System.out.println();
		
		if(isPrint){
			for(int i=0; i<subtree.size(); i++){
				System.out.print(subtree.get(i).itemID + " > ");
			}
			System.out.print("sup=" + subtreeSupport);
		}
		return subtreeSupport;
	}
	
	private List<List<FPNode>> findK_1Permutation(List<FPNode> subtree){
		List<List<FPNode>> subtreeCombination = new ArrayList<List<FPNode>>();		
		List<FPNode> combi;
		
		for(int i=0; i<subtree.size(); i++){
			combi = new ArrayList<>();
			for(FPNode item : subtree){
				//Create a subset of items that excludes current item (i.e. k-1 itemset)
				if(item != subtree.get(i)){
					combi.add(item);
				}
			}
			subtreeCombination.add(combi);
		}
		//System.out.println(subtreeCombination.toString());
		
		return subtreeCombination;
	}
	
	/**
	 * Mine an FP-Tree having more than one path.
	 * @param tree  the FP-tree
	 * @param prefix  the current prefix, named "alpha"
	 * @param mapSupport the frequency of items in the FP-Tree
	 * @throws IOException  exception if error writing the output file
	 */
	private void fpgrowth(FPTree tree, String [] prefix, int prefixLength, int prefixSupport, Map<String, Integer> mapSupport) throws IOException {
//		======= DEBUG ========
//		System.out.print("###### Prefix: ");
//		for(int k=0; k< prefixLength; k++) {
//			System.out.print(prefix[k] + "  ");
//		}
//		System.out.println("\n");
//				========== END DEBUG =======
//		System.out.println(tree);
		
		// We will check if the FPtree contains a single path
		boolean singlePath = true;
		// We will use a variable to keep the support of the single path if there is one
		int singlePathSupport = 0;
		// This variable is used to count the number of items in the single path
		// if there is one
		int position = 0;
		// if the root has more than one child, than it is not a single path
		if(tree.root.childs.size() > 1) {
			singlePath = false;
		}else {
			
			// Otherwise,
			// if the root has exactly one child, we need to recursively check childs
			// of the child to see if they also have one child
			FPNode currentNode = tree.root.childs.get(0);
			while(true){
				// if the current child has more than one child, it is not a single path!
				if(currentNode.childs.size() > 1) {
					singlePath = false;
					break;
				}
				// otherwise, we copy the current item in the buffer and move to the child
				// the buffer will be used to store all items in the path
				fpNodeTempBuffer[position] = currentNode;
				
				position++;
				// if this node has no child, that means that this is the end of this path
				// and it is a single path, so we break
				if(currentNode.childs.size() == 0) {
					break;
				}
				currentNode = currentNode.childs.get(0);
			}
		}
		
		// Case 1: the FPtree contains a single path
		if(singlePath && singlePathSupport >= minSupportRelative){	
			// We save the path, because it is a maximal itemset
			saveAllCombinationsOfPrefixPath(fpNodeTempBuffer, position, prefix, prefixLength);
		}else {
			// For each frequent item in the header table list of the tree in reverse order.
			for(int i = tree.headerList.size()-1; i>=0; i--){
				// get the item
				String item = tree.headerList.get(i);
				
				// get the item support
				int support = mapSupport.get(item);
	
				// Create Beta by concatening prefix Alpha by adding the current item to alpha
				prefix[prefixLength] = item;
				
				// calculate the support of the new prefix beta
				// The support of the suffix will always be the lowest possible support of an item in that suffix
				int betaSupport = (prefixSupport < support) ? prefixSupport: support;
				
				// save beta to the output file
				saveItemset(prefix, prefixLength+1, betaSupport);
				
				// === (A) Construct beta's conditional pattern base ===
				// It is a subdatabase which consists of the set of prefix paths
				// in the FP-tree co-occuring with the prefix pattern.
				List<List<FPNode>> prefixPaths = new ArrayList<List<FPNode>>();
				FPNode path = tree.mapItemNodes.get(item);
				
				// Map to count the support of items in the conditional prefix tree
				// Key: item   Value: support
				Map<String, Integer> mapSupportBeta = new HashMap<String, Integer>();
				
				while(path != null){
					// if the path is not just the root node
					if(path.parent.itemID != "-1"){
						// create the prefixpath
						List<FPNode> prefixPath = new ArrayList<FPNode>();
						// add this node.
						prefixPath.add(path);   // NOTE: we add it just to keep its support,
						// actually it should not be part of the prefixPath
						
						// ####
						int pathCount = path.counter;
						
						//Recursively add all the parents of this node.
						FPNode parent = path.parent;
						while(parent.itemID != "-1"){
							prefixPath.add(parent);
							
							// FOR EACH PATTERN WE ALSO UPDATE THE ITEM SUPPORT AT THE SAME TIME
							// if the first time we see that node id
							if(mapSupportBeta.get(parent.itemID) == null){
								// just add the path count
								mapSupportBeta.put(parent.itemID, pathCount);
							}else{
								// otherwise, make the sum with the value already stored
								mapSupportBeta.put(parent.itemID, mapSupportBeta.get(parent.itemID) + pathCount);
							}
							parent = parent.parent;
						}
						// add the path to the list of prefixpaths
						prefixPaths.add(prefixPath);
					}
					// We will look for the next prefixpath
					path = path.nodeLink;
				}

				// (B) Construct beta's conditional FP-Tree
				// Create the tree.
				FPTree treeBeta = new FPTree();
				// Add each prefixpath in the FP-tree.
				for(List<FPNode> prefixPath : prefixPaths){
					treeBeta.addPrefixPath(prefixPath, mapSupportBeta, minSupportRelative); 
				}  
				
				// Mine recursively the Beta tree if the root has child(s)
				if(treeBeta.root.childs.size() > 0){

					// Create the header list.
					treeBeta.createHeaderList(mapSupportBeta); 
					// recursive call
					fpgrowth(treeBeta, prefix, prefixLength+1, betaSupport, mapSupportBeta);
				}
			}
		}
		
	}


	/**
	 * This method saves all combinations of a prefix path if it has enough support
	 * @param prefix the current prefix
	 * @param prefixLength the current prefix length
	 * @param prefixPath the prefix path
	 * @throws IOException if exception while writting to output file
	 */
	private void saveAllCombinationsOfPrefixPath(FPNode[] fpNodeTempBuffer, int position, 
			String[] prefix, int prefixLength) throws IOException {

		int support = 0;
		// Generate all subsets of the prefixPath except the empty set
		// and output them
		// We use bits to generate all subsets.
		for (long i = 1, max = 1 << position; i < max; i++) {
			
			// we create a new subset
			int newPrefixLength = prefixLength;
			
			// for each bit
			for (int j = 0; j < position; j++) {
				// check if the j bit is set to 1
				int isSet = (int) i & (1 << j);
				// if yes, add the bit position as an item to the new subset
				if (isSet > 0) {
					prefix[newPrefixLength++] = fpNodeTempBuffer[j].itemID;
					if(support == 0) {
						support = fpNodeTempBuffer[j].counter;
					}
				}
			}
			// save the itemset
			saveItemset(prefix, newPrefixLength, support);
		}
	}
	

	/**
	 * This method scans the input database to calculate the support of single items
	 * @param input the path of the input file
	 * @throws IOException  exception if error while writing the file
	 * @return a map for storing the support of each item (key: item, value: support)
	 */
	private  Map<String, Integer> scanDatabaseToDetermineFrequencyOfSingleItems(String input)
			throws FileNotFoundException, IOException {
		// a map for storing the support of each item (key: item, value: support)
		 Map<String, Integer> mapSupport = new HashMap<String, Integer>();
		//Create object for reading the input file
		BufferedReader reader = new BufferedReader(new FileReader(input));
		String line;
		// for each line (transaction) until the end of file
		while( ((line = reader.readLine())!= null)){ 
			// if the line is  a comment, is  empty or is a
			// kind of metadata
			if (line.isEmpty() == true ||  line.charAt(0) == '#' || line.charAt(0) == '%' 	|| line.charAt(0) == '@') {
				continue;
			}
			
			// split the line into items
			String[] lineSplited = line.split(SEPARATOR);
			// for each item
			for(String itemString : lineSplited){  
				// increase the support count of the item
				//Integer item = Integer.parseInt(itemString);
				// increase the support count of the item
				Integer count = mapSupport.get(itemString);
				if(count == null){
					mapSupport.put(itemString, 1);
				}else{
					mapSupport.put(itemString, ++count);
				}
			}
			// increase the transaction count
			transactionCount++;
		}
		// close the input file
		reader.close();
		
		return mapSupport;
	}


	/**
	 * Write a frequent itemset that is found to the output file or
	 * keep into memory if the user prefer that the result be saved into memory.
	 */
	private void saveItemset(String [] itemset, int itemsetLength, int support) throws IOException {
		
		// increase the number of itemsets found for statistics purpose
		itemsetCount++;
		
		// if the result should be saved to a file
		if(writer != null){
			// copy the itemset in the output buffer and sort items
			System.arraycopy(itemset, 0, itemsetOutputBuffer, 0, itemsetLength);
			Arrays.sort(itemsetOutputBuffer, 0, itemsetLength);
			
			// Create a string buffer
			StringBuilder buffer = new StringBuilder();
			// write the items of the itemset
			for(int i=0; i< itemsetLength; i++){
				buffer.append(itemsetOutputBuffer[i]);
				if(i != itemsetLength-1){
					buffer.append(' ');
				}
			}
			// Then, write the support
			buffer.append(" #SUP: ");
			buffer.append(support);
			// write to file and create a new line
			writer.write(buffer.toString());
			writer.newLine();
			
		}// otherwise the result is kept into memory
		else{
			// create an object Itemset and add it to the set of patterns 
			// found.
			String[] itemsetArray = new String[itemsetLength];
			System.arraycopy(itemset, 0, itemsetArray, 0, itemsetLength);
			
			// sort the itemset so that it is sorted according to lexical ordering before we show it to the user
			Arrays.sort(itemsetArray);
			
			Itemset itemsetObj = new Itemset(itemsetArray);
			itemsetObj.setAbsoluteSupport(support);
			patterns.addItemset(itemsetObj, itemsetLength);
		}
	}

	/**
	 * Print statistics about the algorithm execution to System.out.
	 */
	public void printStats() {
		System.out.println("=============  FP-GROWTH 0.96r19 - STATS =============");
		long temps = endTime - startTimestamp;
		System.out.println(" Transactions count from database : " + transactionCount);
		System.out.print(" Max memory usage: " + MemoryLogger.getInstance().getMaxMemory() + " mb \n");
		System.out.println(" Frequent itemsets count : " + itemsetCount); 
		System.out.println(" Total time ~ " + temps + " ms");
		System.out.println("===================================================");
	}

	/**
	 * Get the number of transactions in the last transaction database read.
	 * @return the number of transactions.
	 */
	public int getDatabaseSize() {
		return transactionCount;
	}
	
	/*private void traverseTree(FPNode curNode, List<FPNode> parents, Map<String, Integer> mapSupport){		
		//Traverse recursively to the left-most child node first.
		//This is Pre-order traversal (read the data in the parent first, then traverse from the left-most child to the right-most child)
		
		//Do whatever you want
		if(curNode.itemID != "-1"){
			calculateMaxConf(curNode, parents, mapSupport);
			parents.add(curNode);
		}
		
		for(int i=0; i<curNode.childs.size(); i++){
			traverseTree(curNode.childs.get(i), parents, mapSupport);
		}
	}
	
	private void calculateMaxConf(FPNode curNode, List<FPNode> parents, Map<String, Integer> mapSupport){
		if(parents.size() == 0){
			//Cannot find a max confidence of a single item
			return;
		}
	
		//Min Conf = Sup(parents+currentNode)/Sup(x), where x has 1 item to give max sup
		//Max Conf = Sup(parents+currentNod)/Sup(x), where x has k-1 items to give min sup, k = number of items in (parents+currentNod)
		//We want to find max conf
		//The support of the itemset is equal to the lowest support counter of its member  in the FP-tree
		
		//int support = mapSupport.get(curNode.itemID);
		int support = curNode.counter;
		
		for(FPNode parent: parents){
			//If there's a support that's lower than the current node's support, use that support instead
			
			//support = (support > mapSupport.get(parent.itemID)) ? mapSupport.get(parent.itemID) : support;
			support = (support > parent.counter) ? parent.counter : support;
		}
		System.out.println(curNode.itemID + " minsup=" + support);
		
	}*/
}
