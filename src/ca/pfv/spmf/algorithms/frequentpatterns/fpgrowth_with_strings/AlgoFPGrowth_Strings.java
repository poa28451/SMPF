package ca.pfv.spmf.algorithms.frequentpatterns.fpgrowth_with_strings;

/* This file is copyright (c) 2008-2013 Philippe Fournier-Viger
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

/** 
 * This is an implementation of the FPGROWTH algorithm (Han et al., 2004) that take
 * as input a transaction database where items are represented by strings rather
 * than integers.
 * FPGrowth is described here:
 * <br/><br/>
 * 
 * Han, J., Pei, J., & Yin, Y. (2000, May). Mining frequent patterns without candidate generation. In ACM SIGMOD Record (Vol. 29, No. 2, pp. 1-12). ACM
 * <br/><br/>
 * 
 * This is an optimized version that saves the result to a file.
 *
 * @see FPTree_Strings
 * @author Philippe Fournier-Viger
 */
public class AlgoFPGrowth_Strings {

	
	// for statistics
	private long startTimestamp; // start time of the latest execution
	private long endTime; // end time of the latest execution
	private int transactionCount = 0; // transaction count in the database
	private int itemsetCount; // number of freq. itemsets found
	
	// minimum support threshold
	public int relativeMinsupp;
	
	// object to write the output file
	BufferedWriter writer = null; 
	
	/**
	 * Editted ver.
	 */
	protected Itemsets patterns = null;
	private String SEPARATOR;
	private double pruneConf;
	
	/**
	 * Default constructor
	 */
	public AlgoFPGrowth_Strings() {
		
	}

	/**
	 * Run the algorithm.
	 * @param input the file path of an input transaction database.
	 * @param output the path of the desired output file
	 * @param minsupp minimum support threshold as a percentage (double)
	 * @throws IOException exception if error while writing the file
	 */
	public Itemsets runAlgorithm(String input, String output, double minsupp, double pruneConf, String separator, boolean isPrune) throws FileNotFoundException, IOException {
		// record the start time
		startTimestamp = System.currentTimeMillis();
		// reinitialize the number of itemsets found to 0
		itemsetCount =0;
		
		/**
		 * Editted ver.
		 */
		this.SEPARATOR = separator;
		this.pruneConf = pruneConf;
		if(output == null){
			writer = null;
			patterns =  new Itemsets("FREQUENT ITEMSETS");
	    }else{ // if the user want to save the result to a file
			patterns = new Itemsets("FREQUENT ITEMSETS");
			writer = new BufferedWriter(new FileWriter(output)); 
		}
		// Prepare the output file
		//writer = new BufferedWriter(new FileWriter(output)); 
		
		// (1) PREPROCESSING: Initial database scan to determine the frequency of each item
		// The frequency is store in a map where:
		// key: item   value: support count
		final Map<String, Integer> mapSupport = new HashMap<String, Integer>();
		// call this method  to perform the database scan
		scanDatabaseToDetermineFrequencyOfSingleItems(input, mapSupport);
		
		// convert the absolute minimum support to a relative minimum support
		// by multiplying by the database size.
		this.relativeMinsupp = (int) Math.ceil(minsupp * transactionCount);
		
		/**
		 * Print support counts of every items
		 */
		System.out.println("124:Support Count: (minsup = " + relativeMinsupp + " (" + minsupp + "%))\n " + mapSupport.size());
		
		// (2) Scan the database again to build the initial FP-Tree
		// Before inserting a transaction in the FPTree, we sort the items
		// by descending order of support.  We ignore items that
		// do not have the minimum support.
		
		// create the FPTree
		FPTree_Strings tree = new FPTree_Strings();
		
		
		BufferedReader reader = new BufferedReader(new FileReader(input));
		String line;
		// for each line (transaction) in the input file until the end of file
		while( ((line = reader.readLine())!= null)){ 
			// if the line is  a comment, is  empty or is a
			// kind of metadata
			if (line.isEmpty() == true ||
					line.charAt(0) == '#' || line.charAt(0) == '%'
							|| line.charAt(0) == '@') {
				continue;
			}
			
			// split the transaction into items
			String[] lineSplited = line.split(SEPARATOR);
			// create an array list to store the items
			List<String> transaction = new ArrayList<String>();
			// for each item in the transaction
			for(String itemString : lineSplited){  
				// if it is frequent, add it to the transaction
				// otherwise not because it cannot be part of a frequent itemset.
				if(mapSupport.get(itemString) >= relativeMinsupp){
					transaction.add(itemString);	
				}
			}
			// sort item in the transaction by descending order of support
			Collections.sort(transaction, new Comparator<String>(){
				public int compare(String item1, String item2){
					// compare the support
					int compare = mapSupport.get(item2) - mapSupport.get(item1);
					// if the same support, we check the lexical ordering!
					if(compare == 0){ 
						return item1.compareTo(item2);
					}
					// otherwise use the support
					return compare;
				}
			});
			// add the sorted transaction to the fptree.
			tree.addTransaction(transaction);
		}
		// close the input file
		reader.close();
		
		// We create the header table for the tree
		/**
		 * i.e. this is the list of every items that Support > minsup, sorted in order from highest sup to lowest sup. 
		 */
		tree.createHeaderList(mapSupport);
		
		//System.out.println("219:Ordered Items: (" + tree.headerList.size() + " items)\n " + tree.headerList + "\n");
		System.out.println("221:Ordered Items: (" + tree.headerList.size() + " items)");
		System.out.print("[");
		for(String item : tree.headerList){
			System.out.print(item + "=" + mapSupport.get(item) + ", ");
		}
		System.out.println("]");
		
		/**
		 * At this point, FP-tree has been fully constructed. What's left is to find Conditional Pattern Bases
		 * 	and further find Conditional FP-trees which are used to create Frequent patterns
		 */
		
		/**
		 * Prune the FP-tree
		 */
		//testTraverse(tree.root);
		if(isPrune && tree.headerList.size() > 0) {
			fptreePruning(tree, mapSupport);
		}
				
		// (5) We start to mine the FP-Tree by calling the recursive method.
		// Initially, the prefix alpha is empty.
		String[] prefixAlpha = new String[0];
		if(tree.headerList.size() > 0) {
			fpgrowth(tree, prefixAlpha, transactionCount, mapSupport);
		}
		
		// close the output file if the result was saved to a file		
		//writer.close();
		if(writer != null){
			writer.close();
		}
		// record the end time
		endTime= System.currentTimeMillis();
		
//		print(tree.root, " ");
		return patterns;
	}
	
	private void fptreePruning(FPTree_Strings tree, Map<String, Integer> mapSupport){
		traverseTree(tree.root, tree, mapSupport);
		
		/**
		 * Try to cut tree branches here
		 */
		//cutBranches(tree.root, tree);
		//System.out.println("-------------------------");
		//testTraverse(tree.root);
		System.out.println(testConfCount(tree.root));
	}
	
	private double testConfCount(FPNode_Strings curNode){
		double totalConf;
		if(curNode.maxConf == -1){
			totalConf = 0;
		} else{
			totalConf = curNode.maxConf;
		}
		for(int i=0; i<curNode.childs.size(); i++){
			totalConf += testConfCount(curNode.childs.get(i));
		}
		return totalConf;
	}
	
	private void testTraverse(FPNode_Strings curNode){
		if(curNode.itemID != null){
			FPNode_Strings temp = curNode;
			System.out.print(curNode.itemID);
			while(temp.parent != null && temp.parent.itemID != null){
				temp = temp.parent;
				System.out.print(" > " + temp.itemID);
			}
		}
		System.out.println();
		
		for(int i=0; i<curNode.childs.size(); i++){
			testTraverse(curNode.childs.get(i));
		}
		//System.out.println(curNode.itemID + " " + curNode.maxConf);
	}
	
	private boolean cutBranches(FPNode_Strings curNode, FPTree_Strings tree){
		boolean isNodeDeleted = false;
		
		for(int i=0; i<curNode.childs.size(); i++){
			FPNode_Strings child = curNode.childs.get(i);
			boolean isChildDeleted = cutBranches(child, tree);
			if(isChildDeleted){
				//The reason that I don't directly remove the node from the parent is because
				//	it'll make the pointer (int i) skip the next child since the next child
				//	will be moved to position i instead
				curNode.childs.remove(i);
				child.parent = null;
				i--;
			}
		}
		
		//Cut a branch
		if(curNode.itemID == null){//Ignore the node if it's the root
			return isNodeDeleted;
		}
		
		//If it's a leaf (node without child) and its maxConf is lower than the specified Conf from user)
		//if(curNode.maxConf != -1 && curNode.childs.size() == 0 && curNode.maxConf <= this.pruneConf){
		if(curNode.childs.size() == 0 && curNode.maxConf <= this.pruneConf){
			//System.out.println(curNode.itemID + " " + curNode.counter + " " + curNode.maxConf);
			isNodeDeleted = true;//Flag to indicate that the node has been removed
			cutNode(curNode, tree);//Remove the node from the tree
		}
		
		return isNodeDeleted;
	}
	
	private void cutNode(FPNode_Strings node, FPTree_Strings tree){
		//node.parent.childs.remove(node);//Remove the node from the tree
		FPNode_Strings tempNode = tree.mapItemNodes.get(node.itemID);
		
		/*while(tempNode != null){
			System.out.print(tempNode.itemID + "/" + tempNode.maxConf + " ");
			tempNode = tempNode.nodeLink;
		}
		System.out.println();
		tempNode = tree.mapItemNodes.get(node.itemID);*/
		
		if(tempNode == node){//Also remove the node from the linked list of nodes with the same item ID
			tree.mapItemNodes.replace(tempNode.itemID, tempNode.nodeLink);//In case if it's the head of the list
		} else{
			while(tempNode.nodeLink != null){
				if(tempNode.nodeLink == node){
					tempNode.nodeLink = node.nodeLink;
					break;
				}
				tempNode = tempNode.nodeLink;
			}
		}
		
		/*tempNode = tree.mapItemNodes.get(node.itemID);
		while(tempNode != null){
			System.out.print(tempNode.itemID + "/" + tempNode.maxConf + " ");
			tempNode = tempNode.nodeLink;
		}
		System.out.println();*/
	}
	
	private void traverseTree(FPNode_Strings curNode, FPTree_Strings tree, Map<String, Integer> mapSupport){		
		//Traverse recursively to the left-most child node first.
		//This is Pre-order traversal (read the data in the parent first, then traverse from the left-most child to the right-most child)
		
		//Do whatever you want
		if(curNode.itemID != null){
			double nodeConf = calculateMaxConf(curNode, tree, mapSupport);
			curNode.maxConf = nodeConf;
			System.out.println(" (" + nodeConf + "%)");
		}
		
		for(int i=0; i<curNode.childs.size(); i++){
			traverseTree(curNode.childs.get(i), tree, mapSupport);
		}
	}
	
	private double calculateMaxConf(FPNode_Strings curNode, FPTree_Strings tree, Map<String, Integer> mapSupport){
		boolean isPrint = false;
		
		if(curNode.parent.itemID == null){
			//Cannot find a max confidence of a single item
			if(isPrint) System.out.print(curNode.itemID);
			return -1;
		}

		// Store all the nodes in this path up to the root, make it into current subtree
		List<FPNode_Strings> subtree = new ArrayList<>();
		FPNode_Strings temp = curNode;		
		while (temp.itemID != null) {
			subtree.add(temp);
			temp = temp.parent;
		}
		
		//Min Conf = Sup(parents+currentNode)/Sup(x), where x has 1 item to give max sup
		//Max Conf = Sup(parents+currentNod)/Sup(x), where x has k-1 items to give min sup, k = number of items in (parents+currentNod)
		//We want to find max conf
		
		//Find the whole subtree's Sup
		int subtreeSupport = calculateSubtreeSupport(subtree, tree, mapSupport, isPrint);
		
		//Find the minimum possible Sup of the subtree
		int subtreeMinSup = transactionCount;
		List<List<FPNode_Strings>> subtreeCombination = findK_1Permutation(subtree); //K-1 subset of this subtree will always have the minimum Sup
		for(List<FPNode_Strings> combi : subtreeCombination){
			int combiSup = calculateSubtreeSupport(combi, tree, mapSupport, false);
			subtreeMinSup = (subtreeMinSup < combiSup) ? subtreeMinSup : combiSup;
			
			if(subtreeMinSup == 1){//There's no need to calculate Sup anymore if it reaches 1, because that's the minimum value possible
				break;
			}
		}
		
		if(isPrint){
			System.out.print(" // minsup = " + subtreeMinSup);
			System.out.print(" // Max conf = " + subtreeSupport + "/" + subtreeMinSup);
		}
		
		return (double) subtreeSupport / subtreeMinSup;
	}
	
	private int calculateSubtreeSupport(List<FPNode_Strings> subtree, FPTree_Strings tree, Map<String, Integer> mapSupport, boolean isPrint){
		int subtreeSupport = 0;
		/*List<FPNode> subtree = new ArrayList<>();
		FPNode temp = curNode;
		
		// Store all the nodes in this subtree, make it into current subtree
		while(temp.itemID != null){
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
		
		FPNode_Strings temp = tree.mapItemNodes.get(subtree.get(0).itemID);
		//System.out.println(temp.itemID + ": ");
		int matchedParents = 0;
		while(temp != null){// Do a loop until we have already considered every nodes with this itemID
			int i = 1;
			int nodeSupport = temp.counter; //Remember this FP-tree node's support
			FPNode_Strings parent = temp.parent;
			//System.out.print(parent.itemID + " ");
			
			while(parent.itemID != null && i < subtree.size()){ //If it isn't a root and the whole of current subtree hasn't been found yet
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
	
	private List<List<FPNode_Strings>> findK_1Permutation(List<FPNode_Strings> subtree){
		List<List<FPNode_Strings>> subtreeCombination = new ArrayList<List<FPNode_Strings>>();		
		List<FPNode_Strings> combi;
		
		for(int i=0; i<subtree.size(); i++){
			combi = new ArrayList<>();
			for(FPNode_Strings item : subtree){
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

//	private void print(FPNode node, String indentation) {
//		System.out.println(indentation + "NODE : " + node.itemID + " COUNTER" + node.counter);
//		for(FPNode child : node.childs) {
//			print(child, indentation += "\t");
//		}
//	}

	/**
	 * This method scans the input database to calculate the support of single items
	 * @param input the path of the input file
	 * @param mapSupport a map for storing the support of each item (key: item, value: support)
	 * @throws IOException  exception if error while writing the file
	 */
	private void scanDatabaseToDetermineFrequencyOfSingleItems(String input,
			final Map<String, Integer> mapSupport)
			throws FileNotFoundException, IOException {
		//Create object for reading the input file
		BufferedReader reader = new BufferedReader(new FileReader(input));
		String line;
		// for each line (transaction) until the end of file
		while( ((line = reader.readLine())!= null)){ 
			// if the line is  a comment, is  empty or is a
			// kind of metadata
			if (line.isEmpty() == true ||
					line.charAt(0) == '#' || line.charAt(0) == '%'
							|| line.charAt(0) == '@') {
				continue;
			}
			
			// split the transaction into items
			String[] lineSplited = line.split(SEPARATOR);
			 // for each item in the transaction
			for(String itemString : lineSplited){ 
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
	}


	/**
	 * This method mines pattern from a Prefix-Tree recursively
	 * @param tree  The Prefix Tree
	 * @param prefix  The current prefix "alpha"
	 * @param mapSupport The frequency of each item in the prefix tree.
	 * @throws IOException   exception if error writing the output file
	 */
	private void fpgrowth(FPTree_Strings tree, String[] prefixAlpha, int prefixSupport, Map<String, Integer> mapSupport) throws IOException {
		// We need to check if there is a single path in the prefix tree or not.
		if(tree.hasMoreThanOnePath == false){
			// That means that there is a single path, so we 
			// add all combinations of this path, concatenated with the prefix "alpha", to the set of patterns found.
			addAllCombinationsForPathAndPrefix(tree.root.childs.get(0), prefixAlpha); // CORRECT?
			
		}else{ // There is more than one path
			fpgrowthMoreThanOnePath(tree, prefixAlpha, prefixSupport, mapSupport);
		}
	}
	
	/**
	 * Mine an FP-Tree having more than one path.
	 * @param tree  the FP-tree
	 * @param prefix  the current prefix, named "alpha"
	 * @param mapSupport the frequency of items in the FP-Tree
	 * @throws IOException   exception if error writing the output file
	 */
	private void fpgrowthMoreThanOnePath(FPTree_Strings tree, String [] prefixAlpha, int prefixSupport, Map<String, Integer> mapSupport) throws IOException {
		// We process each frequent item in the header table list of the tree in reverse order.
		for(int i= tree.headerList.size()-1; i>=0; i--){
			String item = tree.headerList.get(i);
			
			int support = mapSupport.get(item);
			// if the item is not frequent, we skip it
			if(support <  relativeMinsupp){
				continue;
			}
			// Create Beta by concatening Alpha with the current item
			// and add it to the list of frequent patterns
			String [] beta = new String[prefixAlpha.length+1];
			System.arraycopy(prefixAlpha, 0, beta, 0, prefixAlpha.length);
			beta[prefixAlpha.length] = item;
			
			/*for(int k=0; k<beta.length; k++){
				System.out.print(beta[k] + " ");
			}System.out.println();*/
			
			// calculate the support of beta
			int betaSupport = (prefixSupport < support) ? prefixSupport: support;
			// save beta to the output file
			writeItemsetToFile(beta, betaSupport);
			
			// === Construct beta's conditional pattern base ===
			// It is a subdatabase which consists of the set of prefix paths
			// in the FP-tree co-occuring with the suffix pattern.
			List<List<FPNode_Strings>> prefixPaths = new ArrayList<List<FPNode_Strings>>();
			FPNode_Strings path = tree.mapItemNodes.get(item);
			while(path != null){
				// if the path is not just the root node
				if(path.parent.itemID != null){
					// create the prefixpath
					List<FPNode_Strings> prefixPath = new ArrayList<FPNode_Strings>();
					// add this node.
					prefixPath.add(path);   // NOTE: we add it just to keep its support,
					// actually it should not be part of the prefixPath
					
					//Recursively add all the parents of this node.
					FPNode_Strings parent = path.parent;
					while(parent.itemID != null){
						prefixPath.add(parent);
						parent = parent.parent;
					}
					// add the path to the list of prefixpaths
					prefixPaths.add(prefixPath);
					
					for(FPNode_Strings node : prefixPath){
						System.out.print(node.itemID + "=" + node.counter + " ");
					}System.out.println();
				}
				// We will look for the next prefixpath
				path = path.nodeLink;
			}
			
			// (A) Calculate the frequency of each item in the prefixpath
			Map<String, Integer> mapSupportBeta = new HashMap<String, Integer>();
			// for each prefixpath
			for(List<FPNode_Strings> prefixPath : prefixPaths){
				// the support of the prefixpath is the support of its first node.
				int pathCount = prefixPath.get(0).counter;  
				 // for each node in the prefixpath,
				// except the first one, we count the frequency
				for(int j=1; j<prefixPath.size(); j++){ 
					FPNode_Strings node = prefixPath.get(j);
					// if the first time we see that node id
					if(mapSupportBeta.get(node.itemID) == null){
						// just add the path count
						mapSupportBeta.put(node.itemID, pathCount);
					}else{
						// otherwise, make the sum with the value already stored
						mapSupportBeta.put(node.itemID, mapSupportBeta.get(node.itemID) + pathCount);
					}
				}
			}
			
			// (B) Construct beta's conditional FP-Tree
			FPTree_Strings treeBeta = new FPTree_Strings();
			// add each prefixpath in the FP-tree
			for(List<FPNode_Strings> prefixPath : prefixPaths){
				treeBeta.addPrefixPath(prefixPath, mapSupportBeta, relativeMinsupp); 
			}  
			// Create the header list.
			treeBeta.createHeaderList(mapSupportBeta); 
			
			// Mine recursively the Beta tree if the root as child(s)
			if(treeBeta.root.childs.size() > 0){
				// recursive call
				fpgrowth(treeBeta, beta, betaSupport, mapSupportBeta);
			}
		}
		
	}

	/**
	 * This method is for adding recursively all combinations of nodes in a path, concatenated with a given prefix,
	 * to the set of patterns found.
	 * @param nodeLink the first node of the path
	 * @param prefix  the prefix
	 * @param minsupportForNode the support of this path.
	 * @throws IOException 
	 */
	private void addAllCombinationsForPathAndPrefix(FPNode_Strings node, String[] prefix) throws IOException {
		// Concatenate the node item to the current prefix
		String [] itemset = new String[prefix.length+1];
		System.arraycopy(prefix, 0, itemset, 0, prefix.length);
		itemset[prefix.length] = node.itemID;

		// save the resulting itemset to the file with its support
		writeItemsetToFile(itemset, node.counter);
			
		if(node.childs.size() != 0) {
			addAllCombinationsForPathAndPrefix(node.childs.get(0), itemset);
			addAllCombinationsForPathAndPrefix(node.childs.get(0), prefix);
		}
	}
	

	/**
	 * Write a frequent itemset that is found to the output file.
	 */
	private void writeItemsetToFile(String [] itemset, int support) throws IOException {
		// increase the number of itemsets found for statistics purpose
		itemsetCount++;
		
		if(writer != null){
			// create a string buffer 
			StringBuilder buffer = new StringBuilder();
			// write items from the itemset to the StringBuilder
			for(int i=0; i< itemset.length; i++){
				buffer.append(itemset[i]);
				if(i != itemset.length-1){
					buffer.append(' ');
				}
			}
			// append the support of the itemset
			buffer.append(" #SUP: ");
			buffer.append(support);
			// write the strinbuffer and create a newline so that we are
			// ready for the next itemset to be written
			writer.write(buffer.toString());
			writer.newLine();
		}
		
		// create an object Itemset and add it to the set of patterns 
		// found.
		int itemsetLength = itemset.length;
		String[] itemsetArray = new String[itemsetLength];
		System.arraycopy(itemset, 0, itemsetArray, 0, itemsetLength);

		// sort the itemset so that it is sorted according to lexical ordering before we show it to the user
		Arrays.sort(itemsetArray);

		Itemset itemsetObj = new Itemset(itemsetArray);
		itemsetObj.setAbsoluteSupport(support);
		patterns.addItemset(itemsetObj, itemsetLength);
	}

	/**
	 * Print statistics about the algorithm execution to System.out.
	 */
	public void printStats() {
		System.out
				.println("=============  FP-GROWTH - STATS =============");
		long temps = endTime - startTimestamp;
		System.out.println(" Transactions count from database : " + transactionCount);
		System.out.println(" Frequent itemsets count : " + itemsetCount); 
		System.out.println(" Total time ~ " + temps + " ms");
		System.out
				.println("===================================================");
	}
	
	/**
	 * Editted ver.
	 */
	public int getDatabaseSize() {
		return transactionCount;
	}
}
