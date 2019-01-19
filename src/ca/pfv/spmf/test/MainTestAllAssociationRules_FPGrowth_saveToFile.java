package ca.pfv.spmf.test;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URL;

import ca.pfv.spmf.algorithms.associationrules.agrawal94_association_rules.AlgoAgrawalFaster94;
import ca.pfv.spmf.algorithms.frequentpatterns.fpgrowth.AlgoFPGrowth;
import ca.pfv.spmf.patterns.itemset_array_integers_with_count.Itemsets;
/**
 * Example of how to mine all association rules with FPGROWTH and save
 * the result to a file, from the source code.
 * 
 * @author Philippe Fournier-Viger (Copyright 2008)
 */
public class MainTestAllAssociationRules_FPGrowth_saveToFile {

	public static void main(String [] arg) throws IOException{
		//String input = "C:\\Users\\Lightning\\Desktop\\contextIGB(comma)test.txt";
		//String input = "C:\\Users\\Lightning\\Desktop\\contextIGB(comma).txt";
		String input = "D:\\Dropbox\\Thesis\\Data\\transaction.txt";
		//String input = "C:\\Users\\Lightning\\Desktop\\randomdata.txt";
		String output = "C:\\Users\\Lightning\\Desktop\\output.txt";


		
		// STEP 1: Applying the FP-GROWTH algorithm to find frequent itemsets
		double minsupp = 0.1;
		//double minsupp = 0.03;
		double pruneConf = 0.6;
		AlgoFPGrowth fpgrowth = new AlgoFPGrowth();
		Itemsets patterns = fpgrowth.runAlgorithm(input, null, minsupp, pruneConf);
		patterns.printItemsets(1);
		fpgrowth.printStats();
		int databaseSize = fpgrowth.getDatabaseSize();
		
		// STEP 2: Generating all rules from the set of frequent itemsets (based on Agrawal & Srikant, 94)
		double  minconf = 0.5;
		AlgoAgrawalFaster94 algoAgrawal = new AlgoAgrawalFaster94();
		algoAgrawal.runAlgorithm(patterns, output, databaseSize, minconf);
		algoAgrawal.printStats();
	}
	
	public static String fileToPath(String filename) throws UnsupportedEncodingException{
		URL url = MainTestAllAssociationRules_FPGrowth_saveToFile.class.getResource(filename);
		 return java.net.URLDecoder.decode(url.getPath(),"UTF-8");
	}
}
