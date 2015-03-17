import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by dashazhernakova on 07.10.14.
 */
public class MetaAnalyser {
	HashMap<String, InteractionTriplet[]> interactionResults;
	int numCohorts;
	String[] cohortNames;
	public MetaAnalyser(int numberCohorts){
		interactionResults = new HashMap<String, InteractionTriplet[]>();
		numCohorts = numberCohorts;
		cohortNames = new String[numCohorts];
	}

	/**
	 * Loads interaction analysis results to the results hashmap
	 * @param fname - path to the interaction analysis results file
	 * @param SNPalleles - previously loaded tested SNP alleles
	 * @param cohortIndex - index of the cohort
	 * @throws IOException
	 */
	private void readInteractionResults(String fname, HashMap<String, char[]> SNPalleles, int cohortIndex) throws IOException {
		//HashMap<String, InteractionTriplet> interactionResults = new HashMap<String, InteractionTriplet>();
		TextFile tf = new TextFile(fname, false);
		String line = tf.readLine();

		String curSnp = "";
		char[] curAlleles = null;

		while ((line = tf.readLine()) != null){
			String[] els = line.split("\t");
			InteractionTriplet triplet = new InteractionTriplet(line);
			String snp = els[0];

			//load allele info
			if (snp.equals(curSnp)){
				triplet.readSNPInfo(curAlleles);
			}
			else{
				curSnp = snp;
				curAlleles = SNPalleles.get(curSnp);
				triplet.readSNPInfo(curAlleles);
			}

			//add to the big results map
			InteractionTriplet[] resultList = interactionResults.get(triplet.id);
			if (resultList == null){
				resultList = new InteractionTriplet[numCohorts];
				interactionResults.put(triplet.id, resultList);
			}
			resultList[cohortIndex] = triplet;

		}
		System.out.println("Read " + interactionResults.size() + " interaction analysis results from " + fname);
		tf.close();


	}

	/**
	 * Loads SNP information: what alleles does it have and which one was tested
	 * @param fname - file path of the SNPSummaryStatistics.txt
	 * @return - hashmap of SNP vs an array containing it's 2 alleles and the tested allele
	 * @throws IOException
	 */
	private HashMap<String, char[]> readSNPstats(String fname) throws IOException {
		HashMap<String, char[]> SNPalleles = new HashMap<String, char[]>();
		TextFile tf = new TextFile(fname, false);
		String[] els = tf.readLineElems(TextFile.tab);

		while ((els = tf.readLineElems(TextFile.tab)) != null){
			String[] strAlleles = els[3].split("/");
			char [] alleles = new char[] {strAlleles[0].charAt(0), strAlleles[1].charAt(0), els[4].charAt(0)};
			SNPalleles.put(els[0], alleles);
		}
		tf.close();

		System.out.println("Read allele info for " + SNPalleles.size() + " SNPs from " + fname);
		return SNPalleles;
	}

	/**
	 * Loads all inteaction results for a cohort
	 * @param interactionFname - the path to the InteractionResults.txt file
	 * @param SnpStatsFname - the path to the SNPSummaryStatistics.txt file
	 * @param cohortIndex - - index of the current cohort
	 */
	public void loadData(String interactionFname, String SnpStatsFname, int cohortIndex){
		System.out.println("Loading data from: " + interactionFname + "\nSNP stats from: " + SnpStatsFname);
		HashMap<String, char[]> SNPalleles = null;

		String[] fname_spl = interactionFname.split("/");
		cohortNames[cohortIndex] = fname_spl[fname_spl.length - 2];

		try {
			SNPalleles = readSNPstats(SnpStatsFname);
		} catch (IOException e) {
			e.printStackTrace();
		}


		try {
			readInteractionResults(interactionFname, SNPalleles, cohortIndex);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Calculates meta z-scores for a triplet of eSNP, cis-gene and covariate gene
	 * @param triplets
	 * @return
	 */
	private String calculateWeightedZscore(InteractionTriplet[] triplets){
		char[] alleles = new char[2];
		char alleleAssesed = 0;
		float interactionZscore = 0, mainZscore = 0, snpZscore = 0, denominator = 0, interactionZscoreFlipped = 0;
		String outputStr = "";

		// assign base alleles as the alleles of the first triplet, write triplet info
		for (int i = 0; i < triplets.length; i++){
			InteractionTriplet triplet = triplets[i];

			if (triplet != null){
				alleles = triplet.alleles;
				alleleAssesed = triplet.alleleAssessed;
				outputStr += triplet.snp + "\t" + triplet.gene + "\t" + triplet.covariate + "\t" + alleles[0] + "/" + alleles[1] + "\t" + alleleAssesed;
				break;
			}
		}

		//calculate z-scores, write results
		for (int i = 0; i < triplets.length; i++){
			InteractionTriplet triplet = triplets[i];

			if (triplet == null){
				outputStr += "\tNA\tNA\tNA\tNA";
				continue;
			}

			// the SNP has same alleles
			if (Arrays.equals(triplet.alleles, alleles)){
				//the assesed allele is different => revert all z-scores signs
				if (triplet.alleleAssessed != alleleAssesed){
					triplet.revertSigns();
				}

				//with reverting the signs relative to main Z score
				//interactionZscore += triplet.getCorrectInteractionValueSign() * triplet.numSamples;
				//snpZscore += triplet.getCorrectSnpZscoreValueSign() * triplet.numSamples;


				//not reverting of the z-scores
				interactionZscore += triplet.interactionZ * triplet.numSamples;
				interactionZscoreFlipped += triplet.interactionZflipped * triplet.numSamples;

				snpZscore += triplet.snpZ * triplet.numSamples;

				mainZscore += triplet.mainZ * triplet.numSamples;

				denominator += Math.pow(triplet.numSamples, 2);
				outputStr += "\t" + triplet.mainZ + "\t" + triplet.interactionZ + "\t" + triplet.snpZ + "\t" + triplet.numSamples;
			}
			else{
				outputStr +="\tNA\tNA\tNA\tNA";
			}

		}

		outputStr += "\t" + mainZscore/Math.sqrt(denominator) + "\t" + interactionZscore/Math.sqrt(denominator) + "\t" + snpZscore/Math.sqrt(denominator) + "\t" + interactionZscoreFlipped/Math.sqrt(denominator);
		return outputStr;
	}

	/**
	 * Runs the meta-analysis
	 * @param outFname
	 * @throws IOException
	 */
	public void runMetaAnalysis(String outFname) throws IOException {

		TextFile out = new TextFile(outFname, true);

		String header = "SNP\tgene\tTF\talleles\talleleAssesed\t";
		for (String cohortName : cohortNames)
			header += cohortName + "_main_zscore\t" + cohortName + "_interaction_zscore\t" + cohortName + "_snp_zscore\t" + cohortName + "_numSamples";

		header += "\tmeta_main_z-score\tmeta_interaction_z-score\tmeta_snp_z-score\tmeta_interaction_z-score_flipped";
		out.writeln(header);
		for (Map.Entry<String, InteractionTriplet[]> entry : interactionResults.entrySet()){
			out.writeln(calculateWeightedZscore(entry.getValue()));
		}
		out.close();
		System.out.println("Finished");
	}

	/**
	 *
	 * @param args - paths to InteractionResults.txt files, followed by the path to output file
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {

		//String[] args = new String[]{"/Users/dashazhernakova/Documents/UMCG/data/BBMRI/interactionWithTFs/szymon+maarten+sasha_TFs/LL/InteractionResults.txt.gz", "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/interactionWithTFs/szymon+maarten+sasha_TFs/LLS/InteractionResults.txt.gz", "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/interactionWithTFs/szymon+maarten+sasha_TFs/RS/InteractionResults.txt.gz", "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/interactionWithTFs/szymon+maarten+sasha_TFs/CODAM/InteractionResults.txt.gz", "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/interactionWithTFs/szymon+maarten+sasha_TFs/meta_tmp.txt"};
		MetaAnalyser meta = new MetaAnalyser(args.length - 1);

		for (int i = 0; i < args.length - 1; i++){
			meta.loadData(args[i], args[i].replace("InteractionResults","SNPSummaryStatistics"), i);
		}
		meta.runMetaAnalysis(args[args.length - 1]);

	}
}
