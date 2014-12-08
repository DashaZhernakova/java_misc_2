import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by dashazhernakova on 07.10.14.
 */
public class MetaAnalyser_synch {
	HashMap<String, InteractionTriplet[]> interactionResults;
	HashMap<String, char[]> SNPalleles;
	int numCohorts;
	String[] cohortNames;
	public MetaAnalyser_synch(int numberCohorts){
		interactionResults = new HashMap<String, InteractionTriplet[]>();
		SNPalleles = new HashMap<String, char[]>();
		numCohorts = numberCohorts;
		cohortNames = new String[numCohorts];
	}

	private TextFile[] initFiles(String[] fnames) throws IOException {
		TextFile[] files = new TextFile[numCohorts];
		for (int i = 0; i < numCohorts; i++){
			TextFile f = new TextFile(fnames[i], false);
			files[i] = f;
		}
		return files;
	}
	private void readInteractionResults(String[] fnames) throws IOException {
		TextFile[] files = initFiles(fnames);
		InteractionTriplet[] curTriplets = new InteractionTriplet[numCohorts];
		String line = "";
		String curSnp = "";
		char[] curAlleles = null;

		//read first lines of each cohort's file
		for (int i = 0; i < numCohorts; i++){
			TextFile f = files[i];
			f.readLine();
			curTriplets[i] = initTriplet(f.readLine(), i);
		}


		while (true){
			// Choose those triplets that have the same id and that are lexicographically first
			ArrayList<Integer> upperLevelOriginalIndices = chooseUpperIndices(curTriplets);
			InteractionTriplet[] upperTriplets = getArrayElementsByIndices(curTriplets, upperLevelOriginalIndices);
			//run meta-analysis on them
			String outString = calculateWeightedZscore(upperTriplets);
			System.out.println(outString);

			// For those cohorts that were used in meta-analsyis read next results line
			for (int idx : upperLevelOriginalIndices){
				TextFile f = files[idx];
				if ((line = f.readLine()) != null)
					curTriplets[idx] = initTriplet(line, idx);
			}
		}
	}

	private InteractionTriplet initTriplet(String line, int i){
		String snp = line.split("\t")[0];
		InteractionTriplet triplet = new InteractionTriplet(line, i);
		triplet.readSNPInfo(SNPalleles.get(snp));
		return triplet;
	}

	/**
	 * Leaves only the triplets that are situated at indices specified by @param upperLevelOriginalIndices. The rest are set to null
	 * @param triplets - original triplets
	 * @param upperLevelOriginalIndices - indices of elements to keep
	 * @return
	 */
	private InteractionTriplet[] getArrayElementsByIndices(InteractionTriplet[] triplets, ArrayList<Integer> upperLevelOriginalIndices){
		InteractionTriplet[] upperTriplets = new InteractionTriplet[triplets.length];

		for (int idx : upperLevelOriginalIndices){
			upperTriplets[idx] = triplets[idx];
		}
		return upperTriplets;
	}

	/**
	 * Chooses the indices of triplets to use in the current meta-analysis: those having the same id and being first alphabetically
	 * @param curTriplets - current triplets coming from all cohprts
	 * @return
	 */
	private ArrayList<Integer> chooseUpperIndices(InteractionTriplet[] curTriplets){
		InteractionTriplet[] sortedCurTriplets = curTriplets.clone();
		Arrays.sort(sortedCurTriplets);

		//Lexicographically top triplet
		InteractionTriplet upperTriplet = sortedCurTriplets[0];
		String upperId = upperTriplet.id;

		// Array of triplets that have the same id and will be used in interaction analysis
		//InteractionTriplet[] upperTriplets = new InteractionTriplet[numCohorts];
		//upperTriplets[0] = upperTriplet;

		// indices of triplets that will be used for meta-analysis in the original cohort array (in order to increment the file counter for these cohorts)
		ArrayList<Integer> upperLevelOriginalIndices = new ArrayList<Integer>(numCohorts);
		upperLevelOriginalIndices.add(upperTriplet.cohortIndex);

		for (int i = 1; i < numCohorts; i++){
			InteractionTriplet triplet = sortedCurTriplets[i];
			if (triplet.id.equals(upperId)){
				//upperTriplets[i] = triplet;
				upperLevelOriginalIndices.add(triplet.cohortIndex);
			}
		}
		return upperLevelOriginalIndices;

	}



	/**
	 * Loads SNP information from cohort number @param cohortIndex: what alleles does it have and which one was tested
	 * @param fname - file path of the SNPSummaryStatistics.txt
	 * @throws java.io.IOException
	 */
	private void readSNPstats(String fname, int cohortIndex) throws IOException {
		System.out.println("Loading SNP stats from: " + fname);
		TextFile tf = new TextFile(fname, false);
		String[] els = tf.readLineElems(TextFile.tab);

		while ((els = tf.readLineElems(TextFile.tab)) != null){
			String[] strAlleles = els[3].split("/");
			String snp = els[0];
			char [] alleles = SNPalleles.get(snp);

			if (alleles == null){
				alleles = new char[numCohorts*3];
			}
			//fill allele info into the array (2 SNP alleles followed by the assessed allele)
			alleles[cohortIndex*3] = strAlleles[0].charAt(0);
			alleles[cohortIndex*3 + 1] = strAlleles[1].charAt(0);
			alleles[cohortIndex*3 + 2] = els[4].charAt(0);

			SNPalleles.put(snp, alleles);
		}
		tf.close();

		System.out.println("Read allele info for " + SNPalleles.size() + " SNPs from " + fname);

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

			//skip empty triplets
			if (triplet == null)
				continue;

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
	 * @param fnames - paths to InteractionResults.txt files
	 * @throws java.io.IOException
	 */
	public void runMetaAnalysis(String[] fnames) throws IOException {

		String outFname = fnames[fnames.length - 1];
		TextFile out = new TextFile(outFname, true);

		for (int idx = 0; idx < numCohorts; idx++){
			String fname = fnames[idx];

			//manage cohort names from the file path (folders containg the InteractionResults.txt file are considered to be cohort names)
			String[] fname_spl = fname.split("/");
			cohortNames[idx] = fname_spl[fname_spl.length - 2];

			//read allele info
			readSNPstats(fname.replace("InteractionResults","SNPSummaryStatistics"), idx);
		}


		// Write header
		String header = "SNP\tgene\tTF\talleles\talleleAssesed\t";
		for (String cohortName : cohortNames){
			header += cohortName + "_main_zscore\t" + cohortName + "_interaction_zscore\t" + cohortName + "_snp_zscore\t" + cohortName + "_numSamples";
		}

		header += "\tmeta_main_z-score\tmeta_interaction_z-score\tmeta_snp_z-score\tmeta_interaction_z-score_flipped";
		out.writeln(header);

		// Run meta-analysis
		readInteractionResults(fnames);

		//TODO: write results

		out.close();
		System.out.println("Finished");
	}

	public static void main(String[] args) throws IOException {

		args = new String[] {"/Users/dashazhernakova/Documents/UMCG/data/BBMRI/interactionWithTFs/szymon+maarten+sasha_TFs/LL/InteractionResults.txt.gz",
				"/Users/dashazhernakova/Documents/UMCG/data/BBMRI/interactionWithTFs/szymon+maarten+sasha_TFs/RS/InteractionResults.txt.gz",
				"/Users/dashazhernakova/Documents/UMCG/data/BBMRI/interactionWithTFs/szymon+maarten+sasha_TFs/LLS/InteractionResults.txt.gz",
				"/Users/dashazhernakova/Documents/UMCG/tmp.txt"};
		MetaAnalyser_synch meta = new MetaAnalyser_synch(args.length - 1);

		meta.runMetaAnalysis(args);
	}
}
