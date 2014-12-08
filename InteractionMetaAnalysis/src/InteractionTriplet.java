import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by dashazhernakova on 07.10.14.
 */
public class InteractionTriplet implements Comparable<InteractionTriplet> {
	String id;
	String covariate;
	String snp;
	String gene;

	float mainZ;
	float snpZ;
	float covariateZ;
	float interactionZ;
	float interactionZflipped;

	int numSamples;

	char[] alleles;
	char alleleAssessed;

	public int cohortIndex;

	public InteractionTriplet(String line){
		String[] els = line.split("\t");
		snp = els[0];
		gene = els[1];
		covariate = els[2];
		id = snp + "_" + gene + "_" + covariate;

		snpZ = Float.parseFloat(els[3]);
		covariateZ = Float.parseFloat(els[4]);
		interactionZ = Float.parseFloat(els[5]);
		mainZ = Float.parseFloat(els[6]);
		interactionZflipped = Float.parseFloat(els[7]); // remove for the old version of Interaction Analysis

		numSamples = Integer.parseInt(els[8]);
	}

	public InteractionTriplet(String line, int cnt){
		cohortIndex = cnt;

		String[] els = line.split("\t");
		snp = els[0];
		gene = els[1];
		covariate = els[2];
		id = snp + "_" + gene + "_" + covariate;

		snpZ = Float.parseFloat(els[3]);
		covariateZ = Float.parseFloat(els[4]);
		interactionZ = Float.parseFloat(els[5]);
		mainZ = Float.parseFloat(els[6]);
		//interactionZflipped = Float.parseFloat(els[7]); // remove for the old version of Interaction Analysis

		//numSamples = Integer.parseInt(els[8]);
		numSamples = Integer.parseInt(els[7]);
	}

	public void readSNPInfo(String line){
		String[] els = line.split("\t");
		/*
		if (!els[0].equals(snp))
			return;
		 */

		String[] strAlleles = els[3].split("/");
		alleles = new char[] {strAlleles[0].charAt(0), strAlleles[1].charAt(0)};

		Arrays.sort(alleles);
		alleleAssessed = els[4].charAt(0);
	}

	// for the paralel reading of the interaction results
	/*public void readSNPInfo(char[] inAlleles){
		alleles = new char[] {inAlleles[cohortIndex*3], inAlleles[cohortIndex*3 + 1]};
		Arrays.sort(alleles);
		alleleAssessed = inAlleles[cohortIndex*3 + 2];
	}*/

	public void readSNPInfo(char[] inAlleles){
		alleles = new char[] {inAlleles[0], inAlleles[1]};
		Arrays.sort(alleles);
		alleleAssessed = inAlleles[2];
	}

	public void revertSigns(){
		mainZ = - mainZ;
		snpZ = - snpZ;
		interactionZ = - interactionZ;
		interactionZflipped = - interactionZflipped;
		//covariate doesn't depend on the genotypes - no need to revert the sign
	}

	/*
	public void correctAllSignsRelativeToMainEffectSign(){
		if (mainZ < 0){
			mainZ = - mainZ;
			snpZ = - snpZ;
			covariateZ = - covariateZ;
			interactionZ = - interactionZ;
		}
	}
	public float getCorrectInteractionValueSign(){
		if (mainZ < 0)
			return - interactionZ;
		return interactionZ;
	}

	public float getCorrectSnpZscoreValueSign(){
		if (mainZ < 0)
			return - snpZ;
		return snpZ;
	}

	public float getCorrectMainZscoreValueSign(){
		if (mainZ < 0)
			return - mainZ;
		return snpZ;
	}*/

	public int compareTo(InteractionTriplet anotherTriplet) {
		return id.compareTo(anotherTriplet.id);
	}
}
