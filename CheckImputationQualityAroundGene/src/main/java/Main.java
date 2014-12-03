import umcg.genetica.io.bed.BedFile;
import umcg.genetica.io.text.TextFile;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by dashazhernakova on 20.11.14.
 */
public class Main {
	HashMap<String, ArrayList<Cohort>> geneToInfos; //gene -> list of cohorts -> list of info scores
	int numCohorts;
	public Main(){

		//geneToInfos = new HashMap<String, ArrayList<Float>>();
		geneToInfos = new HashMap<String, ArrayList<Cohort>>();

	}

	private ArrayList<Cohort> initCohortList(){
		ArrayList<Cohort> cohorts = new ArrayList<Cohort>();
		for (int i = 0; i < numCohorts; i++){
			Cohort c = new Cohort();
			cohorts.add(c);
		}
		return cohorts;

	}
	private void readGeneToInfo(String fname, int datasetNum) throws IOException {
		System.out.println("Reading " + fname);
		TextFile tf = new TextFile(fname, false);
		String[] els;
		while ((els = tf.readLineElems(TextFile.tab)) != null){
			String gene = els[0];
			String[] snp_spl = els[1].split(":");
			String snp = snp_spl[0];
			float info = Float.parseFloat(snp_spl[1]);
			Cohort curCohort = null;

			ArrayList<Cohort> curCohorts = geneToInfos.get(gene);
			if(curCohorts == null){
				curCohorts = initCohortList();
			}
			else{
				curCohort = curCohorts.get(datasetNum);
			}
			if (curCohort == null)
				curCohort = new Cohort();
			curCohort.addInfo(info);
			curCohorts.set(datasetNum, curCohort);
			geneToInfos.put(gene, curCohorts);
		}

		tf.close();
		System.out.println("Finished with loading the file");
	}

	public void readGeneToInfoFiles(String[] fnames) throws IOException {
		numCohorts = fnames.length;
		int cnt = 0;
		for (String fname : fnames){
			readGeneToInfo(fname, cnt);
			cnt += 1;
		}
	}


	private void getStats(){
		int i = 0;
		for (Map.Entry<String,ArrayList<Cohort>> entry: geneToInfos.entrySet()){
			String gene = entry.getKey();
			System.out.print(gene);
			for (Cohort cohort : entry.getValue()){
				System.out.print("\t" + cohort.countMeanInfo() + "\t" + cohort.countNumSNPsInfoLessThanThreshold(0.4f));
			}
			System.out.println();
		}

	}

	public void runAnalysis(String[] fnames) throws IOException {
		readGeneToInfoFiles(fnames);
		getStats();
	}


	public static void main(String[] args) throws IOException {
		Main m = new Main();
		String p = "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/findingCausalVariant/eQTLsFDR0.05-ProbeLevel_genes_pr_noSec_topSnpPvalues_geneCenters_";
		m.runAnalysis(new String[]{p + "LLinfo.txt", p + "RSinfo.txt", p + "LLSinfo.txt", p + "CODAMinfo.txt"});

	}
}
