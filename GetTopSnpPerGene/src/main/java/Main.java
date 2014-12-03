import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

import java.io.IOException;
import java.util.*;

/**
 * Created by dashazhernakova on 05.11.14.
 */
public class Main {
	HashMap<String, ArrayList<EQTL>> eqtlsPerGene;

	public void readEqtls(String fname) throws IOException {
		System.out.println("Reading eQTLs");
		QTLTextFile file = new QTLTextFile(fname, false);
		eqtlsPerGene = new HashMap<String, ArrayList<EQTL>>();

		for (Iterator<EQTL> it = file.getEQtlIterator(); it.hasNext();) {
			EQTL eqtl = it.next();

			ArrayList<EQTL> geneEqtls = eqtlsPerGene.get(eqtl.getProbe());

			if (geneEqtls == null){
				geneEqtls = new ArrayList<EQTL>();
			}

			geneEqtls.add(eqtl);
			eqtlsPerGene.put(eqtl.getProbe(), geneEqtls);

		}

		file.close();
		System.out.println("Finished reading. Loaded " + eqtlsPerGene.size() + " genes");
	}

	public void getTopSnpPvalueCutoff(String outFname) throws IOException {
		System.out.println("Writing eQTLs");
		TextFile out = new TextFile(outFname, true);
		out.writeln("gene\tSNP1_id\tSNP1_chr\tSNP1_pos\tSNP1_pvalue\tSNP2_id\tSNP2_chr\tSNP2_pos\tSNP2_pvalue");
		for (String gene : eqtlsPerGene.keySet()){
			ArrayList<EQTL> eqtls = eqtlsPerGene.get(gene);
			Collections.sort(eqtls, new PvalueComparator());
			if (eqtls.size() > 1){
				EQTL eqtl1 = eqtls.get(0);
				EQTL eqtl2 = eqtls.get(1);
				if (eqtl2.getPvalue() / eqtl1.getPvalue() > 1000)
					out.writeln(gene + "\t" + eqtl1.getRsName() + "\t" + eqtl1.getRsChr() + "\t" + eqtl1.getRsChrPos() + "\t" + eqtl1.getPvalue() + "\t" + eqtl2.getRsName() + "\t" + eqtl2.getRsChr() + "\t" + eqtl2.getRsChrPos() + "\t" + eqtl2.getPvalue());
			}
			else{
				EQTL eqtl1 = eqtls.get(0);
				out.writeln(gene + "\t" + eqtl1.getRsName() + "\t" + eqtl1.getRsChr() + "\t" + eqtl1.getRsChrPos() + "\t" + eqtl1.getPvalue() + "\tNA\tNA\tNA\tNA");
			}
		}
		out.close();
		System.out.println("Finished");
	}

	public class PvalueComparator implements Comparator<EQTL> {

		@Override
		public int compare(EQTL a, EQTL b) {

			Double pval1 = a.getPvalue();
			Double pval2 = b.getPvalue();
			return pval1.compareTo(pval2);

		}
	}

	public static void main(String[] args) throws IOException {
		Main m = new Main();
		m.readEqtls("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/first_run_final/genes/LL+RS+CODAM+LLS_eqtls_genes_snpsInAllDatasets_18072014/eQTLsFDR0.05-ProbeLevel_genes_pr_noSec.txt.gz");
		m.getTopSnpPvalueCutoff("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/first_run_final/genes/LL+RS+CODAM+LLS_eqtls_genes_snpsInAllDatasets_18072014/eQTLsFDR0.05-ProbeLevel_genes_pr_noSec_topSnpPvalues.txt");
	}
}
