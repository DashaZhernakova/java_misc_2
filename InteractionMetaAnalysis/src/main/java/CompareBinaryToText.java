import org.molgenis.genotype.Allele;
import umcg.genetica.io.binInteraction.BinaryInteractionFile;
import umcg.genetica.io.binInteraction.BinaryInteractionFileException;
import umcg.genetica.io.binInteraction.BinaryInteractionQtlZscores;
import umcg.genetica.io.binInteraction.BinaryInteractionZscores;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariant;
import umcg.genetica.io.binInteraction.variant.BinaryInteractionVariantCreator;
import umcg.genetica.io.text.TextFile;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

import static junit.framework.Assert.assertEquals;

/**
 * Created by dashazhernakova on 30.01.15.
 */
public class CompareBinaryToText {

	public void compare(String binaryFname, String txtFname, String snpStatsFname) throws IOException, BinaryInteractionFileException {

		BinaryInteractionFile loadedInteractions = BinaryInteractionFile.load(new File(binaryFname));
		TextFile txtFile = new TextFile(txtFname, false);

		String[] els = txtFile.readLineElems(TextFile.tab);

		String snp, gene, covariate;
		int zsnp, zcov, zint, zmain, zintflip, numSamples;

		Allele [] alleles = new Allele[2];

		HashMap<String, BinaryInteractionVariant> snpStats = getSNPstats(snpStatsFname);

		while ((els = txtFile.readLineElems(TextFile.tab)) != null){
			snp = els[0];
			gene = els[1];
			covariate = els[2];
			BinaryInteractionVariant variant = snpStats.get(snp);
			String id = snp + " " + gene + " " + covariate;
			if (Math.abs(loadedInteractions.readQtlResults(snp, gene).getZscores()[0] - Double.parseDouble(els[6])) > 1e-10)
				System.out.println(id + ": different QTL z-scores: " + loadedInteractions.readQtlResults(snp, gene).getZscores()[0] + " and " + els[6]);
			BinaryInteractionZscores interRes = loadedInteractions.readInteractionResults(snp, gene, covariate);
			if (Math.abs(interRes.getZscoreInteractionCohort()[0] - Double.parseDouble(els[5])) > 1e-10)
				System.out.println(id + ": different interaction z-scores: " + interRes.getZscoreInteractionCohort()[0] + " and " + els[5]);
			if (Math.abs(interRes.getZscoreCovariateCohort()[0] - Double.parseDouble(els[4])) > 1e-10)
				System.out.println(id + ": different covariate z-scores: " + interRes.getZscoreCovariateCohort()[0] + " and " + els[4]);
			if (Math.abs(interRes.getZscoreSnpCohort()[0] - Double.parseDouble(els[3])) > 1e-10)
				System.out.println(id + ": different SNP z-scores: " + interRes.getZscoreSnpCohort()[0] + " and " + els[3]);

			BinaryInteractionVariant loadedVariant = loadedInteractions.getVariant(snp);

			if (! loadedVariant.getName().equals(variant.getName()))
				System.out.println(id + "\t" + loadedVariant.getName() + " " + variant.getName());
			if (! loadedVariant.getChr().equals(variant.getChr()))
				System.out.println(id + "\t" + loadedVariant.getChr() + " " + variant.getChr());
			if (loadedVariant.getPos() != (variant.getPos()))
				System.out.println(id + "\t" + loadedVariant.getPos() + " " + variant.getPos());
			if (! loadedVariant.getRefAllele().equals(variant.getRefAllele()))
				System.out.println(id + " ref allele\t" + loadedVariant.getRefAllele() + " " + variant.getRefAllele());
			if (! loadedVariant.getAltAllele().equals(variant.getAltAllele()))
				System.out.println(id + " alt allele\t" + loadedVariant.getAltAllele() + " " + variant.getAltAllele());
		}

	}

	private Allele[] getAlleles(String[] els){
		Allele minor = Allele.create(els[4]);
		Allele major;
		String[] strAlleles = els[3].split("/");
		if (strAlleles[0] == els[4])
			major = Allele.create(strAlleles[1]);
		else
			major = Allele.create(strAlleles[0]);
		return new Allele[] {major, minor};
	}

	private HashMap<String, BinaryInteractionVariant> getSNPstats(String fname) throws IOException {
		TextFile file = new TextFile(fname, false);
		String[] els = file.readLineElems(TextFile.tab);
		HashMap<String, BinaryInteractionVariant> snpStats = new HashMap<String, BinaryInteractionVariant>();

		while ((els = file.readLineElems(TextFile.tab)) != null){
			Allele[] alleles = getAlleles(els);
			BinaryInteractionVariant variant = new BinaryInteractionVariantCreator(els[0], els[1], Integer.parseInt(els[2]), alleles[0], alleles[1]);
			snpStats.put(els[0], variant);
		}
		file.close();
		return snpStats;
	}


}
