import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by dashazhernakova on 12.02.15.
 */
public class Main {
	DoubleMatrixDataset expressionTable;
	ArrayList<PolyABin> bins;
	public Main(String exprTableFname, String converstionFname) throws IOException {
		expressionTable = new DoubleMatrixDataset(exprTableFname);
		bins = convertBinsToRowNumbers(converstionFname);
	}
	private  ArrayList<PolyABin> convertBinsToRowNumbers(String fname) throws IOException {
		TextFile tagToBin = new TextFile(fname, false);
		HashMap<String, PolyABin> binsToRows = new HashMap<String, PolyABin>();

		String[] els;
		while( (els = tagToBin.readLineElems(TextFile.tab)) != null){
			PolyABin bin = binsToRows.get(els[9]);
			if (bin == null){
				bin = new PolyABin();
				bin.chr = els[6];
				bin.start = Integer.parseInt(els[7]);
				bin.end = Integer.parseInt(els[8]);
				bin.name = els[9];
				bin.strand = els[11].charAt(0);

			}

			String tag = els[3];
			if (expressionTable.hashRows.containsKey(tag)){
				int row = (Integer) expressionTable.hashRows.get(tag);
				bin.tagRows.add(row);
			}
			binsToRows.put(bin.name, bin);
		}
		tagToBin.close();
		ArrayList result = new ArrayList<PolyABin> ();
		result.addAll(binsToRows.values());
		return result;
	}
	public void run(String outPath) throws IOException {
		double [][] rawData = expressionTable.getRawData();
		TextFile out;
		for (int sample = 0; sample < expressionTable.nrCols; sample ++){
			String sampleId = (String) expressionTable.colObjects.get(sample);
			out = new TextFile(outPath + "/" + sampleId + ".polya_bins.count", true);

			for (PolyABin bin : bins){
				double count = 0;
				for (int row : bin.tagRows){
					count += rawData[row][sample];
				}
				out.writeln(bin.chr + "\t" + bin.start + "\t" + bin.end + "\t" + (int) count + "\t" + (int) count + "\t" + bin.strand + "\t" + bin.name);
			}
			out.close();
		}
	}
	public static void main(String[] args) throws IOException {
		Main m = new Main("/Users/dashazhernakova/Documents/UMCG/deepSAGE/new/deepSAGE_tag/tagwise_expression_table_SNP_in_recognition_sequence_tags_excluded.txt.gz",
				"/Users/dashazhernakova/Documents/UMCG/deepSAGE/replicationBBMRI/tag_to_polyA_bin_multipleTags.txt");

		m.run("/Users/dashazhernakova/Documents/UMCG/deepSAGE/replicationBBMRI/counts/");
	}
}
