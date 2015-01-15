/**
 * Created by dashazhernakova on 14.01.15.
 */
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.MatrixHandling;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;


import java.io.IOException;
import java.util.HashSet;

/**
 * Created by dashazhernakova on 14.01.15.
 */
public class CorrelatePerSampleExpressionWithScoresColumn {
	DoubleMatrixDataset<String,String> expression;
	DoubleMatrixDataset<String,String> scoresVector;


	public CorrelatePerSampleExpressionWithScoresColumn(String expressionFname, String scoresColumnFname) throws IOException {
		expression = DoubleMatrixDataset.loadDoubleTextData(expressionFname, "\t");
		scoresVector = DoubleMatrixDataset.loadDoubleTextData(scoresColumnFname, "\t");
	}


	private void reorderVector(){

	}

	private void correlate(String outFname) throws IOException {

		TextFile out = new TextFile(outFname, true);

		DoubleMatrixDataset<String,String> cut_expression = MatrixHandling.CreatSubsetBasedOnRows(expression, new HashSet<String>(scoresVector.getRowObjects()), false);
		expression = null;
		DoubleMatrixDataset<String,String> cut_scoresVector = MatrixHandling.CreatSubsetBasedOnRows(scoresVector, new HashSet<String>(cut_expression.getRowObjects()), false);
		scoresVector = null;
		System.out.println("Number of shared genes: " + cut_scoresVector.rows());
		cut_scoresVector.reorderRows(cut_expression.getHashRows());
		double [] scoresVectorArray = cut_scoresVector.getMatrix().viewColumn(0).toArray();
		for (String sample : cut_expression.getColObjects()){
			out.write("\t" + sample);
		}
		out.writeln();
		out.write("Correlation_public_PC2");

		for (int i = 0; i < cut_expression.columns(); i++ ){
			double cor = new SpearmansCorrelation().correlation(cut_expression.getMatrix().viewColumn(i).toArray(), scoresVectorArray);
			out.write("\t" + cor);
		}
		out.writeln();
		out.close();
	}
	public static void main(String[] args) throws IOException {
		/*CorrelatePerSampleExpressionWithScoresColumn c = new CorrelatePerSampleExpressionWithScoresColumn(
				"/Users/dashazhernakova/Documents/UMCG/data/BBMRI/first_run_final/expression/combined_gene_count_run_1.summed.passQC.TMM.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed_tmp.txt.gz",
				"/Users/dashazhernakova/Documents/UMCG/data/BBMRI/geuvadis_gene_fdr0.05_replication_in_bbmri_fdr0.05/publicRNAseq_PC2_scores/public_RNAseq_data_PC2_scores.txt");

		c.correlate();*/
		CorrelatePerSampleExpressionWithScoresColumn c = new CorrelatePerSampleExpressionWithScoresColumn(args[0], args[1]);
		c.correlate(args[2]);
	}
}
