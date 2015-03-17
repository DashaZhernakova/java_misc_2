import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;

/**
 * Created by dashazhernakova on 13.03.15.
 */
public class Make2DMatrixFrom1DFile {
	DoubleMatrixDataset matrixDataset;


	public void make2DMatrix(String fname, int colNamesColumn, int rowNamesColumn, int valuesColumn, boolean header, String outFname) throws IOException {
		System.out.println("input: " + fname + "\noutoput: " + outFname + "\ncolumn names column: " + colNamesColumn + "\nrow names column: " + rowNamesColumn + "\ncolumn containing values: " + valuesColumn + "\nheader: " + header);
		// initialize the matrix, fill column and row names
		initMatrix(fname, colNamesColumn, rowNamesColumn, header);

		//fill the cell values
		readMatrix(fname, colNamesColumn, rowNamesColumn, valuesColumn, header);

		// write the matrix to file
		matrixDataset.save(outFname);
	}

	private void readMatrix(String fname, int colNamesColumn, int rowNamesColumn, int valuesColumn, boolean header) throws IOException {
		TextFile file = new TextFile(fname, false);
		String[] els;

		if (header){
			els = file.readLineElems(TextFile.tab);
		}
		LinkedHashMap<String, Integer> hashCols = matrixDataset.getHashCols();
		LinkedHashMap<String, Integer> hashRows = matrixDataset.getHashRows();

		DoubleMatrix2D matrix;
		if ((matrixDataset.rows() * (long) matrixDataset.columns()) < (Integer.MAX_VALUE - 2)) {
			matrix = new DenseDoubleMatrix2D(matrixDataset.rows(), matrixDataset.columns());
		} else {
			matrix = new DenseLargeDoubleMatrix2D(matrixDataset.rows(), matrixDataset.columns());
		}

		String colName, rowName;
		int colNum = -9, rowNum = -9;
		while ((els = file.readLineElems(TextFile.tab)) != null){
			colName = els[colNamesColumn];
			rowName = els[rowNamesColumn];
			colNum = hashCols.get(colName);
			rowNum = hashRows.get(rowName);

			matrix.setQuick(rowNum, colNum, Double.parseDouble(els[valuesColumn]));
		}

		matrixDataset.setMatrix(matrix);

		file.close();
	}

	private void initMatrix(String fname, int colNamesColumn, int rowNamesColumn, boolean header) throws IOException {
		TextFile file = new TextFile(fname, false);

		LinkedHashSet<String> colNames = new LinkedHashSet<String>();
		LinkedHashSet<String> rowNames = new LinkedHashSet<String>();

		String[] els;
		if (header){
			els = file.readLineElems(TextFile.tab);
		}

		while ((els = file.readLineElems(TextFile.tab)) != null){
			colNames.add(els[colNamesColumn]);
			rowNames.add(els[rowNamesColumn]);
		}

		file.close();

		matrixDataset = new DoubleMatrixDataset(new ArrayList<String>(rowNames), new ArrayList<String>(colNames));
	}
	public static void main(String[] args) throws IOException {
		Make2DMatrixFrom1DFile m = new Make2DMatrixFrom1DFile();
		m.make2DMatrix("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/interactionWithTFs/allGenes/noZtransform/covariates_interactionResults_ids.txt.gz",
				3,0,44, false, "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/interactionWithTFs/allGenes/noZtransform/covariates_interactionResults_ids_2Dtable.txt.gz");
	}
}
