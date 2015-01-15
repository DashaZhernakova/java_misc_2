import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;


public class TranscriptRate {

	private float[][] table;
	private ArrayList<String> allTranscr;
	private ArrayList<String> sampleNames;
	private ArrayList<String> fNames;
	private boolean[] featureExpressed;
	private int trSize;
	private int samplesSize;
	private HashMap<String,Float> geneCounts;
	private HashMap<String,String> transcrToGene;



	TranscriptRate(String dir, String fname, String annot, String featureType) throws IOException {
		getAllTranscripts(annot, featureType);
		getFileNames(dir, fname);

	}


	TranscriptRate(String fileListPath, String annot, String featureType) throws IOException {
		getAllTranscripts(annot, featureType);
		getFileNames(fileListPath);

	}



	/**
	 * Main function
	 * @throws java.io.IOException
	 */
	public void run(String outFile) throws IOException {

		//all transcripts from ensembl, sorted
		Collections.sort(allTranscr);
		trSize = allTranscr.size(); //number of transcripts

		//make expression table
		makeTable();

		//print the result
		printTable(outFile);

	}

	/**
	 * Checks if file is empty or doesn't exist
	 * @param fName
	 * @return true if file is empty or doesn't exist
	 */
	private boolean fileIsEmpty(String fName) {
		File file = new File(fName);
		if (file.length() == 0)
			return true;
		return false;

	}


	/**
	 * Makes expression table with feature counts.
	 * @throws java.io.IOException
	 */
	public void makeTable() throws IOException{

		TextFile in = null;
		String line = "";
		String[] splLine;

		int index = 0;

		//table with counts
		table = new float[trSize][samplesSize];

		//for each transcript, true if counts > 0 for any of the samples
		featureExpressed = new boolean[trSize];

		transcrToGene = new HashMap<String, String>();
		int curId = 0;
		System.out.println("Started generating the expression table");
		String trId;
		float count = 0;
		for (String fName : fNames){
			int rCount=0;
			in = new TextFile(fName, false);

			System.out.println("\nProcessing file: " + fName);

			//get number of mapped reads from filtered sam file

			geneCounts = new HashMap<String, Float>();
			//transcrToGene = new HashMap<String, String>();

			while ((line = in.readLine()) != null){
				//read the expression from simple txt or from gtf
				count = 0;
				try{

					splLine = readExpressionLine(fName, line);
					count = Float.valueOf(splLine[1]);



					//skip not expressed transcripts
					if (count == 0)
						continue;

					trId = splLine[0];

					if (trId != null){
						//get the transcript position in all transcripts list
						index = Collections.binarySearch(allTranscr, trId);

						if (index >= 0){ //if transcript found
							featureExpressed[index] = true; //set this feature as expressed
							rCount+=count;
							table[index][curId] += count;

						}
					}
				} catch (Exception e){
					System.out.println("Line: " + line);
				}
			}
			convertToFractions(curId);
			curId++;
			System.out.println("overall number of reads=" + rCount);

			if (rCount == 0){
				System.out.println("WARNING: No reads mapping into specified interval type! Check if the feature names in read count files per sample correspond to features specified in the annotation!");
			}
			in.close();

		}

	}

	private void convertToFractions(int sampleId){
		String transcr, gene;
		float geneCount = 0, trCount;
		for (int trId = 0; trId < trSize; trId++){
			trCount = table[trId][sampleId];
			if (trCount != 0){
				transcr = allTranscr.get(trId);
				if ((transcr.equals("ENST00000367771")) || (transcr.equals("ENST00000367771"))){
					System.out.println("ENSG00000000457");
				}
				gene = transcrToGene.get(transcr);
				geneCount = geneCounts.get(gene);

				table[trId][sampleId] = trCount/geneCount;
			}
		}
	}
	/**
	 * Reads an expression line from a per-sample expression file. Supports gtf files from Flux Capacitor and files of "transcript\tcount" format (or transcript<space>count format)
	 * @param fName - file name to determine the format
	 * @param line - expression line
	 * @return String[feature_id, count]
	 */
	private String[] readExpressionLine(String fName, String line){
		String[] spl;
		String count, gene;
		if (fName.replace(".gz", "").endsWith(".gtf")){
			spl = line.split("\"");

			count = spl[6].split(" ")[2].replace(";", "");
			if (count.equals("0.000000"))
				return new String[]{spl[1], count};
			/*if (spl[5].equals("ENSG00000000457")){
				System.out.println("ENSG00000000457");
			}*/
			gene = spl[5];
			transcrToGene.put(spl[1], gene);
			Float geneCount = geneCounts.get(gene);
			if (geneCount == null){
				geneCounts.put(gene, Float.valueOf(count));
			}
			else{
				geneCounts.put(gene, Float.valueOf(count) + geneCount);
			}

			return new String[]{spl[1], count};
		}
		//System.out.println(geneCounts.get("ENSG00000000457"));
		if (! line.contains("\t"))
			return line.split(" ");
		return line.split("\t");

	}

	/**
	 * gets all transcripts from Ensembl. If normalizeByLength = true, then trLengths is filled by transcript/gene lengths
	 * @return
	 * @throws java.io.IOException
	 */
	public void getAllTranscripts(String annotationFname, String feature) throws IOException{
		Annotation annot = new Annotation();
		System.out.println("file " + annotationFname);
		System.out.println("feature " + feature);
		if (feature != null){
			annot.getAllTranscr(annotationFname, feature);
		}
		else{
			annot.getAllTranscr(annotationFname);
		}
		allTranscr = annot.allFeatures;

	}


	/**
	 * Gets file names ending with fnamePattern in dir dirName. Or from the file list, if it is provided. Also gets sample names.
	 * @throws IOException
	 */
	public void getFileNames(String fileList) throws IOException {
		fNames = new ArrayList<String>();
		sampleNames = new ArrayList<String>();

		if (fileList != null){
			System.out.println("Started processing files from file list: " + fileList);
			TextFile flist = new TextFile(fileList, false);
			String[] els;
			while ((els = flist.readLineElems(TextFile.tab)) != null){
				if (! fileIsEmpty(els[1])){
					if (els.length > 1){
						fNames.add(els[1]);
						sampleNames.add(els[0]);
					}
					else{
						fNames.add(els[0]);
						sampleNames.add(els[0].split("/")[els[0].split("/").length - 1]);
					}
				}
			}
			flist.close();
		}

		samplesSize = fNames.size();

	}

	/**
	 *
	 * @param dirName - folder containing sample folders with read counts files
	 * @param fnamePattern - name of the files with counts
	 */
	public void getFileNames(String dirName, String fnamePattern){
		System.out.println("Getting all files with reads per transcripts counts. \n\tFolder name: " + dirName + "\n\tFiles names finish with " + fnamePattern);
		File dir = new File(dirName);
		fNames = new ArrayList<String>();
		sampleNames = new ArrayList<String>();
		for (File ch : dir.listFiles()) {
			if (ch.isDirectory())
				for (File child : ch.listFiles()){
					if (child.getName().equals(fnamePattern)){
						String fName = child.getPath();
						if (! fileIsEmpty(fName)){
							fNames.add(fName);
							String[] splName = fName.split("/");
							sampleNames.add(splName[splName.length - 2]);
						}
					}
				}
		}

		samplesSize = fNames.size();
	}

	/**
	 * Checks if average feature expression is higher than threshold
	 * @param thres
	 * @param values
	 * @return
	 */
	private boolean avgHigherThanThreshold(float thres, float[] values){
		float sum = 0;
		for (float count : values){
			sum+=count;
		}
		if (sum/values.length > thres)
			return true;
		return false;
	}




	/**
	 * Prints the expression table (not normalized by length)
	 * @throws java.io.IOException
	 */
	public void printTable(String outFile) throws IOException{
		System.out.println("Printing the resulting expression table to " + outFile);

		TextFile outReads = new TextFile(outFile, true);
		float reads=0, fraction = 0;
		String transcr, gene;
		outReads.write("probe\t");
		outReads.writelnTabDelimited(sampleNames.toArray());

		DecimalFormat df = new DecimalFormat("0.#####");
		df.setRoundingMode(RoundingMode.FLOOR);
		DecimalFormatSymbols custom = new DecimalFormatSymbols();
		custom.setDecimalSeparator('.');
		df.setDecimalFormatSymbols(custom);

		for (int tr = 0; tr < trSize; tr++){
			//if ( (featureExpressed[tr]) && (avgHigherThanThreshold(0, table[tr])) ){
			if (featureExpressed[tr]){
				transcr = allTranscr.get(tr);
				gene = transcrToGene.get(transcr);
				if (transcr.equals("ENST00000485017"))
					System.out.println("ENST00000485017");
				outReads.write(transcr + "\t" + gene);
				for (int sam = 0; sam < samplesSize; sam++){

					reads = table[tr][sam];
					outReads.write("\t" + df.format(reads));

				}
				//System.out.println(transcr);
				outReads.writeln();
			}
		}
		outReads.close();
	}






	public static void main(String[] args) throws IOException {
		TranscriptRate p = new TranscriptRate("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/mappedData/compare/fileList.txt",
				"/Users/dashazhernakova/Documents/UMCG/hg19/v71/genomicIntervals/gtf/Homo_sapiens.GRCh37.71.cut.sorted.gtf.gz", null
				);

		p.run("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/mappedData/compare/expression_table.txt");
		//p.readCounts();
		// }
		//else
		//  System.out.println("not enough arguments given");

	}
}


