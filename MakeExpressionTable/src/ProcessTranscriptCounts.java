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

/**
 * processes files of the type:
 * trId\tcount
 * trId\tcount
 * or gtf files produced by Flux Capacitor
 *
 * generates a table of ids vs trCounts
 * ...
 * @author dashazhernakova
 */
public class ProcessTranscriptCounts {

	private String alnFnamePattern;
	private boolean normalize;
    private boolean ignoreMissing;
	private boolean normalizeByLength;
	public boolean preprocess;

    private float[][] table;
	private ArrayList<String> allTranscr;
	private ArrayList<String> sampleNames;
	private ArrayList<String> fNames;
	private float[] readsPerSample;
	private boolean[] featureExpressed;
	private int trSize;
	private int samplesSize;
	private HashMap<String,Integer> trLengths;



	ProcessTranscriptCounts(){};

    ProcessTranscriptCounts(String dir, String fname, String annot, String featureType, boolean norm) throws IOException {
        getAllTranscripts(annot, featureType);
        getFileNames(dir, fname);
		normalize = norm;
	}

    ProcessTranscriptCounts(String dir, String fname, String annot, String featureType, String aln, boolean norm) throws IOException {
        this(dir, fname, annot, featureType, norm);
        alnFnamePattern = aln;
    }

    ProcessTranscriptCounts(String fileListPath, String annot, String featureType) throws IOException {
        getAllTranscripts(annot, featureType);
		getFileNames(fileListPath);
	}

	ProcessTranscriptCounts(String fileListPath, String annot, String featureType, boolean norm, boolean lenNorm, boolean skipMissing) throws IOException {
		//this(fileListPath, annot);
		normalizeByLength = lenNorm;
		normalize = norm;
		ignoreMissing = skipMissing;
		getAllTranscripts(annot, featureType);
		getFileNames(fileListPath);
		System.out.println("ignoreMissing=" + ignoreMissing);
	}

	ProcessTranscriptCounts(String dir, String fname, String annot, String featureType, String aln, boolean norm, boolean lenNorm, boolean skipMissing) throws IOException {
        //this(dir, fname, annot, aln, norm);
		normalize = norm;
		alnFnamePattern = aln;
        normalizeByLength = lenNorm;
		ignoreMissing = skipMissing;
		getAllTranscripts(annot, featureType);
		getFileNames(dir, fname);

	}

	public ArrayList<String> outputFnames(){
		return fNames;
	}
    /**
     * Main function
     * @throws java.io.IOException
     */
    public void run(String outFile) throws IOException {

        //all transcripts from ensembl, sorted
        Collections.sort(allTranscr);
        trSize = allTranscr.size(); //number of transcripts

		//preprocess the expression files before running (currently: separate multiple genes per line)
        if (preprocess){
			ExpressionFilesPreprocessor preprocessor = new ExpressionFilesPreprocessor(fNames);
			preprocessor.separateMultiGenesPerLine();
		}

        //make expression table
        makeTable();

		//print the result
        if (! normalizeByLength)
            printTable(outFile);
        else
            printTableNormByLength(outFile );
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
        readsPerSample = new float[samplesSize];

        //for each transcript, true if counts > 0 for any of the samples
        featureExpressed = new boolean[trSize];


        int curId = 0;
        System.out.println("Started generating the expression table");
        String trId;
		float count = 0;
        for (String fName : fNames){
            int rCount=0;
            in = new TextFile(fName, false);

            System.out.println("\nProcessing file: " + fName);

            //get number of mapped reads from filtered sam file
            if (normalize)
                readsPerSample[curId] = getNumMappedReads(fName);

            while ((line = in.readLine()) != null){
                //read the expression from simple txt or from gtf
                count = 0;
				try{
					//check if 2 genes per line (separated by comma) & separate

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
            curId++;
            System.out.println("overall number of reads=" + rCount);
            System.out.println("overall number of mapped reads=" + readsPerSample[curId - 1]);

			if (rCount == 0){
				System.out.println("WARNING: No reads mapping into specified interval type! Check if the feature names in read count files per sample correspond to features specified in the annotation!");
			}
            in.close();

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
        String count;
        if ((fName.replace(".gz", "").endsWith(".gtf")) || (line.split("\t").length == 9)){ //TODO: change the check
            spl = line.split("\"");
            count = spl[6].split(" ")[2].replace(";", "");
            return new String[]{spl[1], count};
        }
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
		System.out.println("normalizeByLength " + normalizeByLength);
		if (feature != null){
			annot.getAllTranscr(annotationFname, feature, normalizeByLength);
		}
		else{
			annot.getAllTranscr(annotationFname);
		}
		allTranscr = annot.allFeatures;
		if (normalizeByLength)
			trLengths = annot.lengths;
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
			if (! ((ignoreMissing) && (fileIsEmpty(els[1])))){
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
						if (! ((ignoreMissing) && (fileIsEmpty(fName)))){
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
     * gets number of mapped reads from a sam file or from an idxstats file
     * @param path - path to a folder with mapped reads or with an idxstats file
     * @return number of mapped reads
     * @throws java.io.IOException
     */
    public int getNumMappedReads(String path){
        int numReads = 0;

		StringBuilder strbuilder = new StringBuilder(path);
		String alnFname = strbuilder.replace(path.lastIndexOf("/"), path.length(), "/").toString();

		try {
            if (alnFnamePattern.endsWith("idxstats")){
                TextFile stats = new TextFile(alnFname + alnFnamePattern, false);
                String[] els;

                while ((els = stats.readLineElems(TextFile.tab)) != null){
                    numReads += Integer.parseInt(els[2]);
                }
                stats.close();
            }
            else{
                TextFile filtSam = null;
                if (alnFnamePattern.isEmpty())
                    filtSam = new TextFile(alnFname + "accepted_hits.filtered.sam", false);
                else if (alnFnamePattern.endsWith(".sam")){
                    filtSam = new TextFile(alnFname + alnFnamePattern, false);
                }
                String line ="";

                while ( (line = filtSam.readLine()) != null){
                    if (! line.startsWith("@")){
                        numReads++;
                    }
                }
                filtSam.close();
            }
            return numReads;
        } catch (IOException ex) {
            Logger.getLogger(ProcessTranscriptCounts.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println("No alignment file or summary file (" + alnFnamePattern + ") for " + path);
        }
        return numReads;
    }

	/**
	 * Prints the expression table (not normalized by length)
	 * @throws java.io.IOException
	 */
	public void printTable(String outFile) throws IOException{
		System.out.println("Printing the resulting expression table to " + outFile);
		System.out.println("Normalize=" + normalize);

		TextFile outReads = new TextFile(outFile, true);
		float reads=0;

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
				outReads.write(allTranscr.get(tr));
				for (int sam = 0; sam < samplesSize; sam++){

					if (normalize){ //normalise by total number of mapped reads
						reads = (table[tr][sam] * 1000000) / readsPerSample[sam];
						outReads.write("\t" + (df.format(reads)));
					}
					else{ //print raw counts
						reads = table[tr][sam];
						outReads.write("\t" + df.format(reads));
					}
				}
				outReads.writeln();
			}
		}
		outReads.close();
	}

	/**
	 * Prints expression table, normalising read counts by the genomic feature length
	 * @throws java.io.IOException
	 */
	public void printTableNormByLength(String outFile) throws IOException{
		System.out.println("Printing the resulting expression table to " + outFile);
		System.out.println("Normalize=" + normalize);

		TextFile outReads = new TextFile(outFile, true);
		String tr_id;
		float reads=0;

		outReads.write("probe\t");
		outReads.writelnTabDelimited(sampleNames.toArray());

		DecimalFormat df = new DecimalFormat("0.#####");
		df.setRoundingMode(RoundingMode.FLOOR);
		DecimalFormatSymbols custom = new DecimalFormatSymbols();
		custom.setDecimalSeparator('.');
		df.setDecimalFormatSymbols(custom);

		for (int tr = 0; tr < trSize; tr++){
			//if ( (featureExpressed[tr]) && (avgHigherThanThreshold(1, table[tr])) ){
			if (featureExpressed[tr]){
				tr_id = allTranscr.get(tr);
				outReads.write(tr_id);
				if (normalize){
					for (int sam = 0; sam < samplesSize; sam++){
						reads = (table[tr][sam] * 1000000) / (readsPerSample[sam] * trLengths.get(tr_id));
						outReads.write("\t" + (df.format(reads)));

					}
				}
				else{
					for (int sam = 0; sam < samplesSize; sam++){
						reads = (table[tr][sam] * 1000000) / trLengths.get(tr_id);
						outReads.write("\t" + (df.format(reads)));
					}
				}
				outReads.writeln();
			}
		}
		outReads.close();
	}




	public static void main(String[] args) throws IOException {
        ProcessTranscriptCounts p = new ProcessTranscriptCounts("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/reseq/fileList_tr.txt",
                "/Users/dashazhernakova/Documents/UMCG/hg19/v71/annotations/annotation_transcr_v71.txt", null,
                false, false, true);

        p.run("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/reseq/expression_table_reseq_transcr.txt");
        //p.readCounts();
       // }
        //else
          //  System.out.println("not enough arguments given");
        
    }
}

