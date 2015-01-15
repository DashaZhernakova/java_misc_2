import java.io.IOException;

/**
 *
 * @author dashazhernakova
 */
public class MakeExpressionTable {

   	public static void usage(){
        System.out.println("Options:");
        System.out.println("--in (Input folder. Should contain folders with alignment files and read counts files for each sample). Required together with \"pattern\" if \"fileList\" is not specified.\n  "
                + "--annot (Annotation file in gtf format or in the following format: Platform\tProbeName(as in readCounts file)\tGeneName\tProbeChr\tProbe start\t Probe end). Required.\n  "
                + "--feature (Feature type to use when getting transcripts/genes from annotation. In case of gtf annotation use feature types used in attributes field) Default: transcript_id\n  "
				+ "--out (Path to output expression table file). Required.\n  "
                + "--pattern (Filename of the file with read counts). Required if \"in\" is specified. \n  "
                + "--alnFnamePattern (Filename of the file with mapped reads or of the idxstats file). Default: accepted_hits.filtered.sam\n  "
                + "--normalize (true if you want the expression values normalized by the total number of mapped reads per sample). Default: false\n  "
                + "--normalizeByLength (true if you want the expression values normalized by the gene length). Default: false\n  "
                + "--fileList (File with sample names and their file paths to process). Required if \"in\" and \"pattern\" are not specified");
    }
	public static void printArgs(String inDir, String pattern, String fileList, String annotation, String feature, String outFile, String alnPattern, boolean normalize, boolean normLen, boolean skipMissing){
		System.out.println("\nInput arguments: " +
				"\nin: " + inDir +
				"\npattern: " + pattern +
				"\nfileList: " + fileList +
				"\nannot: " + annotation +
				"\nfeature: " + feature +
				"\nout: " + outFile +
				"\nnormalize: " + normalize +
				"\nnormalizeByLength: " + normLen +
				"\nalnFnamePattern: " + alnPattern +
				"\nignoreMissing: " + skipMissing +
				"\n"
		);
	}
    public static void main(String[] args) throws IOException {
        
        ProcessTranscriptCounts p;
        String annot = null, pattern = null, dir = null, outFile = null, alnPattern = null, fileList = null, feature = null;
		boolean normalize = false, normLen = false, rpkm = false, skipmissing = false, preprocess = false;
		for (int j = 0; j < args.length; j++) {
			String arg = args[j];
			String val = null;

			if (j + 1 < args.length) {
				val = args[j + 1];
			}

			if (arg.equals("--annot")) {
				annot = val;
			}
			if (arg.equals("--feature")) {
				feature = val;
			}
			if (arg.equals("--in")) {
				dir = val;
			}
			if (arg.equals("--out")) {
				outFile = val;
			}
			if (arg.equals("--pattern")) {
				pattern = val;
			}
			if (arg.equals("--alnFnamePattern")) {
				alnPattern = val;
			}
			if (arg.equals("--fileList")) {
				fileList = val;
			}

			if (arg.equals("--normalize")) {
				normalize = Boolean.valueOf(val);
			}
			if (arg.equals("--normalizeByLength")){
				normLen = Boolean.valueOf(val);
			}
			if (arg.equals("--ignoreMissing")){
				skipmissing = Boolean.valueOf(val);
			}
			if (arg.equals("--preprocess")) {
				preprocess = Boolean.valueOf(val);
			}
		}
		if ((annot == null) || (outFile == null)){
			System.out.println("Not enough arguments: annotation file or output file is missing");
			printArgs(dir, pattern, fileList, annot, feature, outFile, alnPattern, normalize, normLen, skipmissing);
			System.exit(-1);
		}

		printArgs(dir, pattern, fileList, annot, feature, outFile, alnPattern, normalize, normLen, skipmissing);


		if (fileList == null){
			if (((dir == null) || (pattern == null))){
				System.out.println("Folder and/or pattern are missing!");
				printArgs(dir, pattern, fileList, annot, feature, outFile, alnPattern, normalize, normLen, skipmissing);
				System.exit(-1);
			}
			p = new ProcessTranscriptCounts(dir, pattern, annot, feature, alnPattern, normalize, normLen, skipmissing);
		}
		else{
			p = new ProcessTranscriptCounts(fileList, annot, feature, normalize, normLen, skipmissing);
		}

		p.preprocess = preprocess;

		p.run(outFile);

	}
}
