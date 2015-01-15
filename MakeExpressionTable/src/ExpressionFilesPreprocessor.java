import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class ExpressionFilesPreprocessor {
	private ArrayList<String> fNames;

	public ExpressionFilesPreprocessor(String fileList) throws IOException {
		ProcessTranscriptCounts processor = new ProcessTranscriptCounts();
		processor.getFileNames(fileList);
		fNames = processor.outputFnames();
	}

	public ExpressionFilesPreprocessor(String dir, String pattern) throws IOException {
		ProcessTranscriptCounts processor = new ProcessTranscriptCounts();
		processor.getFileNames(dir, pattern);
		fNames = processor.outputFnames();
	}
	public ExpressionFilesPreprocessor(ArrayList<String> fileNames){
		fNames = fileNames;
	}

	public void separateMultiGenesPerLine(){
		for (String fName : fNames){
			moveFile(fName, fName + ".orig");

			try {
				TextFile in = new TextFile(fName + ".orig", false);
				TextFile out = new TextFile(fName, true);
				String[] els;
				while ((els = in.readLineElems(TextFile.tab)) != null){
					String[] genes = els[0].split(",");
					for (String g : genes){
						out.writeln(g + "\t" + els[1]);
					}
				}

				in.close();
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	private void moveFile(String fName, String newFName){
		File f = new File(fName);
		File newFile = new File(newFName);
		f.renameTo(newFile);
	}
}
