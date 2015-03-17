import umcg.genetica.io.binInteraction.BinaryInteractionFileException;

import java.io.IOException;

/**
 * Created by dashazhernakova on 30.01.15.
 */
public class Main {
	public static void usage(){

	}
	public static void main(String[] args) throws IOException, BinaryInteractionFileException {

		String arg, val, in = null, mode = null, out = null, in2 = null;

		int i = 0;
		for (i = 0; i < args.length; i++) {
			arg = args[i];
			val = null;

			if (i + 1 < args.length) {
				val = args[i + 1];
			}

			if (arg.equals("--mode")) {
				mode = val;
				//System.out.println("mode");
				break;
			}

		}
		if (mode == null) {
			System.out.println("ERROR: Please supply --mode");
			usage();
		}
		else if (mode.equals("meta")){
			String[] paths = new String[args.length - i - 2];
			int idx = 0;
			for (int j = i+2; j < args.length; j++){
				paths[idx] = args[j];
				idx++;
			}
			System.out.println("Input paths:");
			for (String p : paths)
				System.out.println(p);
			MetaAnalyser_synch analyser = new MetaAnalyser_synch(paths.length - 1);
			analyser.runMetaAnalysis(paths);
		}
		else if (mode.equals("compare")){
			String txt = null, bin = null, snp = null;
			for (int j = i; j < args.length; j++){
				arg = args[j];
				val = null;

				if (j + 1 < args.length) {
					val = args[j + 1];
				}
				if (arg.equals("--text"))
					txt = val;
				if (arg.equals("--binary"))
					bin = val;
				if (arg.equals("--snpQC"))
					snp = val;

			}
			CompareBinaryToText comp = new CompareBinaryToText();
			comp.compare(bin, txt, snp);
		}
	}
}
