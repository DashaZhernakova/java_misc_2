import java.io.IOException;

/**
 * Created by dashazhernakova on 13.03.15.
 */
public class Main {
	public static void usage(){
		System.out.println("Available modes:");
		System.out.println("Make2DMatrix (--in, --out, --col, --row, --val, --header)");
	}
	public static void main(String[] args) throws IOException {
		String arg, val, in = null, mode = null, out = null;

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
		else if (mode.equals("Make2DMatrix")){
			boolean header = false;
			int col, row, valCol;
			for (int j = i; j < args.length; j++){
				arg = args[j];
				val = null;

				if (j + 1 < args.length) {
					val = args[j + 1];
				}
				if (arg.equals("--in"))
					in = val;
				if (arg.equals("--col"))
					col = Integer.parseInt(val);
				if (arg.equals("--row"))
					row = Integer.parseInt(val);
				if (arg.equals("--value"))
					valCol = Integer.parseInt(val);
				if (arg.equals("--out"))
					out = val;
				if (arg.equals("--header"))
					header = true;
			}
		}
	}
}
