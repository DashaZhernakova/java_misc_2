import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class Annotation {
	public ArrayList<String> allFeatures;
	public HashMap<String, Integer> lengths;

	public ArrayList<String> getAllTranscr(String fName) throws IOException {
		getAllTranscr(fName, "transcript_id", false);
		return allFeatures;
	}

	public ArrayList<String> getAllTranscr(String fName, String feature) throws IOException {
		getAllTranscr(fName, feature, false);
		return allFeatures;
	}

	public ArrayList<String> getAllTranscr(String fName, String feature, boolean withLengths) throws IOException {
		System.out.println("Getting all features of the type: " + feature + " from file: " + fName);
		System.out.println("Getting feature lengths: " + withLengths);
		if (! formatOk(fName, feature)){
			System.out.println("Wrong annotation format! Exiting!");
			System.exit(-1);
		}

		if (fName.replace(".gz","").endsWith(".gtf")){
			getAllFeaturesFromGtf(fName, feature);
		}
		else if (fName.replace(".gz","").endsWith(".txt")){
			getAllFeaturesFromTxt(fName, feature, withLengths);
		}
		else{
			System.out.println("Check the annotation file format. The extension should be either .gtf[.gz] or .txt[.gz]");
			System.exit(-1);
		}
		return allFeatures;
	}
	private void getAllFeaturesFromGtf(String fname, String feature) throws IOException {
		HashSet<String> allFeaturesSet = new HashSet<String>();
		TextFile annot = new TextFile(fname, false);
		String[] els;
		HashMap<String, String> attrs = new HashMap<String, String>();

		while ((els = annot.readLineElems(TextFile.tab)) != null){
			attrs = splitAttrToMap(els[8]);
			allFeaturesSet.add(attrs.get(feature));
		}
		allFeatures = new ArrayList<String>(allFeaturesSet);
		System.out.println("Loaded " + allFeaturesSet.size() + " features of the type " + feature);

	}
	private void getAllFeaturesFromTxt(String fname, String feature, boolean withLengths) throws IOException {
		HashSet<String> allFeaturesSet = new HashSet<String>();
		TextFile annot = new TextFile(fname, false);
		String[] els = annot.readLineElems(TextFile.tab);
		int pos = 1;
		if (feature.startsWith("gene")){
			pos = 2;
			if (withLengths){
				System.out.println("ERROR! Using gene lengths not supported! If the annotation contains one gene per line, then run the command without feature type specified");
				System.exit(-1);
			}
		}
		if (! withLengths){
			while ((els = annot.readLineElems(TextFile.tab)) != null){
				allFeaturesSet.add(els[pos]);
			}
		}
		else{
			lengths = new HashMap<String, Integer>();
			while ((els = annot.readLineElems(TextFile.tab)) != null){

			allFeaturesSet.add(els[1]);
			int len = Integer.parseInt(els[5]) - Integer.parseInt(els[4]);
			lengths.put(els[1], len);

		}
		}
		allFeatures = new ArrayList<String>(allFeaturesSet);
		System.out.println("Loaded " + allFeaturesSet.size() + " features of the type " + feature);
	}

	private boolean formatOk(String fName, String feature) throws IOException {
		TextFile annot = new TextFile(fName, false);
		String[] els = annot.readLineElems(TextFile.tab);
		annot.close();
		if (fName.replace(".gz","").endsWith(".gtf")){
			if ((els.length == 9) && (els[8].contains(feature)))
				return true;
		}
		if (fName.replace(".gz","").endsWith(".txt")){
			if (els.length == 6)
				return true;
		}
		System.out.println("The annotation format is wrong! Please use gtf file or a txt file used for eQTL mapping");
		return false;
	}
	private HashMap<String, String> splitAttrToMap(String in) {
		String[] attrs = in.split("; ");
		HashMap<String, String> map = new HashMap<String, String>();
		for (String attr : attrs){
			String key = attr.split(" ")[0];
			String value = attr.split(" ")[1].replace("\"","");
			map.put(key,value);
		}
		return map;
	}
}
