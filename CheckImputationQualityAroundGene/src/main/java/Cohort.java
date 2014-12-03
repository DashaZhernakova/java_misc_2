import java.util.ArrayList;

/**
 * Created by dashazhernakova on 21.11.14.
 */
public class Cohort {
	private ArrayList<Float> infos;

	public Cohort(){
		infos = new ArrayList<Float>();
	}

	public Cohort(ArrayList<Float> inputInfos){
		infos = inputInfos;
	}

	public void addInfo (float info){
		infos.add(info);

	}
	public int countNumSNPsInfoLessThanThreshold(float threshold){
		int num = 0;
		for (float info : infos){
			if (info < threshold)
				num++;
		}
		return num;
	}

	public float countMeanInfo(){
		int num = 0;
		float sum = 0;
		for (float info : infos){
			sum += info;
			num++;
		}
		return sum/num;
	}
}
