package p2.datagenerator;

import de.lmu.ifi.dbs.elki.application.GeneratorXMLSpec;

public class ElkiDataGenerator {

	public static void main(String[] args) {
		
		
		
		
		//TODO not finished
		
		
		
		
		
		
		String[] generatorArgs = { "-app.out", "data/data_elki.csv",
								   "-bymodel.spec", "clustermodels/mouse.xml",
								   "-verbose", "true"};
		
		GeneratorXMLSpec.main(generatorArgs);
		

		

	}

}
