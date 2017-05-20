package p2.clustering;

public class ExtDensityData {
	
	private double[][] data;
	private double[][] dummy;
	
	public ExtDensityData(double[][] data, double[][] dummy) {
		this.data = data;
		this.dummy = dummy;
	}

	/**
	 * Normalizes the dummy enchoding to the variance of the data, then appends the dummy encoding to the data.
	 * @return Concat between data and dummy encoding.
	 */
	public double[][] getExtData() {
		
		if (data.length != dummy.length)
			throw new IllegalStateException("data.length != dummy.length");
		
		//TODO normalization
		
		int numRows = data.length;
		int numCols = data[0].length + dummy[0].length;
		
		double[][] extData = new double[numRows][numCols];
		
		//TODO concat data and dummy
		
		
		return extData;
	}
}
