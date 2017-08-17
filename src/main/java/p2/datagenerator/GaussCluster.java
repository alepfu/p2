package p2.datagenerator;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class GaussCluster {

	private int numPoints;
	private List<double[]> points;
	private double deviation;
	private double[] mean; 
	private Random rand;
	
	public GaussCluster(int numPoints, double[] mean, double deviation, Random rand) {
		
		this.numPoints = numPoints;
		this.points = new ArrayList<double[]>();
		this.deviation = deviation;
		this.mean = mean;
		this.rand = rand;
		
		for (int i = 0; i < numPoints; i++) {
			double[] point = new double[2];
			point[0] = rand.nextGaussian() * deviation + mean[0]; 
			point[1] = rand.nextGaussian() * deviation + mean[1];
			points.add(point);
		}
	}

	public int getNumPoints() {
		return numPoints;
	}

	public void setNumPoints(int numPoints) {
		this.numPoints = numPoints;
	}

	public List<double[]> getPoints() {
		return points;
	}

	public void setPoints(List<double[]> points) {
		this.points = points;
	}

	public double getDeviation() {
		return deviation;
	}

	public void setDeviation(double deviation) {
		this.deviation = deviation;
	}

	public double[] getMean() {
		return mean;
	}

	public void setMean(double[] mean) {
		this.mean = mean;
	}

	public Random getRand() {
		return rand;
	}

	public void setRand(Random rand) {
		this.rand = rand;
	}
	
}
