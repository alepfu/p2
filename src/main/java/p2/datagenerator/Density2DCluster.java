package p2.datagenerator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class Density2DCluster {
	
	public final static String TYPE_ARC_UP = "arcup";
	public final static String TYPE_ARC_DOWN = "arcdown";
	public final static String TYPE_CIRCLE = "circle";
	
	public final static double SCALING_ARC = 100;
	public final static double SCALING_CIRCLE = 50;
	
	public final static int FRACTION_CORE_POINTS_CIRCLE = 2;
	public final static int FRACTION_CORE_POINTS_ARC = 4;
	
	private int numPoints;
	private String type;
	private List<double[]> points;
	private double spacing;
	private Random rand;
	double[] position;
	
	public Density2DCluster(int numPoints, String type, double[] position, Random rand) {
		
		this.numPoints = numPoints;
		this.type = type;
		this.points = new ArrayList<double[]>();
		this.rand = rand;
		this.position = position;
		
		switch (type) {
		case TYPE_ARC_UP:
			pointsArcUp();
			break;
		case TYPE_ARC_DOWN:
			pointsArcDown();
			break;
		case TYPE_CIRCLE:
			pointsCircle();
			break;
		}
	}
	
	private void pointsArcUp() {
		
		//Core points
		List<double[]> corePoints = new ArrayList<double[]>();
		int numCorePoints = numPoints / FRACTION_CORE_POINTS_ARC;
		double[] l = linspace(0, Math.PI, numCorePoints);
		this.spacing = l[1] * SCALING_ARC;
		for (int i = 0; i < numCorePoints; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = Math.cos(l[i]) * SCALING_ARC; 
			corePoint[1] = Math.sin(l[i]) * SCALING_ARC;
			corePoints.add(corePoint);
		}
		this.points.addAll(corePoints);
		
		//Border Points
		List<double[]> borderPoints = new ArrayList<double[]>();
		int numBorderPoints = numPoints - numCorePoints;
		for (int i = 0; i < numBorderPoints; i++) {
			double[] corePoint = corePoints.get(rand.nextInt(corePoints.size()));
			double[] borderPoint = new double[2];
			borderPoint[0] = corePoint[0] + (spacing * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoint[1] = corePoint[1] + (spacing * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoints.add(borderPoint);
		}
		this.points.addAll(borderPoints);
		
		//Randomize the data
		Collections.shuffle(this.points);
		
		//Move to position
		for (double[] point : this.points) {
			point[0] += position[0];
			point[1] += position[1];
		}
	}
	
	private void pointsArcDown() {
		
		//Core points
		List<double[]> corePoints = new ArrayList<double[]>();
		int numCorePoints = numPoints / FRACTION_CORE_POINTS_ARC;
		double[] l = linspace(0, Math.PI, numCorePoints);
		this.spacing = l[1] * SCALING_ARC;
		for (int i = 0; i < numCorePoints; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = 1 - Math.cos(l[i]) * SCALING_ARC; 
			corePoint[1] = 1 - Math.sin(l[i]) * SCALING_ARC;
			corePoints.add(corePoint);
		}
		this.points.addAll(corePoints);
		
		//Border Points
		List<double[]> borderPoints = new ArrayList<double[]>();
		int numBorderPoints = numPoints - numCorePoints;
		for (int i = 0; i < numBorderPoints; i++) {
			double[] corePoint = corePoints.get(rand.nextInt(corePoints.size()));
			double[] borderPoint = new double[2];
			borderPoint[0] = corePoint[0] + (spacing * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoint[1] = corePoint[1] + (spacing * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoints.add(borderPoint);
		}
		this.points.addAll(borderPoints);
		
		//Randomize the data
		Collections.shuffle(this.points);
		
		//Move to position
		for (double[] point : this.points) {
			point[0] += position[0];
			point[1] += position[1];
		}
	}
	
	private void pointsCircle() {
		
		//Core points
		List<double[]> corePoints = new ArrayList<double[]>();
		int numCorePoints = numPoints / FRACTION_CORE_POINTS_CIRCLE;
		double[] l = linspace(0, 2 * Math.PI, numCorePoints);
		this.spacing = l[1] * SCALING_CIRCLE;
		for (int i = 0; i < numCorePoints; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = Math.cos(l[i]) * SCALING_CIRCLE; 
			corePoint[1] = Math.sin(l[i]) * SCALING_CIRCLE;
			corePoints.add(corePoint);
		}
		this.points.addAll(corePoints);
		
		//Border Points
		List<double[]> borderPoints = new ArrayList<double[]>();
		int numBorderPoints = numPoints - numCorePoints;
		for (int i = 0; i < numBorderPoints; i++) {
			double[] corePoint = corePoints.get(rand.nextInt(corePoints.size()));
			double[] borderPoint = new double[2];
			borderPoint[0] = corePoint[0] + (spacing * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoint[1] = corePoint[1] + (spacing * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoints.add(borderPoint);
		}
		this.points.addAll(borderPoints);
		
		//Randomize the data
		Collections.shuffle(this.points);
		
		//Move to position
		for (double[] point : this.points) {
			point[0] += position[0];
			point[1] += position[1];
		}
	}

	public int getNumPoints() {
		return numPoints;
	}

	public void setNumPoints(int numPoints) {
		this.numPoints = numPoints;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public List<double[]> getPoints() {
		return points;
	}

	public void setPoints(List<double[]> points) {
		this.points = points;
	}
	
	private static double[] linspace(double f, double t, int n) {  
	    double[] d = new double[n];  
	    for (int i = 0; i < n; i++)  
	        d[i] = f + i * (t - f) / (n - 1);  
	    return d;  
	}  
	
}
