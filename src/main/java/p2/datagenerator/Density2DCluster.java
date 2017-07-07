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
	private boolean doLowDensityRegions;
	
	private List<Integer> lowIds;
	private List<Integer> highIds;
	
	public Density2DCluster(int numPoints, String type, double[] position, Random rand, boolean doLowDensityRegions) {
		
		this.numPoints = numPoints;
		this.type = type;
		this.points = new ArrayList<double[]>();
		this.rand = rand;
		this.position = position;
		this.doLowDensityRegions = doLowDensityRegions;
		
		switch (type) {
		
			case TYPE_ARC_UP:
				if (this.doLowDensityRegions)
					pointsArcUpLow();
				else
					pointsArcUp();
				break;
			case TYPE_ARC_DOWN:
				if (this.doLowDensityRegions)
					pointsArcDownLow();
				else
					pointsArcDown();
				break;
			case TYPE_CIRCLE:
				pointsCircle();
				break;
		}
	}
	
	private void pointsArcUp() {
		
		//Add core points
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
		
		//Add border points to core points
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
		
		//Move points to position
		for (double[] point : this.points) {
			point[0] += position[0];
			point[1] += position[1];
		}
	}
	
	private void pointsArcUpLow() {
	
		this.lowIds = new ArrayList<Integer>();
		this.highIds = new ArrayList<Integer>();
		
		int numCorePoints = numPoints / FRACTION_CORE_POINTS_ARC;
		
		int nStart = numCorePoints / 10;
		int nEnd = numCorePoints / 10;
		int nMid = numCorePoints - nStart - nEnd;
		int numBorderPoints = numPoints - nStart - nMid - nEnd;
		int numBorderPointsStart = numBorderPoints / 10;
		int numBorderPointsEnd = numBorderPoints / 10;
		int numBorderPointsMid = numBorderPoints - numBorderPointsStart - numBorderPointsEnd;
		
		//Start
		List<double[]> coreStart = new ArrayList<double[]>();
		double[] linStart = linspace(0, Math.PI / 3, nStart);
		for (int i = 0; i < nStart; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = Math.cos(linStart[i]) * SCALING_ARC; 
			corePoint[1] = Math.sin(linStart[i]) * SCALING_ARC;
			coreStart.add(corePoint);
		}
		this.points.addAll(coreStart);

		//Start border
		List<double[]> borderStart = new ArrayList<double[]>();
		double spacingStart = (linStart[1] -linStart[0]) * SCALING_ARC;
		for (int i = 0; i < numBorderPointsStart; i++) {
			double[] corePoint = coreStart.get(rand.nextInt(coreStart.size()));
			double[] borderPoint = new double[2];
			borderPoint[0] = corePoint[0] + (spacingStart * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoint[1] = corePoint[1] + (spacingStart * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderStart.add(borderPoint);
		}
		this.points.addAll(borderStart);
		
		
		//End
		List<double[]> coreEnd = new ArrayList<double[]>();
		double[] linEnd = linspace(2 * Math.PI / 3, Math.PI, nEnd);
		for (int i = 0; i < nEnd; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = Math.cos(linEnd[i]) * SCALING_ARC; 
			corePoint[1] = Math.sin(linEnd[i]) * SCALING_ARC;
			coreEnd.add(corePoint);
		}
		this.points.addAll(coreEnd);
		
		//End border
		List<double[]> borderEnd = new ArrayList<double[]>();
		double spacingEnd = (linEnd[1] - linEnd[0]) * SCALING_ARC;
		for (int i = 0; i < numBorderPointsEnd; i++) {
			double[] corePoint = coreEnd.get(rand.nextInt(coreEnd.size()));
			double[] borderPoint = new double[2];
			borderPoint[0] = corePoint[0] + (spacingEnd * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoint[1] = corePoint[1] + (spacingEnd * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderEnd.add(borderPoint);
		}
		this.points.addAll(borderEnd);
		
		
		//Set low density ids
		for (int id = 0; id < this.points.size(); id++)
			lowIds.add(id);
		
		
		
		//Mid
		List<double[]> coreMid = new ArrayList<double[]>();
		double[] linMid = linspace(Math.PI / 3, 2 * Math.PI / 3, nMid);
		for (int i = 0; i < nMid; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = Math.cos(linMid[i]) * SCALING_ARC; 
			corePoint[1] = Math.sin(linMid[i]) * SCALING_ARC;
			coreMid.add(corePoint);
		}
		this.points.addAll(coreMid);
		
		//Mid border
		List<double[]> borderMid = new ArrayList<double[]>();
		double spacingMid = (linMid[1] - linMid[0]) * SCALING_ARC;
		for (int i = 0; i < numBorderPointsMid; i++) {
			double[] corePoint = coreMid.get(rand.nextInt(coreMid.size()));
			double[] borderPoint = new double[2];
			borderPoint[0] = corePoint[0] + (spacingMid * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoint[1] = corePoint[1] + (spacingMid * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderMid.add(borderPoint);
		}
		this.points.addAll(borderMid);
		
		
		//Set high density ids
		for (int id = lowIds.size(); id < this.points.size(); id++)
			highIds.add(id);
		
		
		//Move points to position
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
		
		//Move to position
		for (double[] point : this.points) {
			point[0] += position[0];
			point[1] += position[1];
		}
	}
	
	private void pointsArcDownLow() {		
		
		this.lowIds = new ArrayList<Integer>();
		this.highIds = new ArrayList<Integer>();
		
		int numCorePoints = numPoints / FRACTION_CORE_POINTS_ARC;
		
		int nStart = numCorePoints / 10;
		int nEnd = numCorePoints / 10;
		int nMid = numCorePoints - nStart - nEnd;
		int numBorderPoints = numPoints - nStart - nMid - nEnd;
		int numBorderPointsStart = numBorderPoints / 10;
		int numBorderPointsEnd = numBorderPoints / 10;
		int numBorderPointsMid = numBorderPoints - numBorderPointsStart - numBorderPointsEnd;
		
		//Start
		List<double[]> coreStart = new ArrayList<double[]>();
		double[] linStart = linspace(0, Math.PI / 3, nStart);
		for (int i = 0; i < nStart; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = 1 - Math.cos(linStart[i]) * SCALING_ARC; 
			corePoint[1] = 1 - Math.sin(linStart[i]) * SCALING_ARC;
			coreStart.add(corePoint);
		}
		this.points.addAll(coreStart);

		//Start border
		List<double[]> borderStart = new ArrayList<double[]>();
		double spacingStart = (linStart[1] -linStart[0]) * SCALING_ARC;
		for (int i = 0; i < numBorderPointsStart; i++) {
			double[] corePoint = coreStart.get(rand.nextInt(coreStart.size()));
			double[] borderPoint = new double[2];
			borderPoint[0] = corePoint[0] + (spacingStart * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoint[1] = corePoint[1] + (spacingStart * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderStart.add(borderPoint);
		}
		this.points.addAll(borderStart);
		
		
		//End
		List<double[]> coreEnd = new ArrayList<double[]>();
		double[] linEnd = linspace(2 * Math.PI / 3, Math.PI, nEnd);
		for (int i = 0; i < nEnd; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = 1 - Math.cos(linEnd[i]) * SCALING_ARC; 
			corePoint[1] = 1 - Math.sin(linEnd[i]) * SCALING_ARC;
			coreEnd.add(corePoint);
		}
		this.points.addAll(coreEnd);
		
		//End border
		List<double[]> borderEnd = new ArrayList<double[]>();
		double spacingEnd = (linEnd[1] - linEnd[0]) * SCALING_ARC;
		for (int i = 0; i < numBorderPointsEnd; i++) {
			double[] corePoint = coreEnd.get(rand.nextInt(coreEnd.size()));
			double[] borderPoint = new double[2];
			borderPoint[0] = corePoint[0] + (spacingEnd * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoint[1] = corePoint[1] + (spacingEnd * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderEnd.add(borderPoint);
		}
		this.points.addAll(borderEnd);
		
		
		//Set low density ids
		for (int id = 0; id < this.points.size(); id++)
			lowIds.add(id);
			
		
		//Mid
		List<double[]> coreMid = new ArrayList<double[]>();
		double[] linMid = linspace(Math.PI / 3, 2 * Math.PI / 3, nMid);
		for (int i = 0; i < nMid; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = 1 - Math.cos(linMid[i]) * SCALING_ARC; 
			corePoint[1] = 1 - Math.sin(linMid[i]) * SCALING_ARC;
			coreMid.add(corePoint);
		}
		this.points.addAll(coreMid);
		
		//Mid border
		List<double[]> borderMid = new ArrayList<double[]>();
		double spacingMid = (linMid[1] - linMid[0]) * SCALING_ARC;
		for (int i = 0; i < numBorderPointsMid; i++) {
			double[] corePoint = coreMid.get(rand.nextInt(coreMid.size()));
			double[] borderPoint = new double[2];
			borderPoint[0] = corePoint[0] + (spacingMid * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderPoint[1] = corePoint[1] + (spacingMid * rand.nextDouble() * (rand.nextBoolean() ? 1 : -1)); 
			borderMid.add(borderPoint);
		}
		this.points.addAll(borderMid);
		
		
		//Set high density ids
		for (int id = lowIds.size(); id < this.points.size(); id++)
			highIds.add(id);
		
		//Move points to position
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

	public List<Integer> getLowIds() {
		return lowIds;
	}  
	
	public List<Integer> getHighIds() {
		return highIds;
	} 
	
}
