package p2.datagenerator;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class DensityCluster {
	
	public final static String TYPE_ARC_UP_1 = "arcup1";
	public final static String TYPE_ARC_DOWN_1 = "arcdown1";
	
	public final static String TYPE_ARC_UP_2 = "arcup2";
	public final static String TYPE_ARC_DOWN_2 = "arcdown2";
	
	public final static String TYPE_BOX = "box";
	
	
	public final static double SCALING_ARC = 100;
	
	public final static double SCALING_ARC_2 = 33;
	
	public final static int FRACTION_CORE_POINTS_ARC = 4;
	
	public final static int FRACTION_CORE_POINTS_ARC_2 = 2;
	
	private int numPoints;
	private String type;
	private List<double[]> points;
	private double spacing;				
	private Random rand;
	double[] position;
	
	private List<Integer> lowIds;
	private List<Integer> highIds;
	
	public DensityCluster(int numPoints, String type, double[] position, Random rand) {
		
		this.numPoints = numPoints;
		this.type = type;
		this.points = new ArrayList<double[]>();
		this.rand = rand;
		this.position = position;
		
		switch (type) {
		
			case TYPE_ARC_UP_1:
				pointsArcUpLow1();
				break;
			case TYPE_ARC_DOWN_1:
				pointsArcDownLow1();
				break;
				
			case TYPE_ARC_UP_2:
				pointsArcUpLow2();
				break;
			case TYPE_BOX:
				pointsBoxLow();
				break;
		}
	}
	
	private void pointsArcUpLow1() {
	
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
	
	private void pointsArcDownLow1() {		
		
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

	private void pointsArcUpLow2() {
		
		this.lowIds = new ArrayList<Integer>();
		this.highIds = new ArrayList<Integer>();
		
		int numCorePoints = numPoints / FRACTION_CORE_POINTS_ARC_2;
		
		int nStart = numCorePoints / 10;
		int nEnd = numCorePoints / 10;
		int nMid = numCorePoints - nStart - nEnd;
		int numBorderPoints = numPoints - nStart - nMid - nEnd;
		int numBorderPointsStart = numBorderPoints / 10;
		int numBorderPointsEnd = numBorderPoints / 10;
		int numBorderPointsMid = numBorderPoints - numBorderPointsStart - numBorderPointsEnd;
		
		//Start
		List<double[]> coreStart = new ArrayList<double[]>();
		double[] linStart = linspace(0, (Math.PI * 1.2) / 3, nStart);
		for (int i = 0; i < nStart; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = Math.cos(linStart[i]) * SCALING_ARC_2; 
			corePoint[1] = Math.sin(linStart[i]) * SCALING_ARC_2;
			coreStart.add(corePoint);
		}
		this.points.addAll(coreStart);

		//Start border
		List<double[]> borderStart = new ArrayList<double[]>();
		double spacingStart = (linStart[1] -linStart[0]) * SCALING_ARC_2;
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
		double[] linEnd = linspace(2 * (Math.PI * 1.2) / 3, (Math.PI * 1.2), nEnd);
		for (int i = 0; i < nEnd; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = Math.cos(linEnd[i]) * SCALING_ARC_2; 
			corePoint[1] = Math.sin(linEnd[i]) * SCALING_ARC_2;
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
		double[] linMid = linspace((Math.PI * 1.2) / 3, 2 * (Math.PI * 1.2) / 3, nMid);
		for (int i = 0; i < nMid; i++) {
			double[] corePoint = new double[2];
			corePoint[0] = Math.cos(linMid[i]) * SCALING_ARC_2; 
			corePoint[1] = Math.sin(linMid[i]) * SCALING_ARC_2;
			coreMid.add(corePoint);
		}
		this.points.addAll(coreMid);
		
		//Mid border
		List<double[]> borderMid = new ArrayList<double[]>();
		double spacingMid = (linMid[1] - linMid[0]) * SCALING_ARC_2;
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
	
	
	private void pointsBoxLow() {
		
		this.lowIds = new ArrayList<Integer>();
		this.highIds = new ArrayList<Integer>();
		
		int numCorePoints = numPoints / 2;
		
		double coreMinX = -30;
		double coreMaxX = -15;
		double coreMinY = 70;
		double coreMaxY = 80;
		
		List<double[]> corePoints = new ArrayList<double[]>();
		for (int i = 0; i < numCorePoints; i++) {
			
			double[] corePoint = new double[2];
			
			corePoint[0] = coreMinX + (coreMaxX - coreMinX) * rand.nextDouble();
			corePoint[1] = coreMinY + (coreMaxY - coreMinY) * rand.nextDouble();
			
			corePoints.add(corePoint);
		}
		
		double borderMinX = -30;
		double borderMaxX = 25;
		double borderMinY = 0;
		double borderMaxY = 80;
		
		List<double[]> borderPoints = new ArrayList<double[]>();
		for (int i = 0; i < (numPoints - numCorePoints); ) {
			
			double[] borderPoint = new double[2];
			
			borderPoint[0] = borderMinX + (borderMaxX - borderMinX) * rand.nextDouble();
			borderPoint[1] = borderMinY + (borderMaxY - borderMinY) * rand.nextDouble();
			
			if (! (borderPoint[0] >= coreMinX && borderPoint[0] <= coreMaxX && borderPoint[1] >= coreMinY && borderPoint[1] <= coreMaxY)) {
				borderPoints.add(borderPoint);
				++i;
			}
			
		}
		
		
		this.points.addAll(corePoints);
		this.points.addAll(borderPoints);
		
		//Set high density ids
		for (int id = 0; id < corePoints.size(); id++)
			highIds.add(id);
		
		//Set low density ids
		for (int id = corePoints.size(); id < this.points.size(); id++)
			lowIds.add(id);
		
		
		//Move points to position
		for (double[] point : this.points) {
			point[0] += position[0];
			point[1] += position[1];
		}
		
		
	}
	
}
