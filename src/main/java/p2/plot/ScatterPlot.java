package p2.plot;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class ScatterPlot {
	
	public static void plot(List<double[]> dataPoints, String filename) {

		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series = new XYSeries("Cluster");

		for (double[] dataPoint : dataPoints)
			series.add(dataPoint[0], dataPoint[1]);

		dataset.addSeries(series);
		
		JFreeChart chart = ChartFactory.createScatterPlot("", "", "", dataset, PlotOrientation.VERTICAL, false, false, false);
		ChartFrame frame = new ChartFrame("", chart);
		frame.pack();
		frame.setVisible(true);
		
		int width = 660;
		int height = 420;

		try {
			ChartUtilities.saveChartAsJPEG(new File(filename), chart, width, height);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
}
