package p2.plot;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class ScatterPlot {
	
	public static void show(List<double[]> dataPoints) {

		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series = new XYSeries("Cluster");

		for (double[] dataPoint : dataPoints)
			series.add(dataPoint[0], dataPoint[1]);

		dataset.addSeries(series);
		
		JFreeChart chart = ChartFactory.createScatterPlot("", "", "", dataset, PlotOrientation.VERTICAL, false, false, false);
		ChartFrame frame = new ChartFrame("", chart);
		frame.pack();
		frame.setVisible(true);
	}
	
}
