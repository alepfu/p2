package p2;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import p2.datagenerator.DataGenerator;

public class TestDataGenerator {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void test() {
		
		
		
		
		for (long i = 1; i <= 100; i++) {
		
			String[] args = {"--seed", i+"", "--overlap"};
			
			DataGenerator.main(args);
		
		}
		
	}

}
