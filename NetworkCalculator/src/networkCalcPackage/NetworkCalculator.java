package networkCalcPackage;

import java.io.File;
import java.io.PrintWriter;
import java.util.Scanner;
import java.io.FileNotFoundException;
import java.util.ArrayList;

import org.apache.commons.cli.*;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

public class NetworkCalculator {

	public static void main(String[] args) {
		CommandLineParser parser = new BasicParser();
		Options options = buildOptions();
		String pathIn = null;
		String SimOut = null;
		String AdjOut = null;
		try{
			CommandLine cmd = parser.parse(options, args);
			HelpFormatter formatter = new HelpFormatter();
			if(cmd.hasOption("h")){
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("d")){
				pathIn=cmd.getOptionValue("d");
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("a")){
				AdjOut=cmd.getOptionValue("a");
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("s")){
				SimOut=cmd.getOptionValue("s");
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
		}
		catch(ParseException exp){
			System.err.println("Problem parsing arguments:\n" + exp.getMessage());
			System.err.println("Exiting...\n");
			System.exit(0);
		}
		//String pathIn="C:\\Users\\HPriest\\Documents\\Active Projects\\Java Development\\Zm.OrthoMCL.DF.t500.tab";
		//String pathIn = args[0];
		//String SimOut = args[1];
		//String AdjOut = args[2];
		File file = null;
		try{
			file = new File(pathIn);
		}
		catch (NullPointerException e){
			System.err.println("No file found to read.\n");
			System.exit(0);
		}
		ArrayList<String> Loci = loadLoci(file);
		System.err.println("Loading Data File\n");
		ArrayList<ArrayList<Double>> DataFrame = loadData(file);
		System.err.println("Calculating Similarity\n");
		ArrayList<ArrayList<Double>> Similarity = calculateSimilarity(DataFrame);
		System.err.println("Printing similarity to file...\n");
		printMatrixToFile(Similarity,Loci,SimOut);			
		System.err.println("Calculating Adjacency...\n");
		ArrayList<ArrayList<Double>> Adjacency = calculateSigmoidAdjacency(Similarity,0.8,15);
		System.err.println("Printing Adjacency to file...\n");
		printMatrixToFile(Adjacency,Loci,AdjOut);			
	}
	
	private static Options buildOptions (){
		Options options = new Options();
		Option help = new Option( "h", "print this message" );
		Option datafile = OptionBuilder.withArgName("datafile")
				.hasArg()
				.withDescription("Data frame, tab delimited, with header, of per-gene, per-condition expression values")
				.create("d");
		Option similarity = OptionBuilder.withArgName("similarity")
				.hasArg()
				.withDescription("File for output of similarity matrix")
				.create("s");
		Option adjacency = OptionBuilder.withArgName("adjacency")
				.hasArg()
				.withDescription("File for output of adjacency matrix")
				.create("a");
		options.addOption(help);
		options.addOption(datafile);
		options.addOption(similarity);
		options.addOption(adjacency);
		return options;
	}
	
	private static ArrayList<ArrayList<Double>> calculateSigmoidAdjacency (ArrayList<ArrayList<Double>> DataFrame,double mu, double alpha){
		/* Lifted shamelessly from WGCNA:
			function (ss, mu = 0.8, alpha = 20){
				1/(1 + exp(-alpha * (ss - mu)))
			}
		 */
		ArrayList<ArrayList<Double>> Adjacency = new ArrayList<ArrayList<Double>>();
		int S=DataFrame.size();
		for(int i=0;i<S;i++){
			ArrayList<Double> OldRow = DataFrame.get(i);
			ArrayList<Double> NewRow = new ArrayList<Double>();
			for(int j=0;j<S;j++){
				double adjacency=0.0;
				if(i==j){
					adjacency = 1.0;
				}else{
					adjacency = 1/(1+Math.exp(alpha*-1*(OldRow.get(j)-mu)));
				}
				NewRow.add(adjacency);
			}
			Adjacency.add(NewRow);
		}
		return Adjacency;
	}
	
	private static void printMatrixToFile (ArrayList<ArrayList<Double>> DataFrame,ArrayList<String> Loci, String path){
		int S=DataFrame.size();
		int L=Loci.size();
		if(L != S){
			System.out.println("Loci and matrix are not the same size\n");
			System.exit(0);
		}
		try {
			PrintWriter writer = new PrintWriter(path,"UTF-8");
			for(int i=0;i<S;i++){
				ArrayList<Double> Row = new ArrayList<Double>();
				Row=DataFrame.get(i);
				String row = StringUtils.join(Row,",");
				writer.println(row+"\n");
				//System.out.println(ArrayUtils.toString(Row)+"\n");
			}
			writer.close();
		} catch (Exception e){
			
		}
		
	}
	
	private static void printMatrix (ArrayList<ArrayList<Double>> DataFrame){
		int S=DataFrame.size();
		for(int i=0;i<S;i++){
			ArrayList<Double> Row = new ArrayList<Double>();
			Row=DataFrame.get(i);
			System.out.println(ArrayUtils.toString(Row)+"\n");
		}
	}
	
	private static ArrayList<ArrayList<Double>> calculateSimilarity (ArrayList<ArrayList<Double>> DataFrame){
		ArrayList<ArrayList<Double>> Similarity = new ArrayList<ArrayList<Double>>();
		int S=DataFrame.size();
		for(int i=0;i<S;i++){
			ArrayList<Double> Row = new ArrayList<Double>();
			for(int j=0;j<S;j++){
				double correlation=0.0;
				if(i==j){
					correlation = 1.0;
				}else{
					ArrayList<Double> I = DataFrame.get(i);
					ArrayList<Double> J = DataFrame.get(j);
					//Double[] I_data = new Double[I.size()];
					//Double[] J_data = new Double[J.size()];
					Double[] Id = new Double[I.size()];
					Double[] Jd = new Double[J.size()];
					Id = I.toArray(Id);
					Jd = J.toArray(Jd);
					double[] I_data = ArrayUtils.toPrimitive(Id);
					double[] J_data = ArrayUtils.toPrimitive(Jd);
					PearsonsCorrelation corr = new PearsonsCorrelation(); 
					correlation = corr.correlation(I_data,J_data);
				}
				Row.add(correlation);
			}
			Similarity.add(Row);
		}
		return Similarity;
	}
	
	private static ArrayList<String> loadLoci (File file) {
		ArrayList<String> Loci = new ArrayList<String>();
		try {
			Scanner scanner = new Scanner(file);
			String header[] = scanner.nextLine().split("\t");
			while(scanner.hasNextLine()){
				String line=scanner.nextLine();
				String[] Line = line.split("\t");
				Loci.add(Line[0]);
			}
			scanner.close();		
		} catch (FileNotFoundException e){
			e.printStackTrace();
		}
		return Loci;
	}
	
	private static ArrayList<ArrayList<Double>> loadData (File file) {
		ArrayList<ArrayList<Double>> DataFrame = new ArrayList<ArrayList<Double>>();
		try {
			Scanner scanner = new Scanner(file);
			String header[] = scanner.nextLine().split("\t");
			while(scanner.hasNextLine()){
				String line=scanner.nextLine();
				String[] Line = line.split("\t");
				ArrayList<Double> Data = new ArrayList<Double>();
				for(int i=1;i<Line.length;i++){
					try {
						double value = Double.parseDouble(Line[i]);
						Data.add(value);
					}catch(NumberFormatException e){
						e.printStackTrace();
					}
				}
				DataFrame.add(Data);
			}
			scanner.close();
		} catch (FileNotFoundException e){
			e.printStackTrace();
		}
		return DataFrame;
	}

}
