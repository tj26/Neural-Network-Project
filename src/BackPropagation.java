
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;

public class BackPropagation {
	/*
	 * nH = number of neurons in hidden layer
	 * nI = number of neurons in input layer
	 * zH = output at hidden neurons before applying transfer function
	 * aH = output of hidden neurons after applying transfer function
	 * zO = output of neuron in the output layer before transfer function
	 * aO - output from neuron in the output layer after transfer function
	 * WeightsToHidden = weights from input layer to hidden layer
	 * WeightsToOutput = weights from hidden layer to output layer
	 * sH = sensitivity of the hidden layer neurons
	 * sO = sensitivity of the output layer neurons
	 */
	static int nH,nI;		
	static int numInputs;
	
	static double [][] Input ;
	double [][] WeightsToHidden;
	double [] HiddenBias;
	double [][] zH;
	double [][] aH;
	double [] sH;	
	double learningRate;
	
	double [][]WeightsToOutput;
	double OutputBias;
	
	static double[] Output;
	double zO[][] , aO[][];
	double sO;
	
	double T = 0.05;
	double error;


/* * Uncomment the following and comment the above to test
 * against weights and bias as defined in the book	
 
	int nH=4,nI=2,numInputs=4;
	double [][] Input = {{1,1}, {1,-1}, {-1,1}, {-1,-1}};
	double [] Output =  {-1, 1, 1, -1};
	
	double learningRate = 0.2;
	
	double OutputBias = -0.1401;
	double WeightsToHidden[][] = {{0.1970, 0.3191, -0.1448, 0.3594}, {0.3099,0.1904, -0.0347, -0.4861}};
	double HiddenBias[] = {-0.3378, 0.2771, 0.2859, -0.3329};
	double[][] WeightsToOutput = {{0.4919}, {-0.2913}, {-0.3979}, {0.3581}};
	
	double [][] zH = new double[1][nH];
	double [][] aH = new double[1][nH];
	double [] sH = new double[nH];
	double sO;
	double[][] zO = new double[1][1];
	double[][] aO = new double[1][1];
	double error;
	double T = 0.05;
	
*/	// Generate random numbers 
	public double generateRandom(double t){
	
		double leftLimit = (double)(-t) ;
		double rightLimit = (double) t;
		return (leftLimit+new Random().nextDouble() * (rightLimit - leftLimit));
	}
	
	// This function uses randomly generated weights and bias but variable number 
	// of neurons in input and hidden layer
	public void getIntputs(double d) {
	
		WeightsToHidden = new double[nI][nH];
		HiddenBias = new double[nH];
		WeightsToOutput = new double[nH][1];
		zH = new double[1][nH];
		aH = new double[1][nH];
		zO = new double[1][1];
		aO = new double[1][1];
		sH = new double[nH];
		
	
		// Initialize weights from input to hidden layer
		for(int i = 0; i < nI; i++){
			for(int j = 0 ; j < nH; j++){
				WeightsToHidden[i][j] = (double) generateRandom(d);
			
			}
		}
		
		// Initialize bias for hidden layer
		for(int i =0 ; i < nH;i++){
			HiddenBias[i] =  (double)  generateRandom(d);
			
		}
		
		// Initialize weights to output layer
		for(int i = 0 ; i < nH; i++){
			WeightsToOutput[i][0] =  (double) generateRandom(d);
		
		}
		
		// Initialize bias to output layer
		OutputBias =  (double) generateRandom(d);
		
	}
	
	// This function takes user input for weights and bias
	// and variable number of neurons in input and hidden layer
	public void getIntputsUser() {
		Scanner s = new Scanner(System.in);
		
		System.out.println("Enter the number of input neurons");
		nI = s.nextInt();
	
		
		System.out.println("Enter the number of neurons in hidden layer");
		nH = s.nextInt();
		
		
		System.out.println("Enter the number of inputs");
		numInputs = s.nextInt();
		
		
		Input = new double[numInputs][nI];
		Output = new double[numInputs];
		WeightsToHidden = new double[nI][nH];
		HiddenBias = new double[nH];
		WeightsToOutput = new double[nH][1];
		zH = new double[1][nH];
		aH = new double[1][nH];
		zO = new double[1][1];
		aO = new double[1][1];
		sH = new double[nH];
		
		for(int i = 0; i < numInputs; i++){
			for(int j = 0 ; j < nI; j++){
				System.out.println("Enter input "+ (i+1) + " for neuron "+ (j+1));
				Input[i][j] = s.nextInt();
			}
		}
		
		
		for(int i = 0; i < numInputs; i++){
			System.out.println("Enter the output " + (i+1));
			Output[i] = s.nextInt();
		}
		
		System.out.println("Enter Weights for hidden layer");
		// Initialize weights from input to hidden layer
		for(int i = 0; i < nI; i++){
			for(int j = 0 ; j < nH; j++){
				System.out.println("Enter Weight from input neuron" + (i+1)+ " to hidden neuron "+ (j+1));
				WeightsToHidden[i][j] = s.nextDouble();
			
			}
		}
		
		System.out.println("Enter Bias for hidden layer");
		// Initialize bias for hidden layer
		for(int i =0 ; i < nH;i++){
			System.out.println("Enter bias for neuron" + (i+1));
			HiddenBias[i] = s.nextDouble();
			
		}
		
		// Initialize weights to output layer
		System.out.println("Enter Weights for output layer");
		for(int i = 0 ; i < nH; i++){
			System.out.println("Enter weight for neuron" + (i+1));
			WeightsToOutput[i][0] = s.nextDouble();
		
		}
		
		// Initialize bias to output layer
		System.out.println("Enter bias for output neuron");
		OutputBias = s.nextDouble();
		
	}

	// Bipolar sigmoid transfer function. 
	public double transfer(double x, double xO){
		return (1-Math.exp(-x/xO))/(1+Math.exp(-x/xO));
	}
	
	// Differentiated Bipolar sigmoid function
	public double difftransfer(double x,double xO){
		return ((0.5/xO) *(1+(1-Math.exp(-x/xO))/(1+Math.exp(-x/xO)))* (1-(1-Math.exp(-x/xO))/(1+Math.exp(-x/xO))));
	}
	
	// Function to calculate matrix multiplication
	 public static double[][] multMatrix(double a[][], double b[][]){//a[m][n], b[n][p]
		   if(a.length == 0) return new double[0][0];
		   if(a[0].length != b.length) return null; 

		   int n = a[0].length;
		   int m = a.length;
		   int p = b[0].length;
		   double ans[][] = new double[m][p];

		   for(int i = 0;i < m;i++){
		      for(int j = 0;j < p;j++){
		         for(int k = 0;k < n;k++){
		            ans[i][j] += a[i][k] * b[k][j];
		         }
		      }
		   }
		   return ans;
	 }
	 
	public Integer train(double xO,double lrate, int choice){
		
		double[][] temp = new double[1][nI];
		Integer epoch = 20;
		this.learningRate = lrate;
		
		for(Integer k=0; k<epoch; k++){
			error =0;
			for(int i=0; i < numInputs; i++){
				temp[0] = Input[i];
				
				zH = multMatrix(temp,WeightsToHidden);
				for(int j =0; j < nH;j++){
					zH[0][j] += HiddenBias[j];
					aH[0][j] = transfer(zH[0][j],xO);
					System.out.println("Activation of hidden layer " + aH[0][j]);
				}
				
				zO = multMatrix(aH,WeightsToOutput) ;
				
				zO[0][0] += OutputBias;
				aO[0][0] = transfer(zO[0][0],xO);		// Feed forward output
				
				System.out.println("a0 "+aO[0][0]);
				
				//Sensitivity of output layer
				if(choice==1){
					sO = (aO[0][0] - Output[i]) * difftransfer(zO[0][0],xO);
				}
				else{
					sO = -2/((aO[0][0]) + Output[i]) * difftransfer(zO[0][0],xO);
				}
				
				System.out.println("Sensitivity of output layer "+ sO);
				
				// Sensitivity of hidden layer
				for(int j=0; j < nH; j++){
					sH[j] = difftransfer(zH[0][j],xO) * WeightsToOutput[j][0] *sO;
					//System.out.println("Sensitivity of outer layer "+ sH[j]);
					
				}
				
				// Update hidden layer weights
				for(int l= 0;l<nI; l++){
					for(int m=0; m <nH;m++){
						WeightsToHidden[l][m] = WeightsToHidden[l][m] - (learningRate *Input[i][l] * sH[m]);
					}
				}
				
				//Update Hidden Layer Bias
				for(int l =0 ; l < nH; l++){
					HiddenBias[l] = HiddenBias[l] - (learningRate * sH[l]);
				}
				
				// Update Output layer weights
				for(int l=0; l < nH; l++){
					WeightsToOutput[l][0] = WeightsToOutput[l][0] - (learningRate * aH[0][l] * sO);
				}
				
				// Update Output layer Bias
				OutputBias = OutputBias - (learningRate * sO);
				
				error += Math.pow(aO[0][0]- Output[i], 2);
				//System.out.println("K "+ k);
				//printWeightBias();
			}
			
			double round = Math.round(error*100.0)/100.0;
			if(round <= T){
				return k;
			}
		}
		return 0;
	}
	
	// Print weights and bias
	public void printWeightBias(){
		
		System.out.println("Hidden layer Weights");
		
		for(int i =0; i < nI; i++){
			for(int j=0; j< nH; j++){
				System.out.print(WeightsToHidden[i][j] + " ");
			}
			System.out.println();
		}
		
		System.out.println("Hidden layer Bias");
		
		for(int i=0; i < nH;i++){
			System.out.print(HiddenBias[i] + " ");
		}
		System.out.println();
		
		System.out.println("Output layer weights");
		
		for(int i = 0; i < nH; i++){
			for(int j=0;j <1 ; j++){
				System.out.print(WeightsToOutput[i][j] + " ");
			}
			System.out.println();
		}
		
		System.out.println("Output neuron  Bias: "+ OutputBias);
		
	}
	
	@SuppressWarnings("null")
	public static void getInfo(int [] array){
		int converged=0;
		int index=0;
		for(int i=0; i < array.length; i++){
			if(array[i] !=0){
				converged++;
			}
		}
		int []epochs = new int[converged];
		for(int i=0; i < array.length; i++){
			if(array[i] !=0){
				epochs[index] =array[i];
				index++;
			}
		}
		
		Arrays.sort(epochs);
		// average
		int sum=0;
		for(int i=0; i < epochs.length; i++){
			 sum = sum+epochs[i];
		}
		
		System.out.println(" Total number of times convergence occurred = " + converged);
		System.out.println(" Total number of times convergence did not occur =" + (array.length-converged) );
		
		if(converged !=0){
			System.out.println(" Average number of epochs needed: " + sum/(epochs.length));
			
			// max 
			System.out.println(" Maximum number of epochs needed: " + epochs[epochs.length -1]);
			
			
			// min 
			System.out.println(" Minimum number of epochs needed: " + epochs[0]);
			
			// median
			if(epochs.length %2 !=0){
				System.out.println(" Median of the epochs: " + epochs[(epochs.length)/2]);
			}
			else{
				System.out.println(" Median of the epochs: " + (epochs[(epochs.length)/2]+ epochs[(epochs.length-1)/2])/2);
			}
		}
		
	}
	
	public static void main(String args[]){
		Scanner s = new Scanner(System.in);
		BackPropagation bp = new BackPropagation();	
		
		System.out.println("Enter the number of input neurons");
		nI = s.nextInt();
		
		System.out.println("Enter the number of inputs");
		numInputs = s.nextInt();
		
		Input = new double[numInputs][nI];
		Output = new double[numInputs];
		
		for(int i = 0; i < numInputs; i++){
			for(int j = 0 ; j < nI; j++){
				System.out.println("Enter input "+ (i+1) + " for neuron "+ (j+1));
				Input[i][j] = s.nextInt();
			}
		}
			
		for(int i = 0; i < numInputs; i++){
			System.out.println("Enter the output " + (i+1));
			Output[i] = s.nextInt();
		}
				
		System.out.println("Enter 1 for Quadratic cost function 2 for Cross Entropy");
		int choice = s.nextInt();
		
		double[] xO = {0.5,1.0,1.5};
		
		ArrayList<Integer> arl = new ArrayList<Integer>();
				for(int n=4; n<=4; n= n+2)
				{	
					nH = n;
					int Case =1;	// 27 cases for each N1, total=135
					for(int random =5; random <=5; random =random+5)
					{		
						for(int lRate = 3; lRate<=3; lRate++ )
							{
							for(int i =1; i < 2; i++)
								{
								arl = new ArrayList<Integer>();
								for(int k=0; k < 100; k++){
									bp.getIntputs((0.1 * random));
									int e;
									//bp.printWeightBias();
									if(choice ==1){
										 e = bp.train(xO[i],(0.1*lRate), 1);
										 //System.exit(0);
									}
									else{
										 e = bp.train(xO[i],(0.1*lRate), 2);
										 
									}
									
									arl.add(e);
								}
								
								int[] intArray = new int[arl.size()];
								for (int j = 0; j < intArray.length; j++) {
								    intArray[j] = arl.get(j); 
								}
								System.out.println();
								System.out.println("######### Case " + Case + " #########");
								System.out.println(" N1: " + n + " Learning Rate: "+ (0.1* lRate) + " Random: " + (0.1*random) + " XO: " +xO[i]);
								getInfo(intArray);
								Case++;
								
						}
					}
				  }
				}
							
		
	}// main
} //class
