#include<iostream>
#include<fstream>
#include<cmath>
#include <cstdlib>
#include <ctime>
using namespace std;

#define RAND ((long double)(rand())/(RAND_MAX+1.0))

/*-----------------------start of heap sort-----------------------------*/
// Function to heapify a subtree rooted at index i
void heapify(int arr[], double pb[], int** matrix, int n, int i) {
    int largest = i; // Initialize largest as root
    int left = 2 * i + 1; // Left child
    int right = 2 * i + 2; // Right child

    // If left child is larger than root
    if (left < n && arr[left] > arr[largest])
        largest = left;

    // If right child is larger than largest so far
    if (right < n && arr[right] > arr[largest])
        largest = right;

    // If largest is not root
    if (largest != i) {
        // Swap the fitness values
        int tempFit = arr[i];
        arr[i] = arr[largest];
        arr[largest] = tempFit;

        // Swap the values in pb array
        double tempPb = pb[i];
        pb[i] = pb[largest];
        pb[largest] = tempPb;

        // Swap the rows in the matrix
        for (int j = 0; j < 3; j++) {
            int tempVal = matrix[i][j];
            matrix[i][j] = matrix[largest][j];
            matrix[largest][j] = tempVal;
        }

        // Recursively heapify the affected sub-tree
        heapify(arr, pb, matrix, n, largest);
    }
}

// Main function to perform heap sort
void heapSort(int arr[], double pb[], int** matrix, int n) {
    // Build heap (rearrange array)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, pb, matrix, n, i);

    // Extract elements one by one from the heap
    for (int i = n - 1; i >= 0; i--) {
        // Move current root to the end
        int tempFit = arr[0];
        arr[0] = arr[i];
        arr[i] = tempFit;

        // Swap the values in pb array
        double tempPb = pb[0];
        pb[0] = pb[i];
        pb[i] = tempPb;

        // Swap the rows in the matrix
        for (int j = 0; j < 3; j++) {
            int tempVal = matrix[0][j];
            matrix[0][j] = matrix[i][j];
            matrix[i][j] = tempVal;
        }

        // call max heapify on the reduced heap
        heapify(arr, pb, matrix, i, 0);
    }
}

/*---------------------------end of heap sort--------------------*/


/*---------------------start of Fitness function----------------------------*/
int fitness_function(int **F_LIST, int *pop, int *demand, int n, int p){
	int temp = 0;
	int i;
	int j;
	int sum = 0;	
	for(i=0; i<n; i++){
		for(int j=0; j<p;j++){
			if(i<pop[j]){
				if(F_LIST[i][pop[j]-i] == 1){
					sum += demand[i];
					break;
				}
			}
			else{
				if(F_LIST[pop[j]][i-pop[j]] == 1){
					sum += demand[i];
					break;
				}
			}
		}
	}
	return sum;
}
/*-------------------------end of fitness function--------------------------*/


/*------------------------start of main function----------------------------*/
main(){
	/*looping variables*/
	int i,j,k;
	
	/*file management variables*/
	fstream infile,demand_file;
	string tp;
	/*end of file management variables*/
	
	
	/*static variables*/
	int GEN = 1000;	
	int D = 10;
	int gen = 0;
	int n, ct, x1, y1, x2, y2, temp, flag;
	int p, s;
	int x_max, x_min=0;
	double Nc;
	double lambda = 0.2;
	double beta = 0.5;
	/*end of static variables*/
	
	/*Random variables*/
	int *r = new int[p];
	/*end of random variables*/
	
	
	/*reading input files*/
	infile.open("DataSet/sample.txt",ios::in);
	demand_file.open("DataSet/demand-sample.txt",ios::in);
	/*end of redinginput files*/
	
	/*extracting data*/
	infile>>n;
	infile>>ct;
	infile>>p;
	infile>>s;
	cout<<"No of nodes: "<<n<<endl;
	cout<<"no of facilities: "<<p<<endl;
	cout<<"Service distance: "<<s<<endl;
	/*end of data extraction*/
	
	/*preprocessing variables*/
	int **adj = new int*[n];
	int **F_LIST = new int*[n];
	int **coor = new int*[n];
	int *demand = new int[n];
	/*end of preprocessing variables*/
	
	/*populatin variables*/
	int **pop = new int*[D];	//population
	double *Pb = new double[D];	//probability vector
	int *FITNESS = new int[D];	//fitness vector
	int *FSort = new int[D];
	int *Ftemp = new int[D];
	/*end of population variables*/
	
	/*memory allocation*/
	for(i=0; i<D; i++){
		pop[i] = new int[p];
		adj[i] = new int[n-i];
		coor[i] = new int[2];
		F_LIST[i] = new int[n-i];
	}
	/*end of memory allocation*/
	
	/*extracting demand value */
	for(int i=0; i<n; i++){			
		demand_file>>demand[i];	//extracting damand values
	}
	/*end of extracting demand values*/
	
	x_max = n;
	
	/*------------------------------pre-processing---------------------*/	
	for(i=0; i<n; i++){		
		infile>>coor[i][0]>>coor[i][1];	//reading co-ordinate values
	}
	
	//creating adjacancy list
	for(i=0; i<n; i++){
		k=i;
		for(j=0; j<n-i; j++){				
			adj[i][j] = pow( pow(coor[k][0]-coor[i][0], 2) + pow(coor[k][1]-coor[i][1], 2), 0.5);
			k++;
		}
	}
	//F_LIST-processing
	for(i=0; i<n; i++){  
		k=i;
		for(j=0; j<n-i; j++){
			if(i<k){                      
				if(adj[i][k] <= s){
					F_LIST[i][j] = 1;
				}
				else{
					F_LIST[i][j] = 0;
				}
			}
			else{
				if(adj[j][k] <= s){
					F_LIST[j][i] = 1;
				}
				else{
					F_LIST[j][i] = 0;
				}
			}
		}
	}	
//	cout<<"F_LIST"<<endl;	//displaying F_LIST
//	for(i=0; i<n; i++){
//		for(j=0; j<n-i; j++){
//			cout<<F_LIST[i][j]<<"\t";
//		}
//		cout<<endl;
//	}
	/*----------------------end of pre-processing-----------------------*/
	
	/*----------------------initial solution generation-----------------*/
	srand(time(nullptr));
	for(i=0; i<n; i++){
		for(j=0; j<p; j++){
			do{
				temp = rand() % n ;
				flag = 0;
				for(k=0; k<j; k++){
					if(r[k] == temp){
						flag = 1;
						break;
					}
				}	
			}while(flag == 1);
			r[j] = temp;
			pop[i][j] = r[j];
		}
	}	
//	cout<<"Initial population"<<endl;	//displaying the initial population
//	for(i=0; i<D; i++){
//		for(j=0; j<p; j++){	
//			cout<<pop[i][j]<<"\t";
//		}
//		cout<<endl;
//	}
	/*-------------------end of initial solution generation-----------------*/
	
	/*initial fitness calculation*/
	for(i=0; i<D; i++){
		FITNESS[i] = fitness_function(F_LIST,pop[i],demand,n,p);
	}
	cout<<"FITNESS"<<endl;	//displaying fitness values
	for(i=0; i<D; i++){
		cout<<FITNESS[i]<<"\t";	
	}
	cout<<endl;
	/*end of fitness calculation*/   
	
	/*-----initial probability vector--------------*/
	for(i=0;i<n;i++){
		int sum = 0;
		for(j=0;j<D;j++){
			for(k=0;k<p;k++){
				if(i == pop[j][k]){
					sum++;
					break;
				}
			}
		}
		Pb[i] = (double)(sum)/(10);
	}
	cout<<"initial Probability vector"<<endl;
	for(int i = 0; i< 10; i++){
		cout<<Pb[i]<<"\t";	//displaying intial probability vector
	}
	cout<<endl; 
	/*----------end of initial probability vector------------------*/
	
	/*------------------------------start of EDA algorithm----------------------------*/
	while(gen < GEN){
		
		heapSort(FITNESS, Pb, pop, n); //sorting using heapsort
		
		/*-------------------------updating the probability vector-------------------*/
		for(i=0;i<n;i++){
			int sum = 0;
			for(j=D/2;j<D;j++){
				for(k=0;k<p;k++){
					if(i == pop[j][k]){
						sum++;
						break;
					}
				}
			}
			Nc = (double)(sum)/(D/2);
			Pb[i] = (1-lambda) * Pb[i] + lambda * Nc ;
		}
		/*-------------------------end of probability vector updation-----------------*/
		
		/*------------------------generating new solution ----------------------------*/
		//for(i=0;)
		/*------------------------end of new solution generation-----------------------*/
		
		break;
		
	}
	/*------------------------------end of EDA algorithm------------------------*/
}
/*--------------------------end of main function----------------------------*/