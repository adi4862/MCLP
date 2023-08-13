#include<iostream>
#include<fstream>
#include<cmath>
#include <cstdlib>
#include <ctime>
#include<vector>
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

        // call m1 heapify on the reduced heap
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
	int GEN = 10;	
	int D = 10;
	int gen = 0;
	int n, ct, x1, y1, x2, y2, temp, flag;
	int p, s;
	int x_max, x_min=0;
	double Nc;
	double lambda = 0.2;
	double beta = 0.9;
	/*end of static variables*/
	
	/*Random variables*/
	int *r = new int[p];
	double r1,r2;
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
	vector<vector<int>> NEW_POP;
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
	for(i=0; i<n; i++){			
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
	
	/*----------------------initial probability vector-----------------------*/
	double temp1 = 0;
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
		Pb[i] = (double)(sum)/(30);
		temp1 += Pb[i];
	}
	cout<<"initial Probability vector"<<endl;
	for(int i = 0; i< 10; i++){
		cout<<Pb[i]<<"\t";	//displaying intial probability vector
	}
//	cout<<endl; 
//	cout<<"Sum of intial probability vector "<<temp1<<endl;
//	cout<<endl<<endl;
	cout<<"Global solution :"<<endl;
	/*--------------------end of initial probability vector--------------------------*/
	
	/*------------------------------start of EDA algorithm----------------------------*/
	while(gen < GEN){
		gen++;
		heapSort(FITNESS, Pb, pop, n); //sorting using heapsort
		
		/*-------------------------updating the probability vector-------------------*/
		temp1 = 0;
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
			Nc = (double)(sum)/(30/2);
			Pb[i] = (1-lambda) * Pb[i] + lambda * Nc ;
			temp1 += Pb[i];
		}
//		cout<<endl<<endl;
//		cout<<"Updated Probability vector"<<endl;
//	for(int i = 0; i< 10; i++){
//		cout<<Pb[i]<<"\t";	//displaying intial probability vector
//	}
//	cout<<endl; 
//	cout<<"Sum of updated probability vector "<<temp1<<endl;		
		/*-------------------------end of probability vector updation-----------------*/
		
		/*------------------------generating new solution ----------------------------*/
		int *gbest = new int[p];
		int best_fitness = 0;
		for(int i=0;i<D;i++){
			if(best_fitness < FITNESS[i] ){
				best_fitness = FITNESS[i];
				for(j=0; j<p; j++){
					gbest[j] = pop[i][j];
				}
			}			
		}
		cout<<"Global best Fitness "<<best_fitness<<endl;
		cout<<"Solution :"<<endl;
		for(j=0; j<p; j++){
					cout<<gbest[j]<<" ";
				}
		cout<<endl;
//		cout<<"percentage of coverage :"<<best_fitness/312*100<<endl;
		
		//int** NEW_POP =  new int *[D/2];
		
//		for(i=0; i<D/2; i++){
//			NEW_POP[i] = new int [n];
//		}
		int t = 0;
		for(i=0; i<D/2; i++){
			t = 0;
			vector<int> vec;
			for(j=0; j<n; j++){
				
				r1 = (double)rand() / RAND_MAX;
				r2 = (double)rand() / RAND_MAX;
				if(r1 < beta){
					if(r2 < Pb[j]){
						//NEW_POP[i][t++] = j;
						vec.push_back(j);
					}
				}
				else{
					for(k=0; k<p; k++){
						if( j == gbest[k]){
							//NEW_POP[i][t++] = j;
							vec.push_back(j);
							break;
						}
					}
				}
				
			}
			NEW_POP.push_back(vec);
		}
//		cout<<"Solution from EA/G "<<endl;
//		for(i=0; i<NEW_POP.size(); i++){
//			for(j=0; j<NEW_POP[i].size(); j++){
//				cout<<NEW_POP[i][j]<<" ";
//			}
//			cout<<endl;
//		}
		/*------------------------end of new solution generation-----------------------*/
		
		
		/*-------------------------Repairing----------------------------------------*/
		double p1 = 0;
		flag = 0;
		int l;
		double m1;
		int max_index;
		temp = 0;
		double sum =0;
		vector<int> covered_node;
		for(k=0; k<NEW_POP.size(); k++){
			covered_node.clear();
			for(i=0; i<n; i++){
				flag =0;
				for(j=0; j<NEW_POP[k].size(); j++){
					if(i<NEW_POP[k][j]){
						if(F_LIST[i][NEW_POP[k][j]-i] == 1){
							sum += demand[i];
							covered_node.push_back(i);
							break;
						}
					}
					else{
						if(F_LIST[NEW_POP[k][j]][i-NEW_POP[k][j]] == 1){
							sum += demand[i];
							covered_node.push_back(i);
							break;
						}
					}
				}
			}
			
			int s;
			if(NEW_POP[k].size() < 3){
				s = NEW_POP[k].size();
				for(i=0; i<(3-s); i++){
					m1 = 0;
					for(j=0; j<n; j++){
						//checking if node is already in covered node
						flag = 0;
						for(l=0; l<NEW_POP[k].size(); l++){
							if(j == NEW_POP[k][l]){
								flag =1;
								break;
							}
						}
						if(flag == 1)
							continue;
						else{
							sum = 0;
							for(l=0; l<n; l++){
								flag = 0;
								temp = 0;
								for(int m=0; m<covered_node.size(); m++){
									if( l == covered_node[m]){
										flag =1;
										break;
									}
								}
								if(flag == 1)
									continue;
								else{
									if(l<j){
										if(F_LIST[l][j-l] == 1){
											sum += demand[l];
										}
									}
									else{
										if(F_LIST[j][l-j] == 1){
											sum += demand[l];
										}
									}									
								}
							}							
						}
						if(sum > m1){
							m1 = sum;
							max_index = j;
						}
					}
					//add the node that covers maximum demand
//					cout<<max_index<<endl;
					NEW_POP[k].push_back(max_index);
					for(l=0; l<n; l++){
						flag = 0;
						temp = 0;
						for(int m=0; m<covered_node.size(); m++){
							if( l == covered_node[m]){
								flag =1;
								break;
							}
						}
						if(flag == 1)
							continue;
						else{
							if(l<max_index){
								if(F_LIST[l][max_index-l] == 1){
									covered_node.push_back(l);
								}
							}
							else{
								if(F_LIST[max_index][l-max_index] == 1){
									covered_node.push_back(l);
								}
							}									
						}
					}
				}
			}
			
			
			
		}
		
//		if(NEW_POP[i].size() < 3){
//			for(i=0; i<(3-NEW_POP.size()); i++){
//				
//			}
//		}
//		
//		
//		for(i=0; i<NEW_POP.size(); i++){
//			if(NEW_POP[i].size() < 3){
//				for(i=0; i<n; i++){
//					flag =0;
//					for(j=0; j<NEW_POP[i].size(); j++){
//						if(i<NEW_POP[i][j]){
//							if(F_LIST[i][NEW_POP[i][j]-i] == 1){
//								sum += demand[i];
//								flag = 1;
//								break;
//							}
//						}
//						else{
//							if(F_LIST[pop[j]][i-pop[j]] == 1){
//								sum += demand[i];
//								flag = 1;
//								break;
//							}
//						}
//					}
//				}
//			}
//			
//			
//			
//			
//			if(NEW_POP[i].size() < 3){
//				for(k=NEW_POP[i].size(); k<3; k++){
//					p1 = 0;
//					for(j=0; j<n; j++){
//						flag = 0;
//						for(int l=0; l<NEW_POP[i].size(); l++){
//							if(NEW_POP[i][l] == j){
//								flag =1;
//								break;
//							}
//						}
//						if(flag == 1)
//							break;
//						else{
//							temp = 
//							if(Pb[j] > p1){
//								p1 = Pb[j];
//								temp = j;
//							}
//						}
//					}
//					NEW_POP[i].push_back(temp);
//				}
//			}
//			else{
//				
//			}
//		}
//		cout<<"Solution after applying repairing operator :"<<endl;
//		for(i=0; i<NEW_POP.size(); i++){
//			for(j=0; j<3; j++){
//				cout<<NEW_POP[i][j]<<" ";
//			}
//			cout<<endl;
//		}
		/*------------------------end of repairing----------------------------------*/
		
		/*-------------------------generating new population------------------------*/
		for(i=0;i<D/2; i++){
			for(j=0; j<p; j++){
				pop[i][j] = NEW_POP[i][j];
			}
			FITNESS[i] = fitness_function(F_LIST,pop[i],demand,n,p);
		}
		/*------------------------end of new solution generation-------------------*/
		
		//break;
	}
	/*------------------------------end of EDA algorithm------------------------*/
}
/*--------------------------end of main function----------------------------*/