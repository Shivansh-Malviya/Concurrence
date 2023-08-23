#include<iostream>
#include<armadillo>
#include<cmath>
#include<complex>
#include<string>
#include<chrono>

#define QICLIB_DONT_USE_NLOPT
#include "/home/frostt/Desktop/QIClib/include/QIClib"

using namespace std;
using namespace arma;
using namespace qic;

void hz(int n, sp_mat &Hz)
{
/* For Sigma_Z elements in H */

	int i, j;
	sp_mat I(2, 2), Z(2, 2), temp(2, 2);
	 	
	I(0, 0) = 1;								
	I(1, 1) = 1;
	Z(0, 0) = 1;
	Z(1, 1) = -1;
	temp = I;
	
	for(i = 1; i <= n; i++)					
	{	
		temp = I;

		for(j = 1; j <= n; j++)										// For 1 entire tensor product in the summation
		{
			if(i == j && j == 1)									// 1st element of the 1st tensor product
				temp = Z;

		  	else if(j > 1)
		     {
				if(i == j)										// Position of the Z matrix
					temp = kron(Z, temp);

				else												// Filling the rest with I	
					temp = kron(I, temp);
			}

			else
					temp = I;
					
			//cout<< "i, j : " << i << " " << j << endl;
		}

		Hz = Hz + temp;											// Adding each final tensor product
		
	}
     	//cout<<"\n Hz = \n"<<Hz;
}

void hx(int n, sp_mat &Hx)
{
/* For sigma_X */    
	
	int i, j;
	sp_mat I(2,2), X(2, 2), temp, S;
	
	I(0, 0) = 1;								
	I(1, 1) = 1;
	X(0, 1) = 1;
	X(1, 0) = 1;
	temp = I;
	
	S = kron(X, X);
	
	for(i = 1; i < n; i++)				
	{
		temp = I;

		for(j = 1; j < n; j++)										// different tensor-product elements of the summation
		{
			if(i == j && j == 1)									// For the first 2 elements of the first tensor product
				temp = S;

			else if(j > 1)
			{
				if(i == j)
					temp = kron(S, temp);
			
				else
					temp = kron(I, temp);
			}

			else
				temp = I;
		}

		Hx = Hx + temp;
	}
	
		//cout<<"\n Hx = \n"<<Hx;


     //PBC
     temp = X;

     for(j = 2; j <= n - 1; j++)
		temp = kron(temp, I);

	temp = kron(temp, X);

	Hx += temp;

     //cout<<"\n Hx = \n"<<Hx;
}

float Conc(sp_mat &H, int n)
{
	unsigned int min;
	cx_mat rho12(4, 4), rhot(4, 4), Y(2, 2);
	cx_mat con, rho;
	cx_vec psi;
	vec eigenvalues;
	mat eigen_vectors;
	uvec subsystemInd, sortedIndices;
	
	vector<unsigned int> subsystemind;
	
	Y(0, 0) = 0;	Y(1, 1) = 0;
	Y(0, 1) = -complex<double>(0, 1);
	Y(1, 0) = complex<double>(0, 1);
	
	// cout<<"\n Y = " << Y;					// ??
	
	eigs_sym(eigenvalues, eigen_vectors, H, 1, "sa");						// Getting Eigenvectors and Eigenvalues for the corresponding value of h
		
	//min = eigenvalues.index_min();									// to get the index of the ground state in Hamiltonian
	cout<< "\n Min = " << min;
	psi = conv_to<cx_vec>::from(eigen_vectors.col(0));							
	rho = psi * psi.t();
	
	for(int num = 3; num <= n; num++)
		subsystemind.push_back(num);									// To append the indices of the subsystem to vec subsystemind; starting from 3
		
	subsystemInd = conv_to<uvec>::from(subsystemind);	
	
	rho12 = TrX(rho, subsystemInd);									// Traces rho over the subsystem subsystemInd. A function from QIClib.
	rhot = kron(Y, Y) * rho12 * kron(Y, Y);								// Used in the formula
	con = sqrt(sqrt(rho12) * rhot * sqrt(rho12));						// The formula for concurrence
	
	cx_mat eigenvectors;
	
	eig_sym(eigenvalues, eigenvectors, con);							// Calculating the eigenvectors and eigenvalues of "con" matrix
	
	/*sortedIndices = arma::sort_index(eigenvalues, "descend");				// Sorting eigenvalues in descending order
	eigenvalues = eigenvalues(sortedIndices);
	eigenvectors = eigenvectors.cols(sortedIndices);*/
  			
	float concurrence = eigenvalues(3)-eigenvalues(2)-eigenvalues(1)-eigenvalues(0);	// because eigenvalues in desceding order			
  						
	if(concurrence < 0)
    		concurrence = 0;											// By definition			
   		
	return concurrence;
}

int main()
{	
	int n, k;
	cout<<"\n Input the dimension of the spin chain : ";
	cin>>n;
	
	//cout << "\n Input the number of iterations : ";
	//cin >> N;
	
L :	k = pow(2, n);
	cout<< "\n k = " << k;
	
	int N, l = 0;
	float J = 1.0, h = 1.0, min_ev, concurrence, lambda;
	
	vec eigenvalues;
	mat eigenvectors;	
	
	system("clear");
	
	string filename = "sp_con" + to_string(n) + ".dat";
	ofstream fspcon;
	fspcon.open(filename, ios::trunc);
	cout<<"\n Filename : "<<filename;
	
	sp_mat H(k, k), Hx(k, k), Hz(k, k);
	sp_vec v0(k);
	
	//v0(0) = 1;
	//cout<< "\n initial vec = " << v0;
	
	auto start = std::chrono::high_resolution_clock::now();
	
	hx(n, Hx);
	hz(n, Hz);
	
	auto end = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration = end - start;
	cout << "\n Execution time: " << duration.count() << " seconds" << endl;

// Concurrence	
	for (h = 0.3; h <= 3; h = h + 0.05)								// Running h from 0.1 to 3 with intervals of 0.05
	{
    		H = h*Hz - J*Hx;											// Defining H with J = 1 and h as parameter
			
		concurrence = Conc(H, n);
		
		lambda = J/h;
		fspcon<<lambda<<"   "<<concurrence<<endl;						// Store values in the file
	
		//cout<<"\n Flag j"<<i<<endl;
		//cout<<"\n Flag k"<<i<<endl;
		cout<< "itr : "<<20*h - 1 << endl;
    			    			
	}
    
	cout<<"\n Concurrence data for "<<n<<" spin chain written to : "<<filename;

   	fspcon.close();
   	
	n = n-2;														// Automation
	//n = 2;
	if(n >= 4)
	{
		cout << "\n n = " << n;	
		goto L;
	}
	
	return 0;
}

// p 'sp_con4.dat' w l, 'sp_con6.dat' w l, 'sp_con8.dat' w l, 'sp_con10.dat' w l, 'sp_con12.dat' w l, 'sp_con14.dat' w l, 'sp_con16.dat' w l, 'sp_con18.dat' w l, 'sp_con20.dat' w l 


