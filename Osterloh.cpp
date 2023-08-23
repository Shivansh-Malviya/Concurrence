
#include<iostream>
#include<armadillo>
#include<cmath>
#include<complex>
#include<string>

#define QICLIB_DONT_USE_NLOPT
#include "/home/frostt/Desktop/QIClib/include/QIClib"

using namespace std;
using namespace arma;
using namespace qic;

#define cout std::cout


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
		//cout<< "\n Hz, first loop " << i;
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
					
			//cout<< "\n j : " << j << endl;
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
		//cout<< "\n Hx, first loop " << i;
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


cx_mat trace(int k, cx_vec psi)
{
	int r = 4;
	cx_mat rho12;
	rho12.zeros(4, 4);
	
	for(int i = 0; i < (k/r); i++)
	{
		for(int j = 0; j < r; j++)
		{                
			for(int l = 0; l < r; l++)
				rho12((i*r + j)%r , (i*r + l)%r) += psi(i*r + j, 0)*psi(i*r + l, 0);                               
		}
	}
	
	return rho12; 

}





cx_vec lanczos(sp_mat &H, const sp_vec v0, int N, int k)
{
	int i, j = 0, l = 0;
	float min_ev, eigtemp = 0;
	//string filename;	
	//ofstream fgev;
	cx_vec psi;
	vec eigenvalues;
	mat eigenvectors;
	//filename = "LG_sp" + to_string(int(log2(k))) + "_" + to_string(N) + ".dat";
	//fgev.open(filename, ios::trunc);												
	
	//cout<<"\n Flag : " <<l++;

	for(j = 2; j <= N ; j++)											// j acts as a proxy for N
	{
		sp_mat H_f(j, j);									
		sp_vec a(j), b(j), v(j), u(j), t(j);
			
		cout<<"\n Flag : " <<l++;	
		//cout<<"\n Flag : " <<l++;	
		
		for(i = 1; i < j; i++)		//j-1 cuz we already have the data prior to that		// for getting the N dimensional Krylov space.
		{
			if(i == 1)
				u = v0;											// Initial vector	//cout<<"\n Flag : " <<j++;

			v = H*u; 
			a(i-1) = dot(u, v);	
		
			cout<<" " <<i;		
			//cout<<" " <<i;
		
			v -= a(i-1) * u;
			if(i != 1)
				  v -= b(i-1)*t;
			
			//cout<<"\n Flag : " <<l++;
			
			b(i) = sqrt(dot(v, v)) ;
			v = v/b(i);											// b_n = <v_n|Hv_(n-1)>, n >= 1. But the indexing starts from 0
		
			t = u;
			u = v;
		}		
		a(i-1) = dot(u, H*u);
		
		cout<<" Hf ";
		
		
		H_f(0, 0) = a(0);
		H_f(0, 1) = b(1);
		
		for(i = 1; i < j-1; i++)										// j(N) different equations relating H|v_n> to a(n) and b(n)
		{
			H_f(i, i-1) = b(i);
			H_f(i, i) = a(i);
			H_f(i, i+1) = b(i+1);
		}
			
		H_f(j-1, j-2) = b(j-1);
		H_f(j-1, j-1) = a(j-1);
		
		//cout << "\n Tridiagonal sp_mat : \n\n" << H_f;
		
		eigs_sym(eigenvalues, eigenvectors, H_f, 1, "sa");
		//H_f = diagmat(eigenvalues);
		//min_ev = ( eigenvalues(0) < eigenvalues(j-2) )? eigenvalues(0) : eigenvalues(j-2);
		//cout << "\n Diagonal sp_mat : \n\n" << H_f;

		min_ev = eigenvalues(0);
		
		if( abs(eigtemp - min_ev) > 0.001 )
		{
			eigtemp = min_ev;
			//fgev<<j<<"	"<<eigtemp<<endl;
		}
		 
		else
		{
			psi = conv_to<cx_vec>::from(eigenvectors.col(0));
			break;
		}
		//cout << "\n N, min gev : " << j << " " << min_ev << endl;	
		//fgev << j << " " << min_ev << endl;
	}
	
	if(j == N+1)
	{
		cout << "\n Not converged for the given N!";
	}
	cout<<"\n Converged for N = " << j;
	
	irlm(H, H_f, min_ev, v0, psi, N);
	
	return psi;
	//fgev.close();
}





int main()
{

	system("clear");
	
	unsigned int i, j, n, k, p, N;
	unsigned int min, num;
	float J = 1.0, h = 1.0;
	float concurrence, derivative;
	char choice;
	string filename;
	
	cout<<"\n\n\n Enter number of spins(>= 3) : ";
	cin>>n;
	cout << "\n Input the limit of iterations : ";
	cin >> N;
	
L :	k = pow(2, n);													// Dimension of Hilbert space	
	cout<<"\n\n k = " << k << endl;
	
	ofstream fsp;
	
	vec eigenvalues;
	sp_vec v0(k);
     uvec subsystemInd, sortedIndices;
     vector<unsigned int> subsystemind;
     vector<float> conc, dconc, lambda;
	
	cx_mat rho12, rhot, con, eigenvectors, rho, Y;
	//sp_vec psi(k);
	cx_vec psi;
	sp_mat H(k, k), rho12_sp(4, 4), Hz(k, k), Hx(k, k);	
	mat eigen_vectors;
	
	filename = "stl_con" + to_string(n) + ".dat";
	fsp.open(filename, ios::trunc);
	//cout<< "\n Filename : "<<filename;

	hx(n, Hx);
	hz(n, Hz);									
    
	auto start = std::chrono::high_resolution_clock::now();
	
/* Concurrence */
	
	Y = spm.S(2);	
	
	for (h = 0.1; h <= 3; h = h + 0.05)								// Running h from 0.1 to 3 with intervals of 0.05
	{	
		v0(0) = 1;	
    		H = h*Hz - J*Hx;											// Defining H with J = 1 and h as parameter
			
		psi = lanczos(H, v0, N, k);
		cout<<"\n Psi recieved "<<psi;
		//cout<<"\n Psi recieved ";
    		//eigs_sym(eigenvalues, eigen_vectors, H, 1, "sa");				// Getting GS Eigenvector and Eigenvalue for the corresponding value of h
		//psi = conv_to<cx_vec>::from(eigen_vectors.col(0));				
		
		rho12 = trace(k, psi);	
		cout<<"\n Rho_12 calculated ";	
		cout<<"\n Rho_12 calculated ";						
		rhot = kron(Y, Y) * rho12 * kron(Y, Y);							// Used in the formula
		con = sqrtmat(sqrtmat(rho12) * rhot * sqrtmat(rho12));				// The formula for concurrence
		
		eig_sym(eigenvalues, eigenvectors, con);						// Calculating the eigenvectors and eigenvalues of "con" matrix
   			
  		concurrence = eigenvalues(3) - eigenvalues(2) - eigenvalues(1) - eigenvalues(0);	// because eigenvalues in ascending order			
   						
		if(concurrence < 0)
	    		concurrence = 0;										// By definition			

		fsp<<J/h<<"   "<<concurrence<<endl;							// Store values in the file
	
		//cout<<"\n Flag j"<<i<<endl;
		//cout<<"\n Flag k"<<i<<endl;
		
		conc.push_back(concurrence);
      	lambda.push_back(J/h); 
    			    			
	}
    
   	fsp.close();
   	
   	cout<<"\n Concurrence data for "<<n<<" spin chain written to : "<<filename;


/* Derivative Block */
	filename = "stl_dcon" + to_string(n) + ".dat";
	fsp.open(filename, ios::trunc);
	
	for(i = 1; i < conc.size(); i++)									// 59 elements in conc which will yield 58 averages. But the index goes from 0 to 58
	{
		derivative = (conc[i] - conc[i-1]) / (lambda[i] - lambda[i-1]);
		dconc.push_back(derivative);
	}
	 
	i = 0; 
	for(h = 0.1; h < 2.95 ; h = h + 0.05)								
		fsp << (J/h + J/(h+0.05)) / 2 << "   "<<dconc[i++]<<endl;			// Taking average of lambda because we calculated the derivative at the midpoint
	  
	cout<< "\n Data for derivative of concurrence written to : "<<filename; 
	fsp.close();
	
	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> duration = end - start;
	
	cout << "\n Execution time: " << duration.count() << " seconds" << endl;

/* Choice Block 														
M:   cout<<"\n\n Do you want to continue? [y/n] : ";
   	cin>>choice;
   	
   	if(choice == 'y')
   		goto L;
   		
   	else if(choice != 'n' && choice != 'y')
   	{
   		cout<<"\n Invalid Input! Choose again. ";
   		goto M;
   	}
*/

	n = n + 2;
	goto L;
   	
}
