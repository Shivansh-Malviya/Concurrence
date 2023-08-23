
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

void hz(int n, mat &Hz)
{
/* For Sigma_Z elements in H */

	int i, j;
	mat I, Z, temp;
	 	
	I = {{1, 0}, {0, 1}};											
	Z = {{1, 0}, {0, -1}};;
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
		}

		Hz = Hz + temp;											// Adding each final tensor product
		
	}
     	//cout<<"\n Hz = \n"<<Hz;
}

void hx(int n, mat &Hx)
{
/* For sigma_X */    
	
	int i, j;
	mat I, X, temp, S;
	
	I = {{1, 0}, {0, 1}};											
	X = {{0, 1}, {1, 0}};;
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


void lanczos(mat &H, vec v0, int N, int k)
{
	vec a, b, eigenvalues;
	mat H_f, K, eigenvectors;
	
	int i, j = 0, l = 0;
	float min_ev;
	string filename;	
	ofstream fgev;
	
	filename = "Lan_gev_" + to_string(int(log2(k))) + "_" + to_string(N) + ".dat";
	fgev.open(filename, ios::trunc);
	
	a.zeros(N);
	b.zeros(N);
	H_f.zeros(N, N);
	K.zeros(k, N);													// K_l matrix - N dimensional matrix
	
	//cout<<"\n Flag : " <<j++;
		
	K.col(0) = v0;													// Initial vector	//cout<<"\n Flag : " <<j++;
	
	for(j = 2; j <= N ; j++)											// j acts as a proxy for N
	{
		for(i = j-1; i < j; i++)		//j-1 cuz we already have the data prior to that		// for getting the N dimensional Krylov space.
		{
			K.col(i) = H*K.col(i-1); 
			a(i-1) = dot(K.col(i-1), K.col(i));	
		
			K.col(i) -= a(i-1)* K.col(i-1);
			if(i != 1)
				  K.col(i) -= b(i-1)*K.col(i-2);
			
			b(i) = sqrt(dot(K.col(i), K.col(i))) ;
			K.col(i) = K.col(i)/b(i);					// b_n = <v_n|Hv_(n-1)>, n >= 1. But the indexing starts from 0
		}		
		a(i-1) = dot(K.col(i-1), H*K.col(i-1));
		
		H_f(0, 0) = a(0);
		H_f(0, 1) = b(1);
		
		for(i = 1; i < j-1; i++)										// j(N) different equations relating H|v_n> to a(n) and b(n)
		{
			H_f(i, i) = a(i);
			H_f(i, i-1) = b(i);
			H_f(i, i+1) = b(i+1);
		}
			
		H_f(j-1, j-1) = a(j-1);
		H_f(j-1, j-2) = b(j-2);
		
		//cout << "\n K_l : \n " << K;
		//cout << "\n Tridiagonal Mat : \n\n" << H_f;
		
		eig_sym(eigenvalues, eigenvectors, H_f);
		H_f = diagmat(eigenvalues);
		min_ev = ( eigenvalues(0) < eigenvalues(j-1) )? eigenvalues(0) : eigenvalues(j-1);
		//cout << "\n Diagonal Mat : \n\n" << H_f;
		
		//cout << "\n N, min gev : " << j << " " << min_ev << endl;	
		fgev << j << " " << min_ev << endl;
	}
	
	fgev.close();
	
	/* 
     //cout<<v0;
     v.col(0)=v0;
     N[0]= sqrt(dot(v0,v0));
     vec v1= H*v0;
     a[0] =dot(v0,v1);
     v1=v1-a[0]*v0;
     //cout<<v1;
     N[1]=dot(v1,v1);
     //cout<<N[1]<<"Hello";
     N[1]=sqrt(N[1]);
     v1=v1/(N[1]);
     v.col(1)=v1;
     vec vtemp;
     double q;
    // cout<<N[1];
     for(i=2;i<=k-1;i++)
     {
       v.col(i)=H*v.col(i-1);
       a[i-1]=dot(v.col(i-1),v.col(i));
       v.col(i)=v.col(i)-a[i-1]*v.col(i-1)-N[i-1]*v.col(i-2);
       N[i]=sqrt(dot(v.col(i),v.col(i)));
       
       v.col(i)=v.col(i)/(N[i]);
      // vtemp.zeros(k);
      
      }
      a[k-1]=dot(v.col(k-1),v.col(k-2));
       for(i=0;i<k;i++)
     {
      for(j=0;j<k;j++)
      {
        if(i==j)
         H_diag(i,i)=a[i];
        else if(i==j-1)
         H_diag(i,j)=N[j];
        else if(i==j+1)
         H_diag(i,j)=N[i];
       }
      }
      vec eigval;
      mat eigvec;
      mat diagonal;//(k,k,fill::zeros);
      eig_sym(eigval,eigvec,H_diag);
      diagonal = diagmat(eigval);
      */
	
}


int main()
{	
	mat H, Hx, Hz;
	vec v0;
	int N, n, k;
	float J = 1.0, h = 1.0;
	
	system("clear");
	
	cout<<"\n Input the dimension of the spin chain : ";
	cin>>n;
	
	cout << "\n Input the number of iterations : ";
	cin >> N;
	
L :	k = pow(2, n);
	cout<< "\n k = " << k;
	
	H.zeros(k, k);
	Hx.zeros(k, k);
	Hz.zeros(k, k);	
	
	v0.zeros(k);
	v0(0) = 1;
	
	auto start = std::chrono::high_resolution_clock::now();
	
	hx(n, Hx);
	hz(n, Hz);
	
	auto end = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration = end - start;
	cout << "\n Execution time: " << duration.count() << " seconds" << endl;
	
	H = J*Hx - h*Hz;
	cout<<"\n Hamiltonian : \n"<<H;
	
	lanczos(H, v0, N, k);
	
	n = n-2;
	
	if(n >= 4)
	{
		cout << "\n n = " << n;	
		goto L;
	}
	
	return 0;
}































