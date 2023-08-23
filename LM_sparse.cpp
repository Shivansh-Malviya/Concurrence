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
	sp_mat I(2,2), Z(2, 2), temp(2, 2);
	 	
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

void lanczos(sp_mat &H, sp_vec v0, int N, int k)
{
	int i, j = 0, l = 0;
	float min_ev, eigtemp = 0;
	string filename;	
	ofstream fgev;
	
	filename = "LG_sp" + to_string(int(log2(k))) + "_" + to_string(N) + ".dat";
	fgev.open(filename, ios::trunc);												
	
	//cout<<"\n Flag : " <<l++;

	for(j = 2; j <= N ; j++)											// j acts as a proxy for N
	{
		sp_mat H_f(j, j);									
		sp_vec a(j), b(j), v(j), u(j), t(j);
		vec eigenvalues;
		mat eigenvectors;	
		
		//cout<<"\n Flag : " <<l++;	
		//cout<<"\n Flag : " <<l++;	
		
		for(i = 1; i < j; i++)		//j-1 cuz we already have the data prior to that		// for getting the N dimensional Krylov space.
		{
			if(i == 1)
				u = v0;											// Initial vector	//cout<<"\n Flag : " <<j++;

			v = H*u; 
			a(i-1) = dot(u, v);	
		
			//cout<<"\n Flag : " <<l++;
		
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
		
		//cout<<"\n Flag : " <<l++;
		
		
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
		
		
		//cout << "\n K_l : \n " << K;
		//cout << "\n Tridiagonal sp_mat : \n\n" << H_f;
		
		eigs_sym(eigenvalues, eigenvectors, H_f, j-1, "sa");
		//H_f = diagmat(eigenvalues);
		min_ev = ( eigenvalues(0) < eigenvalues(j-2) )? eigenvalues(0) : eigenvalues(j-2);
		//cout << "\n Diagonal sp_mat : \n\n" << H_f;
		
		//min_ev = eigenvalues(0);
		
		if( abs(eigtemp - min_ev) > 0.001 )
		{
			eigtemp = min_ev;
			fgev<<j<<"	"<<eigtemp<<endl;
		}
		 
		else
			break;
		
		//cout << "\n N, min gev : " << j << " " << min_ev << endl;	
		//fgev << j << " " << min_ev << endl;
	}
	
	cout<<"\n Converged for N = " << j;
	
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
      sp_mat eigvec;
      sp_mat diagonal;//(k,k,fill::zeros);
      eig_sym(eigval,eigvec,H_diag);
      diagonal = diagmat(eigval);
      */
	
}


int main()
{	
	int N, n, k, l = 0;
	float J = 1.0, h = 1.0;
	
	system("clear");
	
	cout<<"\n Input the dimension of the spin chain : ";
	cin>>n;
	
	cout << "\n Input the number of iterations : ";
	cin >> N;
	
L :	k = pow(2, n);
	cout<< "\n k = " << k;
	
	sp_mat H(k, k), Hx(k, k), Hz(k, k);
	sp_vec v0(k);
	
	v0(0) = 1;
	cout<< "\n initial vec = " << v0;
	
	auto start = std::chrono::high_resolution_clock::now();
	
	hx(n, Hx);
	hz(n, Hz);
	
	auto end = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration = end - start;
	cout << "\n Execution time: " << duration.count() << " seconds" << endl;
	
	H = J*Hx - h*Hz;
	//cout<<"\n Hamiltonian : \n"<<H;
	
	lanczos(H, v0, N, k);
	
	n = n-2;														// Automation
	//n = 2;
	if(n >= 4)
	{
		cout << "\n n = " << n;	
		goto L;
	}
	
	return 0;
}

//p 'LG_sp4_50.dat' w line, 'LG_sp6_50.dat' w line, 'LG_sp8_50.dat' w line, 'LG_sp10_50.dat' w line, 'LG_sp12_50.dat' w line, 'LG_sp14_50.dat' w line, 'LG_sp16_50.dat' w line,  'LG_sp18_50.dat' w line, 'LG_sp20_50.dat' w line






























