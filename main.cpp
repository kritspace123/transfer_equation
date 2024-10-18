#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

double u(double x,  double t){
	return pow(x+0.1, 2) - (sin(2*M_PI*t))/2+x - 3.5*t;
}

double f(double x, double t){
	return -M_PI*cos(2*M_PI*t) - 3.164+0.56*x;
}
double fi(double x){
	return pow(x+0.1, 2) + x;//u(x, 0)
}
double g1(double t){
	return pow(0.1, 2) - 0.5 * (sin(2*M_PI*t)) - 3.5*t ; //  u(0, t)
}

int main(){
	int N, J;
	ofstream file("value.txt");
	ofstream fout("value_function.txt");


	double tao = 0.01;
	double h;
	double a = 0.28;

	cout << "h>"<< a * tao << endl;
	cin >> h;

	N = 1 / h;
	J = 1 / tao;

	double A = (a*tao) / h;
	double B  = 1 - A;	

	double* old_layer = new double[N+1];
	double* new_layer = new double[N+1];

	double fun, raz;
	double maxim = 0;

	//Вычисляем значения функции на НУЛЕВОМ слое
	for(int i = 0; i < N+1; i++){
		old_layer[i] = fi(h*i);
		file << old_layer[i] << " ";
		fout << u(h*i, 0) << " ";
	}
	file << endl;
	fout << endl;
	
	//Вычисляем значения функции на остальных слоях...
	for(int j = 1; j < J+1; j++){
		//Вычисляем значение функции при x = 0
		new_layer[0] = g1(j*tao);
		file << new_layer[0] << " ";
		fout << u(0, tao*j) << " ";
		
		//вычисляем остальное
		for (int k = 1; k < N+1; k++){
			new_layer[k] = B * old_layer[k] + A * old_layer[k-1] + f(k*h, (j-1)*tao)*tao;
			fun = u(h*k, tao*j);

			raz = abs(fun - new_layer[k]);
			if(raz > maxim){ 
				maxim = raz;
			}

			file << new_layer[k] << " ";
			fout << fun << " ";
		}
		file << endl;
		fout << endl;
		for(int t = 0; t < N+1; t++){
			old_layer[t] = new_layer[t];
		}
	}

	cout << maxim << endl;
	cout << tao + h << endl;

	// double old, now1, fun, s;
	// double maxim = 0;
	
	// for(int j = 0; j <= N; j++){
	// 	if(j == 0){
	// 		for(int k = 0; k <= N; k++){ 
	// 			array[k] = fi(k*h);
	// 			file << array[k] << " ";
	// 			fout << u(h*k, j*tao)<< " ";
	// 		}
	// 		file << endl;
	// 		fout << endl;
	// 		continue;
	// 	}
	// 	for(int k = 0; k <= N; k++){
	// 		if(k == 0){
	// 			old = array[0];
	// 			array[0] = g1(j*tao);
	// 			file << array[0] << " ";
	// 			fout << u(0, j*tao)<< " ";
	// 			continue;
	// 		}
	// 		now1 = B * array[k] + A * old + f(h*k, tao*(j-1));
	// 		old = array[k];
	// 		array[k] = now1;
	// 		fun = u(h*k, j*tao);

	// 		fout << fun << " ";
	// 		file << now1 << " ";

	// 		s = abs(fun - now1);
	// 		if (s > maxim){
	// 			maxim = s;
	// 		}
	// 	}
	// 	file << endl;
	// 	fout << endl;
	// }
	// cout << maxim<< endl;
	delete[] old_layer;
	delete[] new_layer;
}