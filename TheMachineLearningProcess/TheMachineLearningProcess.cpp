#include "pch.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <time.h>

using namespace std;

void Read(double** data, string name);

void kNN(double** data, double** data2, const int columns, const int rows, const int rows2, const int k);

vector <vector<int> > Save_classes(double** data, const int columns, const int rows);

vector <double> Distance(string heuristics, double ** data, double ** data2, const int columns, int i, const int rows2);

void Bayes(double** data, double** data2, const int columns, const int rows, const int rows2);

int main()
{
	//Read data from files
	const int columns = 15;
	const int rows = 345;
	const int rows2 = 345;

	double **data;

	data = new double*[columns];

	for (int i = 0; i < columns; i++)
	{
		data[i] = new double[rows];
	}

	Read(data, "australian_TST.txt");

	double **data2;

	data2 = new double*[columns];

	for (int i = 0; i < columns; i++)
	{
		data2[i] = new double[rows2];
	}

	Read(data2, "australian_TRN.txt");


	//Determining the data
	const int k = 2;
	
	//Calling the algorithm
	//double sorted_distances[rows2]; //required to use a constant to allocate an array
	kNN(data, data2, columns, rows, rows2, k);

	//Calling the algorithm
	Bayes(data, data2, columns, rows, rows2);


	for (int i = 0; i < columns; i++)
	{
		delete[] data[i];
		delete[] data2[i];
	}
	delete[] data;
	delete[] data2;
}

void Read(double** data, string name)
{
	fstream file;
	file.open(name, ios::in);

	if (file.good())
	{
		string line; //line from file
		string digit; //single value from file

		int number_of_row = 0;

		while (!file.eof())
		{
			getline(file, line);

			int counter = 0;

			for (int i = 0; i <= line.length(); i++)
			{
				if (line[i] == ' ' || line[i] == '\0')
				{
					double number = stod(digit);

					data[counter][number_of_row] = number;

					digit = "";
					counter++;
				}
				else
				{
					digit += line[i];
				}
			}

			number_of_row++;
		}
	}
	else
	{
		cout << "File error. Check the catalog with the project." << endl;
		system("PAUSE");
		exit(0);
	}

	file.close();
}

void kNN(double** data, double** data2, const int columns, const int rows, const int rows2, const int k)
{
	cout << "kNN" << endl << endl;

	vector <double> distances;

	//Save available classes into vector
	vector< vector<int> > classes = Save_classes(data, columns, rows); //{class, number of occurrences}
	
	//Finding the shortest distance
	//string heuristics = "Euclid";
	vector <string> heuristics = { "Euclid", "Canberr", "Chebyshev", "Manhattan", "Pearson" };

	for (int a = 0; a < heuristics.size(); a++)
	{
		if(heuristics[a] != "Pearson")
			cout << "Heuristics: " << heuristics[a] << endl << endl;
		else cout << "Heuristics: Absolute Pearson correlation coefficient" << endl << endl;

		vector<double> acc;
		vector<double> cov;
		vector<double> tpr;

		vector <double> sorted_distances;
		vector <double> sum;
		double partial_sum = 0;
		vector <int> partial_classification; //part of classification vector { decision(from system)/not catch(number -1), (if catch) true/false classification }
		vector <vector <int> > classification; //classification for objects

		for (int i = 0; i < rows; i++)
		{
			sum.clear();
			distances.clear();

			int decision = data[columns - 1][i];

			distances = Distance(heuristics[a], data, data2, columns, i, rows2);

			for (int j = 0; j < classes.size(); j++)
			{
				sorted_distances.clear();

				for (int l = 0; l < rows2; l++)
				{
					if (data2[columns - 1][l] == classes[j][0])
					{
						sorted_distances.push_back(distances[l]);
					}
				}

				sort(sorted_distances.begin(), sorted_distances.end());

				partial_sum = 0;

				for (int m = 0; m < k; m++)
				{
					partial_sum += sorted_distances[m];
				}

				sum.push_back(partial_sum);
			}

			//cout << endl;
			//for (int j = 0; j < sum.size(); j++)
			//	cout << sum[j] << endl;
			//cout << endl;

			int number_of_index = 0;
			double check_distance = 1000;
			partial_classification.clear();

			for (int j = 0; j < sum.size(); j++)
			{
				if (sum[j] < check_distance)
				{
					check_distance = sum[j];
					number_of_index = j;
				}
				if (sum[0] == sum[1] && partial_classification.size() == 0)
					partial_classification.push_back(-1);
			}

			if (partial_classification.size() > 0)
			{
				classification.push_back(partial_classification);
			}
			else
			{
				partial_classification.push_back(classes[number_of_index][0]);
			}

			if (partial_classification[0] != -1 && partial_classification[0] == decision)
			{
				partial_classification.push_back(1);
				classification.push_back(partial_classification);
			}
			else if (partial_classification[0] != -1 && partial_classification[0] != decision)
			{

				partial_classification.push_back(0);
				classification.push_back(partial_classification);
			}
		}

		//cout << endl;
		//for (int i = 0; i < classification.size(); i++)
		//{
		//	for (int j = 0; j < classification[i].size(); j++)
		//		cout << classification[i][j] << " ";
		//	cout << endl;
		//}

		//--------------Create raport
		int number_of_true_classificated_objects;
		int number_of_catch_objects;
		int number_of_objets_wrongly_ending_in_our_class;

		for (int j = 0; j < classes.size(); j++)
		{
			number_of_true_classificated_objects = 0;
			number_of_catch_objects = 0;
			number_of_objets_wrongly_ending_in_our_class = 0;

			for (int i = 0; i < classification.size(); i++)
			{
				if (classification[i].size() > 1)
				{
					if (classification[i][1] == 1 && classification[i][0] == classes[j][0])
						number_of_true_classificated_objects++;

					if (classification[i][1] != 1 && classification[i][0] == classes[j][0])
						number_of_objets_wrongly_ending_in_our_class++;
				}
			}

			for (int i = 0; i < rows; i++)
			{
				if (data[columns - 1][i] == classes[j][0] && classification[i][0] != -1)
					number_of_catch_objects++;
			}

			//cout << "number_of_true_classificated_objects: " << number_of_true_classificated_objects << endl;
			//cout << "number_of_catch_objects: " << number_of_catch_objects << endl;
			//cout << "number_of_objets_wrongly_ending_in_our_class: " << number_of_objets_wrongly_ending_in_our_class << endl;

			acc.push_back((double)number_of_true_classificated_objects / (double)number_of_catch_objects);
			cov.push_back((double)number_of_catch_objects / (double)classes[j][1]);
			tpr.push_back((double)number_of_true_classificated_objects / ((double)number_of_true_classificated_objects + (double)number_of_objets_wrongly_ending_in_our_class));

			cout << "acc_" << classes[j][0] << ": " << acc[j] << endl;
			cout << "cov_" << classes[j][0] << ": " << cov[j] << endl;
			cout << "TPR_" << classes[j][0] << ": " << tpr[j] << endl << endl;
		}

		number_of_true_classificated_objects = 0;
		number_of_catch_objects = 0;
		double acc_global = 0;
		double cov_global = 0;

		for (int i = 0; i < classification.size(); i++)
		{
			if (classification[i].size() > 1)
			{
				number_of_catch_objects++;

				if (classification[i][1] == 1)
					number_of_true_classificated_objects++;
			}
		}

		acc_global = (double)number_of_true_classificated_objects / (double)number_of_catch_objects;
		cov_global = (double)number_of_catch_objects / (double)rows;

		cout << "acc_global: " << acc_global << endl;
		cout << "cov_global: " << cov_global << endl;
		cout << endl << endl;
	}
}

vector <vector<int> > Save_classes(double** data, const int columns, const int rows)
{
	vector< vector<int> > classes;

	bool flag = false;
	int counter = 0;
	for (int i = 0; i < rows; i++)
	{
		flag = false;
		for (int j = 0; j < classes.size(); j++)
		{
			if (classes[j][0] == (int)data[columns - 1][i])
				flag = true;
		}

		if (flag == false)
		{
			counter = 0;
			for (int j = 0; j < rows; j++)
			{
				if (data[columns - 1][j] == data[columns - 1][i])
					counter++;
			}
			classes.push_back({ (int)data[columns - 1][i], counter });
		}
	}

	sort(classes.begin(), classes.end());

	return classes;
}

vector <double> Distance(string heuristics, double ** data, double ** data2, const int columns, int i, const int rows2)
{
	vector <double> distances;

	//We count the distance for each object
	if (heuristics == "Euclid")
	{
		double sum = 0;

		for (int j = 0; j < rows2; j++)
		{
			sum = 0;

			for (int k = 0; k < columns-1; k++) //data[k][i] -> data2[k][j]
			{
				sum += pow(data[k][i] - data2[k][j], 2);
			}

			distances.push_back(sqrt(sum));
		}
	}
	else if (heuristics == "Canberr")
	{
		double sum = 0;

		for (int j = 0; j < rows2; j++)
		{
			sum = 0;

			for (int k = 0; k < columns - 1; k++)
			{
				sum += fabs((data[k][i] - data2[k][j]) / (data[k][i] + data2[k][j]));
			}

			distances.push_back(sum);
		}
	}
	else if (heuristics == "Chebyshev")
	{
		vector<double> sum;
		double maximum = 0;

		for (int j = 0; j < rows2; j++)
		{
			for (int k = 0; k < columns - 1; k++)
			{
				sum.push_back(fabs(data[k][i] - data2[k][j]));
			}

			distances.push_back(*max_element(sum.begin(), sum.end()));

			sum.clear();
		}
	}
	else if (heuristics == "Manhattan")
	{
		double sum = 0;

		for (int j = 0; j < rows2; j++)
		{
			sum = 0;

			for (int k = 0; k < columns - 1; k++)
			{
				sum += fabs((data[k][i] - data2[k][j]));
			}

			distances.push_back(sum);
		}
	}
	else if (heuristics == "Pearson")
	{
		double average_x;
		double average_y;
		double denominator_x;
		double denominator_y;
		double r_xy;
		double d_xy;

		for (int j = 0; j < rows2; j++)
		{
			average_x = 0;
			average_y = 0;
			denominator_x = 0;
			denominator_y = 0;
			r_xy = 0;
			d_xy = 0;

			for (int k = 0; k < columns - 1; k++)
			{
				average_x += data[k][i];
				average_y += data2[k][j];
			}
			average_x /= columns - 1;
			average_y /= columns - 1;

			for (int k = 0; k < columns - 1; k++)
			{
				denominator_x += pow(data[k][i] - average_x, 2);
				denominator_y += pow(data2[k][j] - average_y, 2);
			}
			
			denominator_x /= (double)columns - 1;
			denominator_y /= (double)columns - 1;

			denominator_x = sqrt(denominator_x);
			denominator_y = sqrt(denominator_y);

			//cout << denominator_x;

			for (int k = 0; k < columns - 1; k++)
			{
				r_xy += ((data[k][i] - average_x) / denominator_x) * ((data2[k][j] - average_y) / denominator_y);
			}

			r_xy /= columns - 1;

			d_xy = fabs(r_xy);

			distances.push_back(1 - d_xy);
		}
	}

	//cout << "Distances: " << distances[0] << " " << distances[1] << " " << distances[2] << " " << distances[3] << " " << distances[4] << " " << distances[5];

	return distances;
}

void Bayes(double** data, double** data2, const int columns, const int rows, const int rows2)
{
	cout << "Naive Bayes classifier" << endl;
	//Save available classes into vector
	vector <vector<int> > classes = Save_classes(data2, columns, rows2); //{class, number of occurrences}

	//Determinig the data
	int counter;
	double sum;
	double class_size;
	vector<double> param_c;
	int decision = -1;
	vector <vector <int> > classification; //classification for objects { decision(from system), true/false classification }

	for (int i = 0; i < rows; i++)
	{
		param_c.clear();
		for (int j = 0; j < classes.size(); j++)
		{
			class_size = (double)classes[j][1] / (double)rows2;
			sum = 0;

			for (int k = 0; k < columns - 1; k++)
			{
				counter = 0;
				for (int l = 0; l < rows2; l++)
				{
					if (data2[k][l] == data[k][i] && (int)data2[columns - 1][l] == classes[j][0])
						counter++;
				}
				sum += (double)counter / (double)classes[j][1];
			}
			
			param_c.push_back(sum * class_size);
		}
		
		//counting param_c is correct
		if (param_c[0] > param_c[1])
			decision = classes[0][0];
		else if (param_c[0] < param_c[1])
			decision = classes[1][0];
		else
		{
			srand(time(NULL));
			int number_of_index = rand() % classes.size();
			decision = classes[number_of_index][0];
		}

		if (decision == data[columns - 1][i])
			classification.push_back({decision, 1});
		else classification.push_back({ decision, 0 });
	}

	//Save decisions into file
	fstream file;
	file.open("dec_bayes.txt", ios::out);

	for (int i = 0; i < classification.size(); i++)
	{
		file << classification[i][0];
		if (i < rows-1)
			file << endl;
	}

	file.close();

	int number_of_true_classificated_objects = 0;

	for (int i = 0; i < classification.size(); i++)
	{
		if (classification[i][1] == 1)
			number_of_true_classificated_objects++;

		//cout << classification[i][0] << " " << classification[i][1] << endl;
	}

	double global_accuracy = (double)number_of_true_classificated_objects / (double)classification.size();

	int number_of_classificated_objects;
	double sum_accuracy = 0;

	for (int i = 0; i < classes.size(); i++)
	{
		number_of_true_classificated_objects = 0;
		number_of_classificated_objects = 0;

		for (int j = 0; j < classification.size(); j++)
		{
			if (classification[j][0] == classes[i][0] && classification[j][1] == 1)
				number_of_true_classificated_objects++;
			if (data[columns-1][j] == classes[i][0])
				number_of_classificated_objects++;
		}
		
		if (number_of_classificated_objects == 0)
			sum_accuracy += 0;
		else sum_accuracy += (double)number_of_true_classificated_objects / (double)number_of_classificated_objects;
	}

	double balanced_accuracy = sum_accuracy / (double)classes.size();

	cout << endl << "Global accuracy: " << global_accuracy << endl << endl;
	cout << "Balanced accuracy: " << balanced_accuracy << endl;

	fstream file2;
	file2.open("acc_bayes.txt", ios::out);

	file2 << "Global accuracy: " << global_accuracy << endl;
	file2 << "Balanced accuracy: " << balanced_accuracy;

	file2.close();
}