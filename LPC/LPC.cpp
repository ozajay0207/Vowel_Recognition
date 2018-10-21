// LPC.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include "iostream"
#include "fstream"
#include "string"
#define FRAME_SIZE 320
#define FRAMES 5

const int total_window = (320 + (FRAMES - 1) * 80) + 1;
int frame_cutting_index = 1;
int start_marker = 0;
int end_marker = 0;

using namespace std;

//Used Files
char* input_file = "input.txt";
char* reference_file[5] = { "ref_file_prime_a.txt", "ref_file_prime_e.txt", "ref_file_prime_i.txt", "ref_file_prime_o.txt", "ref_file_prime_u.txt" };
char* normalized_file = "Normalized.txt";
char* silence_file = "silence.txt";
char* trimmed_file = "trim.txt";
char* ri_file = "ri_file.txt";
char* ai_file = "ai_file.txt";
char* ci_file = "ci_file.txt";
char* c_prime_file = "c_prime.txt";
char* universe_file = "universe_vowel.txt";
char* hamming_file = "Hamming_window.txt";

//I/O Pointer
ifstream in, in1;
ofstream out,out1;

float tokhura_weight[12] = { 1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61 };
float tokhura_dist[5] = { 0 };
const int p = 12;
int flag = 0;
int max_sample_index = 0;
int remove_header_counter = 5;
double hamming_window[320] = { 0 };
long double current_value;
long double sum_samples = 0;
double dc_shift_value = 0;
double normalization_ratio = 1;
long int zcr_count = 0;
long int no_of_samples = 0;
long int max_sample_value = 1;
long double sample_array[total_window] = { 0 };
long double r[12] = { 0 }, k[12] = { 0 }, alpha[13][13] = { 0 }, E[13] = { 0 }, a[13] = { 0 }, c[13] = { 0 }, c_prime[13] = { 0 }, w[12] = { 0 };


//To Remove header from the input files
void remove_header(ifstream& in){
	if (flag){
		remove_header_counter = 5;
		string temp = "";
		while (getline(in, temp) && remove_header_counter){
			remove_header_counter--;
		}
		flag = 0;
	}
	else{
		remove_header_counter = 5;
		string temp = "";
		while (getline(in, temp) && remove_header_counter){
			remove_header_counter--;
		}
	}
}


//To calculate DC Shift
void calculate_dc_shift(){
	in.open(silence_file);
	flag = 1;
	cout << "\n................Calculating DC shift..................." << endl;
	remove_header(in);
	string temp = "";

	while (in >> temp && !remove_header_counter){
		current_value = stod(temp);
		sum_samples += current_value;
		no_of_samples++;
	}

	dc_shift_value = sum_samples / no_of_samples;
	cout << "DC Shift value:" << dc_shift_value << endl;
	in.close();
}

//To calculate normalization ratio
void calculate_normalization_ratio(){
	int index_count = 0;
	cout << "\n............Reading from " << input_file << "................" << endl;
	in.open(input_file);
	remove_header(in);
	string temp = "";
	cout << "\n............Calculating Normalization Ratio................" << endl;
	while (in >> temp && !remove_header_counter){
		index_count++;
		current_value = stod(temp);

		//Saving maximum index value and maximum sample value
		if (current_value > max_sample_value){
			max_sample_value = current_value;
			max_sample_index = index_count;
		}
	}
	cout << "Max Sample value:" << max_sample_value << endl;

	//Calculating normalization ratio
	normalization_ratio = 5000.0 / max_sample_value;
	cout << "Normalization ratio:" << normalization_ratio << endl;
	in.close();
}

//To get Hamming Windows values to array from the pre-calculated File
void get_hamming_window(){
	int index_count = 0;
	in.open(hamming_file);
	string temp = "";
	while (in >> temp){
		hamming_window[index_count++] = stod(temp);
	}
	in.close();
}

//To remove DC Shift and normalize
void dc_normalize(){

	//Initializing Markers as + or - 640 from the max sample value
	start_marker = max_sample_index - (total_window / 2);
	end_marker = max_sample_index + (total_window / 2);
	cout << "Start Marker:" << start_marker << endl;
	cout << "End Marker:" << end_marker << endl;
	cout << "Total Window:" << total_window << endl;

	int index_count = 0;
	int arr_index = 0;
	in.open(input_file);
	remove_header(in);
	string temp = "";
	cout << "\n........Removing DC shift and Normalizing File.........." << endl;
	out.open(normalized_file);

	//Subtracting DC shift and Multiplying Normalization ratio
	while (in >> temp && !remove_header_counter){
		index_count++;
		current_value = stod(temp);
		current_value = current_value - dc_shift_value;
		current_value = current_value * normalization_ratio;

		//Making the Array of 5 Frames using start and end markers
		if (index_count >= start_marker && index_count < end_marker)
			sample_array[arr_index++] = current_value;

		//Writing the Normalized values to file "normalized.txt"
		out << to_string(current_value) << endl;
	}
	out.close();
	in.close();

}

//Calculating Ri's
void calculate_Ris(){
	cout << "\n........Writing Ri's for 5 Frames to file.........." << endl;
	int count = 0;
	long double first_value = 0;
	long double second_value = 0;
	string temp;

	//To test file for single frame
	/*in.open("sample_output.txt");
	while (in >> temp){
	if (count<total_window)
	sample_array[count++] = stod(temp);
	}

	in.close();	*/
	out.open(ri_file);

	//Calculating Ri and also multiplying with hamming window appropriately
	for (int k = 0; k < FRAMES; k++){
		for (int j = 0; j <= p; j++){
			for (int i = 0; i < FRAME_SIZE - j; i++){
				first_value = sample_array[i + (80 * k)];
				second_value = sample_array[(i + j) + (80 * k)];
				r[j] = r[j] + ((first_value*hamming_window[i])*(second_value*hamming_window[i + j]));
			}
		}		
		for (int j = 0; j <= p; j++){
			//printf("R[%d] : %Lf \n", j, r[j]);
			out << fixed << r[j] << " ";
			r[j] = 0;
		}
		out << endl;
		//cout << endl;
	}
	out.close();
	cout << endl;
}


//Calculate Ai's using Levenson-Durbin
void levenson_durbin(long double r[]){
	int i, j;
	long double summation = 0;

	E[0] = r[0];
	for (i = 1; i <= p; i++){
		summation = 0.0;
		for (j = 1; j <= i - 1; j++){
			summation += alpha[j][i - 1] * r[i - j];
		}
		k[i] = (r[i] - summation) / E[i - 1];
		alpha[i][i] = k[i];
		for (j = 1; j <= i - 1; j++){
			alpha[j][i] = alpha[j][i - 1] - (k[i] * alpha[i - j][i - 1]);
		}
		E[i] = (1 - (k[i] * k[i]))*E[i - 1];
	}
	//cout << endl;

	//Calculating the Ai's as the last last column of index i.e. 12th column of matrix alpha
	out.open(ai_file, std::ios_base::app);
	a[0] = 0.0;
	out << a[0] << " ";
	for (int i = 1; i <= p; i++){
		a[i] = alpha[i][12];
		//printf("\nA[%d] : %Lf ", i, a[i]);
		out << fixed << a[i] << " ";
	}
	out << endl;
	out.close();
}

//calling levenson Durbin for 5 frames one after the another
void calculate_Ais(){
	cout << "\n........Writing Ai's for 5 Frames to file.........." << endl;
	int abc = 1;
	string temp;
	int count = 0;
	in.open(ri_file);
	out.open(ai_file);
	out.close();
	for (int k = 0; k < FRAMES; k++){
		while (count <= 12){
			in >> temp;
			r[count] = stod(temp);
			count++;
		}
		levenson_durbin(r);
		count = 0;
	}
	in.close();
}

//Calculate Ci's
void calculate_Cis(){
	cout << "\n........Writing Ci's for 5 Frames to files.........." << endl;
	string temp;
	long double summation = 0;
	int count = 0;
	long double a1 = 0.0, a2 = 0.0;
	in.open(ai_file);
	in1.open(ri_file);
	out.open(ci_file);
	out.close();

	//Calculating Cis by reading files of Ais and Ris
	out.open(ci_file, std::ios_base::app);
	for (int k = 0; k < FRAMES; k++){
		while (count <= 12){
			in >> temp;
			//cout << temp << endl;
			a[count] = stod(temp);
			in1 >> temp;
			r[count] = stod(temp);
			count++;
		}
		c[0] = 2 * log(r[0]) / log(2.0);
		for (int m = 1; m <= p; m++){
			summation = 0;
			for (int k = 1; k <= m - 1; k++)
			{
				summation += ((double)k / m)*c[k] * a[m - k];
			}
			c[m] = a[m] + summation;
		}
		for (int m = 1; m <= p; m++){
			out << fixed << c[m] << " ";
			//printf("\nC[%d] : %Lf ", m, c[m]);
			c[m] = 0;
		}
		//cout << endl;
		out << endl;
		count = 0;
	}
	out.close();
	in1.close();
	in.close();
}

//To calculate C Prime(Raised Sine wave)
void calculate_c_prime(){
	cout << "\n........Writing Ci Prime for 5 Frames to file.........." << endl;
	string temp;
	int count = 0;
	in.open(ci_file);
	out.open(c_prime_file);	
	out.close();

	//out1.open(universe_file, std::ios_base::app);

	//Calculting CiPrimes by reading Ci values from File
	out.open(c_prime_file, std::ios_base::app);
	for (int k = 0; k < FRAMES; k++){
		while (count < 12){
			in >> temp;
			c[count] = stod(temp);
			count++;
		}
		for (int m = 1; m <= p; m++){
			w[m] = 1.0 + 6.0 * sin((22.0 / 7.0)*m / 12.0);
			c_prime[m] = c[m - 1] * w[m];
		}
		for (int m = 1; m <= p; m++){
			//	printf("\nC_Prime[%d] : %Lf ", m, c_prime[m]);
			out << c_prime[m] << " ";
			//out1 << c_prime[m] << " ";
			c[m] = 0;
		}
		//cout << endl;
		out << endl;
		//out1 << endl;
		count = 0;
	}
	in.close();
	out.close();

	//out1.close();
}

//Calculating Tokhura's Distance Using Reference Files for vowels
void calculate_tokhura_distance(){
	int count = 0, j = 0, min_index = 0;
	float ref_sample, input_sample, min = 9999;
	float sum[5] = { 0 }, grand_sum = 0;
	string temp, temp1;
	for (int i = 0; i < 5; i++){
		in.open(reference_file[i]);
		in1.open(c_prime_file);
		//cout << "----New File----" << endl;

		while (in >> temp && in1 >> temp1){
			count++;
			ref_sample = stof(temp);
			input_sample = stof(temp1);
			sum[j] += ((ref_sample - input_sample)*(ref_sample - input_sample))*tokhura_weight[count - 1];
			if (count == 12){
				grand_sum += sum[j];
				//cout << "Grand_Sum : " << grand_sum << endl;
				sum[j] = 0;
				j++;
				count = 0;
			}
			j = 0;
		}
		tokhura_dist[i] = grand_sum / 5;
		grand_sum = 0;
		in1.close();
		in.close();
	}
	cout << "\n\nTokhura Distances are in order of A E I O U" << endl;
	for (int i = 0; i < 5; i++){
		printf("\nTokhura :%.6f", tokhura_dist[i]);
		if (tokhura_dist[i] < min){
			min = tokhura_dist[i];
			min_index = i;
		}
	}
	if (min_index == 0)
		cout << "\n\nVowel is A" << endl;
	else if (min_index == 1)
		cout << "\n\nVowel is E" << endl;
	else if (min_index == 2)
		cout << "\n\nVowel is I" << endl;
	else if (min_index == 3)
		cout << "\n\nVowel is O" << endl;
	else
		cout << "\n\nVowel is U" << endl;

	cout << endl;
}

int _tmain(int argc, _TCHAR* argv[])
{
	//init_file_names();

	cout << ".....Recording will be for 3 seconds......" << endl;
	cout << "....... Recording Silence......" << endl;
	system("Recording_Module.exe 3 silence.wav silence.txt");
	cout << "\nSilence recorded. **Press Enter to record your VOWEL**" << endl;
	system("Recording_Module.exe 3 input.wav input.txt");
	cout << "\nRecording successfull **Press ENTER to proceed with program**" << endl;

	cout << "\nReading Input from : " << input_file << endl;

	
	get_hamming_window();
	calculate_dc_shift();
	calculate_normalization_ratio();
	dc_normalize();


	calculate_Ris();
	calculate_Ais();
	calculate_Cis();
	calculate_c_prime();

	calculate_tokhura_distance();

	return 0;
}

