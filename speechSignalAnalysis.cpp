//Author: Emily Garcia, Steven Groves, Israel Andrade
#include <iostream>
#include <vector>
#include <fstream> 
#include <cmath>
#include <cstdlib>
#include <iomanip>
using namespace std;

//Description: The values with the smallest percent difference are the standard deviations of the two data files. The values with the largest percent difference are the zeros crossing.
//Are there other statistical measures that you could suggest?
/* We can detect the voice frequency from both files to detect the if the voices match*/
//Can you determine if these sound recordings are from the same person?
/* We can determine if the sound recording came from the same person if we have a standard percent limit, but viewing the persent difference we may assume that
the vioce is not the same given that the percent difference is too great to determine*/

double mean(vector<double> v);
// Summary: Obtain the mean of a vector 
// Precondition: Mean function will add all elements in the array then divides by the number of elements in that vector
// Postcondition: Calculate the average of both files.

void avgPower(vector<double> &a, vector<double> &b, double &avePowerA, double &avePowerB);
// Summary: Obtain all the average powers 
// Precondition: Square each element in our vector then sum up all values. After obtaining the sum we divided by the number of elements in our vector.
// Postcondition: Obtain the average power.

double difference(double value_1, double value_2);
// Summary: Finds the difference of the arguments.
// Precondition: Function will be have two parameter 
// Postcondition: difference will output the percent difference between two values

double variance(vector<double> v);
// Summary: Determines the variance value in a given vector.
// Precondition: Takes in a vector.
// Postcondition: Returns the variance value in the vector.

double standardDeviation(double var);
// Summary: Determines the standard deviation value in a given vector.
// Precondition: Takes in a vector.
// Postcondition: Returns the standard deviation value in the vector
double zeroCross(vector<double> var);
// Summary: Function will count the number of times where values change from positive to negative or negative or positive
// Precondition: A vector will be passed in to the function
// Postcondition: Vector will be checked to view the zero crossings

double avgMagnitude (vector<double> v);
// Summary: Adds together the absolute values of the elements in the vector.
// Precondition: 
// Postcondition: Obtain the average magnitude.

int main()
{
    ifstream fin_1, fin_2;// input stream
    ofstream fout;// output stream
    fin_1.open("two_a.txt");// input files
    fin_2.open("two_b.txt");
    fout.open("comparison.txt");// output file
    
    double percentDifference;// declare variables
    double valueA, valueB, avePowerA, avePowerB;
    vector<double> a;// declare vectors
    vector<double> b;
    
    while(fin_1 >> valueA)
    {
        a.push_back(valueA);
    }
    while(fin_2 >> valueB)
    {
        b.push_back(valueB);
    }
    
    mean(a);// call mean()
    avgPower(a, b, avePowerA, avePowerB);// call avgPower()
    
    // output to file
    fout << "Team members: Emily Garcia, Steven Groves, Israel Andrade " << endl;
    fout << endl << endl;
    fout << setw(40) << "two_a.txt" << setw(20) << "two_b.txt" << setw(20) << "%" << "difference" << endl;
    fout << "Mean" << setw(37) << mean(a) << setw(20)  << mean(b) << setw(25)  << difference(mean(a), mean(b)) << endl;
    fout << "Standard"  << setw(31) << standardDeviation(variance(a))  << setw(19) << standardDeviation(variance(b))  << setw(28) << difference(standardDeviation(variance(a)), standardDeviation(variance(b)))<< endl;
    fout << "Variance" << setw(31) <<  variance(a) << setw(20) <<  variance(b) << setw(27) << difference(variance(a),variance(b)) << endl;
    fout << "Average Power" << setw(26) << avePowerA << setw(19)  << avePowerB << setw(27)  << difference(avePowerA, avePowerB) << endl;
    fout << "Average Magnitude" << setw(22) << avgMagnitude(a) << setw(20) << avgMagnitude(b) << setw(27) << difference(avgMagnitude(a), avgMagnitude(b)) << endl;
    fout << "Number of Zero Crossings" << setw(10) << zeroCross(a) << setw(19) << zeroCross(b) << setw(33) << difference(zeroCross(a), zeroCross(b)) << endl;
    

    
    
    fin_1.close();// close input stream
    fin_2.close(); //close second output file
    fout.close(); //close output file
    
}
double mean(vector<double> v)
{
    double sumOfA = 0, mean = 0;
    for(int i = 0; i < v.size(); i++)
    {
        sumOfA += v.at(i);
    }
    
    mean = sumOfA/v.size();
    return mean;
}

void avgPower(vector<double> &a, vector<double> &b, double &avePowerA, double &avePowerB)
{
    double sumPowerOfA = 0, sumPowerOfB = 0;
    
    for(int i = 0; i < a.size(); i++)
    {
        sumPowerOfA += pow(a.at(i), 2); 
    }
    
    for(int i = 0; i < b.size(); i++)
    {
        sumPowerOfB += pow(b.at(i), 2); 
    }
    avePowerA = sumPowerOfA / a.size();
    avePowerB = sumPowerOfB / b.size();
}

double difference(double value_1, double value_2)
{
    //% differece is mean1 - mean2 absoluta value divided my mean1+ mean2 divided total by 2 then multiply by 100
    //     |mean1 - mean2| / ((mean1 + mean2)/ 2) * 100
    double percentDifference = (fabs(value_1 - value_2) / ((value_1 + value_2) / 2.0)) * 100;
    return percentDifference;
}
double variance(vector<double> v)
{
   double mean1, mean2, var, size;
   mean1 = mean(v);
   size = v.size();
   for(int ix = 0; ix < size; ix++)
   {
      v[ix] = mean1 - v[ix];
      v[ix] = pow(v[ix],2);
   }
   
   var = mean(v);
   return var;
}

double standardDeviation(double var)
{
   return sqrt(var);  
}
double zeroCross(vector<double> v)
{
    double zeroCrossing = 0;
    for(int i = 0; i < v.size() - 1; i++)
    {
        if((v[i] > 0 && v[i+1] < 0) || (v[i] < 0 && v[i+1] > 0))
        {
            zeroCrossing++;   
        }   
    }    
    return zeroCrossing;
}

double avgMagnitude (vector<double> v)
{
   double total = 0, avg;
   int count = 0;
   for(int i = 0; i < v.size(); i++)
   {
      total += fabs(v[i]);
      count++;
   }
   avg = total / count;
   return avg;
}
