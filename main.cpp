#include <iostream>
#include <map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <utility>
#include <string>
using namespace std;

// Define global variable
map<int, double> qvaluehat;
map<int, long double> bias;
map<int, long double> absolute_bias;
int alpha = 1, beta = 6;
double convergence = 0.1, Gamma = 0.1, budget = 50;

// Count total lines in a file
int CountLines(string filename)
{
    ifstream ReadFile;
    int n = 0;
    string tmp;
    ReadFile.open(filename.c_str()); //ios::in read oly
    if (ReadFile.fail())             // Fail to open the file
    {
        return 0;
    }
    else
    {
        while (getline(ReadFile, tmp, '\n'))
        {
            n++;
        }
        ReadFile.close();
        return n;
    }
}

string ReadLine(string filename, int line, int lines)
{
    int i = 0;
    string temp;
    fstream file;
    file.open(filename.c_str());
    //lines = CountLines(filename);

    if (line <= 0)
    {
        return "Error 1: line should be positive";
    }
    if (file.fail())
    {
        return "Error 2: file none exist";
    }
    if (line > lines)
    {
        return "Error 3: line out of limit";
    }
    while (getline(file, temp) && i < line - 1)
    {
        i++;
    }
    file.close();
    return temp;
}

int timetoseconds(string Time) // convert the time in the dataset to seconds to compare
{
    stringstream ss(Time);
    string tmp;
    getline(ss, tmp, ':');
    string hour = tmp; //cout<<hour<<endl;
    int h = atoi(hour.c_str());
    getline(ss, tmp, ':');
    string minute = tmp; //cout<<minute<<endl;
    int m = atoi(minute.c_str());
    getline(ss, tmp, '.');
    string second = tmp; //cout<<second<<endl;
    int s = atoi(second.c_str());
    return h * 3600 + m * 60 + s;
}

int getsubdata(string data, int &taxi_ID, string &Date, int &time, long double &latitude, long double &longitude, long double &temperature)
{
    stringstream ss(data);
    string tmp, Taxi_ID, Time, Latitude, Longitude, Temperature;
    getline(ss, tmp, ',');
    Taxi_ID = tmp;
    //cout << Taxi_ID << endl;
    taxi_ID = atoi(Taxi_ID.c_str());
    //cout << taxi_ID << endl;
    getline(ss, tmp, ',');
    Date = tmp;
    //cout << Date << endl;
    getline(ss, tmp, ',');
    Time = tmp;
    //cout << Time << endl;
    time = timetoseconds(Time);
    //cout << time << endl;
    getline(ss, tmp, ',');
    Latitude = tmp;
    //cout << Latitude << endl;
    latitude = atof(Latitude.c_str());
    //cout << setprecision(15) << latitude << endl;
    getline(ss, tmp, ',');
    Longitude = tmp;
    //cout << Longitude << endl;
    longitude = atof(Longitude.c_str());
    // cout << setprecision(15) << longitude << endl;
    getline(ss, tmp, ',');
    Temperature = tmp;
    //cout << Temperature << endl;
    temperature = atof(Temperature.c_str());
    //cout << setprecision(15) << temperature << endl;
    return 1;
}

// calculate the latitude and longitude range according to the region radius in km
int getradius(long double Latitude, long double radiusinkm, long double &lati, long double &longi)
{
    lati = radiusinkm / 111.574;
    long double x = Latitude * 3.14159 / 180;
    longi = radiusinkm / (111.574 * cos(x));
    // cout << setprecision(15) << lati << endl;
    // cout << setprecision(15) << longi << endl;
    return 1;
}

long double truthdiscovery(vector<pair<int, long double> > observations, vector<pair<int, long double> > &qvalue_true, vector<pair<int, long double> > &delta)
{
    /* Initialize */
    long double o_tao = 0, o_tao_2 = 0;
    for (auto it = observations.begin(); it != observations.end(); it++)
    {
        qvalue_true.push_back(make_pair(it->first, 0)); //initial is 0
        delta.push_back(make_pair(it->first, 0));       //initial is 0
    }
    int E = qvalue_true.size();
    do
    {
        long double sum2 = 0;
        long double sum3 = 0;
        for (int i = 0; i < qvalue_true.size(); i++)
        {
            long double sum1 = 0;
            for (auto it = observations.begin(); it != observations.end(); it++)
            {
                sum1 += pow((it->second - o_tao), 2);
            }
            qvalue_true[i].second = log(sum1 / pow(observations[i].second - o_tao, 2));
            sum2 += qvalue_true[i].second * observations[i].second;
            sum3 += qvalue_true[i].second;
        }
        o_tao_2 = o_tao;
        o_tao = sum2 / sum3;
    } while (fabs(o_tao - o_tao_2) > convergence);

    long double sum4 = 0;
    for (int i = 0; i < qvalue_true.size(); i++)
    {
        sum4 += qvalue_true[i].second;
    }
    for (int i = 0; i < qvalue_true.size(); i++)
    {
        delta[i].second = (1 - exp((-beta) * E * (qvalue_true[i].second / sum4 - 1 / E))) / (1 + exp((-beta) * E * (qvalue_true[i].second / sum4 - 1 / E)));
        qvalue_true[i].second = exp(-alpha * fabs(observations[i].second - o_tao));
    }
    return o_tao;
}

long double getDeviation(vector<pair<int, long double> > observations)
{
    long double sum = 0, mean = 0, dev = 0;
    int n = observations.size();
    for (auto it = observations.begin(); it != observations.end(); it++)
    {
        sum += it->second;
    }
    mean = sum / n;
    for (auto it = observations.begin(); it != observations.end(); it++)
    {
        dev += pow(it->second - mean, 2);
    }
    return sqrt(dev / n);
}
int qualityupdate(long double observation_truth, vector<pair<int, long double> > observations, vector<pair<int, long double> > &qvalue_true, vector<pair<int, long double> > &delta)
{

    for (int i = 0; i < delta.size(); i++)
    {
        map<int, long double>::iterator iter_map;
        iter_map = bias.find(delta[i].first);
        long double bias_2 = iter_map->second;
        if (delta[i].second >= 0)
        {
            long double temp = fabs(delta[i].second) * bias_2 + (1 - fabs(delta[i].second)) * (observations[i].second - observation_truth);
            bias[iter_map->first] = temp;
        }
        else
        {
            long double temp = fabs(delta[i].second) * bias_2 + (1 - fabs(delta[i].second)) * (observations[i].second - observation_truth);
            bias[iter_map->first] = temp;
        }
        long double difference;
        iter_map = absolute_bias.find(delta[i].first);
        difference = iter_map->second;
        absolute_bias[iter_map->first] = observations[i].second - observation_truth;
        map<int, double>::iterator iter_map2;
        long double temp = exp((-alpha) * fabs((1 - Gamma) * difference + Gamma * (observations[i].second - observation_truth)));
        iter_map2 = qvaluehat.find(delta[i].first);
        qvaluehat[iter_map2->first] = temp;
    }
    return 1;
}

template <typename T>
void printVectorElements(vector<T> &vec)
{
    for (auto i = 0; i < vec.size(); ++i)
    {
        cout << "(" << vec.at(i).first << ","
             << vec.at(i).second << ")"
             << "; ";
    }
    cout << endl;
}

/*****************************************************************************
 * Function Name: allocation
 * Description: Algorithm 1, quality-driven sensing duty allocation
 * @ Q[]:  data quality indications
 * @ budget:  the budget of current task
 * @ employee_id[]: the set of employees recruited by task
 * @ duty[]: duty allocations of
 *
 * Return
 * @ vector<pair<int, double>>
 */

bool sortbysec(const pair<int, double> &a,
               const pair<int, double> &b)
{
    return (a.second > b.second);
}

vector<pair<int, double> > allocation(vector<pair<int, long double> > temperature, vector<pair<int, long double> > &observations)
{
    vector<pair<int, long double> >::iterator iter_temp;
    map<int, double>::iterator iter_map;
    vector<pair<int, double> > vtMap;
    for (auto it = temperature.begin(); it != temperature.end(); it++)
    {
        iter_map = qvaluehat.find(it->first);
        vtMap.push_back(make_pair(it->first, iter_map->second));
    }
    sort(vtMap.begin(), vtMap.end(), sortbysec);

    for (auto it = vtMap.begin(); it != vtMap.end(); it++)
        cout << it->first << ':' << it->second << '\n';

    vector<pair<int, double> >::iterator iter;
    iter = vtMap.begin();

    vector<pair<int, double> > employee_id;
    employee_id.push_back(make_pair(iter->first, iter->second));
    cout << "employee_id:";
    printVectorElements(employee_id);
    cout << endl;
    vector<pair<int, double> >::iterator _iter;
    for (iter = vtMap.begin() + 1; iter != vtMap.end(); iter++)
    {
        double _q = 0;
        for (_iter = employee_id.begin(); _iter != employee_id.end(); _iter++)
        {
            _q = _q + pow(_iter->second, -2);
        }
        _q = _q + pow(iter->second, -2);
        double div = 0;
        div = (double)employee_id.size() / _q;
        div = sqrt(div);

        if (iter->second > div)
        {
            employee_id.push_back(make_pair(iter->first, iter->second));
            for (auto it = temperature.begin(); it != temperature.end(); it++)
            {
                if (it->first == iter->first)
                {
                    map<int, long double>::iterator iter_map;
                    iter_map = bias.find(it->first);
                    observations.push_back(make_pair(it->first, (it->second) + (iter_map->second)));
                }
            }
        }
        cout << "q:" << iter->second << ";  "
             << "div:" << div << endl;
    }
    cout << "employee_id:";
    printVectorElements(employee_id);
    cout << endl;

    // duty allocation
    vector<pair<int, double> > duty_allocation;
    // vector<pair<int, double>>::iterator iter_q;

    vector<pair<int, double> >::iterator iter_duty;
    for (_iter = employee_id.begin(); _iter != employee_id.end(); _iter++)
    {
        double sum = 0;
        for (iter_duty = employee_id.begin(); iter_duty != employee_id.end(); iter_duty++)
        {
            sum = sum + pow(iter_duty->second, -2);
        }
        double left = 0;
        left = ((employee_id.size() - 1) * budget) / (_iter->second * sum);
        double right = 0;
        right = (employee_id.size() - 1) / (_iter->second * _iter->second * sum);
        right = 1 - right;
        double duty = 0;
        duty = left * right;
        cout << duty << endl;
        duty_allocation.push_back(make_pair(_iter->first, duty));
    }

    return duty_allocation;
}

int main()
{

    /* Read data from file*/
    string filename = "crowd_temperature.txt";
    int lines = CountLines(filename);
    int task = 0;
    ofstream result;
    result.open("result.csv", ios::out);
    result << "Task"
           << ","
           << "Observation truth"
           << ","
           << "standard deviation" << endl;

    ofstream Agent135;
    Agent135.open("Agent135.csv", ios::out);
    Agent135 << "AgebtID"
             << ","
             << "qvalue"
             << ","
             << "profit" << endl;
    ofstream Agent352;
    Agent352.open("Agent352.csv", ios::out);
    Agent352 << "AgebtID"
             << ","
             << "qvalue"
             << ","
             << "profit" << endl;
    ofstream Agent8;
    Agent8.open("Agent8.csv", ios::out);
    Agent8 << "Agent8"
           << ","
           << "qvalue"
           << ","
           << "profit" << endl;
    for (int i = 0; i < 368; i++) //368 agents in total
    {
        qvaluehat.insert(pair<int, double>(i + 1, 1.0)); //initial q value is 1
        bias.insert(pair<int, long double>(i + 1, 0));   //initial q value is 1
        absolute_bias.insert(pair<int, long double>(i + 1, 0));
    }

    map<int, double>::iterator itr;
    /*
    for (itr = qvaluehat.begin(); itr != qvaluehat.end(); ++itr)
    {
        cout << '\t' << itr->first
             << '\t' << itr->second << '\n';
    }
    */
    int count = 1000; // to control the number of tasks
    while (count-- > 0)
    {
        int line = (rand() % (lines - 7)) + 2; // randomly pick an entry in the dataset as a task location
        //cout<<lines<<endl;

        // get the data of the picked entry
        string data = ReadLine(filename, line, lines);
        //cout<<"attempt:";
        //cout << data << endl;
        int Taxi_ID;
        string Date;
        int Time;
        long double pivot_Latitude, lati;
        long double pivot_Longitude, longi;
        long double Temperature;

        // save the subdata in the according variable
        getsubdata(data, Taxi_ID, Date, Time, pivot_Latitude, pivot_Longitude, Temperature);

        // pair the agentID with temperature for future use
        vector<pair<int, long double> > vtMap;
        vtMap.push_back(make_pair(Taxi_ID, Temperature));
        //cout << Time << endl;

        /* calculate the range of latitude and longtitude according to the task's location*/
        getradius(pivot_Latitude, 2.0, lati, longi);

        /* Pick qualified users*/
        int l = 1, r = 1, users = 0, time = Time;
        while ((line - l) > 1) //dataset is sorted by time, increasing
        {
            string data = ReadLine(filename, line - l, lines);
            int Taxi_ID;
            string Date;
            int Time;
            long double Latitude;
            long double Longitude;
            long double Temperature;
            getsubdata(data, Taxi_ID, Date, Time, Latitude, Longitude, Temperature);

            if (Time < (time - 60)) // time is unqualified
            {
                break;
            }
            else
            {
                l++;
                if ((pivot_Latitude + lati) >= Latitude && (pivot_Latitude - lati) <= Latitude)
                {
                    if ((pivot_Longitude + longi) >= Longitude && (pivot_Longitude - longi) <= Longitude)
                    {
                        users++;
                        vtMap.push_back(make_pair(Taxi_ID, Temperature)); // save qualified users' ID and temperature for future use
                        //cout << "Got one!  ";
                        //cout << data << endl;
                    }
                }
            }
        }
        while ((line + r) < lines) //dataset is sorted by time, increasing
        {
            string data = ReadLine(filename, line + r, lines);
            int Taxi_ID;
            string Date;
            int Time;
            long double Latitude;
            long double Longitude;
            long double Temperature;

            getsubdata(data, Taxi_ID, Date, Time, Latitude, Longitude, Temperature);
            if (Time > (time + 60))
            {
                break;
            }
            else
            {
                r++;
                if ((pivot_Latitude + lati) <= Latitude && (pivot_Latitude - lati) >= Latitude)
                {
                    if ((pivot_Longitude + longi) <= Longitude && (pivot_Longitude - longi) >= Longitude)
                    {
                        users++;
                        vtMap.push_back(make_pair(Taxi_ID, Temperature));
                        //cout << "Got one!  ";
                        //cout << data << endl;
                    }
                }
            }
        }
        if (users >= 4)
        {
            cout << users << endl;
            cout << "users:    ";
            for (auto it = vtMap.begin(); it != vtMap.end(); it++)
                cout << it->first << ':' << it->second << '\n';
            vector<pair<int, long double> > observations;
            vector<pair<int, long double> > qvalue_true;
            vector<pair<int, long double> > delta;
            vector<pair<int, double> > picked_agents = allocation(vtMap, observations);
            if (picked_agents.size() < 4)
            {
                continue;
            }
            task++;
            cout << "picked users:    ";
            for (auto it = observations.begin(); it != observations.end(); it++)
            {
                cout << it->first << ':' << it->second << '\n';
            }
            long double observation_truth = truthdiscovery(observations, qvalue_true, delta);
            qualityupdate(observation_truth, observations, qvalue_true, delta);
            long double deviation = getDeviation(observations);
            result << task << "," << observation_truth << "," << deviation << endl;
            for (auto it = picked_agents.begin(); it != picked_agents.end(); it++)
            {
                map<int, double>::iterator iter_map;
                map<int, long double>::iterator iter_map2;
                if (it->first == 135)
                {
                    iter_map = qvaluehat.find(it->first);
                    Agent135 << iter_map->first << "," << iter_map->second << endl;
                }
                if (it->first == 352)
                {
                    iter_map = qvaluehat.find(it->first);
                    Agent352 << iter_map->first << "," << iter_map->second << endl;
                }
                if (it->first == 8)
                {
                    iter_map = qvaluehat.find(it->first);
                    Agent8 << iter_map->first << "," << iter_map->second << endl;
                }
                iter_map = qvaluehat.find(it->first);
                cout << "qvalue:" << endl;
                cout << iter_map->first << ":" << iter_map->second;
                iter_map2 = bias.find(it->first);
                cout << "bias:" << endl;
                cout << iter_map2->first << ":" << iter_map2->second;
                cout << endl;
            }
        }
    }

    result.close();
    Agent135.close();
    Agent352.close();
    Agent8.close();
    return 1;
}