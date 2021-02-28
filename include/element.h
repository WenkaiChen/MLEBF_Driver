#ifndef MBE_DRIVER_ELEMENT_H
#define MBE_DRIVER_ELEMENT_H

#include <map>
#include <string>

using namespace std;
static map<int, string> atmnb2ele={{1,"H"},    {2,"He"},
                            {3,"Li"},   {4,"Be"},  {5,"B"},   {6,"C"},   {7,"N"},   {8,"O"},   {9,"F"},   {10,"Ne"},
                            {11,"Na"},  {12,"Mg"}, {13,"Al"}, {14,"Si"}, {15,"P"},  {16,"S"},  {17,"Cl"}, {18,"Ar"},
                            {19,"K"},   {20,"Ca"}};
static map<string, int> ele2atmnb={{"H",1},    {"He",2},
                            {"Li",3},   {"Be",4},  {"B",5},   {"C",6},   {"N",7},   {"O",8},   {"F",9},   {"Ne",10},
                            {"Na",11},  {"Mg",12}, {"Al",13}, {"Si",14}, {"P",15},  {"S",16},  {"Cl",17}, {"Ar",18},
                            {"K",19},   {"Ca",20}};
// ele2mass in AMU
static map<int, double> ele2mass={{1,1.008},   {2,4.003},
                           {3,6.938},   {4,9.012},   {5,10.806},  {6,12.010},  {7,14.006},  {8,15.999},  {9,18.998},  {10,20.180},
                           {11,22.990}, {12,24.305}, {13,26.981}, {14,28.084}, {15,30.974}, {16,32.059}, {17,35.446}, {18,39.948},
                           {19,39.098}, {20,40.078}};




#endif //MBE_DRIVER_ELEMENT_H
