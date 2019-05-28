//random number functions form project 3

#include "randomNumberFunctions.hpp"

float Ranf(float low, float high) {
    float r = (float)rand();
    float t = r / (float)RAND_MAX;

    return low + t*(high-low);
}

//random int function
int Ranf(int ilow, int ihigh) {
    float low = (float)ilow;
    float high = ceil((float)ihigh);

    return (int)Ranf(low, high);
}

//time of day seed function
void TimeOfDaySeed() {
    struct tm y2k = { 0 };
    y2k.tm_hour = 0;   
    y2k.tm_min = 0; 
    y2k.tm_sec = 0;
    y2k.tm_year = 100; 
    y2k.tm_mon = 0; 
    y2k.tm_mday = 1;

    time_t  timer;
    time(&timer);
    double seconds = difftime(timer, mktime(&y2k));
    unsigned int seed = (unsigned int)(1000.*seconds); //msec
    srand(seed);
}
