#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <fstream>
#include <string>
#include <cmath>
#include <queue>
#include <numeric>
#include <chrono>
#include <utility>
#include <unordered_set>
#include <iomanip>
#include <algorithm>
#include <iterator>

#include <thread>
#include "randchoice.hpp"

double randomFloat(){return rand() / (RAND_MAX+1.);}

std::ostream& operator<< (std::ostream &out, std::vector<int> const& data) {
    for(auto x : data)  out << x << " ";
    return out;
}

int sample(std::vector<double>cdf){
    int n = cdf.size();
    double x = randomFloat();
    if(x < cdf[0])      return 1;
    int low = 0, high = n-1;
    while(high - low > 1){
        int mid = (low + high) >> 1;
        if(x > cdf[mid])    low = mid;
        else                high = mid;
    }
    return high + 1;
}

std::unordered_set<int> pickSet(int N, int k, std::mt19937& gen)
{
    std::uniform_int_distribution<> dis(1, N);
    std::unordered_set<int> elems;

    while (elems.size() < k) {
        elems.insert(dis(gen));
    }

    return elems;
}

std::vector<int> pick(int N, int k) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::unordered_set<int> elems = pickSet(N, k, gen);

    // ok, now we have a set of k elements. but now
    // it's in a [unknown] deterministic order.
    // so we have to shuffle it:

    std::vector<int> result(elems.begin(), elems.end());
    std::shuffle(result.begin(), result.end(), gen);
    return result;
}