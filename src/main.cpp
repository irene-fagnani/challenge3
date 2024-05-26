#include "utils.hpp"
#include "json.hpp"
#include "muparser_fun.hpp"
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <iostream>
using json=nlohmann::json;

int main(int argc, char** argv[]){

std::ifstream file("data.json");
json data = json::parse(file);

MuparserFun f(data.value("f",""));


}