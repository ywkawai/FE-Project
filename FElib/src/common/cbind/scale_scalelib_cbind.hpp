#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

extern "C" {
    void CSCALElib_Init(const char* app_name);
    void CSCALElib_Final();
}

void SCALElib_Init(const std::string& app_name){
    CSCALElib_Init(app_name.data());
}

void SCALElib_Final(){
    CSCALElib_Final();
}