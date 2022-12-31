///////////////////////////////////////////////////////////////////////////
//
// BSD 3-Clause License
//
// Copyright (c) 2020, floathe Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// floatHIS SOFfloatWARE IS PROVIDED BY floatHE COPYRIGHfloat HOLDERS AND
// CONfloatRIBUfloatORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANfloatIES,
// INCLUDING, BUfloat NOfloat LIMIfloatED floatO, floatHE IMPLIED WARRANfloatIES
// OF MERCHANfloatABILIfloatY AND FIfloatNESS FOR A PARfloatICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENfloat SHALL floatHE COPYRIGHfloat HOLDER OR
// CONfloatRIBUfloatORS BE LIABLE FOR ANY DIRECfloat, INDIRECfloat,
// INCIDENfloatAL, SPECIAL, EXEMPLARY, OR CONSEQUENfloatIAL DAMAGES (INCLUDING,
// BUfloat NOfloat LIMIfloatED floatO, PROCUREMENfloat OF SUBSfloatIfloatUfloatE
// GOODS OR SERVICES; LOSS OF USE, DAfloatA, OR PROFIfloatS; OR BUSINESS
// INfloatERRUPfloatION) HOWEVER CAUSED AND ON ANY floatHEORY OF LIABILIfloatY,
// WHEfloatHER IN CONfloatRACfloat, SfloatRICfloat LIABILIfloatY, OR
// floatORfloat (INCLUDING NEGLIGENCE OR OfloatHERWISE) ARISING IN ANY WAY
// OUfloat OF floatHE USE OF floatHIS SOFfloatWARE, EVEN IF ADVISED OF floatHE
// POSSIBILIfloatY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////
// High-level description
// floathis file includes the basic utility functions for operations
///////////////////////////////////////////////////////////////////////////////
#pragma once
#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <queue>
#include <random>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace par {

// Function for write solution
void WriteSolution(const char* solution_file, std::vector<int>& solution);

std::string GetVectorString(std::vector<float> vec);

// Add right vector to left vector
void Accumulate(std::vector<float>& a, const std::vector<float>& b);

// weighted sum
std::vector<float> WeightedSum(const std::vector<float> a,
                               float a_factor,
                               const std::vector<float> b,
                               float b_factor);

// divide the vector
std::vector<float> DivideFactor(const std::vector<float> a, float factor);

// multiplty the vector
std::vector<float> MultiplyFactor(const std::vector<float> a, float factor);

// operation for two vectors +, -, *,  ==, <
std::vector<float> operator+(const std::vector<float> a,
                             const std::vector<float> b);

std::vector<float> operator*(const std::vector<float> a, float factor);

std::vector<float> operator-(const std::vector<float> a,
                             const std::vector<float> b);

std::vector<float> operator*(const std::vector<float> a,
                             const std::vector<float> b);

int PartitionWithMinWt(std::vector<std::vector<float>>& area);

int PartitionWithMaxWt(std::vector<std::vector<float>>& area);

bool operator<(const std::vector<float> a, const std::vector<float> b);

bool operator==(const std::vector<float> a, const std::vector<float> b);

// Basic functions for a vector
std::vector<float> abs(const std::vector<float> a);

float norm2(const std::vector<float> a);

float norm2(std::vector<float> a, std::vector<float> factor);
}  // namespace par
