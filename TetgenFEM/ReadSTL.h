#pragma once
#include "tetgen.h"
#include <iostream>
#include <fstream>
void readSTL(const std::string& filename, tetgenio& in);

void readOBJ(const std::string& filename, tetgenio& in);