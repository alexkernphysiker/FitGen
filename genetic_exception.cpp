// this file is distributed under 
// GPL v 3.0 license
#include <iostream>
#include "genetic_exception.h"
using namespace std;
GeneticException::GeneticException(std::string msg):exception(){m_msg=msg;}
const char* GeneticException::what() const throw(){
	return m_msg.c_str();
}
