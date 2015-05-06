#include <iostream>
#include "genetic_exception.h"
using namespace std;
function<void(string)> output=[](string){};
void SetGeneticErrorOutput(function<void(string)> func){
	output=func;
}
GeneticException::GeneticException(std::string msg):exception(){
	m_msg=msg;
	output(msg);
}
const char* GeneticException::what() const throw(){
	return m_msg.c_str();
}
