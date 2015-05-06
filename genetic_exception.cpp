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
	printf("\nThrowing GeneticException: %s\n",msg.c_str());
}
const char* GeneticException::what() const throw(){
	return m_msg.c_str();
}
