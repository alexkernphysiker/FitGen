#include <iostream>
#include "fitexception.h"
FitException::FitException(std::string msg):exception(){
	m_msg="libapro_gen: "+msg;
	printf("\nThrowing FitException: %s\n",msg.c_str());
}
const char* FitException::what() const throw(){
	return m_msg.c_str();
}
