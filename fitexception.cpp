#include "fitexception.h"
FitException::FitException(std::string msg):exception(){
	m_msg="libapro_gen: "+msg;
}
const char* FitException::what() const throw(){
	return m_msg.c_str();
}
